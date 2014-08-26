(function(){

  var graphSpecDefaults = {
	// TODO: determine me dynamically
	leftPadding: 120,
	
	axisLabelFontsize: 11,
	axisFontcolor: "#444",
	axisTitleFontsize: 20,
	borderColor: "#444",
	rightPadding: 25,    
	profileHeight: 100,
	profileTicks: 5,
	
	// header labels and arrows
	headerSizes:{'extract': 30, 'hits': 15},
	headerLabelPadding: 10,    
	headerLabelFontcolor: "#444",
	headerLabelFontsize: 11,
	headerArrowHeight: 12,
	headerArrowWidth: 6,
	headerArrowFontsize: 9,
	axisTitle: '% GC'
	
      },
      lineColors = {
	r: "red",
	g: "blue",
	b: "green"
      },
      colorBlindLineColors = {
	r: "rgb(213, 94, 0)",
	g: "rgb(204, 121, 167)",
	b: "rgb(0, 114, 178)"
      };
  
  // the `GraphDealer` hands out graph data and processes events
  function GraphDealer($log, Utils, $q, $rootScope){

    var hasProfileData = false,
	_onProfileData = $q.defer(),
	onProfileData = _onProfileData.promise,
	_onWidth = $q.defer(),
	onWidth = _onWidth.promise,
	pendingRedraws = 0,
	opts = {
	  basesPerGraph: 10000,
	  colorBlindFriendly: false,
	  page: 0,
	  graphsPerPage: 5,
	  length: 0,
	  profile: null,
	  extracts: {},
	  graphSpecs: []
	};

    // once we have data, load it and rebuild the graphs
    onProfileData.then(loadProfileData).then(rebuildGraphs);


    // public interface
    return {
      opts:opts,
      setProfile:function(profile){
	$log.log('setProfile', profile.length);
	_onProfileData.resolve(profile);
	return onProfileData;
      },

      addExtract:function(name, data){
	$log.log('setExtract', name,  data.length);
	// save it for later
	opts.extracts[name] = data;
	// if we've got profile data already, rebuild
	if(hasProfileData){
	  // TODO: maybe a simpler way to add extracts to existing
	  // graphs?
	  rebuildGraphs();
	}
	// if we're still waiting on profile data, then when it hits
	// this extract will be processed.
      },

      setColors:function(colorBlindFriendly){
	$log.log('setColors', arguments);
	opts.colorBlindFriendly = colorBlindFriendly;
	// TODO: lighter change here; see if we can update colors
	// in-place without a full rebuild
	rebuildGraphs();
      },
      setZoom:function(basesPerGraph){
	$log.log('setZoom', arguments);
	opts.basesPerGraph = basesPerGraph;
	rebuildGraphs();
      },
      zoomTo: function(){
	// TODO: figure out a good API for this
      },
      panTo: function(oldStartBase, newStartBase){
	// TODO: find the graph with the old start base
	// TODO: recalculate graph start/ends
	// TODO: repartition data
	// TODO: redraw graphs with new data
      },
      hasNextPage: function(){
	return opts.page < maxPages();
      },
      nextPage: function(){
	$log.log('nextPage', arguments);
	opts.page++;
	rebuildGraphs();
      },
      hasPreviousPage: function(){ return opts.page > 0;},
      previousPage: function(){
	$log.log('previousPage', arguments);
	opts.page--;
	rebuildGraphs();
      },
      setWidth: function(w){
	$log.log('setting graph width to', w);
	// let folks know it's ready
	_onWidth.resolve(w);

	// TODO: handle changing this after the fact
      }
    };

    function makeGraphSpec(startBase, width){      
      var endBase = startBase + opts.basesPerGraph,
	  // TODO: reduce duplication between here and
	  // `GraphingCalculator.stops`
	  stopInterval = Utils.orderOfMagnitude(opts.basesPerGraph, -1),
	  spec = {
	    range: [startBase, endBase],
	    // get some margin in here
	    startBase: startBase - stopInterval,
	    endBase: endBase + stopInterval,
	    extracts: {},
	    profile: [],
	    headers: [],
	    width: width,
	    colors: opts.colorBlindFriendly ? colorBlindLineColors : lineColors
	  };
      return angular.extend(spec, graphSpecDefaults);
    }
    
    function makeGraphSpecs(width){
      var base = opts.page * opts.basesPerGraph * opts.graphsPerPage;
      return _.range(0, opts.graphsPerPage)
	.map(function(n){
	  return makeGraphSpec(base + n*opts.basesPerGraph, width);
	});
    }


    function loadProfileData(profile){
      hasProfileData = true;
      opts.profile = profile;
      var start = profile[0],
	  end = _.last(profile);

      opts.startBase = start.coordinate;
      opts.endBase = end.coordinate;
      opts.length = opts.endBase - opts.startBase;
      // find a sensible zoom level
      var basesPerGraph = opts.length / opts.graphsPerPage;
      // if we're really short, reset out bases per graph
      if (basesPerGraph < opts.basesPerGraph) {
	opts.basesPerGraph = Utils.orderOfMagnitude(basesPerGraph);
      }
      return profile; // enable chaining
    }

    /**
     * distribute profile data to the graph specs
     *
     * @param {Array} graphSpecs - list of graph specifications
     * @returns {Promise} list of modified graph specifications
     */
    function attachProfileData(graphSpecs){     
      // assumes the profile is ordered by coordinate from low to high      
      // eg `[{coordinate: 0}, {coordinate: 10}]`

      // reset the profiles
      graphSpecs.forEach(function(gs){ gs.profile = []; });

      return onProfileData.then(function(profile){
	if (profile == null) throw 'Need profile data';

	// search for the index in our profile where we would insert `c`
	// and have it still be sorted
	function sortedIdx(c){
	  return _.sortedIndex(profile, {'coordinate': c}, 'coordinate');
	}

	graphSpecs.forEach(function(gs){
	  var startIdx = Math.max(0, sortedIdx(gs.startBase) - 1),
	      endIdx = Math.min(sortedIdx(gs.endBase) + 1, profile.length - 1);

	  // shallow copy of the relevant data
	  gs.profile = profile.slice(startIdx, endIdx);
	});

	return graphSpecs;
      });
    }

    function attachExtractData(name, graphSpecs){
      // reset the extracts
      graphSpecs.forEach(function(gs){
	gs.extracts[name] = [];
	gs.headers.push({title: name, lineType:'extract'});
      });

      // TODO: use _.sortedIndex to binary search and array.slice to
      // make shallow copies onto the graph specs
      return Utils.forEachAsync(opts.extracts[name], function(dataPoint){
	graphSpecs.forEach(function(gs){
	  // extract starts in this range?
	  var startsInRange = dataPoint.start >= gs.startBase
		&& dataPoint.start <= gs.endBase,
	      // extract ends in this range?
	      endsInRange = dataPoint.end >= gs.startBase
		&& dataPoint.end <= gs.endBase;
	  if(startsInRange || endsInRange){
	    gs.extracts[name].push(dataPoint);
	  }
	});
      })
      // promise resolves as the graphspecs
	.then(function(){ return graphSpecs;});
    }

    function attachAllExtracts(graphSpecs){
      var promises = [];
      // iterate over the named extracts
      _.forOwn(opts.extracts, function(value, key, obj){
	promises.push(attachExtractData(key, graphSpecs));
      });
      return $q.all(promises).then(function(){ return graphSpecs;});
    }

    function maxPages(){
      return Math.ceil(opts.length / opts.basesPerGraph*opts.graphsPerPage);
    }

    
    function redrawRequest(){
      pendingRedraws++;
      return function(graphSpecs){
	pendingRedraws--;
	if(pendingRedraws == 0){
	  return redraw(graphSpecs);
	}else{
	  $log.log('skipping redraw, still have pending work');
	}

	return graphSpecs;
      };
    }

    function redraw(graphSpecs){
      // TODO: tell everyone to redraw
      $log.log('I should redraw', graphSpecs);
      $rootScope.graphSpecs = graphSpecs;
      return graphSpecs;
    }

    function rebuildGraphs(){
      var t1 = new Date();
      return onWidth.then(makeGraphSpecs)
	.then(attachProfileData)
	.then(attachAllExtracts)
	.then(function(graphSpecs){
	  return opts.graphSpecs = graphSpecs;
	})
	.then(redrawRequest())
	.then(function(graphSpecs){
	  $log.log('rebuild:', new Date() - t1, 'ms');
	  return graphSpecs;
	}).catch(function(){
	  $log.log('failed to rebuild, resetting state', arguments);
	  pendingRedraws = 0;
	});
    }
  }

  angular.module('npact')
    .factory('GraphDealer', GraphDealer);
}());
