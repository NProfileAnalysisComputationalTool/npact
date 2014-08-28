(function(){

  var graphSpecDefaults = {
    // TODO: determine me dynamically
    leftPadding: 120,

    axisLabelFontsize: 11,
    axisFontcolor: "#444",
    axisTitleFontsize: 20,
    borderColor: "#444",
    rightPadding: 25,
    // TODO: reduce duplication between here and `page.html`
    height: 200,
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
  function GraphDealer($log, Utils, $q, $rootScope, GraphingCalculator){

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
	  offset: 0,
	  graphsPerPage: 5,
	  length: 0,
	  profile: null,
	  extracts: {},
	  graphSpecs: [],
	  get maxPages (){
	    return Math.ceil(this.length / this.basesPerGraph*this.graphsPerPage);
	  },
	  get hasPreviousPage(){ return this.page > 0; },
	  get hasNextPage(){ return this.page < this.maxPages; }
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
      zoomTo: zoomTo,
      panTo: panTo,
      nextPage: function(){
	$log.log('nextPage', arguments);
	opts.page++;
	rebuildGraphs();
      },
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

    function zoomTo(startBase, zoomPct, zoomingOut){
      // TODO: make this a options object
      var res = GraphingCalculator.zoom(startBase, zoomPct, opts.basesPerGraph, opts.offset, zoomingOut);
      opts.offset = res.offset;
      opts.basesPerGraph = res.basesPerGraph;
      rebuildGraphs();
    }

    /**
     * pan all the graphs
     *
     * @param {Number} oldStartBase - the start of the pan, in gene space
     * @param {Number} newStartBase - the end of the pan, in gene space
     */
    function panTo(oldStartBase, newStartBase){
      var offset = Math.floor(newStartBase - oldStartBase);
      opts.offset += offset;
      $log.log('panTo', opts.offset);
      rebuildGraphs();
    }

    function makeGraphSpec(startBase, width){
      var endBase = startBase + opts.basesPerGraph,
	  spec = angular.extend({
	    range: [startBase, endBase],
	    // get some margin in here
	    startBase: startBase,
	    endBase: endBase,
	    extracts: {},
	    profile: [],
	    headers: [],
	    width: width,
	    colors: opts.colorBlindFriendly ? colorBlindLineColors : lineColors
	  }, graphSpecDefaults);

      spec.chart = GraphingCalculator.chart(spec);
      spec.xaxis = GraphingCalculator.xaxis(spec);
      return spec;
    }

    function makeGraphSpecs(width){
      var base = opts.page * opts.basesPerGraph * opts.graphsPerPage
	    + opts.offset;
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

      // want to grab data a little out of range so we can scroll a
      // little bit
      var margin = Utils.orderOfMagnitude(opts.basesPerGraph, -1);
      // TODO: reduce duplication between here and
      // `GraphingCalculator.stops`
      return onProfileData.then(function(profile){

	graphSpecs.forEach(function(gs){
	  var startIdx = sortedIdx(gs.startBase - margin) - 1,
	      endIdx = sortedIdx(gs.endBase + margin) + 1;

	  // shallow copy of the relevant data, guarding against
	  // out-of-bounds indexes
	  gs.profile = profile.slice(Math.max(startIdx,0),
				     Math.min(endIdx, profile.length - 1));
	});

	return graphSpecs;

	// search for the index in our profile where we would insert `c`
	// and have it still be sorted
	function sortedIdx(c){
	  return _.sortedIndex(profile, {'coordinate': c}, 'coordinate');
	}
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

    /**
     * tell the world to redraw everything
     */
    function redraw(graphSpecs){
      $rootScope.graphSpecs = graphSpecs;

      opts.visible = {
	startBase : _(graphSpecs).pluck('startBase').min().value(),
	endBase : _(graphSpecs).pluck('endBase').max().value()
      };
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
