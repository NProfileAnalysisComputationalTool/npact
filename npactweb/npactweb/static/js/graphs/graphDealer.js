angular.module('npact')
  .constant('npactConstants', {
    graphSpecDefaults: {
      // TODO: determine me dynamically
      leftPadding: 120,

      axisLabelFontsize: 11,
      axisFontcolor: '#444',
      axisTitleFontsize: 20,
      borderColor: '#444',
      rightPadding: 25,
      height: 200,
      profileTicks: 5,


      // header labels and arrows
      headerY: 5,
      headerLabelPadding: 10,
      headerLabelFontcolor: '#444',
      headerLabelFontsize: 11,
      headerArrowHeight: 12,
      headerArrowWidth: 6,
      headerArrowFontsize: 9,
      axisTitle: '% GC'

    },
    lineColors : {
      r: 'red',
      g: 'blue',
      b: 'green'
    },
    colorBlindLineColors : {
      r: 'rgb(213, 94, 0)',
      g: 'rgb(204, 121, 167)',
      b: 'rgb(0, 114, 178)'
    },
    // how much vertical space to leave for different kinds of headers
    headerSizes:{'extracts': 30, 'hits': 20}
  })

  .factory('GraphDealer', function($log, Utils, $q, $rootScope, GraphingCalculator, npactConstants, ExtractParser, ProfileReader) {
    'use strict';
    // the `GraphDealer` hands out graph data and processes events
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
          startBase: 0,
          endBase: 0,
          get length(){ return this.endBase - this.startBase; },
          profile: null,
          extracts: {},
          hits: {},
          graphSpecs: [],
          get maxPages (){
            return Math.ceil(this.length / this.basesPerGraph*this.graphsPerPage);
          },
          get hasPreviousPage(){ return this.page > 0; },
          get hasNextPage(){ return this.page < this.maxPages; }
        };

    // once we have data, load it and rebuild the graphs
    onProfileData.then(rebuildGraphs);


    // public interface
    return {
      opts:opts,
      setProfile:function(profile){
        $log.log('setProfile:', profile.length, 'genes');
        return ProfileReader.load(profile)
          .then(function(summary) {
            opts.profileSummary = summary;
            hasProfileData = true;
            opts.profile = profile;

            // TODO: not need these
            opts.startBase = summary.startBase;
            opts.endBase = summary.endBase;

            // find a sensible zoom level
            var basesPerGraph = opts.profileSummary.length / opts.graphsPerPage;
            // if we're really short, reset out bases per graph
            if (basesPerGraph < opts.basesPerGraph) {
              opts.basesPerGraph = Utils.orderOfMagnitude(basesPerGraph);
            }
            _onProfileData.resolve(profile);
          });
      },

      addExtract:function(exopts){
        return addHeader(opts.extracts, ExtractParser, exopts.name, exopts.data);
      },
      addHits: function(exopts) {
        return addHeader(opts.hits, ExtractParser, exopts.name, exopts.data);
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

        if(_onWidth){
          // we're waiting for the initial width, let folks know it's
          // ready
          _onWidth.resolve(w);
          _onWidth = null;
        }else{
          // setup a new resolved promise, since we don't have to wait
          // for the value this time
          var q = $q.defer();
          onWidth = q.promise;
          q.resolve(w);
          rebuildGraphs();
        }
      },
      showMore: function(){
        opts.graphsPerPage += 5;
        rebuildGraphs();
      }
    };
    function addHeader(collection, parser, name, data) {
        $log.log('addHeader', name);

        return parser.parseAsync(data)
          .then(function(data){
            $log.log('parsed ', name, hasProfileData);
            // save it for later
            collection[name] = data;
            // if we've got profile data already, rebuild
            if(hasProfileData){
              // TODO: maybe a simpler way to add extracts to existing
              // graphs?
              rebuildGraphs();
            }
            // if we're still waiting on profile data, then when it hits
            // this extract will be processed.
          }, function(){
            $log.log('failed to parse', name, arguments);
          });

    }

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

    function makeGraphSpec(range, width){
      var spec = angular.extend({
            startBase: range.startBase,
            endBase: range.endBase,
            extracts: {},
            hits: {},
            headers: [],
            width: width,
            colors: opts.colorBlindFriendly ? npactConstants.colorBlindLineColors
              : npactConstants.lineColors,

            // how much space should be taken by the line graph alone
            get profileHeight(){
              // TODO: reduce duplication between here and GraphingCalculator
              return this.height - this.headerY - this.axisLabelFontsize -
                3*this.profileTicks;
            },
            get range() {return [this.startBase, this.endBase];}
          }, npactConstants.graphSpecDefaults);

      return spec;
    }

    function makeGraphSpecs(width){
      return onProfileData.then(function() {
        var base = opts.page * opts.basesPerGraph * opts.graphsPerPage +
              opts.offset,

            partitions = ProfileReader.partition({
              basesPerGraph: opts.basesPerGraph,
              summary: opts.profileSummary,
              margin: Utils.orderOfMagnitude(opts.basesPerGraph, -1)
            });

        return _.take(partitions, opts.graphsPerPage)
          .map(function(p){
            return makeGraphSpec(p, width);
          });
      });
    }

    function attachAllData(graphSpecs, key) {
      var promises = [];
      // iterate over the named extracts
      _.forOwn(opts[key], function(data, name){
        promises.push(attachData(graphSpecs, name, key));
      });
      return $q.all(promises).then(function(){ return graphSpecs;});
    }

    function attachData(graphSpecs, name, key) {
      var headerHeight = npactConstants.headerSizes[key];

      graphSpecs.forEach(function(gs){
        gs[key][name] = [];
        gs.headers.push({
          text: name,
          lineType:key,
          y: gs.headerY,
          height: headerHeight
        });
        // increment where the next header will start
        gs.headerY += headerHeight;
      });

      // TODO: use _.sortedIndex to binary search and array.slice to
      // make shallow copies onto the graph specs
      return Utils.forEachAsync(opts[key][name], function(dataPoint){
        graphSpecs.forEach(function(gs){
          // extract starts in this range?
          var startsInRange = dataPoint.start >= gs.startBase &&
                dataPoint.start <= gs.endBase,
              // extract ends in this range?
              endsInRange = dataPoint.end >= gs.startBase &&
                dataPoint.end <= gs.endBase;
          if(startsInRange || endsInRange){
            gs[key][name].push(dataPoint);
          }
        });
      })
      // promise resolves as the graphspecs
        .then(function(){ return graphSpecs;});
    }

    function redrawRequest(){
      pendingRedraws++;
      return function(graphSpecs){
        pendingRedraws--;
        if(pendingRedraws === 0){
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
      graphSpecs.forEach(function(spec){
        spec.chart = GraphingCalculator.chart(spec);
        spec.xaxis = GraphingCalculator.xaxis(spec);
      });

      opts.visible = {
        startBase : _(graphSpecs).pluck('startBase').min().value(),
        endBase : _(graphSpecs).pluck('endBase').max().value()
      };
      return graphSpecs;
    }

    function rebuildGraphs(){
      var t1 = new Date();
      return onWidth.then(makeGraphSpecs)
        .then(_.partialRight(attachAllData, 'extracts'))
        .then(_.partialRight(attachAllData, 'hits'))
        .then(function(graphSpecs){
          return (opts.graphSpecs = graphSpecs);
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
  });
