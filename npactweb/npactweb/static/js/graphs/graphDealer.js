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

  .factory('GraphDealer', function($log, Utils, $q, $rootScope, GraphingCalculator, npactConstants, ExtractParser, ProfileReader, TrackReader, GraphConfig) {
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
          offset: 0,
          graphsPerPage: 5,
          startBase: 0,
          endBase: 0,
          get length(){ return this.endBase - this.startBase; },
          profile: null,
          extracts: {},
          hits: {},
          graphSpecs: []
        };

    // once we have data, load it and rebuild the graphs
    onProfileData.then(rebuildGraphs);

    var maybeDrawTrack = function(key, name, data) {
      return TrackReader.load(name, data)
        .then(function() {
          GraphConfig.loadTrack(name, key);
          if(hasProfileData){
            // TODO: maybe a simpler way to add extracts to existing
            // graphs?
            rebuildGraphs();
          }
        });
    };

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
        return maybeDrawTrack('extracts', exopts.name, exopts.data);
      },
      addHits: function(exopts) {
        return maybeDrawTrack('hits', exopts.name, exopts.data);
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
      setWidth: function(w){
        $log.log('setting graph width to', w);

        if(_onWidth){
          // we're waiting for the initial width, let folks know it's
          // ready
          _onWidth.resolve(w);
          _onWidth = null;
        }else{
          onWidth = $q.when(w);
          rebuildGraphs();
        }
      },
      showMore: function(){
        opts.graphsPerPage += 5;
        rebuildGraphs();
      },
      redraw: rebuildGraphs
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

    function makeGraphSpec(range, width){
      return {
        startBase: range.startBase,
        endBase: range.endBase,
        width: width,
        colors: opts.colorBlindFriendly ? npactConstants.colorBlindLineColors
          : npactConstants.lineColors
      };
    }

    function makeGraphSpecs(width){
      return onProfileData.then(function() {
        var partitions = ProfileReader.partition({
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

      opts.visible = {
        startBase : _(graphSpecs).pluck('startBase').min().value(),
        endBase : _(graphSpecs).pluck('endBase').max().value()
      };
      return graphSpecs;
    }

    function rebuildGraphs(){
      var t1 = new Date();
      return onWidth.then(makeGraphSpecs)
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
