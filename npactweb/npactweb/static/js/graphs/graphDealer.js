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
    var opts = {
          graphsPerPage: 5,
          startBase: 0,
          endBase: 0,
          get length(){ return this.endBase - this.startBase; },
          graphSpecs: []
        };

    // public interface
    return {
      opts:opts,
      setProfile:function(profile){
        $log.log('setProfile:', profile.length, 'genes');
        return ProfileReader.load(profile)
          .then(function(summary) {
            $log.log('setProfile: loaded', summary);
            // TODO: not need these
            opts.startBase = summary.startBase;
            opts.endBase = summary.endBase;

            // find a sensible zoom level
            var basesPerGraph = summary.length / opts.graphsPerPage;
            // if we're really short, reset out bases per graph
            if (basesPerGraph < GraphConfig.basesPerGraph) {
              GraphConfig.basesPerGraph = Utils.orderOfMagnitude(basesPerGraph);
            }
            GraphConfig.profileSummary = summary;
          });
      },

      zoomTo: zoomTo,
      showMore: function(){
        opts.graphsPerPage += 5;
        //rebuildGraphs();
      },
      rebuildGraphs:rebuildGraphs
    };

    function zoomTo(startBase, zoomPct, zoomingOut){
      // TODO: make this a options object
      throw new Error('Not implemented');
      var res = GraphingCalculator.zoom(startBase, zoomPct, opts.basesPerGraph, opts.offset, zoomingOut);
      opts.offset = res.offset;
      opts.basesPerGraph = res.basesPerGraph;
      rebuildGraphs();
    }

    function makeGraphSpecs(width){

      var partitions = ProfileReader.partition({
        offset: GraphConfig.offset,
        basesPerGraph: GraphConfig.basesPerGraph
      });

      return _.take(partitions, opts.graphsPerPage)
        .map(function(p){
          return {
            startBase: p.startBase,
            endBase: p.endBase
          };
        });
    }

    function rebuildGraphs(){
      $rootScope.graphSpecs = opts.graphSpecs = makeGraphSpecs();
      opts.visible = {
        startBase : _(opts.graphSpecs).pluck('startBase').min().value(),
        endBase : _(opts.graphSpecs).pluck('endBase').max().value()
      };
    }
  });
