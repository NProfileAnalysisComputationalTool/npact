angular.module('npact')

  .controller('npactGraphContainerCtrl', function($scope, npactConstants, Utils, GraphConfig, ProfileReader, Evt, $log) {
    var self = this,
        graphUpdateStart = null, graphsDrawn = 0,
        visibleGraphs = 5,
        graphSpecs = [],
        getGraphConfig = function() { return GraphConfig; }
    ;

    $scope.graphHeight = npactConstants.graphSpecDefaults.height;
    /**
     * add more visible entries to $scope
     */
    self.addMore = function(){
      if($scope.graphSpecs){
        $log.log('scrolling down via infinite scroller');
        graphUpdateStart = new Date();
        Utils.extendByPage(graphSpecs, $scope.graphSpecs, visibleGraphs);
      }
    };

    $scope.$on(Evt.PAN, function(evt, opts) {
      var offset = Math.floor(opts.newStartBase - opts.oldStartBase);
      GraphConfig.offset += offset;
      // TODO: Event originated from outside ng, but why doesn't
      // `$watch` pick up the `offset` change?
      $scope.$apply();
    });

    $scope.$on(Evt.ZOOM, function(evt, opts) {
      var res = GraphingCalculator.zoom(angular.extend({}, opts, GraphConfig));
      GraphConfig.offset = res.offset;
      GraphConfig.basesPerGraph = res.basesPerGraph;
      // TODO: Event originated from outside ng, but why doesn't
      // `$watch` pick up the `offset` change?
      $scope.$apply();
    });

    $scope.$on(Evt.GRAPH_REDRAW_COMPLETE, function(evt) {
      graphsDrawn++;
      if(graphsDrawn === $scope.graphSpecs.length){
        $log.log('graphs done', new Date() - graphUpdateStart, 'ms');
      }
    });

    $scope.$watch(getGraphConfig, function(newValue, oldValue){
      var cmd = newValue.refreshCommand(oldValue);
      $log.log('graph config changed:', cmd);
      graphUpdateStart = new Date();
      graphsDrawn = 0;
      switch(cmd){
      case Evt.REBUILD:
        graphSpecs = ProfileReader.partition(GraphConfig);
        $scope.graphSpecs = _.take(graphSpecs, visibleGraphs);
        break;
      case Evt.REDRAW:
        $scope.$broadcast(cmd);
        break;
      }
    }, true); // deep-equality


    $scope.$watchCollection(
      function() {
        try {return ProfileReader.summary();} catch(e) { }
        return null;
      },
      function(summary, oldSummary) {
        if(summary){
          // find a sensible zoom level
          var basesPerGraph = summary.length / 5;
          // if we're really short, reset out bases per graph
          if (basesPerGraph < GraphConfig.basesPerGraph) {
            GraphConfig.basesPerGraph = Utils.orderOfMagnitude(basesPerGraph);
          }
          GraphConfig.profileSummary = summary;
        }
      });
  })
  .directive('npactGraphContainer', function(STATIC_BASE_URL) {
    return {
      restrict: 'A',
      scope: {},
      templateUrl: STATIC_BASE_URL + 'js/graphs/graph-container.html',
      controller:'npactGraphContainerCtrl',
      controllerAs: 'ctrl'
    };
  })

  .directive('npactGraph', function npactGraph($log, Grapher, Evt, GraphConfig, Utils, npactConstants, GraphingCalculator){
    return {
      restrict: 'A',
      scope:{
        spec: '=npactGraph'
      },
      controller: function($scope, $element, $attrs){
        var g = null, buildOptions = function(range) {
          var opts = angular.extend(
            {},
            npactConstants.graphSpecDefaults,
            range,
            {
              $scope: $scope,
              element: $element[0],
              headerSpec: GraphConfig.headerSpec(),
              margin: Utils.orderOfMagnitude(GraphConfig.basesPerGraph, -1),
              width: GraphConfig.width,
              colors: GraphConfig.colorBlindFriendly ?
                npactConstants.colorBlindLineColors : npactConstants.lineColors
            });

          opts.headerY = opts.headerSpec.headerY;
          opts.m = GraphingCalculator.chart(opts);
          opts.xaxis = GraphingCalculator.xaxis(opts);
          return opts;
        },
            redraw = function() {
              if(!$scope.spec) { return; }
              if(g !== null) {g.stage.destroy();}
              g = new Grapher(buildOptions($scope.spec));
              g.redraw().then($scope.$emit(Evt.GRAPH_REDRAW_COMPLETE));
            };

        $scope.$on(Evt.REDRAW, redraw);
        redraw();
      }
    };
  })
  .directive('npactExtract', function(STATIC_BASE_URL){
    return {
      restrict: 'A',
      scope: { extract:'=npactExtract'},
      templateUrl:STATIC_BASE_URL+'js/graphs/extract.html'
    };
  })
;
