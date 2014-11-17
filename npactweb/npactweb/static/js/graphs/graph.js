angular.module('npact')

  .controller('npactGraphContainerCtrl', function($scope, npactConstants, Utils, GraphConfig, ProfileReader, Evt, $log, GraphingCalculator) {
    var self = this,
        graphUpdateStart = null, graphsDrawn = 0,
        visibleGraphs = 5,
        graphSpecs = [],
        getGraphConfig = function() { return GraphConfig; },
        getProfileSummary = function() {
          try {return ProfileReader.summary();} catch(e) { }
          return null;
        },
        makeGraphOptions = function() {
          var opts = angular.extend(
            {},
            npactConstants.graphSpecDefaults,
            {
              headerSpec: GraphConfig.headerSpec(),
              margin: Utils.orderOfMagnitude(GraphConfig.basesPerGraph, -1),
              width: GraphConfig.width,
              colors: GraphConfig.colorBlindFriendly ?
                npactConstants.colorBlindLineColors : npactConstants.lineColors
            });
          opts.headerY = opts.headerSpec.headerY;
          return opts;
        },
        onPan = function(evt, opts) {
          var offset = Math.floor(opts.newStartBase - opts.oldStartBase);
          GraphConfig.offset += offset;
          // TODO: Event originated from outside ng, but why doesn't
          // `$watch` pick up the `offset` change?
          $scope.$apply();
        },
        onZoom = function(evt, opts) {
          var zoomOpts = angular.extend({}, opts, GraphConfig),
              res = GraphingCalculator.zoom(zoomOpts);
          GraphConfig.offset = res.offset;
          GraphConfig.basesPerGraph = res.basesPerGraph;
          // TODO: Event originated from outside ng, but why doesn't
          // `$watch` pick up the `offset` change?
          $scope.$apply();
        },
        onGraphRedrawComplete = function(evt) {
          graphsDrawn++;
          if(graphsDrawn === $scope.graphSpecs.length){
            $log.log('graphs done', new Date() - graphUpdateStart, 'ms');
          }
        },
        onGraphConfigChanged = function(newValue, oldValue){
          var cmd = newValue.refreshCommand(oldValue);
          $log.log('graph config changed:', cmd);
          graphUpdateStart = new Date();
          graphsDrawn = 0;
          self.graphOptions = makeGraphOptions();
          switch(cmd){
          case Evt.REBUILD:
            graphSpecs = ProfileReader.partition(GraphConfig);
            $scope.graphSpecs = _.take(graphSpecs, visibleGraphs);
            break;
          case Evt.REDRAW:
            $scope.$broadcast(cmd);
            break;
          }
        },
        onProfileSummaryChanged = function(summary, oldSummary) {
          if(summary){
            // find a sensible zoom level
            var basesPerGraph = summary.length / 5;
            // if we're really short, reset out bases per graph
            if (basesPerGraph < GraphConfig.basesPerGraph) {
              GraphConfig.basesPerGraph = Utils.orderOfMagnitude(basesPerGraph);
            }
            GraphConfig.profileSummary = summary;
          }
        }
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

    // listen for graph events
    $scope.$on(Evt.PAN, onPan);
    $scope.$on(Evt.ZOOM, onZoom);
    $scope.$on(Evt.GRAPH_REDRAW_COMPLETE, onGraphRedrawComplete);

    // watch the environment for changes we care about
    $scope.$watch(getGraphConfig, onGraphConfigChanged, true); // deep-equality
    $scope.$watchCollection(getProfileSummary, onProfileSummaryChanged);
  })
  .directive('npactGraphContainer', function(STATIC_BASE_URL) {
    return {
      restrict: 'A',
      scope: {},
      templateUrl: STATIC_BASE_URL + 'js/graphs/graph-container.html',
      controller: 'npactGraphContainerCtrl as ctrl'
    };
  })

  .directive('npactGraph', function npactGraph(Grapher, Evt, GraphingCalculator) {
    return {
      restrict: 'A',
      scope: { spec: '=npactGraph' },
      require: '^npactGraphContainer',
      link: function($scope, $element, $attrs, ctrl) {
        var g = null,
            buildOptions = function(range) {
              var opts = angular.extend(
                {}, ctrl.graphOptions, range,
                {$scope: $scope, element: $element[0]});

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
  .directive('npactExtract', function(STATIC_BASE_URL) {
    return {
      restrict: 'A',
      scope: { extract: '=npactExtract'},
      templateUrl: STATIC_BASE_URL + 'js/graphs/extract.html'
    };
  })
;
