angular.module('npact')

  .controller('npactGraphContainerCtrl', function($scope, $element, $window, npactConstants, Utils, GraphConfig, ProfileReader, Evt, $log, GraphingCalculator) {

    var baseOpts = angular.extend({ width: $element.width() },
                                  npactConstants.graphSpecDefaults),
        getGraphConfig = function() { return GraphConfig; },
        getProfileSummary = function() {
          try {return ProfileReader.summary();} catch(e) { }
          return null;
        },
        updateBaseOptions = function() {
          baseOpts.headerSpec = GraphConfig.headerSpec();
          baseOpts.margin = Utils.orderOfMagnitude(GraphConfig.basesPerGraph, -1);
          baseOpts.colors = GraphConfig.colorBlindFriendly ?
            npactConstants.colorBlindLineColors : npactConstants.lineColors;
          baseOpts.headerY = baseOpts.headerSpec.headerY;
          return baseOpts;
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
        redraw = function() {
          updateBaseOptions();
          $scope.$broadcast(Evt.REDRAW);
        },
        rebuild = function() {
          $scope.graphSpecs = ProfileReader.partition(GraphConfig);
          updateVisibility();
          redraw();
        },
        onGraphConfigChanged = function(newValue, oldValue){
          var cmd = newValue.refreshCommand(oldValue);
          $log.log('graph config changed:', cmd);
          switch(cmd) {
          case Evt.REBUILD:
            rebuild();
            break;
          case Evt.REDRAW:
            redraw();
            break;
          }
        },
        onProfileSummaryChanged = function(summary, oldSummary) {
          if(summary) {
            // find a sensible zoom level
            var basesPerGraph = summary.length / 5;
            // if we're really short, reset out bases per graph
            if (basesPerGraph < GraphConfig.basesPerGraph) {
              GraphConfig.basesPerGraph = Utils.orderOfMagnitude(basesPerGraph);
            }
            GraphConfig.profileSummary = summary;
          }
        };

    this.graphOptions = function(idx, element) {
      range = $scope.graphSpecs[idx];
      var opts = angular.extend(
        {$scope: $scope, element: element},
        baseOpts, range);
      opts.m = GraphingCalculator.chart(opts);
      opts.xaxis = GraphingCalculator.xaxis(opts);
      return opts;
    };

    // listen for graph events
    $scope.$on(Evt.PAN, onPan);
    $scope.$on(Evt.ZOOM, onZoom);

    // watch the environment for changes we care about
    $scope.$watch(getGraphConfig, onGraphConfigChanged, true); // deep-equality
    $scope.$watchCollection(getProfileSummary, onProfileSummaryChanged);


    /***  Scrolling and GRaph Visibility management ***/
    var $win = angular.element($window),
        winHeight = $win.height(),
        graphHeight = npactConstants.graphSpecDefaults.height,
        topOffset = $element.offset().top,
        topIdx = 0, bottomIdx = 0,
        slack = 50, // how many pixels outside of viewport to render
        updateVisibility = function() {
              recttop = topOffset - $window.scrollY,
          topIdx = Math.floor(-(recttop + slack) / graphHeight);
          bottomIdx = topIdx + Math.ceil((winHeight + slack) / graphHeight);
        },
        onScroll = function() {
          updateVisibility();
          $scope.$apply();
        },
        onResize = function() {
          winHeight = $win.height();
          topOffset = $element.offset().top;
          updateVisibility();
          if($element.width() !== baseOpts.width) {
            baseOpts.width = $element.width();
            redraw();
          }
          $scope.$apply();
        };
    this.visible = function(idx) {
      return idx >= topIdx && idx <= bottomIdx;
    };
    $scope.graphHeight = graphHeight;
    $win.on('resize', onResize);
    $win.on('scroll', onScroll);
  })
  .directive('npactGraphContainer', function(STATIC_BASE_URL) {
    return {
      restrict: 'A',
      scope: {},
      templateUrl: STATIC_BASE_URL + 'js/graphs/graph-container.html',
      controller: 'npactGraphContainerCtrl as ctrl'
    };
  })

  .directive('npactGraph', function npactGraph(Grapher, Evt, GraphingCalculator, $log) {
    return {
      restrict: 'A',
      require: '^npactGraphContainer',
      link: function($scope, $element, $attrs, ctrl) {
        var g = null,
            visible = ctrl.visible,
            idx = $attrs.idx,
            el = $element[0];
            draw = function() {
              if(g || !visible(idx)) return;
              var graphUpdateStart = new Date();
              g = new Grapher(ctrl.graphOptions(idx, el));
              g.redraw().then(function() {
                $log.log('Draw of starBase:',
                         $scope.spec.startBase, 'took',
                         new Date() - graphUpdateStart, 'ms');
              });
            },
            redraw = function() {
              if(g !== null) { g.stage.destroy(); g = null; }
              //Don't need to do anything: `$scope.$watch(draw);`
              //handles it
            };
        $scope.$on(Evt.REDRAW, redraw);
        //Just call draw every time
        $scope.$watch(draw);
        $scope.$on('$destroy', function() {
          if(g) { g.stage.destroy(); }
        });
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
