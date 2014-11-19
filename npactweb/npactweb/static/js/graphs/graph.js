angular.module('npact')

  .controller('npactGraphContainerCtrl', function($scope, $element, $window, npactConstants, Utils, GraphConfig, ProfileReader, Evt, $log, GraphingCalculator) {
    var self = this,
        width = $element.width(),
        $win = angular.element($window),
        winHeight = $win.height(),
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
              width: width,
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
        redraw = function() {
          self.graphOptions = makeGraphOptions();
          $scope.$broadcast(Evt.REDRAW);
        },
        rebuild = function() {
          $scope.graphSpecs = ProfileReader.partition(GraphConfig);
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
          if(summary){
            // find a sensible zoom level
            var basesPerGraph = summary.length / 5;
            // if we're really short, reset out bases per graph
            if (basesPerGraph < GraphConfig.basesPerGraph) {
              GraphConfig.basesPerGraph = Utils.orderOfMagnitude(basesPerGraph);
            }
            GraphConfig.profileSummary = summary;
          }
        },
        onResize = function() {
          winHeight = $win.height();
          if($element.width() !== width) {
            width = $element.width();
            redraw();
          }
          $scope.$apply();
        }
    ;
    self.visible = function(el) {
      // The rect.top and rect.bottom are relative to the viewport
      var slack = 50;
      var rect = el.getBoundingClientRect();
      return (rect.top <= 0 && rect.bottom > -slack) ||
        (rect.top >=0 && rect.top < (winHeight + slack));
    };

    $scope.graphHeight = npactConstants.graphSpecDefaults.height;
    $win.on('resize', onResize);
    $win.on('scroll', $scope.$apply.bind($scope));

    // listen for graph events
    $scope.$on(Evt.PAN, onPan);
    $scope.$on(Evt.ZOOM, onZoom);

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

  .directive('npactGraph', function npactGraph(Grapher, Evt, GraphingCalculator, $log) {
    return {
      restrict: 'A',
      scope: { spec: '=npactGraph' },
      require: '^npactGraphContainer',
      link: function($scope, $element, $attrs, ctrl) {
        var g = null,
            visible = ctrl.visible,
            buildOptions = function(range) {
              var opts = angular.extend(
                {}, ctrl.graphOptions, range,
                {$scope: $scope, element: $element[0]});

              opts.m = GraphingCalculator.chart(opts);
              opts.xaxis = GraphingCalculator.xaxis(opts);
              return opts;
            },
            draw = function() {
              if(g || !$scope.spec || !visible($element[0])) return;
              var graphUpdateStart = new Date();
              g = new Grapher(buildOptions($scope.spec));
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
