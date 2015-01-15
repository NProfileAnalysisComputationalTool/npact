angular.module('npact')
  .controller('npactGraphContainerCtrl', function($scope, $element, $window, $log, $timeout,
                                           npactConstants, Utils, GraphConfig, Evt,
                                           GraphingCalculator) {
    'use strict';

    var getWidth = function() {
      return $element.width() - 15; //from style.css `.graph`
    };

    //The mouse event responders
    var onPan = function(offset) {
          GraphConfig.offset += Math.round(offset);
          $scope.$apply();
        },
        onZoom = function(opts) {
          $log.log('Zoom event:', opts);
          var zoomOpts = angular.extend({}, GraphConfig, opts);
          //updates `offset`, and `basesPerGraph`
          angular.extend(GraphConfig, GraphingCalculator.zoom(zoomOpts));
          $scope.$apply();
        };

    //The baseOpts are the graph options that are the same for every graph
    var baseOpts = { width: getWidth(), onPan: onPan, onZoom: onZoom };
    this.graphOptions = function(idx) {
      // This function builds the specific options for a graph; as
      // many graph rows will never be drawn this only generates the
      // object for a row when needed.
      var start = $scope.graphSpecs[idx];
      var opts = { startBase: start,
                   endBase: start + GraphConfig.basesPerGraph};
      opts = angular.extend(opts, baseOpts);
      return opts;
    };

    var draw = function() { $scope.$broadcast(Evt.DRAW); },
        redraw = function() {
          if (!baseOpts.axisTitle || !baseOpts.m) return;  //too early to do anything
          $scope.$broadcast(Evt.REDRAW);
        },
        rebuild = function() {
          if (!baseOpts.axisTitle || !baseOpts.m) return;  //too early to do anything
          $scope.$broadcast(Evt.REBUILD);
        };


    /*** Watch the config for changes we care about ***/
    $scope.$watchCollection(
      function() {
        return [GraphConfig.length, GraphConfig.basesPerGraph,
                GraphConfig.offset, GraphConfig.startBase, GraphConfig.endBase];
      },
      function() {
        // basic row geometry changed, repartition and rebuild
        if(GraphConfig.startBase === undefined ||
           GraphConfig.endBase === undefined ||
           GraphConfig.basesPerGraph === undefined) { return; }
        $scope.graphSpecs = GraphingCalculator.partition(GraphConfig);
        $log.log('Partitioned into', $scope.graphSpecs.length, 'rows.');
        updateVisibility(); //number of rows might have changed.
        $timeout(rebuild);
      });

    $scope.$watch(GraphConfig.profileTitle, function(title) {
      //A profileTitle change indicates different nucleotides: rebuild
      baseOpts.axisTitle = title;
      rebuild();
    });

    $scope.$watch(GraphConfig.activeTracks, function(val) {
      //Find headers and headerY
      angular.extend(baseOpts, GraphingCalculator.trackSizeCalc(val));
      baseOpts.m = GraphingCalculator.chart(baseOpts);
      updateRowHeight(baseOpts.m.height);
      redraw();
    }, true);

    $scope.$watch(function() { return GraphConfig.colorBlindFriendly; },
                  function(val) {
                    baseOpts.colors = val ?
                      npactConstants.colorBlindLineColors :
                      npactConstants.lineColors;
                    redraw();
                  });


    /***  Scrolling and Graph Visibility management ***/
    var $win = angular.element($window),
        winHeight = $win.height(),
        slack = 50, // how many pixels outside of viewport to render
        topOffset = $element.offset().top,
        topIdx = 0, bottomIdx = 0,
        graphRowHeight = 0,
        updateRowHeight = function(height) {
          $scope.graphHeight = height;
          //add padding and border from style.css `.graph`
          graphRowHeight = height + 8 + 8 + 1;
          updateVisibility();
        },
        updateVisibility = function() {
          if(!baseOpts.m) return;
          var scrollDist = $window.scrollY - topOffset - slack;
          topIdx = Math.floor(scrollDist / graphRowHeight);
          bottomIdx = topIdx + Math.ceil(winHeight / graphRowHeight);
        },
        onScroll = function() {
          updateVisibility();
          draw();
        },
        onResize = function() {
          winHeight = $win.height();
          if(getWidth() !== baseOpts.width) {
            topOffset = $element.offset().top;
            baseOpts.width = getWidth();
            baseOpts.m = GraphingCalculator.chart(baseOpts);
            updateVisibility();
            redraw();
          }
          else {
            //If the width didn't change then its the same as scrolling
            onScroll();
          }
        };

    this.visible = function(idx) { return idx >= topIdx && idx <= bottomIdx; };
    $win.on('resize', onResize);
    $win.on('scroll', onScroll);
    $scope.$on('$destroy', function() {
      $win.off('resize', onResize);
      $win.off('scroll', onScroll);
    });
  })
  .directive('npactGraphContainer', function(STATIC_BASE_URL) {
    return {
      restrict: 'A',
      scope: {},
      templateUrl: STATIC_BASE_URL + 'js/graphs/graph-container.html',
      controller: 'npactGraphContainerCtrl as ctrl'
    };
  })

  .directive('npactGraph', function (Grapher, Evt, GraphingCalculator, $log, $timeout) {
    'use strict';
    return {
      restrict: 'A',
      require: '^npactGraphContainer',
      link: function($scope, $element, $attrs, ctrl) {
        var g = null,
            visible = ctrl.visible,
            idx = $attrs.idx,
            el = $element[0],
            // redraw gets set for all graphs once (e.g. a new track
            // triggers broadcasts redraw), but only gets cleared as
            // the currently visible ones are drawn
            redraw = false,
            draw = function() {
              if(!redraw || !visible(idx)) { return; }
              var opts = ctrl.graphOptions(idx);
              //However long it actually takes to draw, we have the
              //latest options as of this point
              redraw = false;
              (g || (g = new Grapher(el, opts)))
                .redraw(opts)
                .catch(function() {
                  //something went wrong, we will still need to redraw this
                  redraw = true;
                })
              ;
            },
            schedule = function() {
              if(!redraw || !visible(idx)) { return; }
              $timeout(draw, 0, false);
            },
            discard = function() {
              if(g) { g.destroy(); g = null; }
            };
        $scope.$on(Evt.DRAW, schedule);
        $scope.$on(Evt.REDRAW, function() { redraw = true; schedule();});
        $scope.$on(Evt.REBUILD, function() { discard(); redraw = true; schedule(); });
        $scope.$on('$destroy', discard);
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
