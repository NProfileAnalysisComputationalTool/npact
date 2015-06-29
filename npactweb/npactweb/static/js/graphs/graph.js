angular.module('npact')
  .controller('npactGraphContainerCtrl', function($scope, $element, $window, $log, $timeout,
                                           npactConstants, Utils, GraphConfig, Evt,
                                           GraphingCalculator) {
    'use strict';

    var getWidth = function() {
      return $element.width(); //from style.css `.graph`
    };

    //The baseOpts are the graph options that are the same for every graph
    var baseOpts = { width: getWidth() },
        updateMetrics = function() {
          baseOpts.m = GraphingCalculator.chart(baseOpts);
        };
    this.graphOptions = function(idx) {
      // This function builds the specific options for a graph; as
      // many graph rows will never be drawn this only generates the
      // object for a row when needed.
      var start = $scope.graphSpecs[idx];
      return angular.extend({ startBase: start }, baseOpts);
    };

    var ready = function() {
          return GraphConfig.nucleotides && GraphConfig.nucleotides.length &&
            baseOpts.m;
        },
        draw = function() { $scope.$broadcast(Evt.DRAW); },
        redraw = function() {
          if (!ready()) return;  //too early to do anything
          $scope.$broadcast(Evt.REDRAW);
        },
        rebuild = function() {
          if (!ready()) return;  //too early to do anything
          $scope.$broadcast(Evt.REBUILD);
        };


    /*** Watch the config for changes we care about ***/
    $scope.$watchCollection(
      function() { return [GraphConfig.basesPerGraph,
                    GraphConfig.offset, GraphConfig.startBase, GraphConfig.endBase]; },
      function() {
        // basic row geometry changed, repartition and rebuild
        if(isNaN(GraphConfig.startBase) ||
           isNaN(GraphConfig.endBase) ||
           isNaN(GraphConfig.basesPerGraph)) { return; }
        baseOpts.basesPerGraph = GraphConfig.basesPerGraph;
        $scope.graphSpecs = GraphingCalculator.partition(GraphConfig);
        updateMetrics();
        $log.log('Partitioned into', $scope.graphSpecs.length, 'rows.');
        updateVisibility(); //number of rows might have changed.
        $timeout(rebuild);
      });

    $scope.$watch(function() { return GraphConfig.nucleotides; }, function() {
      $scope.$broadcast('n-profile');
    }, true);
    $scope.$watch(function() { return GraphConfig.colorBlindFriendly; }, redraw);

    $scope.$watch(GraphConfig.activeTracks, function(val) {
      //Find headers and headerY
      baseOpts.tracks = val;
      updateMetrics();
      updateRowHeight(baseOpts.m.height);
      redraw();
    }, true);



    /***  Scrolling and Graph Visibility management ***/
    var $win = angular.element($window),
        winHeight = $win.height(),
        slack = 20, // how many pixels outside of viewport to render
        topOffset = 0,
        topIdx = 0, bottomIdx = 0,
        graphRowHeight = 0,
        updateRowHeight = function(height) {
          topOffset = Math.floor($element.offset().top + _.parseInt($element.css('padding-top')));
          $log.log('topOffset:', topOffset);
          $scope.graphHeight = height;
          try {
             graphRowHeight = angular.element('#graph_0', $element).outerHeight();
          }
          catch(e) {
            graphRowHeight = 0; // There are no rows
          }
          updateVisibility();
        },
        scrollToBase = function(base) {
          if(isNaN(base)) return;
          var idx = Math.floor(base / GraphConfig.basesPerGraph);
          var id = '#graph_' + idx;
          $timeout(function() {
            $log.log('scrolling to', id);
            $window.scrollTo(0, $(id).offset().top);
          }, 40);
        },
        updateVisibility = function() {
          if(!baseOpts.m) return;
          var scrollDist = $window.scrollY - topOffset - slack;
          topIdx = Math.floor(scrollDist / graphRowHeight);
          bottomIdx = topIdx + Math.ceil((winHeight + 2* slack) / graphRowHeight);
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
            updateMetrics();
            updateVisibility();
            redraw();
          }
          else {
            //If the width didn't change then its the same as scrolling
            onScroll();
          }
        },
        onKeyDown = _.throttle(function(event) {
          var keyCode = event.which, delta;
          switch(keyCode) {
          case 37: // left key
            delta = -GraphConfig.basesPerGraph / 100;
            GraphConfig.offset += delta;
            $scope.$broadcast('offset', delta);
            break;
          case 39: // right key
            delta = GraphConfig.basesPerGraph / 100;
            GraphConfig.offset += delta;
            $scope.$broadcast('offset', delta);
            break;
          }
        }, 40, {leading:true}),
        onKeyUp = _.debounce(function(event) {
          var keyCode = event.which;
          switch(keyCode) {
          case 37: // left key
          case 39: // right key
            $scope.$apply();
            break;
          }
        }, 800);
    $scope.$watch(function() { return GraphConfig.gotoBase; }, scrollToBase);

    this.visible = function(idx) { return idx >= topIdx && idx <= bottomIdx; };
    $win.on('resize', onResize);
    $win.on('scroll', onScroll);
    $win.on('keydown', onKeyDown);
    $win.on('keyup', onKeyUp);
    $scope.$on('$destroy', function() {
      $win.off('resize', onResize);
      $win.off('scroll', onScroll);
      $win.off('keydown', onKeyDown);
      $win.off('keyup', onKeyUp);
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

  .directive('npactGraph', function (Grapher, Evt, GraphingCalculator, GraphConfig,
                              $log, $timeout) {
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
            draw = function(force) {
              if(!redraw || (!force && !visible(idx))) { return null; }
              var opts = ctrl.graphOptions(idx);
              //However long it actually takes to draw, we have the
              //latest options as of this point
              redraw = false;
              return (g || (g = new Grapher(el, opts)))
                    .redraw(opts)
                    .catch(function() {
                      //something went wrong, we will still need to redraw this
                      redraw = true;
                    });
            },
            schedule = function(force) {
              if(!redraw || (!force && !visible(idx))) { return null; }
              return $timeout(_.partial(draw, force), 0, false);
            },
            discard = function() {
              if(g) { g.destroy(); g = null; }
            };
        $scope.$on(Evt.DRAW, _.partial(schedule, false));
        $scope.$on(Evt.REDRAW, function() { redraw = true; schedule();});
        $scope.$on(Evt.REBUILD, function() { discard(); redraw = true; schedule(); });
        $scope.$on('n-profile', function() {
          redraw = true;
          if(g) { g.clearProfilePoints(); }
          schedule();
        });
        $scope.$on('offset', function(event, dx) {
          if(g && g.offset && visible(idx)) {
            g.offset(dx);
          }
        });
        $scope.$on('$destroy', discard);

        $scope.$on('print', function(evt, callback) {
          var p = schedule(true);
          if(p) {
            callback(p.then(function() { return g.replaceWithImage(); }));
          }
          else { callback(g.replaceWithImage()); }
        });
        $scope.$watch(function() { return GraphConfig.gotoBase; }, function(gotoBase, fromBase) {
          var start = idx * GraphConfig.basesPerGraph,
              end = start + GraphConfig.basesPerGraph;
          if (gotoBase && start <= gotoBase && gotoBase <= end) {
            $log.log('gotoBase triggered redraw:', gotoBase);
            redraw = true;
            schedule(true);
          }
          else if(fromBase && start <= fromBase && fromBase <= end) {
            redraw = true;
          }
        });
      }
    };
  })
  .directive('npactExtract', function(STATIC_BASE_URL, GraphConfig, NProfiler,
                               $log, Translater) {
    return {
      restrict: 'A',
      scope: { extract: '=npactExtract'},
      templateUrl: STATIC_BASE_URL + 'js/graphs/extract.html',
      link: function($scope, $element, $attrs, ctrl) {
        $log.log(NProfiler);
        window.currentScope = $scope;
        $scope.NProfiler = NProfiler;
        // TODO: if extract.complement, reverse the string
        // ?Maybe we get results either way?

        // I *think* we need to subtract 1 because coming from C/NCBI
        // the DNA is treated as 1 indexed
        // at least we got no results till we did this.
        $scope.ddna = NProfiler.ddna.slice($scope.extract.start - 1, $scope.extract.end - 1);
        Translater($scope.ddna, GraphConfig.mycoplasma, $scope.extract.complement)
          .then(function(res) {
            $scope.ddnaP = res.data && res.data.seq;
            $scope.ddnaRC = res.data && res.data.complement;
            $scope.query = $scope.extract.complement ? $scope.ddnaRC : $scope.ddna;
          });
        //window.NProfiler = NProfiler; window.GraphConfig = GraphConfig;
        $scope.config = GraphConfig;
      }
    };
  })
;
