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
        draw = function() {
          updateVisibility();
          $scope.$broadcast(Evt.DRAW);
        },
        redraw = function() {
          if (!ready()) return;  //too early to do anything
          updateVisibility();
          $scope.$broadcast(Evt.REDRAW);
        },
        rebuild = function() {
          if (!ready()) return;  //too early to do anything
          updateVisibility();
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
        $timeout(rebuild);
      });

    $scope.$watch(function() { return GraphConfig.nucleotides; }, function() {
      $scope.$broadcast('n-profile');
    }, true);
    $scope.$watch(function() { return GraphConfig.colorBlindFriendly; }, redraw);

    $scope.$watchCollection(GraphConfig.activeTracks, function(val, old) {
      //Find headers and headerY
      baseOpts.tracks = val;
      updateMetrics();
      updateRowHeight(baseOpts.m.height);
    });

    /***  Scrolling and Graph Visibility management ***/
    var $win = angular.element($window),
        winHeight = $win.height(),
        slack = 20, // how many pixels outside of viewport to render
        topOffset = 0,
        topIdx = 0, bottomIdx = 0,
        graphRowHeight = 0,
        updateRowHeight = function(height) {
          topOffset = $element.offset().top;
          topOffset += _.parseInt($element.css('padding-top'));
          topOffset = Math.floor(topOffset);

          // Keep the topIdx at the top through the height change
          var delta = Math.max(0, topIdx) * (height - $scope.graphHeight);
          if (delta) {
            $window.scrollBy(0, delta);
          }
          $scope.graphHeight = height;  //set the inner height of the canvas container
          $timeout(function() {
            //This code won't work until after the `$scope.graphHeight`
            //above has a chance to take effect. Hence the $timeout.
            try {
              // Find the height including the padding+border for
              // graphRowHeight visibility calculations
              graphRowHeight = angular.element('#graph_0', $element).outerHeight();
            }
            catch(e) {
              graphRowHeight = 0; // There are no rows
            }
            redraw();
          });
        },
        updateVisibility = function() {
          if(!graphRowHeight) return;
          var scrollDist = $window.scrollY - topOffset - slack;
          topIdx = Math.floor(scrollDist / graphRowHeight);
          bottomIdx = topIdx + Math.ceil((winHeight + 2* slack) / graphRowHeight);
        },
        onScroll = draw,
        onResize = function() {
          winHeight = $win.height();
          if(getWidth() !== baseOpts.width) {
            topOffset = $element.offset().top;
            baseOpts.width = getWidth();
            updateMetrics();
            redraw();
          }
          else {
            //If the width didn't change then its the same as scrolling
            onScroll();
          }
        },
        onKeyDown = _.throttle(function(event) {
          //Ignore keys while the user is typing in input controls
          if(event.target.nodeName == "INPUT") return;
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

    $scope.$on('printresize', function(event, printing) {
      if(printing) {
        $element.css({width: '7in'});
        onResize();
      }
      else {
        $element.css({width: 'auto'});
        onResize();
      }
    });

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
    'use strict';
    return {
      restrict: 'A',
      scope: {},
      templateUrl: STATIC_BASE_URL + 'js/graphs/graph-container.html',
      controller: 'npactGraphContainerCtrl as ctrl'
    };
  })

  .directive('npactGraph', function (Grapher, Evt, GraphingCalculator, GraphConfig,
                              $log, $timeout, $window) {
    'use strict';
    return {
      restrict: 'A',
      require: '^npactGraphContainer',
      link: function($scope, $element, $attrs, ctrl) {
        var g = null,
            visible = ctrl.visible,
            idx = $attrs.idx, id = '#graph_' + idx,
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
            discard = function() { if(g) { g.destroy(); g = null; } },
            scrollToHere = function() {
              $log.log('scrolling to', id);
              $window.scrollTo(0, $(id).offset().top);
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
          var startBase = Number($attrs.startBase);
          var endBase = startBase + GraphConfig.basesPerGraph;
          if (_.isFinite(gotoBase) && startBase <= gotoBase && gotoBase <= endBase) {
            $log.log('gotoBase triggered redraw:', startBase, id);
            redraw = true;
            schedule();
            $timeout(scrollToHere());
          }
          else if(_.isFinite(fromBase) && startBase <= fromBase && fromBase <= endBase) {
            redraw = true;
          }
        });
      }
    };
  })
  .directive('npactExtract', function(STATIC_BASE_URL, GraphConfig, DDNA,
                               $log, TranslatePath) {
    'use strict';
    return {
      restrict: 'A',
      scope: { extract: '=npactExtract'},
      templateUrl: STATIC_BASE_URL + 'js/graphs/extract.html',
      link: function($scope, $element, $attrs, ctrl) {
        var e = $scope.extract;
        TranslatePath(e.start, e.end, e.complement).then(function (data) {
          $scope.ddnaP = data.trans;
          $scope.ddna = data.seq;
        });
      }
    };
  })
;
