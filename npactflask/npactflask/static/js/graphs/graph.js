angular.module('npact')
  .directive('npactGraphContainer', function() {
    'use strict';
    return {
      restrict: 'A',
      scope: true,
      controller: 'npactGraphContainerCtrl as ctrl'
    };
  })
  .controller('npactGraphContainerCtrl', function($scope, $element, $window, $log, $timeout,
                                           npactConstants, Utils, GraphConfig, Evt,
                                           GraphingCalculator) {
    'use strict';
    var $win = angular.element($window),
        getWidth =  function() { return $element.width(); },
        baseOpts = {
          width: null,
          m: null,
          tracks: null,
          onDragMove: function (dx) {
          $scope.$broadcast('offset', dx);
          },
          onDragEnd: function (dx) {
            $scope.$evalAsync(function () {

              //$log.log("Finished dragging, new offset is", GraphConfig.offset, 'upper',(GraphConfig.endbase|| 1e15), 'dx',Math.round(dx), 'new', GraphConfig.offset+Math.round(dx));
              GraphConfig.offset =
                Math.min(
                  (GraphConfig.endbase|| 1e15),
                  Math.max(0,GraphConfig.offset+Math.round(dx)));
              redraw();
            });
          },
          onRegionSelected: function (data) {
            $log.log('region-selected', data);
            $scope.$emit('region-selected', data); },
          onOrfSelected: function (data) {
            $log.log('ORF-selected', data);
            $scope.$emit('ORF-selected', data); },
          onHitSelected: function (data) {
            $log.log('hit-selected', data);
            $scope.$emit('hit-selected', data); }
        },
        updateMetrics = function () {
          baseOpts.width = getWidth();
          baseOpts.m = GraphingCalculator.chart(baseOpts);
          $scope.graphHeight = baseOpts.m.height;
        },
        redraw = function (clearing) { $scope.$broadcast(Evt.REDRAW, clearing); };

    $scope.baseOpts = baseOpts;
    this.winHeight = $win.height();

    var onResize = _.throttle(_.bind(function () {
      this.winHeight = $win.height();
      if(getWidth() !== baseOpts.width) {
        updateMetrics();
        redraw();
      }
      $scope.$applyAsync();
    }, this), 1000/30);
    $win.on('resize', onResize);
    $scope.$on('$destroy', function () { $win.off('resize', onResize); });
    $timeout(onResize, 100, false);

    $scope.$watch(function () { return GraphConfig.basesPerGraph; }, updateMetrics);
    $scope.$watchCollection(GraphConfig.activeTracks, function (val, old) {
      //Find headers and headerY
      baseOpts.tracks = val;
      updateMetrics();
      $scope.$broadcast('updateRowHeight', baseOpts.m.height);
      redraw({headers: true});
    });

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

    $scope.$watch(function() { return GraphConfig.nucleotides; },
                  function() { redraw({nProfiles: true, headers: true}); },
                  true);
    $scope.$watch(function() { return GraphConfig.colorBlindFriendly; }, redraw);
  })

  .directive('npactGraphScroller', function ($log, $window, Evt) {
    'use strict';
    return {
      restrict: 'A',
      require: 'npactGraphContainer',
      link: function ($scope, $element, $attrs, npactGraphContainerCtrl) {
        var $win = angular.element($window),
        graphHeight = 0, graphRowHeight = null,
        slack = 20, // how many pixels outside of viewport to render
        topOffset = 0,
        topIdx = 0, bottomIdx = 0;

        var updateVisibility = function () {
          if(!graphRowHeight) return;
          slack = graphRowHeight / 4;
          var scrollDist = $window.scrollY - topOffset - slack;
          topIdx = Math.floor(scrollDist / graphRowHeight);
          bottomIdx = topIdx + Math.ceil(
            (npactGraphContainerCtrl.winHeight + slack) / graphRowHeight);
        };

        $scope.visible = function(idx) { return idx >= topIdx && idx <= bottomIdx; };

        var draw = _.throttle(function () {
          updateVisibility();
          $scope.$broadcast(Evt.DRAW);
        }, 1000/20);

        $scope.$on('updateRowHeight', function ($evt, height) {
          $log.debug("Updating graphHeight to", height);
          topOffset = $element.offset().top;
          topOffset += _.parseInt($element.css('padding-top'));
          topOffset = Math.floor(topOffset);

          // Keep the topIdx at the top through the height change
          var delta = Math.max(0, topIdx) * (height - graphHeight);
          if (delta) { $window.scrollBy(0, delta); }
          graphHeight = height;
          $scope.$evalAsync(function () {
            //This code won't work until after the `$scope.graphHeight`
            //above has a chance to take effect. Hence the async.
            try {
              // Find the height including the padding+border for
              // graphRowHeight visibility calculations
              graphRowHeight = angular.element('.graph', $element).outerHeight();
            }
            catch(e) {
              graphRowHeight = 0; // There are no rows
            }
            draw();
          });
        });

        $win.on('scroll', draw);
        $scope.$on('$destroy', function () { $win.off('scroll', draw); });
      }
    };
  })


  .directive('npactKeyHandlers', function($window, $log, GraphConfig) {
    'use strict';
    return {
      restrict: 'A',
      link: function ($scope, $element, $attrs) {
        var $win = angular.element($window),
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
        $win.on('keydown', onKeyDown);
        $win.on('keyup', onKeyUp);
        $scope.$on('$destroy', function() {
          $win.off('keydown', onKeyDown);
          $win.off('keyup', onKeyUp);
        });
      }
    };
  })

  .directive('npactGraph', function (Grapher, Evt, GraphingCalculator, GraphConfig,
                              $log, $timeout, $window) {
    'use strict';
    return {
      restrict: 'A',
      scope: {
        startBase: '=',
        endBase: '=',
        visible: '&',
        graphOptions: '&'
      },
      link: function($scope, $element, $attrs) {
        var g = null, scheduledOn = null,
            visible = $scope.visible,
            id = $attrs.id,
            // redraw gets set for all graphs once (e.g. a new track
            // triggers broadcasts redraw), but only gets cleared as
            // the currently visible ones are drawn
            redraw = false;



        var draw = function(force) {
          var t1 = scheduledOn;
          scheduledOn = null;
          if(!redraw || (!force && !visible())) { return null; }
          var it;
          var opts = _.clone($scope.graphOptions());
          opts.startBase = $scope.startBase;
          opts.endBase = $scope.endBase;

          //However long it actually takes to draw, we have the
          //latest options as of this point
          redraw = false;
          return ((g = new Grapher($element, $scope, opts)))
            .draw(opts)
            .then(function () {
              $log.log("Finished draw at", opts.startBase, "in", new Date() - t1);
            })
            .catch(function() {
              //something went wrong, we will still need to redraw this
              $log.log('Need to redraw');
              redraw = true;
            });
        },
        schedule = function(force) {
          if(redraw && (force || (visible() && scheduledOn === null))) {
            scheduledOn = new Date();
            return $timeout(_.partial(draw, force), 0, false);
          }
          return null;
        },
        discard = function() { if(g) { g.destroy(); g = null; } },
        scrollToHere = function() {
          $log.log('scrolling to', $element);
          $window.scrollTo(0, $element.offset().top);
        };
        $scope.$on(Evt.DRAW, _.partial(schedule, false));
        $scope.$on(Evt.REDRAW, function(evt, cachesToClear) {
          redraw = true;
          if(g && cachesToClear) {
            if(cachesToClear.nProfiles) { g.clearProfilePoints(); }
            if(cachesToClear.headers) { g.clearLeftLayerCache(); }
          }
          schedule();
        });
        $scope.$on(Evt.REBUILD, function() {
          discard();
          startBase = $scope.startBase;
          endBase = $scope.endBase;
          redraw = true;
          schedule();
        });
        $scope.$on('offset', function(event, dx) {
          if(g && g.offset && visible()) { g.offset(dx); }
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
          if (_.isFinite(gotoBase) && $scope.startBase <= gotoBase && gotoBase <= $scope.endBase) {
            $log.log('gotoBase triggered redraw:', $scope.startBase, id);
            redraw = true;
            schedule();
            $timeout(scrollToHere, 0, false);
          }
          else if(_.isFinite(fromBase) && $scope.startBase <= fromBase && fromBase <= $scope.endBase) {
            redraw = true;
          }
        });
      }
    };
  })
;
