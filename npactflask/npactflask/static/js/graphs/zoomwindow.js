angular.module('npact')
  .service('FocusData', function ($log, $location, GraphConfig) {
    /**
     * focusData interface:
     *
     *  type       : "(hit|orf|region)"
     *  track      : Reference to a Track object
     *  item       : Data Object inside a track
     *  start      : start base of zoom region
     *  end        : end base of zoom region
     *  complement : (0|1), indicates strand
     *  phase      : (0|1|2)
     *  name       : name of item
     *
     */
    this.deserialize = function () {
      if(!(GraphConfig.zoomTrack && GraphConfig.zoomIdx >= 0))
        return null;
      var track = GraphConfig.findTrack(GraphConfig.zoomTrack);
      var cds = track && track.data[GraphConfig.zoomIdx];
      //$log.log('loading', GraphConfig.zoomTrack, GraphConfig.zoomIdx,track, cds);
      var focusData = track && cds && {
        track: track,
        item: cds
      };
      return focusData;
    };
    this.serialize = function (focusData) {
      $log.log('Serializing: ', focusData,
               focusData && focusData.track.filename,
               focusData && focusData.item.cdsidx);
      if(!focusData) return $location.absUrl();
      if(focusData.track)
        GraphConfig.zoomTrack = focusData.track.filename;
      if(focusData.item && focusData.item.cdsidx >= 0)
        GraphConfig.zoomIdx = focusData.item.cdsidx;
      return $location.absUrl();
    };
    this.clearQuery = function () {
      $log.debug("Clearing FocusData from querystring");
      GraphConfig.zoomTrack=null;
      GraphConfig.zoomIdx = null;
    };
  })
  .service('ZoomWindowHandler', function ($log, $uibModal, FocusData, STATIC_BASE_URL, Evt, Utils, $timeout) {
    var self = this;
    self.register = function ($scope) {
      $scope.$on('ORF-selected',
                 function (evt, data) { self.popup(data, null, $scope); });
      $scope.$on('hit-selected',
                 function (evt, data) { self.popup(data, null, $scope); });
      $scope.$on('region-selected', function (evt, data) {
        var start = Utils.floor3(data.start),
            //end is inclusive so need to subtract one for the length to be mod3
            end = Utils.ceil3(data.end) -1,
            item = {
              name: "ORF from region",
              start: start,
              end: end,
              complement: 0,
              selected: true
            },
            track = _.find(GraphConfig.activeTracks(), {name: 'New ORFs'}) ||
                _.first(GraphConfig.activeTracks());
        track.add(item);
        self.popup({ item: item, track: track }, null, $scope);
      });
    };

    self.maybePopup = function ($scope) {
      var fd = FocusData.deserialize();
      if(fd) return self.popup(fd, null, $scope);
      return null;
    };
    self.popup = function (focusData, modalopts, $scope) {
      $log.log("popping!", focusData, modalopts);
      if(focusData.item == null) return null;
      GraphConfig.clearORFSelection();
      focusData.item.selected = true;
      FocusData.serialize(focusData);
      if(GraphConfig.zoomwindow){
        GraphConfig.zoomwindow.$scope.data = focusData;
        $scope.$broadcast(Evt.REDRAW);
        return GraphConfig.zoomwindow.result;
      }

      var child = $scope.$new();
      child.data = focusData;
      child.extendedWindow=false;
      var modalDefaults = {
        templateUrl: STATIC_BASE_URL + 'js/graphs/zoomwindow.html',
        controller: 'ZoomWindowCtrl',
        controllerAs: 'zw',
        bindToController: true,
        scope: child,
        animation: true,
        size: 'lg'
      };
      GraphConfig.zoomwindow = $uibModal.open(_.assign(modalDefaults, modalopts));
      $timeout(function() {
        $('.modal-content').resizable({
          alsoResize: ".modal-header, .modal-body, .modal-footer"
        });
      }, 100);
      GraphConfig.zoomwindow.result.catch(function() {
        // when exiting the window, data may have changed that requires a redraw

        $scope.$broadcast(Evt.REDRAW);
      }).finally( function() {
        FocusData.clearQuery();
        GraphConfig.zoomwindow = null;
      });
      return GraphConfig.zoomwindow.result;
    };
  })
  .controller('ZoomWindowCtrl', function ($scope, $log, FocusData, getPhase,
                                   GraphConfig, Utils, $location, $uibModalInstance,
                                   $window, $timeout ) {
    //$log.log('ZoomWindowCtrl', focusData);
    GraphConfig.stack=[];
    GraphConfig.zoomwindow.$scope = $scope;
    GraphConfig.Utils = Utils;
    $scope.GraphConfig = GraphConfig;
    var type = $scope.data.type;
    var self = this;
    this.cancel = function() {
      $log.log('Cancelling track edit');
      $window.location.reload();
    };
    this.extendWindowClick = function(e) {
      var $el = angular.element(e.target);
      var dir = $el.is('.left')?'start':'end',
          comp = $scope.data.item.complement;
      $log.log('zoom-btn clicked:', dir, comp, $el);
      $scope.extendedWindow = true;
      if((dir === 'start' && !comp)
         || (dir === 'end' && comp))
        $scope.startBase = Math.max(0, $scope.startBase-GraphConfig.graphMargin);
      else
        $scope.endBase = Math.min(
          GraphConfig.endBase, $scope.endBase+GraphConfig.graphMargin);
    };

    this.save = function(track) {
      $log.log('saving track', track);
      $uibModalInstance.dismiss();
      return track.save();
    };
    this.delete = function(orf, track) {
      $log.log('deleting orf', orf, track);
      self.track.remove(orf);
      return self.save(self.track);
    };

    // TODO: I think this can be removed
    if(type === 'region') {
      $scope.data = {
        name: "New",
        start: Utils.floor3($scope.data.start),
        end: Utils.ceil3($scope.data.end) -1, //the end is an inclusive range
        complement: 0
      };
    }
    else {
      this.item = $scope.data.item;
    }
    $scope.startBaseInp = $scope.data.item.start;
    $scope.endBaseInp = $scope.data.item.end;
    this.track = $scope.data.track || _.first(GraphConfig.activeTracks);
    $scope.setGraphBounds = function() {
      var len = $scope.data.item.end - $scope.data.item.start;
      if(!$scope.extendedWindow){
        $scope.startBase = Math.max(
          0, // from 1 behind the start back 1 margin
          (($scope.data.item.start-1) - GraphConfig.graphMargin ));
        $scope.endBase = Math.min(
          GraphConfig.endBase,// from 1 past the end + a margin
          $scope.data.item.end + 1 + GraphConfig.graphMargin);
      }
      //$log.log('setGraphBounds', $scope.data.item, len, GraphConfig.graphMargin, $scope.startBase, $scope.endBase, $scope.extendedWindow);
    };
    $scope.setGraphBounds();

    $scope.$watch(_.partial(getPhase, $scope.data.item),
                  _.bind(function (phase) { $scope.data.item.phase = phase; }, this));

    //$log.log("Focusing on", focusData);
    $scope.$on('offset', _.bind(function ($evt, dx) {
      $evt.stopPropagation();
      //$log.log("Updating offsets", dx);
      dx = Utils.round3(dx);
      $scope.startBase = $scope.startBase + dx;
      $scope.endBase = $scope.endBase + dx;
      $scope.$apply();
    }, this));
    $scope.redraw = _.debounce(function () {
      //GraphConfig.stack.push('redraw'); $log.log('stack', GraphConfig.stack);
      $scope.setGraphBounds();
      $scope.$broadcast('redraw');
    },50);

    $scope.inbounds = function(inp) {
      if(!inp || isNaN(inp)) return false;
      inp = Number(inp);
      if(inp < 0 || inp > GraphConfig.endBase) return false;
      return true;
    };
    $scope._lastGoodStart=null;
    $scope.startChange = function(newStart, oldStart) {
      $log.log('startChange', arguments);
      // start validation
      if(!$scope.inbounds(newStart) || newStart == $scope.data.item.end){
        $log.log('resetting invalid start input');
        $scope.startBaseInp = oldStart || $scope._lastGoodStart;
        return;
      }
      else if(newStart > $scope.data.item.end){
        $scope.startBaseInp = $scope.endBaseInp;
        $scope.endBaseInp = newStart;
        return;
      }
      else if(newStart == $scope.data.item.start){
        // nothing to do
        return;
      }
      // End Validation
      $scope._lastGoodStart = newStart;
      var diff = newStart-oldStart,
          m = (diff)%3;
      //$log.log('Mod',newStart, ' newend', $scope.data.item.end+m, diff, m);
      $scope.data.item.start = newStart;
      // keep triplets in sync
      $scope.endBaseInp = $scope.data.item.end = $scope.data.item.end+m;
      $scope.extendedWindow = false;
      $scope.redraw();
      return;
    };
    $scope._lastGoodEnd=null;
    $scope.endChange = function(newEnd, oldEnd) {
      $log.log('endChange', arguments);
      if(!$scope.inbounds(newEnd) || $scope.data.item.start == newEnd){
        $log.log('resetting invalid end input');
        $scope.endBaseInp = oldEnd || $scope._lastGoodEnd;
        return;
      }else if(newEnd < $scope.data.item.start){
        $scope.endBaseInp = $scope.startBaseInp;
        $scope.startBaseInp = newEnd;
        return;
      }
      else if(newEnd == $scope.data.item.end){
        // nothing to do
        return;
      }
      // End Validation
      var diff = newEnd-oldEnd,
          m = (diff)%3;
      $scope.extendedWindow = false;
      //$log.log('Mod',newEnd, ' newend', newEnd+m, diff, m);
      $scope.startBaseInp = $scope.data.item.start = $scope.data.item.start+m;
      $scope.data.item.end = newEnd;
      $scope.redraw();
      return;
    };
    $scope.orfChange = function() {
      //$log.log("orfChange");
      $scope.startBaseInp = $scope.data.item.start;
      $scope.endBaseInp = $scope.data.item.end;
      $scope.redraw();
    };

    $scope.$watch('data.item.start', $scope.orfChange);
    $scope.$watch('data.item.end', $scope.orfChange);
    $scope.$watch('startBaseInp', $scope.startChange);
    $scope.$watch('endBaseInp', $scope.endChange);
    $scope.$watch('data.item.complement', $scope.redraw);
    $scope.$watch('zw.track', $scope.redraw);
    //Keep phase up to date with the end

    this.buildPermalink = function () {
      //TODO: Params to this function
      return FocusData.serialize();
    };

  })

  .directive('npactSingleGraph', function($log, $timeout,
                                   GraphConfig, GraphingCalculator, Grapher, Evt) {
    'use strict';
    return {
      restrict: 'A',
      scope: {
        startBase: '=',
        endBase: '=',
        extendedWindow: '='
      },
      link: function($scope, $element, $attrs) {
        GraphConfig.graphScope = $scope;
        $scope.extendedWindow = false;
        var g = null;
        var draw = function() {
          //GraphConfig.stack.push('draw'); $log.log('stack', GraphConfig.stack);
          $log.debug("Drawing single graph", $scope.startBase, $scope.endBase);
          var opts = {
            width: $element.width(),
            m: null,
            tracks: GraphConfig.activeTracks(),
            startBase: $scope.startBase,
            endBase: $scope.endBase,
            onDragEnd: function (dx) { $scope.$emit('offset', dx); },
            onRegionSelected: function (data) { $scope.$emit('region-selected', data); },
            onOrfSelected: function (data) { $scope.$emit('ORF-selected', data); },
            onHitSelected: function (data) { $scope.$emit('hit-selected', data); }
          };
          opts.m = GraphingCalculator.chart(opts);
          $scope.graphHeight = opts.m.height;
          $log.debug("Redrawing with params:", opts);
          if(g) { g.destroy(); }
          g = new Grapher($element, $scope, opts);
          return g.draw();
        };
        var schedule = _.debounce(function () {$timeout(draw);}, 500);
        $scope.$on('redraw', schedule);
        $scope.$watchGroup(['startBase', 'endBase'], schedule);
        $scope.$watch(function() {return $element.width();}, schedule);
      }
    };
  })

  .directive('npactProteinTranslation', function (
    $log, TranslatePath, GraphConfig, Evt, Utils, CodonFinder) {
    'use strict';
    return {
      restrict: 'E',
      scope: {
        item: '=',
        start: '=',
        end: '=',
        complement: '='
      },
      link: function ($scope, $element, attrs, ctrl) {
        $element.addClass('npact-protein-translation');

        $scope._debounced_applied = function(fn) {
          return _.debounce(function($el){
            var rtn = null;
            $scope.$apply(function() {
              rtn = fn($el);
            });
            return rtn;
          },25);
        };

        $scope._stopClicked = $scope._debounced_applied(function($el){
          var clicked = $el && Number($el.data('index')),
              prevStop = CodonFinder.startOfNextCodon(
                CodonFinder.findPrevStopCodon(clicked-1, $scope.complement)+1),
              prevStart = CodonFinder.findPrevStartCodon(clicked-1, $scope.complement, prevStop)+1,
              end, cur, start;
          cur = $scope.item.start;
          if(!$scope.item.complement){
            start = cur;
            end = CodonFinder.endOfThisCodon(clicked);
            // our start should never be after our end
            if(!start || start > end) start = prevStart;
            if(!start || start > end || (prevStop && prevStop > start))
              start = prevStop;
            if(!start || start > end) start = Math.max(end - 99, 0);
          }else{
            // shift the window for complements
            //$log.log('!!!=',$scope.item.start,'+ (',$scope.item.end,' - ', clicked,'= ',($scope.item.end - clicked),') = ', ($scope.item.start + ($scope.item.end - clicked))-2);
            start = ($scope.item.start + ($scope.item.end - clicked))-2;
            end =  $scope.item.end;

          }

          $log.log('stop clicked', clicked, 'prevstop',prevStop, 'prevstart', prevStart, 'cur: ',cur, 'start', start, 'end', end);
          if(((end - start) % 3) !== 2){
            $log.log('out of phase shift!  stop clicked',
                     clicked, 'prevstop',prevStop, 'prevstart', prevStart,
                     'cur: ',cur, 'start', start, 'end', end, 'mod',((end - start) % 3));
            return;
          }
          $scope.item.end = end;
          if(start) $scope.item.start = start;
        });
        $scope._startClicked = $scope._debounced_applied(function($el){
          var end=$scope.item.end, cur = $scope.item.start,
              start=cur,
              clicked = $el.data('index'),
              eidx = $scope.item.end,
              nextStop = CodonFinder.findNextStopCodon(clicked-1, $scope.complement)+1;
          if(!$scope.item.complement){
            start = clicked;
            end = (nextStop && CodonFinder.endOfThisCodon(nextStop)) ||
                Math.min(start + GraphConfig.graphMargin,
                         GraphConfig.endBase);
          }else{
            end = clicked;
            //$log.log('start:',$scope.item.start, 'clicked:', clicked, 'd:',d, 'end:', end);
          }
          if(((end - start) % 3) !== 2){
            $log.log('out of phase shift!  start clicked',
                     clicked, 'start', start, 'end', end, 'curstart', cur);
            return;
          }
          if(start) $scope.item.start = start;
          if(end) $scope.item.end = end;
          $log.log('start clicked', clicked, 'cur: ',cur, 'start', start, 'end', end);
        });
        $scope._highlightDnaP= function(ddnaP) {
          if(!_.isArray(ddnaP)) ddnaP = ddnaP.split('');
          var comp = $scope.item.complement;
          return _.map(ddnaP, function(v, k){
            // the DNA indexes run from 1 but our DNA string is indexed from 0
            // to make the display line up, just add one when displaying the inde
            var idx = ((k*3) + $scope.start), c;
            if(idx >= 9250 && idx <= 9260)
              $log.log(idx, comp?($scope.end - k*3):idx);
            if(comp) idx = ($scope.end - k*3);
            // CodonFinder deals in string indexes (0 based)

            if((c = CodonFinder.isStop(idx-1, comp))){
              return '<span class="codon stop"  title="stop: '+c+
                ' ('+idx+')'+ '" data-index="'+idx+'">'+v+'</span>';
            }else if((c = CodonFinder.isStart(idx-1, comp))){
              return '<span class="codon standard start" title="start: '+ c +
                ' ('+idx+')'+ '" data-index="'+idx+'">'+v+'</span>';
            }else if((c = CodonFinder.isAltStart(idx-1, comp))){
              return '<span class="codon rare start" title="alt start:'+c+
                ' ('+idx+')'+ '" data-index="'+idx+'">'+v+'</span>';
            }
            else return '<span title="'+idx+'">'+v+'</span>';
          });
        };

        $scope.$watchGroup(
          ['start', 'end', 'complement'],
          function () {
            $element.text('');
            TranslatePath($scope.start, $scope.end, $scope.complement)
              .then(function (data) {
                if($scope.item && $scope.item.start == $scope.start
                   && $scope.item.end == $scope.end){
                  $scope.item.ddnaP = data.trans;
                  $scope.item.ddna = data.seq;
                }
                var ddnaP = $scope._highlightDnaP(data.trans);

                var html = ['<div>('+$scope.start+"-<wbr>"+$scope.end+")</div>"]
                  .concat(ddnaP);
                $element.html(html);
                $element.bind('click', function(e){
                  var $el = angular.element(e.target);
                  // $log.log('translation clicked', e, e.target, $scope.item, $el.is('.codon'));
                  GraphConfig.zoomwindow.$evttarget = $el;
                  GraphConfig.zoomwindow.$evtscope = $scope;
                  if($el.is('.codon') && $scope.item){
                    // $log.log('codon clicked', $el, $el.data('index'), $scope.item);
                    // GraphConfig.stack.push('codon clicked'); $log.log('stack', GraphConfig.stack);
                    if ($el.is('.stop'))
                      $scope._stopClicked($el);
                    else if ($el.is('.start'))
                      $scope._startClicked($el);
                  }
                });
              });
          }
        );
      }
    };
  })
;
