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
      $log.log('Serializing: ', focusData, focusData.track.filename, focusData.item.cdsidx);
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
  .service('ZoomWindowHandler', function ($log, $uibModal, FocusData, STATIC_BASE_URL, Evt, Utils) {
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
                complement: 0
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
      FocusData.serialize(focusData);
      var modalDefaults = {
        templateUrl: STATIC_BASE_URL + 'js/graphs/zoomwindow.html',
        controller: 'ZoomWindowCtrl',
        controllerAs: 'zw',
        bindToController: true,
        resolve: {
          focusData: function () { return focusData; }
        },
        animation: true,
        size: 'lg'
      };
      GraphConfig.zoomwindow = $uibModal.open(_.assign(modalDefaults, modalopts));
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
  .controller('ZoomWindowCtrl', function ($scope, $log, FocusData, focusData, getPhase,
                                   GraphConfig, Utils, $location, $uibModalInstance,
                                   $window ) {
    //$log.log('ZoomWindowCtrl', focusData);
    var self = this;
    this.item = focusData.item;
    this.track = focusData.track;
    this.cancel = function() {
      console.log('Cancelling track edit');
      $window.location.reload();
    };
    this.save = function(track) {
      console.log('saving track', track);
      $uibModalInstance.dismiss();
      return track.save();
    };
    this.delete = function(orf, track) {
      console.log('deleting orf', orf, track);
      self.track.remove(orf);
      return self.save(self.track);
    };

    var length = this.item.end - this.item.start;
    var margin = Utils.ceil3(Utils.orderOfMagnitude(length, -1));
    this.startBase = Math.max(this.item.start - margin, 0);
    this.endBase = Math.min(this.item.end + margin, GraphConfig.endBase);

    $scope.$watch(_.partial(getPhase, this.item),
                  _.bind(function (phase) { this.item.phase = phase; }, this));

    //$log.log("Focusing on", focusData);
    $scope.$on('offset', _.bind(function ($evt, dx) {
      $evt.stopPropagation();
      $log.log("Updating offsets", dx);
      dx = Utils.round3(dx);
      this.startBase = this.startBase + dx;
      this.endBase = this.endBase + dx;
      $scope.$apply();
    }, this));
    var redraw = function () { $scope.$broadcast('redraw'); };
    $scope.$watchCollection("zw.item", redraw);
    $scope.$watch('zw.track', redraw);
    //Keep phase up to date with the end

    this.buildPermalink = function () {
      //TODO: Params to this function
      window.focusData = focusData;
      return FocusData.serialize(focusData);
    };

  })

  .directive('npactSingleGraph', function($log, $timeout,
                                   GraphConfig, GraphingCalculator, Grapher, Evt) {
    'use strict';
    return {
      restrict: 'A',
      scope: {
        startBase: '=',
        endBase: '='
      },
      link: function($scope, $element, $attrs) {
        var g = null;
        var draw = function() {
          $log.debug("Redrawing single graph");
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
        var schedule = _.debounce(function () {$timeout(draw);}, 400);
        $scope.$on('redraw', schedule);
        $scope.$watchGroup(['startBase', 'endBase'], schedule);
        $scope.$watch(function() {return $element.width();}, schedule);
      }
    };
  })

  .directive('npactProteinTranslation', function ($log, TranslatePath) {
    return {
      restrict: 'E',
      scope: {
        start: '=',
        end: '=',
        complement: '='
      },
      link: function ($scope, $element) {
        $scope.$watchGroup(
          ['start', 'end', 'complement'],
          function () {
            $element.text('');
            TranslatePath($scope.start, $scope.end, $scope.complement)
              .then(function (data) {
                ddnaP = data.trans;
                ddna = data.seq;
                $element.text(ddnaP);
              });
          }
        );
      }
    };
  })
;
