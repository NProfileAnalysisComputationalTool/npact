angular.module('npact')

// {type: 'region', start: 0, end: 10}
// {type: 'hit}

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
      var cds = track.data[GraphConfig.zoomIdx];
      //$log.log('loading', GraphConfig.zoomTrack, GraphConfig.zoomIdx,track, cds);
      var focusData = {
        type: 'orf',
        track: track,
        start: cds.start,
        end: cds.end,
        complement: cds.complement,
        phase: cds.phase,
        name: cds.name,
        item: cds
      };
      return focusData;
    };
    this.serialize = function (focusData) {
      $log.log('Serializing: ', focusData, focusData.track.filename, focusData.item.cdsidx);
      if(!focusData) return $location.absUrl();
      if(focusData.track)
        GraphConfig.zoomTrack = focusData.track.filename;
      if(focusData.item && focusData.item.cdsidx>=0)
        GraphConfig.zoomIdx = focusData.item.cdsidx;
      return $location.absUrl();
    };
    this.clearQuery = function () {
      //$log.log("Clearing FocusData from querystring");
      GraphConfig.zoomTrack=null;
      GraphConfig.zoomIdx = null;
    };
  })
  .service('ZoomWindowHandler', function ($log, $uibModal, FocusData, STATIC_BASE_URL) {

    this.register = function ($scope) {
      $scope.$on('ORF-selected',
                 _.bind(function (evt, data) { this.popup(data); }, this));
      $scope.$on('hit-selected',
                 _.bind(function (evt, data) { this.popup(data); }, this));
      $scope.$on('region-selected',
                 _.bind(function (evt, data) { this.popup(data); }, this));
    };

    this.maybePopup = function () {
      var fd = FocusData.deserialize();
      //console.log('maybe popup', fd);
      if(fd) this.popup(fd);
    };
    this.popup = function (focusData, modalopts) {
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
      GraphConfig.zoomwindow.result.finally( function() {
        FocusData.clearQuery();
        GraphConfig.zoomwindow = null;
      });
    };
  })
  .controller('ZoomWindowCtrl', function ($scope, $log, FocusData, focusData, getPhase,
                                   GraphConfig, Utils, $location, $uibModalInstance,
                                   $window, Evt
                                  ) {
    //$log.log('ZoomWindowCtrl', focusData);
    var type = focusData.type;
    var self = this;
    this.data = focusData;
    this.cancel = function() {
      console.log('Cancelling track edit');
      $window.location.reload();
    };
    this.save = function(track) {
      console.log('saving track', track);
      return track.save().then(function(newtr){
        self.data.track = newtr;
        FocusData.serialize(self.data);
        var idx = null;
        _.each(GraphConfig.tracks,function(t, i) {
          if(t.filename==track.filename) idx=i;
        });
        if(!idx) idx=GraphConfig.tracks.length;
        GraphConfig.tracks.splice(idx, 1, newtr);
        $scope.$broadcast(Evt.REDRAW);
        $uibModalInstance.dismiss();
      });
    };
    this.delete = function(orf, track) {
      console.log('deleting orf', orf, track);
      self.track.remove(orf);
      return self.save(self.track);
    };
    if(type === 'region') {
      this.item = {
        name: "New",
        start: Utils.floor3(focusData.start),
        end: Utils.ceil3(focusData.end) -1, //the end is an inclusive range
        complement: 0
      };
    }
    else {
      this.item = focusData.item;
    }
    this.track = focusData.track || _.first(GraphConfig.activeTracks);

    var margin = Utils.ceil3(Utils.orderOfMagnitude(this.item.end - this.item.start, -1));
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
        var scheduled = false,
            schedule = function () {
              if(scheduled) return;
              scheduled = true;
              $timeout(function () {
                draw();
                scheduled=false;
              });
            };
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
