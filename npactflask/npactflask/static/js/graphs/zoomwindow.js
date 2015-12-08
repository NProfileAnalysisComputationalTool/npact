angular.module('npact')
  .service('FocusData', function ($log, $location, GraphConfig) {
    /**
     * focusData interface:
     *
     *  type       : "(hit|orf|region)"
     *  track      : Reference to a Track object
     *  item       : Data Object  inside a track
     *  start      : start base of zoom region
     *  end        : end base of zoom region
     *  complement : (0|1), indicates strand
     *  phase      : (0|1|2)
     *  name       : name of item
     *
     */
    this.deserialize = function () {
      var search = $location.search();
      if(search['zoom-start'] === undefined) return null;
      var focusData = {
        type: search['zoom-type'],
        track: GraphConfig.findTrack(search['zoom-trackName']),
        start: Number(search['zoom-start']),
        end: Number(search['zoom-end']),
        complement: Number(search['zoom-complement']),
        phase: Number(search['zoom-phase']),
        name: search['zoom-name']
      };
      if(focusData.track && focusData.type === 'orf') {
        focusData.track.findByName(focusData.name);
      }
      return focusData;
    };
    this.serialize = function (focusData) {
      $location.search('zoom-type', focusData.type);
      if(focusData.track) $location.search('zoom-trackName', focusData.track.name);
      $location.search('zoom-start', focusData.start);
      $location.search('zoom-end', focusData.end);
      $location.search('zoom-complement', focusData.complement);
      $location.search('zoom-phase', focusData.phase);
      $location.search('zoom-name', focusData.name);
      $log.debug("Finished serializing focusdata, final url", $location.absUrl());
      return $location.absUrl();
    };
    this.clearQuery = function () {
      $log.log("Clearing FocusData from querystring");
      this.serialize({});
      $location.search('zoom-trackName', undefined);
    };
  })
  .service('ZoomWindowHandler', function ($log, $uibModal, FocusData, STATIC_BASE_URL) {

    this.register = function ($scope) {
      $scope.$on('ORF-selected',
                 _.bind(function (evt, data) { this.popup(data); }, this));
      $scope.$on('hit-selected',
                 _.bind(function (evt, data) { this.popup(data); }, this));
      $scope.$on('region-selected');
    };

    this.maybePopup = function () {
      var fd = FocusData.deserialize();
      if(fd) this.popup(fd);
    };
    this.popup = function (focusData, modalopts) {
      $log.log("popping!", focusData, modalopts);
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
      $uibModal.open(_.assign(modalDefaults, modalopts)).result
        .finally( _.bind(FocusData.clearQuery, FocusData));
    };
  })
  .controller('ZoomWindowCtrl', function ($scope, $log, FocusData, focusData,
                                   GraphConfig, Utils, TranslatePath) {
    var type = focusData.type;
    this.data = focusData;
    var length = this.data.end - this.data.start;
    var margin = Utils.orderOfMagnitude(length, -1);
    this.startBase = Math.max(this.data.start - margin, 0);
    this.endBase = Math.min(this.data.end + margin, GraphConfig.endBase);

    $log.log("Focusing on", focusData);
    $scope.$on('offset', _.bind(function ($evt, dx) {
      $evt.stopPropagation();
      dx = Math.round(dx);
      this.startBase += dx;
      this.endBase += dx;
      $log.log("Updating offsets", dx);
      $scope.$apply();
    }, this));

    TranslatePath(this.data.start, this.data.end, this.data.complement)
      .then(_.bind(function (data) {
        this.ddnaP = data.trans;
        this.ddna = data.seq;
      }, this));
    TranslatePath(this.startBase, this.data.start, this.data.complement)
      .then(_.bind(function (data) {
        this.leftDdnaP = data.trans;
      }, this));
    TranslatePath(this.data.end, this.endBase, this.data.complement)
      .then(_.bind(function (data) {
        this.rightDdnaP = data.trans;
      }, this));
  })

  .directive('zwPermalink', function ($log, $location, FocusData) {
    return {
      restrict: 'A',
      scope: {
        data: '&zwPermalink'
      },
      link: function ($scope, $element, attrs) {
        $element.on('mouseenter focus', function () {
          $element.attr('href', FocusData.serialize($scope.data()));
        });
      }
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
          var opts = {
            width: $element.width(),
            m: null,
            tracks: GraphConfig.activeTracks(),
            startBase: $scope.startBase,
            endBase: $scope.endBase,
            onDragEnd: function (dx) {
              $scope.$emit('offset', dx);
            }
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
              draw();
            };
        $scope.$watchGroup(['startBase', 'endBase', 'zw.data.complement'], schedule);
        $scope.$watch(function() {return $element.width();}, schedule);
      }
    };
  })

;
