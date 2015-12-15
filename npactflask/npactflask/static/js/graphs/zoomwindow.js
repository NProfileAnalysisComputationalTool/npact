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
      $scope.$on('region-selected',
                 _.bind(function (evt, data) { this.popup(data); }, this));
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
  .controller('ZoomWindowCtrl', function ($scope, $log, FocusData, focusData, getPhase,
                                   GraphConfig, Utils) {
    var type = focusData.type;
    this.data = focusData;
    if(type === 'region') {
      this.item = {
        name: "New",
        start: Utils.floor3(focusData.start),
        end: Utils.ceil3(focusData.end),
        complement: 0
      };
      this.item.phase = getPhase(this.item);
    }
    else {
      this.item = focusData.item;
    }
    this.track = focusData.track || _.first(GraphConfig.activeTracks);

    var margin = Utils.orderOfMagnitude(this.item.end - this.item.start, -1);
    this.startBase = Utils.floor3(Math.max(this.item.start - margin, 0));
    this.endBase = Utils.ceil3(Math.min(this.item.end + margin, GraphConfig.endBase));


    $log.log("Focusing on", focusData);
    $scope.$on('offset', _.bind(function ($evt, dx) {
      $evt.stopPropagation();
      $log.log("Updating offsets", dx);
      dx = Math.round(dx);
      this.startBase = Utils.floor3(this.startBase + dx);
      this.endBase = Utils.ceil3(this.endBase + dx);
      $scope.$apply();
    }, this));
    var redraw = function () { $scope.$broadcast('redraw'); };
    $scope.$watchCollection("zw.item", redraw);
    $scope.$watch('zw.track', redraw);


    this.buildPermalink = function () {
      //TODO: Params to this function
      return FocusData.serialize();
    };
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
