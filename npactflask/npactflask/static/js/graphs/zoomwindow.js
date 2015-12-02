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
      $uibModal.open(_.assign(modalDefaults, modalopts));
    };
  })
  .controller('ZoomWindowCtrl', function ($log, FocusData, focusData, GraphConfig) {
    var type = focusData.type;
    this.data = focusData;
    $log.log("Focusing on", focusData);
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
;
