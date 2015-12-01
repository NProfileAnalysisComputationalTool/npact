angular.module('npact')
  .service('ZoomWindowHandler', function ($log, $uibModal, STATIC_BASE_URL) {

    this.register = function ($scope) {
      $scope.$on('ORF-selected',
                 _.bind(function (evt, data) { this.popup(data); }, this));
      $scope.$on('region-selected');
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
  .controller('ZoomWindowCtrl', function ($log, focusData, GraphConfig) {
    $log.log("Focusing on", focusData);

  })
;
