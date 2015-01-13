angular.module('npact', ['ngMessages', 'sticky'])
  .value('K', Kinetic)
  .config(function($locationProvider) {
    $locationProvider.html5Mode({requireBase: false, enabled: true});
  })
  .run(function($rootScope, STATIC_BASE_URL) {
    $rootScope.STATIC_BASE_URL = STATIC_BASE_URL;
  })
;
