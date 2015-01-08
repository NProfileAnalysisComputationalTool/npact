angular.module('npact', ['ngMessages', 'sticky'])
  .value('K', Kinetic)
  .run(function($rootScope, STATIC_BASE_URL) {
    $rootScope.STATIC_BASE_URL = STATIC_BASE_URL;
  })
;
