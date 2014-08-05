angular.module('npact')
  .directive('npactGraphConfig', function(STATIC_BASE_URL){
    return {
      restrict: 'A',
      templateUrl:STATIC_BASE_URL+'js/graphs/config.html'

    };
  });
