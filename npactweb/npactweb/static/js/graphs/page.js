angular.module('npact')
  .directive('npactGraphPage', function(STATIC_BASE_URL){
    return {
      restrict:'A',
      scope:{
	
      },
      templateUrl:STATIC_BASE_URL+'js/graphs/page.html'
    };
  });
