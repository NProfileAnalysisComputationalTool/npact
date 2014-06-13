angular.module('npact')
  .directive('npactGraphPage', function(STATIC_BASE_URL, $http){
    return {
      restrict:'A',
      templateUrl:STATIC_BASE_URL+'js/graphs/page.html',
      link:function(scope, element, attrs){
	scope.title = 'Moorella thermoacetica Y72, whole genome shotgun sequence';
	scope.range = [0, 10000];
	scope.d = {};

	$http.get(STATIC_BASE_URL+'js/nprofile.json').then(function(res){
	  scope.d.profile = res.data;
	});
      }
    };
  });
