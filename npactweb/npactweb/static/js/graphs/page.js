angular.module('npact')
  .directive('npactGraphPage', function(STATIC_BASE_URL, $http){
    return {
      restrict:'A',
      templateUrl:STATIC_BASE_URL+'js/graphs/page.html',
      link:function(scope, element, attrs){
	scope.title = 'Moorella thermoacetica Y72, whole genome shotgun sequence';
	scope.ranges = [
	  [0, 10000],
	  [10000, 20000],
	  [20000, 30000],
	  [30000, 40000],
	  [40000, 50000]
	];
	scope.d = {};

	$http.get(STATIC_BASE_URL+'js/nprofile.json').then(function(res){
	  scope.d.profile = res.data;
	});
	$http.get(STATIC_BASE_URL+'js/extract.json').then(function(res){
	  scope.d.cds = res.data;
	});
      }
    };
  });
