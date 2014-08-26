angular.module('npact')
  .directive('npactGraphPage', function(STATIC_BASE_URL, $http, GraphDealer){
    return {
      restrict:'A',
      templateUrl:STATIC_BASE_URL+'js/graphs/page.html',
      link:function(scope, element, attrs){
	scope.title = 'Moorella thermoacetica Y72, whole genome shotgun sequence';

	$http.get(STATIC_BASE_URL+'js/nprofile.json').then(function(res){
	  GraphDealer.setProfile(res.data);
	});
	$http.get(STATIC_BASE_URL+'js/extract.json').then(function(res){
	  GraphDealer.addExtract('Input file CDS', res.data);
	});
      }
    };
  });
