angular.module('npact')
  .directive('npactGraphPage', function(STATIC_BASE_URL, $http, GraphDealer, Utils){
    return {
      restrict:'A',
      templateUrl:STATIC_BASE_URL+'js/graphs/page.html',
      link:function($scope, $element, $attrs){
	$scope.title = 'Moorella thermoacetica Y72, whole genome shotgun sequence';

	// TODO: watch for changes in width
	// http://stackoverflow.com/questions/23044338/window-resize-directive
	Utils.widthAvailable($element).then(GraphDealer.setWidth);

	$http.get(STATIC_BASE_URL+'js/nprofile.json').then(function(res){
	  GraphDealer.setProfile(res.data);
	});
	$http.get(STATIC_BASE_URL+'js/extract.json').then(function(res){
	  GraphDealer.addExtract('Input file CDS', res.data);
	});
      }
    };
  });
