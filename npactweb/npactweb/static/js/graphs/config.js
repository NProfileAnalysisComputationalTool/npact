angular.module('npact')
  .directive('npactGraphConfig', function(STATIC_BASE_URL, GraphDealer, $log){
    return {
      restrict: 'A',
      templateUrl:STATIC_BASE_URL+'js/graphs/config.html',
      link:function($scope, $element, $attrs){
	$scope.basesPerGraph = GraphDealer.opts.basesPerGraph;
	$scope.changeZoom = _.debounce(GraphDealer.setZoom, 500);
	$scope.GraphDealer = GraphDealer;
      }
    };
  });
