angular.module('npact')
  .directive('npactGraphConfig', function npactGraphConfig(STATIC_BASE_URL, GraphDealer, $log) {
    return {
      restrict: 'A',
      templateUrl:STATIC_BASE_URL+'js/graphs/config.html',
      controller: function(){
        this.nextPage = _.debounce(GraphDealer.nextPage, 250);
        this.previousPage = _.debounce(GraphDealer.previousPage, 250);
        this.setZoom = _.debounce(GraphDealer.setZoom, 500);
        this.setColors = _.debounce(GraphDealer.setColors, 250);
      },
      controllerAs: 'ctrl',
      link:function($scope, $element, $attrs){
        $scope.graphConfig = GraphDealer.opts;
      }
    };
  });
