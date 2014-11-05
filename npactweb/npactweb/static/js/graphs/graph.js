angular.module('npact')
  .directive('npactGraph', function npactGraph($log, Grapher, Evt){
    return {
      restrict: 'A',
      scope:{
        spec: '=npactGraph'
      },
      controller: function($scope, $element, $attrs){
        var opts = { element: $element[0] },
            redraw = function() {
              if(!$scope.spec) { return; }
              // TODO: clean up a pre-existing graph?
              var g = new Grapher(angular.extend({}, $scope.spec, opts), $scope);
              g.redraw();
            };

        $scope.$on(Evt.REDRAW, redraw);
        redraw();
      }
    };
  })
  .directive('npactExtract', function(STATIC_BASE_URL){
    return {
      restrict: 'A',
      scope: { extract:'=npactExtract'},
      templateUrl:STATIC_BASE_URL+'js/graphs/extract.html'
    };
  })
;
