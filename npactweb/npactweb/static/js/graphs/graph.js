angular.module('npact')
  .directive('npactGraph', function npactGraph($log, Grapher, Evt){
    return {
      restrict: 'A',
      scope:{
        spec: '=npactGraph'
      },
      controller: function($scope, $element, $attrs){
        var opts = {
          element: $element[0]
        };
        $scope.$watch('spec', onSpec);

        $scope.$on(Evt.REDRAW, function() {
          onSpec($scope.spec);
        });

        function onSpec(spec, oldval){
          if(!spec) return;
          // TODO: clean up a pre-existing graph?
          var g = new Grapher(angular.extend({}, spec, opts));
          g.redraw();
        }
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
