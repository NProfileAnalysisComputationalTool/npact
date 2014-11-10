angular.module('npact')
  .directive('npactGraph', function npactGraph($log, Grapher, Evt, GraphConfig, Utils, npactConstants, GraphingCalculator){
    return {
      restrict: 'A',
      scope:{
        spec: '=npactGraph'
      },
      controller: function($scope, $element, $attrs){
        var g = null, buildOptions = function(range) {
          var opts = angular.extend(
            {},
            npactConstants.graphSpecDefaults,
            range,
            {
              $scope: $scope,
              element: $element[0],
              headerSpec: GraphConfig.headerSpec(),
              margin: Utils.orderOfMagnitude(GraphConfig.basesPerGraph, -1),
              height: npactConstants.graphSpecDefaults.height,
              width: GraphConfig.width,
              colors: GraphConfig.colorBlindFriendly ?
                npactConstants.colorBlindLineColors : npactConstants.lineColors
            });

          opts.headerY = opts.headerSpec.headerY;
          opts.m = GraphingCalculator.chart(opts);
          opts.xaxis = GraphingCalculator.xaxis(opts);
          return opts;
        },
            redraw = function() {
              if(!$scope.spec) { return; }
              if(g !== null) {g.stage.destroy();}
              g = new Grapher(buildOptions($scope.spec));
              g.redraw().then($scope.$emit(Evt.GRAPH_REDRAW_COMPLETE));
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
