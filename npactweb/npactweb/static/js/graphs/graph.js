(function(){
  function npactGraph($log, Grapher){
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

	function onSpec(spec, oldval){
	  if(!spec) return;
	  // TODO: clean up a pre-existing graph?
	  var g = new Grapher(angular.extend({}, spec, opts));
	  // TODO: make redraw handle extracts
	  g.redraw();
	  _.forOwn(spec.extracts, function(value, key){
	    g.drawCDS(value);
	  });
	}
      }
    };
  }

  // register everything with angular
  angular.module('npact')
    .directive('npactGraph', npactGraph)
    .directive('npactExtract', function(STATIC_BASE_URL){
      return {
	restrict: 'A',
	scope: { extract:'=npactExtract'},
	templateUrl:STATIC_BASE_URL+'js/graphs/extract.html'
      };
    })
  ;
}());
