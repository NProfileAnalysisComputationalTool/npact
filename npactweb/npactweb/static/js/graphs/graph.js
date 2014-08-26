(function(){
  function NpactGraph($log, Grapher){
    return {
      restrict: 'A',
      scope:{
	spec: '=npactGraph'
      },
      link: function($scope, $element, $attrs){
	var opts = {
	  element: $element[0]
	};
	$scope.$watch('spec', onSpec);
	
	function onSpec(spec, oldval){
	  if(!spec) return;
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
    .directive('npactGraph', NpactGraph)
    .directive('npactExtract', function(STATIC_BASE_URL){
      return {
	restrict: 'A',
	scope: { extract:'=npactExtract'},
	templateUrl:STATIC_BASE_URL+'js/graphs/extract.html'
      };
    })
  ;
}());
