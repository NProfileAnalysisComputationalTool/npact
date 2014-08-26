(function(){
  function WidthChecker($q, $log, $interval){
    /**
     * polls for when the element has a width
     *
     * @returns {Promise} the width of the element
     */
    function pollForWidth(element, frequency){
      var d = $q.defer(),
	  p = d.promise,
	  t1 = new Date(),
	  task = $interval(checkForWidth, frequency || 100);
      // stop polling once we've found a value
      p.then(function(){$interval.cancel(task);});
      return p;
      
      function checkForWidth(){
	var w = element.width();
	// indicate progress via the promise
	d.notify(new Date() - t1);
	if(w > 0){
	  d.resolve(w);
	}
      }
    };

    return pollForWidth;
  }

  function NpactGraph(WidthChecker, $log, Grapher){
    return {
      restrict: 'A',
      scope:{
	spec: '=npactGraph'
      },
      link: function($scope, $element, $attrs){
	var onOpts = WidthChecker($element).then(graphOpts);
	$scope.$watch('spec', onSpec);

	function graphOpts(width){
	  return {
	    element: $element[0],
	    width: width,
	    height: $element.height()
	  };
	}

	// TODO: watch for changes in width
	// http://stackoverflow.com/questions/23044338/window-resize-directive
	
	function onSpec(newval, oldval){
	  if(!newval) return;
	  onOpts.then(_.partial(redraw, newval));
	}

	function redraw(spec, opts){
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
    .factory('WidthChecker', WidthChecker)
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
