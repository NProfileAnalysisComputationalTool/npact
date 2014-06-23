angular.module('npact')
  .factory('WidthChecker', function($q, $log, $interval){
    return function(element, frequency){
      var hasWidth = $q.defer(),
	  d1 = new Date(),
	  widthChecker;
      
      function checkForWidth(){
	var w = element.width();
	$log.log('T[', new Date() - d1, '] width:', w);
	if(w > 0){
	  $interval.cancel(widthChecker);
	  hasWidth.resolve(w);
	}
      }      
      widthChecker = $interval(checkForWidth, frequency || 100);
      return hasWidth.promise;
    };
  })
  .constant('GraphSettings', {
    // TODO: determine me dynamically
    leftPadding:120,

    rightPadding:25,
    headerSizes:{'extract':30, 'hits':15},
    profileHeight:100,
    profileTicks:5,
    headerLabelPadding:10,
    
    headerLabelFontcolor:"#444",
    headerLabelFontsize:11,
    axisLabelFontsize:11,
    axisFontcolor:"#444",
    axisTitleFontsize:20,
    borderColor:"#444",
    // profile line colors
    graphRedColor:"red",
    graphBlueColor:"blue",
    graphGreenColor:"green",
    graphRedColorblind:"rgb(213, 94, 0)",
    graphBlueColorblind:"rgb(204, 121, 167)",
    graphGreenColorblind:"rgb(0, 114, 178)"
  })
  .directive('npactGraph', function($log, Grapher, WidthChecker, GraphSettings){

    
    return {
      restrict: 'A',
      scope:{
	data:'=',
	range:'=',
	headers:'=',
	colorblindFriendlyColors:'=',
	axisTitle:'='
      },
      link:function(scope, element, attrs){

	// width takes a minute to get sorted, we might be waiting on jquery
	var p = WidthChecker(element)
	  .then(function(w){
	    var opts = angular.extend({
	      element: element[0],
	      width: w,
	      height: element.height()
	    }, GraphSettings, scope);
	    return new Grapher(opts);
	  });
	
	scope.$watch('data.profile', function(newval, oldval){
	  if(!newval) return;

	  // when we have a graph ready, redraw it.
	  p.then(function(g){
	    g.redraw();
	  });
	});
      }
    };
  })
;
