angular.module('npact')
  .factory('GraphingCalculator', function(K, $log){
    return {
      chart:function(opts){
	return {
	  // the top of the graph
	  // need to leave some room at the bottom
	  y: opts.stageHeight - opts.profileHeight
	    - opts.axisLabelFontsize // a-axis labels
	    - opts.profileTicks, // tick marks
	  // how frequently to draw an Y-axis label
	  tickY: opts.profileHeight / 5,
	  // account for font size
	  labelY: (opts.profileHeight / 5)
	};
      },
      // return a top left x/y position that aligns
      // the centers of the two rectangles
      alignRectangles: function(rect, toAlign){
	// get the top left coordinate of the `toAlign` rect
	return {x: rect.x + (rect.w/2) - (toAlign.w/2),
		y: rect.y + (rect.h/2) - (toAlign.h/2)};
      }
    };
  })

  .factory('Grapher', function($log, GraphingCalculator, K){

    function drawLabels(layer, opts){
      // from top to bottom
      var labels = [
	{text: 'Newly identified ORFs', height: opts.orfHeight},
	{text: 'Input file CDS', height: opts.cdsHeight},
	{text: 'Hits', height: opts.hitsHeight}
      ];

      var defaultTextOpts = {
	  align: 'right', x: 0, 
	  fontSize: opts.headerLabelFontsize,
	  width: opts.leftPadding - opts.headerLabelPadding,
	  fill:'black'
      };
      var y = 0;
      labels.map(function(lbl){
	var txtOpts = angular.extend({
	  y: y
	}, defaultTextOpts, lbl);
	var txt = new K.Text(txtOpts);
	y += parseInt(lbl.height);
	layer.add(txt);
      });
    }

    function drawYAxis(layer, opts){

      var m = GraphingCalculator.chart(opts);

      // draw labels at the right spacing
      var yLabels = [100, 80, 60, 40, 20, 0],
	  y = m.y;

      var defaultTextOpts = {
	  align: 'right', x: 0, 
	  fontSize: opts.axisLabelFontsize,
	  width: opts.leftPadding,
	  fill:'black'
      };

      yLabels.map(function(lbl){
	var txtOpts = angular.extend({
	  y: y, text:lbl
	}, defaultTextOpts);
	var txt = new K.Text(txtOpts);
	y += m.labelY;
	layer.add(txt);
      });

      // the title
      var title = new K.Text({
	y:0, x:0, // reposition this below
	rotation:-90,
	fill:'black',
	fontSize: opts.axisTitleFontsize,
	text: opts.axisTitle	
      });
      // center it in the space left of the axes
      var pos = GraphingCalculator.alignRectangles(
	// define the space we want to center inside
	{x:0, y:m.y,
	 w: opts.leftPadding - opts.headerLabelPadding,
	 h: opts.profileHeight
	},
	// bounding box of the text, pre-rotated so need to swap W and H
	{w: title.getHeight(), h: title.getWidth()});

      // translate to the bottom-left of the new rectangle, so after we rotate
      // we'll be centered
      pos.y += title.getWidth();
      title.position(pos);
      
      layer.add(title);      
      
    }
    
    function createLeftArea(opts){
      var layer = new K.Layer({
	width: opts.leftPadding,
	x:0, y:0,
	height: opts.stageHeight
      });

      // TODO: return something more useful here and change the
      // rendering of the left labels as we get data
      drawLabels(layer, opts);

      drawYAxis(layer, opts);

      return { layer: layer };
    }
    
    function createChartChrome(opts){

      var m = GraphingCalculator.chart(opts);

      var layer = new K.Layer({
	width: opts.stageWidth - opts.leftPadding,
	height: opts.stageHeight - m.y,
	x: opts.leftPadding,
	y: m.y
      });

      var frame = new K.Rect({
	x: opts.profileTicks, y:0,
	width: layer.getWidth() - 2*opts.profileTicks, 
	height: opts.profileHeight,
	stroke: 'black'
      });
      layer.add(frame);

      return { layer: layer };
    }

    function make(stage, opts){

      // how much gene-space are we going to show
      var gLength = opts.range[1] - opts.range[0];
      opts.stageHeight = stage.getHeight();
      opts.stageWidth = stage.getWidth();
      // setup Y axis and static left area
      
      var leftArea = createLeftArea(opts);
      stage.add(leftArea.layer);

      // setup chart lines and background
      var chartChrome = createChartChrome(opts);
      stage.add(chartChrome.layer);
      
      stage.batchDraw();

      return {
	redraw:function(data){
	  $log.log('redraw', arguments);
	  stage.batchDraw();
	}
      };
    }
    
    return function(opts){
      $log.log('new graph', opts);
      var stage = new K.Stage({
	container:opts.element,
	height: opts.height,
	width: opts.width
      });
      return make(stage, opts);
    };
  })
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
  .directive('npactGraph', function($log, Grapher, WidthChecker){

    var defaultGraphOpts = {
      leftPadding: 80,
      profileHeight: 100,
      hitsHeight:15, cdsHeight:15, orfHeight:30,
      headerLabelPadding: 10, headerLabelFontsize:10,
      axisLabelFontsize: 10, axisTitleFontsize:14,
      width:600, height:200, profileTicks:5
    };
    
    return {
      restrict: 'A',
      scope:{
	data:'=',
	range:'='
      },
      link:function(scope, element, attrs){
	// our graph object
	var g;
	// width takes a minute to get sorted, we might be waiting on jquery
	WidthChecker(element)
	  .then(function(w){
	    var opts = angular.extend({}, defaultGraphOpts, attrs, {
	      element: element[0],
	      width: w,
	      height: element.height(),
	      data: scope.data,
	      range: scope.range
	    });
	    // convert `string`s to `int`s
	    angular.forEach(defaultGraphOpts, function(val, key){
	      opts[key] = parseInt(opts[key]);
	    });
	    g = Grapher(opts);
	  });
	
	scope.$watch('data', function(newval, oldval){
	  if(!newval) return;
	  $log.log('data changed:', newval);
	  if(g) g.redraw(newval);
	});
      }
    };
  })
;
