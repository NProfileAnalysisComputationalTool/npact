angular.module('npact')
  .factory('GraphingCalculator', function(K, $log){
    return {
      chart:function(opts){

	var yStops = [100, 80, 60, 40, 20, 0],
	    yAxisTicks = yStops.length - 1,
	    // the line graph, excluding tick marks and axis labels
	    g = {
	      x: opts.leftPadding,
	      y: opts.stageHeight - opts.profileHeight
		- opts.axisLabelFontsize // x-axis labels
		- opts.profileTicks, // tick marks
	      h: opts.profileHeight,
	      w: opts.stageWidth - opts.leftPadding
		- opts.profileTicks	  
	    },
	    yAxisTickSpacing = g.h / yAxisTicks,
	    yAxisTickX = g.x - opts.profileTicks,
	    yticks = yStops.map(function(lbl, n){
		var y = g.y + n*yAxisTickSpacing;
		return {
		  x: yAxisTickX, y: parseInt(y),
		  x2: g.x, y2: parseInt(y)
		};
	    }),
	    // space between the label and the tick
	    yAxisLabelRightPadding = opts.profileTicks*2,
	    yAxisLabelWidth = opts.leftPadding - yAxisLabelRightPadding,
	    ylabels = yStops.map(function(lbl, n){
		// center on the tick
		var y = yticks[n].y - (opts.axisLabelFontsize/2);
		return {
		  x: 0, y: parseInt(y),
		  width: yAxisLabelWidth,
		  text: lbl
		};
	    }),
	    // box to center the title inside of
	    ytitle={
	      x: 0, y: g.y, w: yAxisLabelWidth, h: g.h
	    }
	    ;

	return {
	  graph:g,
	  yaxis:{
	    ticks: yticks,
	    labels: ylabels,
	    titleBox: ytitle
	    }
	};
	
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
	  fill:opts.headerLabelFontcolor
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

      var tickOpts = {x: 0, y:0, stroke: opts.borderColor};
      m.yaxis.ticks
	.map(function(t){
	  tickOpts.points = [t.x, t.y, t.x2, t.y2];
	  layer.add(new K.Line(tickOpts));
	});

      // draw labels at the right spacing
      var defaultTextOpts = {
	  align: 'right',
	  fontSize: opts.axisLabelFontsize,
	  fill:opts.axisFontcolor
      };
      
      m.yaxis.labels.map(function(lbl){
	var txtOpts = angular.extend({}, lbl, defaultTextOpts);
	var txt = new K.Text(txtOpts);
	layer.add(txt);	
      });

      // the title
      var title = new K.Text({
	y:0, x:0, // reposition this below
	rotation:-90,
	fill:opts.axisFontcolor,
	fontSize: opts.axisTitleFontsize,
	text: opts.axisTitle	
      });
      // center it in the space left of the axes
      var pos = GraphingCalculator.alignRectangles(
	// define the space we want to center inside
	m.yaxis.titleBox,
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

      // graph layer, translates origin to graph
      var layer = new K.Layer({
	width: m.graph.w,
	height: m.graph.h,
	x: m.graph.x,
	y: m.graph.y
      });

      // frame around the graph
      var frame = new K.Rect({
	x: 0, y:0,
	width: m.graph.w, 
	height: m.graph.h,
	stroke: opts.borderColor
      });
      layer.add(frame);

      return { layer: layer };
    }

    function drawProfile(layer, opts){
      var profile = opts.data.profile,
	  start = opts.range[0],
	  end = opts.range[1],
	  // total length of profile to draw
	  length = end - start,
	  m = GraphingCalculator.chart(opts),
	  lines = ['r', 'g', 'b'],
	  colors = {
	    r: opts.graphRedColor,
	    g: opts.graphGreenColor,
	    b: opts.graphBlueColor
	  }
      ;
      $log.log(start, end, length, m.graph.w / length);
      var frame = new K.Group({
	x:0, y:0,
	height:100, width:length,
	scaleX: m.graph.w / length
      });
      layer.add(frame);

      var dataToDraw = _(profile)
	    .filter(function(g){
	      return g.coordinate >= start
		&& g.coordinate <= end;
	    })
	    .reduce(function(acc, g){
	      lines.map(function(x){
		acc[x].push(g.coordinate);
		acc[x].push(100 - g[x]);
	      });
	      return acc;
	    }, {r:[], g:[], b:[]});
      $log.log(dataToDraw);
      lines.map(function(x){
	var l = new K.Line({
	  points:dataToDraw[x],
	  stroke:colors[x],
	  strokeWidth:1,
	  strokeScaleEnabled: false
	});
	frame.add(l);
      });
      $log.log(frame);
    }
    function drawXAxis(layer, opts){
      
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
	redraw:function(newOpts){
	  drawProfile(chartChrome.layer, opts);
	  // drawXAxis(?, opts);
	  stage.batchDraw();
	}
      };
    }
    
    return function(opts){
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

	// width takes a minute to get sorted, we might be waiting on jquery
	var p = WidthChecker(element)
	  .then(function(w){
	    var opts = angular.extend({}, defaultGraphOpts, attrs, {
	      element: element[0],
	      width: w,
	      data: scope.data,
	      range: scope.range,
	      height: element.height()
	    });
	    // convert `string`s to `int`s
	    angular.forEach(defaultGraphOpts, function(val, key){
	      opts[key] = parseInt(opts[key]);
	    });
	    return Grapher(opts);
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
