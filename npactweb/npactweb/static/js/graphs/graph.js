angular.module('npact')
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

    function drawAxisTicks(layer, ticks, opts){
      var tickOpts = {x: 0, y:0, 
		      stroke: opts.borderColor,
		      strokeScaleEnabled: false};
      ticks
	.map(function(t){
	  tickOpts.points = [t.x, t.y, t.x2, t.y2];
	  layer.add(new K.Line(tickOpts));
	});
    }

    function drawAxisLabels(layer, labels, opts, textOpts){
      // draw labels at the right spacing
      var defaultTextOpts = {
	  fontSize: opts.axisLabelFontsize,
	  fill:opts.axisFontcolor
      };
      
      return labels.map(function(lbl){
	var txtOpts = angular.extend({}, lbl, defaultTextOpts, textOpts);
	var txt = new K.Text(txtOpts);
	layer.add(txt);	
	return txt;
      });
    }
    
    function drawYAxis(layer, opts){

      var m = GraphingCalculator.chart(opts);
      drawAxisTicks(layer, m.yaxis.ticks, opts);
      // draw labels at the right aligned
      drawAxisLabels(layer, m.yaxis.labels, opts, {align:'right'});

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
      var layer = new K.Group({
	x: m.graph.x,
	y: m.graph.y,
	draggable:true,
	dragBoundFunc:function(pos, evt){
	  $log.log('dragging', pos);
	  // TODO: stop if we reach gene-space 0
	  return {x: pos.x, y: this.y() };
	}
      });

      layer.on('mouseover', function() {
        document.body.style.cursor = 'pointer';
      });
      layer.on('mouseout', function() {
        document.body.style.cursor = 'default';
      });
      // frame around the graph
      var frame = new K.Rect({
	x: m.graph.x,
	y: m.graph.y,
	width: m.graph.w, 
	height: m.graph.h,
	stroke: opts.borderColor
      });
      var frameLayer = new K.Layer({
	clip: {x: m.graph.x-1,
	y: m.graph.y-1,
	width: m.graph.w+2, 
	height: 1000}
      });
      frameLayer.add(frame, layer);
      

      return { layer: layer, frame: frameLayer };
    }

    function drawProfile(layer, opts){
      var profile = opts.data.profile,
	  m = GraphingCalculator.chart(opts),
	  xaxis = GraphingCalculator.xaxis(opts),
	  lines = ['r', 'g', 'b'],
	  colors = {
	    r: opts.graphRedColor,
	    g: opts.graphGreenColor,
	    b: opts.graphBlueColor
	  }
      ;

      var frame = new K.Group({
	x: 0, y: 0,
	height: 100, width: xaxis.length,
	scaleX: xaxis.scaleX
      });
      layer.add(frame);

      var dataToDraw = _(profile)
	    .filter(function(g){
	      return g.coordinate >= xaxis.start
		&& g.coordinate <= xaxis.end;
	    })
	    .reduce(function(acc, g){
	      lines.map(function(x){
		acc[x].push(g.coordinate);
		acc[x].push(100 - g[x]);
	      });
	      return acc;
	    }, {r:[], g:[], b:[]});
      lines.map(function(x){
	var l = new K.Line({
	  points:dataToDraw[x],
	  stroke:colors[x],
	  strokeWidth:1,
	  strokeScaleEnabled: false
	});
	frame.add(l);
      });
    }
    function drawXAxis(layer, opts){
      var profile = opts.data.profile,
	  m = GraphingCalculator.chart(opts),
	  xaxis = GraphingCalculator.xaxis(opts)
      ;

      var frame = new K.Group({
	x: 0, y: m.graph.h,
	width: xaxis.length,
	scaleX: xaxis.scaleX
      });
      layer.add(frame);

      // draw ticks
      drawAxisTicks(frame, xaxis.ticks, opts);
      drawAxisLabels(frame, xaxis.labels, opts)
	.map(function(txt){
	  // center on the tick marks, at this point we know how wide
	  // this text is
	  var w = txt.getWidth() / xaxis.scaleX,
	      x = txt.x();
	  txt.x(x - w/2);
	});

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
      stage.add(chartChrome.frame);
      // need a shape that can be clicked on to allow dragging the entire canvas
      var r = new K.Rect({
	x:0, y:0,
	width:1000, height:1000	
      });
      chartChrome.layer.add(r);
      
      stage.batchDraw();

      return {
	redraw:function(newOpts){
	  drawProfile(chartChrome.layer, opts);
	  drawXAxis(chartChrome.layer, opts);
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
      leftPadding: 80, rightPadding:20,
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
