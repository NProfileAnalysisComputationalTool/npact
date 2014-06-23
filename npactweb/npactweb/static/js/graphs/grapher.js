angular.module('npact')
  .factory('Grapher', function(K, $log, GraphingCalculator){

    function addMany(container, children){
      return container.add.apply(container, children);
    }
    
    function Grapher(opts){
      $log.log('new Grapher:', opts);
      this.opts = opts;
      this.stage = new K.Stage({
	container:opts.element,
	height: opts.height,
	width: opts.width
      });
      // TODO: height should be dynamic
      opts.stageHeight = this.stage.getHeight();
      opts.stageWidth = this.stage.getWidth();

      this.colors = opts.colorblindFriendlyColors ? {
	r: opts.graphRedColorblind,
	g: opts.graphGreenColorblind,
	b: opts.graphBlueColorblind
      } : {
	r: opts.graphRedColor,
	g: opts.graphGreenColor,
	b: opts.graphBlueColor
      };
      
      this.m = GraphingCalculator.chart(opts);
      this.xaxis = GraphingCalculator.xaxis(opts);
      $log.log('chart measurements:', this.m, this.xaxis);
      this.stage.add(this.leftLayer(), this.chartLayer());
      this.stage.batchDraw();

      
    }
    var GP = Grapher.prototype;
    
    GP.drawAxisTicks = function(ticks){
      var tickOpts = {x: 0, y:0, 
		      stroke: this.opts.borderColor,
		      strokeScaleEnabled: false};
      return ticks
	.map(function(t){
	  tickOpts.points = [t.x, t.y, t.x2, t.y2];
	  return new K.Line(tickOpts);
	});
    };

    GP.drawAxisLabels = function(labels, textOpts){
      // draw labels at the right spacing
      var defaultTextOpts = {
	fontSize: this.opts.axisLabelFontsize,
	fill:this.opts.axisFontcolor
      };

      return labels.map(function(lbl){
	var txtOpts = angular.extend({}, lbl, defaultTextOpts, textOpts);
	return new K.Text(txtOpts);
      });
    };

    GP.headerGroup = function(){
	// from top to bottom
      var opts = this.opts,
	  g = new K.Group(),
	  labels = opts.headers.map(function(hdr){
	    return {
	      text:hdr.title,
	      height:opts.headerSizes[hdr.lineType]};
	  }),
	// TODO: derive opts.leftPadding
	// * add all labels
	// * loop over K.Text objects, find max width
	// * leftPadding = maxWidth + headerLabelPadding
	// * loop over K.Text objects, set width to leftPadding

	defaultTextOpts = {
	  align: 'right', x: 0, 
	  fontSize: opts.headerLabelFontsize,
	  width: opts.leftPadding - opts.headerLabelPadding,
	  fill:opts.headerLabelFontcolor
	},
	  y = 0,
	  lbls =  labels.map(function(lbl){
	    var txtOpts = angular.extend({
	      y: y
	    }, defaultTextOpts, lbl);
	    var txt = new K.Text(txtOpts);
	    y += parseInt(lbl.height);
	    return txt;
	  });

      addMany(g, lbls);
      return g;
    };

    GP.yAxisGroup = function(){
      var g = new K.Group(),
	  opts = this.opts,
	  m = this.m;

      addMany(g, this.drawAxisTicks(m.yaxis.ticks));
      	    // draw labels at the right aligned
      addMany(g, this.drawAxisLabels(m.yaxis.labels, {align:'right'}));

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
      
      g.add(title);
      return g;
    };

    GP.leftLayer = function(){
      
      var opts = this.opts,
	  layer = new K.Layer({
	    width: opts.leftPadding,
	    x:0, y:0,
	    height: opts.stageHeight
	  });

      layer.add(this.headerGroup(), this.yAxisGroup());
      
      return layer;
    };

    GP.genomeGroup = function(){
      var g = new K.Group({
	id:'genome',
	x: this.m.graph.x,
	draggable: true,
	dragBoundFunc:function(pos, evt){
	  // TODO: stop if we reach limits of gene-space	  
	  return {x: pos.x, y: this.y() };
	}
      });
      g.on('mouseover', function() {
        document.body.style.cursor = 'pointer';
      });
      g.on('mouseout', function() {
        document.body.style.cursor = 'default';
      });
      // need a shape that can be clicked on to allow dragging the
      // entire canvas
      g.add(new K.Rect({width:1000, height:1000}));
      return g;
    };

    GP.chartLayer = function(){
      var m = this.m,
	  xaxis = this.xaxis,
	  opts = this.opts,
	  l = new K.Layer({
	    clip:{
	      x: m.graph.x-1, y:0,
	      width:m.graph.w+2,
	      height:1000
	    }
	  }),
	  // frame around the graph
	  border = new K.Rect({
	    x: m.graph.x,
	    y: m.graph.y,
	    width: m.graph.w, 
	    height: m.graph.h,
	    stroke: opts.borderColor
	  });
      
      l.add(border, this.genomeGroup());
      return l;
    };

    GP.xAxisGroup = function(){
      var profile = this.opts.data.profile,
	  m = this.m,
	  xaxis = this.xaxis,
	  g = new K.Group({
	    x: 0, y: m.graph.h + m.graph.y,
	    width: xaxis.length,
	    scaleX: xaxis.scaleX
	  });

      addMany(g, this.drawAxisTicks(xaxis.ticks));
      var labels = this.drawAxisLabels(xaxis.labels);
      addMany(g, labels);

      labels .map(function(txt){
	  // center on the tick marks, at this point we know how wide
	  // this text is
	  var w = txt.getWidth() / xaxis.scaleX,
	      x = txt.x();
	  txt.x(x - w/2);
	});

      return g;
    };

    GP.profileGroup = function(){
      var m = this.m, opts = this.opts, xaxis = this.xaxis,
	  colors = this.colors,
	  lines = _.keys(colors),
	  g = new K.Group({
	    x: 0, y:m.graph.y,
	    height:opts.profileHeight, width:xaxis.length,
	    scaleX: xaxis.scaleX
	  }),
	  dataToDraw = _(this.opts.data.profile)
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
	    }, {r:[], g:[], b:[]}),
	  profiles = lines.map(function(x){
	    return new K.Line({
	      points:dataToDraw[x],
	      stroke:colors[x],
	      strokeWidth:1,
	      strokeScaleEnabled: false
	    });})
      ;
      addMany(g, profiles);

      return g;
    };
    
    GP.redraw = function(){
      $log.log('redraw');
      var genomeGroup = this.stage.find('#genome');
      genomeGroup.add(this.xAxisGroup(), this.profileGroup());
      this.stage.batchDraw();
    };
    return Grapher;
  });
