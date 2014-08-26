angular.module('npact')
  .factory('Grapher', function(K, $log, GraphingCalculator, $rootScope, $compile, GraphDealer){

    function addMany(container, children){
      return container.add.apply(container, children);
    }

    function Grapher(opts){
      $log.log('new Grapher:', opts.range);
      this.opts = opts;
      this.zoomLevel = 1;
      this.$element = jQuery(opts.element);
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
      this.baseOffsetX = 0;
      this.stage.add(this.leftLayer(), this.chartLayer());
      this.stage.batchDraw();

//      this.onDblClick = _.bind(this.onDblClick, this);
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
	{width: title.getHeight(), height: title.getWidth()});

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

    function clamp(lower, value, upper){
      if(value <= lower) return lower;
      if(upper && value >= upper) return upper;
      return value;
    }
    
    GP.genomeGroup = function(){
      var self = this,
	  gx = this.m.graph.x,	  
	  dragRes = {x: gx, y:0},
	  g = new K.Group({
	    x: gx,
	    draggable: true,
	    dragBoundFunc:function(pos, evt){
	      // pos is the distance dragged, except `pos.x` always
	      // starts at `gx`. Something due to container coordinate
	      // system vs stage coordinate system. `pos.y` is not
	      // similarly affected.
	      var nextOffset = self.baseOffsetX - (pos.x - gx);
	      // change position via offsetX
	      this.offsetX(nextOffset);
	      // never change position via this return value
	      return dragRes;
	    }
	  });
      // remember where we are at the start/end of a drag, so dragging can
      // match our "scroll"
      g.on('dragstart dragend', function(obj){
	self.$element.trigger('tooltipHide.npact');
	self.baseOffsetX = this.offsetX();
	$log.log('baseOffsetX:', self.baseOffsetX, self.baseOffsetX/(self.zoomLevel * self.xaxis.scaleX) + self.opts.range[0]);
      });

      g.on('dragend', function(evt){
	// TODO: tell the graph dealer
	var oldStartBase = self.opts.range[0],
	    newStartBase = (this.offsetX() / self.xaxis.scaleX) + oldStartBase;
	GraphDealer.panTo(oldStartBase, newStartBase);
      });
      
      g.on('mouseover', function() {
        document.body.style.cursor = 'pointer';
      });
      g.on('mouseout', function() {
        document.body.style.cursor = 'default';
      });
      g.on('dblclick', _.bind(this.onDblClick, this));
      // need a shape that can be clicked on to allow dragging the
      // entire canvas
      // TODO: make this width always match the group width, in genespace
      g.add(new K.Rect({x: -100000, width:1000000, height:1000}));
      return this._genomeGroup = g;
    };

    GP.onDblClick = function(evt){
      // TODO: tell the GraphDealer about the zoom
      this.$element.trigger('tooltipHide.npact');
      var oldZoom = this.zoomLevel;
      $log.log('onDblClick', evt.evt.layerX - this.m.graph.x);
      
      if(evt.evt.shiftKey){
	// can't zoom out more
	if(oldZoom == 1) return;
	this.zoomLevel--;
      }else{
	this.zoomLevel++;
      }
      
      var zoom = GraphingCalculator.zoom({
	oldZoom: oldZoom,
	zoom:this.zoomLevel,
	offsetX: this.baseOffsetX,
	baseScaleX: this.xaxis.scaleX,
	start_gn: this.opts.range[0],
	length_gn: this.xaxis.length,
	// `layerX` is in stage coords, convert to graph coords
	centerOn_px: evt.evt.layerX - this.m.graph.x
      }),

	  xAxisGroup = this._xAxisGroup,
	  xAxisLabels = xAxisGroup.find('Text'),
	  cdsGroup = this._cdsGroup,
	  cdsLabels = cdsGroup.find('Text'),
	  scaleX = zoom.scaleX,
	  textScaleX = zoom.textScaleX
      ;


      // redraw axis with new scale factor
      xAxisGroup.scaleX(scaleX);
      // update labels
      xAxisLabels.map(function(txt){
	txt.scaleX(textScaleX);
	centerXLabel(txt, scaleX);
      });
      
      // rescale extracts
      cdsGroup.scaleX(scaleX);
      // update labels
      cdsLabels.map(function(txt){
	txt.scaleX(textScaleX);
	centerExtractLabel(txt, scaleX);
      });
      
      // rescale profiles
      this._profileGroup.scaleX(scaleX);
      // reset baseOffsetX, center on clicked spot
      this._genomeGroup.offsetX((this.baseOffsetX = zoom.offsetX));

      // TODO: resize the dragging Rect
      
      this.stage.batchDraw();
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

    function centerXLabel(txt, scaleX){
      // center on the tick marks, at this point we know how wide
      // this text is
      var w = txt.getWidth() / scaleX,
	  x = txt.getAttr('coord'),
	  newX = x - w/2;
      txt.x(newX);
    };
    
    GP.xAxisGroup = function(){
      var m = this.m,
	  xaxis = this.xaxis,
	  g = new K.Group({
	    x: 0, y: m.graph.h + m.graph.y,
	    width: xaxis.length,
	    scaleX: xaxis.scaleX,
	    offsetX: this.opts.range[0]
	  });

      addMany(g, this.drawAxisTicks(xaxis.ticks));
      var labels = this.drawAxisLabels(xaxis.labels);
      addMany(g, labels);

      // call `centerXLabel(txt, xaxis.scaleX)`
      labels.map(function(txt){
	centerXLabel(txt, xaxis.scaleX);
      });

      return this._xAxisGroup = g;
    };

    GP.profileGroup = function(){
      var m = this.m, opts = this.opts, xaxis = this.xaxis,
	  colors = this.colors,
	  lines = _.keys(colors),
	  g = new K.Group({
	    x: 0, y:m.graph.y,
	    height:opts.profileHeight, width:xaxis.length,
	    scaleX: xaxis.scaleX,
	    offsetX: this.opts.range[0]
	  }),
	  dataToDraw = _(this.opts.profile)
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

      return this._profileGroup = g;
    };
    
    GP.redraw = function(){
      this._genomeGroup.add(this.xAxisGroup(), this.profileGroup());
      this.stage.batchDraw();
    };

    function centerExtractLabel(txt, scaleX){
      // now that lbl is on the canvas, we can see what it's
      // height/width is
      var arrowBounds = txt.getAttr('arrowBounds'),
	  pos = GraphingCalculator.alignRectangles(
	    arrowBounds,
	    {
	      // convert to gene space
	      width:txt.getWidth()/ scaleX,
	      height:txt.getHeight()
	    });
      pos.x = Math.max(pos.x, arrowBounds.x+1);
      txt.position(pos);
    }
    
    GP.cdsGroup = function(cds){
      var xaxis = this.xaxis, opts = this.opts,
	  colors = this.colors,
	  g = new K.Group({
	    x: 0, y: 0,
	    scaleX: xaxis.scaleX,
	    offsetX: this.opts.range[0]
	  }),
	  colorNames = 'rgb',
	  // TODO: lookup where this should live in the header space
	  y = 30,
	  ahHalfHeight = opts.headerArrowHeight/2,
	  compY = y + opts.headerArrowHeight,
	  ahw = opts.headerArrowWidth/xaxis.scaleX,
	  textOpts = {
	    fontSize: opts.headerArrowFontsize,
	    fill:opts.axisFontcolor,
	    scaleX:1/xaxis.scaleX,
	    strokeScaleEnabled: false
	  },
	  arrowOpts = {
	    x: 0, y:0,
	    closed: true,
	    strokeWidth:1,
	    strokeScaleEnabled: false
	  },
	  
	  arrows = _(cds)
	    .filter(function(x){
	      return x.start >= xaxis.start
		&& x.end <= xaxis.end;
	    }).forEach(function(x){
	      // color determined by `x.end - 1 % 3`
	      var isComplement = x.complement == 0,
		  colorCoord = isComplement ? x.start : x.end,
		  c = colors[colorNames[(colorCoord - 1) % 3]],
		  baseY = isComplement ? y + opts.headerArrowHeight : y,
		  arrowPointY = baseY + ahHalfHeight,
		  arrowMaxY = baseY + opts.headerArrowHeight,
		  shape = isComplement
		    ? [
		      x.start, arrowPointY,
		      x.start + ahw, baseY,
		      x.end, baseY,
		      x.end, arrowMaxY,
		      x.start + ahw, arrowMaxY,
		      x.start, arrowPointY
		    ] : [
		      x.start, baseY,
		      x.end - ahw, baseY,
		      x.end, arrowPointY,
		      x.end - ahw, arrowMaxY,
		      x.start, arrowMaxY,
		      x.start, y
		    ],
		  arrowBounds = isComplement
		    ? {
		      x: x.start+ahw, y: baseY,
		      width: x.end-x.start-ahw, height: opts.headerArrowHeight
		    } : {
		      x: x.start, y: baseY,
		      width: x.end-x.start-ahw, height: opts.headerArrowHeight
		    },
		  line = new K.Line(angular.extend({
		    extract: x,
		    points: shape,
		    stroke:c
		  }, arrowOpts)),
		  // render the name, too
		  lbl = new K.Text(angular.extend({
		    extract: x,
		    arrowBounds: arrowBounds,
		    text: x.name
		  }, textOpts)),
		  // need a dummy group for clipping, `Text` doesn't
		  // support clip directly
		  lblGroup = new K.Group({clip:arrowBounds});
		  ;
	      
	      lblGroup.add(lbl);
	      g.add(line, lblGroup);
	      // now that lbl is on the canvas, we can see what it's
	      // height/width is
	      centerExtractLabel(lbl, xaxis.scaleX);
	    }),
	  $el = this.$element;

      // TODO: move tooltip control to the `npactGraph`, just throw an
      // event or callback here with the extract, let `npactGraph`
      // handle this stuff
      g.on('click', function(evt){
	var scope = $rootScope.$new(),
	    tpl = '<div npact-extract="extract"></div>';
	scope.extract = evt.target.getAttrs().extract;
	$el.qtip({	  
	  content: {text: $compile(tpl)(scope)},
	  position:{
	    my: scope.extract.complement == 0 ? 'top center' : 'bottom center',
	    target:[evt.evt.pageX, evt.evt.pageY]
	  },
	  show: {event:'tooltipShow.npact'},
	  hide: {event:'tooltipHide.npact'}
	});
	$el.trigger('tooltipShow.npact');
      });
      return this._cdsGroup = g;
    };
    
    GP.drawCDS = function(cds){
      $log.log(this.opts.range, 'drawCDS');
      this._genomeGroup.add(this.cdsGroup(cds));
    };
    return Grapher;
  });
