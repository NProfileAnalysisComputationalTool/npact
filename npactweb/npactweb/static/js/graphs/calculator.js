(function(){

  function GraphingCalculator(K, $log){

    var log10 = Math.log(10);
    
    function stops(a,b){
      var length = b - a,
	  // get even stops; look for one lower order of magnitude
	  // than our length - TODO: to many stops for lower ranges
	  interval = Math.pow(10, Math.round(Math.log(length) / log10) - 1),
	  // want to capture some margin
	  start = Math.max(a - interval,0),
	  end = b+interval;
      // want to be start/end INCLUSIVE so pad the end
      return {
	interval: interval,
	stops: _.range(start, end+interval, interval)
      };
    }
    
    function chart(opts){
      var yStops = [100, 80, 60, 40, 20, 0],
	  yAxisTicks = yStops.length - 1,
	  // the line graph, excluding tick marks and axis labels
	  g = {
	    x: opts.leftPadding,
	    y: opts.stageHeight - opts.profileHeight
	      - opts.axisLabelFontsize // x-axis labels
	      - 2*opts.profileTicks, // tick marks
	    h: opts.profileHeight,
	    w: opts.stageWidth - opts.leftPadding
	      - opts.rightPadding
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
	// TODO: return xaxis here
	yaxis:{
	  ticks: yticks,
	  labels: ylabels,
	  titleBox: ytitle
	}
      };
    }


    // return a top left x/y position that aligns
    // the centers of the two rectangles
    function alignRectangles(rect, toAlign){
      // get the top left coordinate of the `toAlign` rect
      return {x: rect.x + (rect.w/2) - (toAlign.w/2),
	      y: rect.y + (rect.h/2) - (toAlign.h/2)};
    }

    function xaxis(opts){
      var m = chart(opts),
	  s = stops(opts.range[0], opts.range[1]),
	  length = opts.range[1] - opts.range[0],
	  interval = s.interval, ss = s.stops,
	  // want to capture some margin
	  start = ss[0],
	  end = _.last(ss),
	  xAxisTicks = ss.length,
	  xAxisTickY = m.graph.y + m.graph.h,
	  ticks = ss.map(function(lbl, n){
	    var x = parseInt(start + n*interval);
	    return {
	      x: x, y: 0,
	      x2: x, y2: opts.profileTicks
	    };
	  }),
	  scaleX = m.graph.w / length,
	  // blow the text back up, we don't want it scaled
	  labelScaleX = 1/scaleX,
	  labels = ss.map(function(lbl, n){
	    return {
	      x: ticks[n].x, y: opts.profileTicks,
	      text: lbl,
	      scaleX: labelScaleX
	    };
	  })
      ;

      return {
	start: start,
	end: end,
	// total length of profile to draw
	length: length,
	interval: interval,
	x: m.graph.x,
	y: m.graph.y + m.graph.h,
	scaleX: scaleX,
	ticks: ticks,
	labels:labels
      };
    }
    return {
      chart:chart,
      xaxis:xaxis,
      alignRectangles:alignRectangles,
      stops:stops
    };
  }


  angular.module('npact')
    .factory('GraphingCalculator', GraphingCalculator)
  ;

}());

