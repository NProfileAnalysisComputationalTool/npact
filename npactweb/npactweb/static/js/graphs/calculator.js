(function(){

  function GraphingCalculator(K, $log){
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
	  start = opts.range[0],
	  end = opts.range[1],
	  length = end - start,
	  // get even stops; look for one lower order of magnitude
	  // than our length 
	  interval = Math.pow(10, Math.floor(Math.log(length) / Math.log(10)) - 1),
	  // want to be start/end INCLUSIVE so pad the end
	  stops = _.range(start, end+interval, interval),
	  xAxisTicks = stops.length,
	  xAxisTickY = m.graph.y + m.graph.h,
	  ticks = stops.map(function(lbl, n){
	    var x = parseInt(start + n*interval);
	    return {
	      x: x, y: 0,
	      x2: x, y2: opts.profileTicks
	    };
	  }),
	  scaleX = m.graph.w / length,
	  // blow the text back up, we don't want it scaled
	  labelScaleX = 1/scaleX,
	  labels = stops.map(function(lbl, n){
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
      alignRectangles:alignRectangles
    };
  }


  angular.module('npact')
    .factory('GraphingCalculator', GraphingCalculator)
  ;

}());

