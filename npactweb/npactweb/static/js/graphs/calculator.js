(function(){

  /**
   * calculate measurements about the chart
   */
  function chartMetrics(opts){
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
	yticks = yStops.map(makeYTick),
	// space between the label and the tick
	yAxisLabelRightPadding = opts.profileTicks*2,
	yAxisLabelWidth = opts.leftPadding - yAxisLabelRightPadding,
	ylabels = yStops.map(makeYLabel),
	// box to center the title inside of
	ytitle={
	  x: 0, y: g.y, width: yAxisLabelWidth, height: g.h
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

    function makeYTick(lbl, n){
      var y = g.y + n*yAxisTickSpacing;
      return {
	x: yAxisTickX, y: parseInt(y),
	x2: g.x, y2: parseInt(y)
      };
    }
    
    function makeYLabel(lbl, n){
      // center on the tick
      var y = yticks[n].y - (opts.axisLabelFontsize/2);
      return {
	x: 0, y: parseInt(y),
	width: yAxisLabelWidth,
	text: lbl
      }; 
    }
  }

  /**
   * find a position that will align the two rectangles on their center
   *
   * @param {Object} rect - the rectangle to align to {width,height}
   * @param {Object} toAlign - the rectangle you want to move {width,height}
   * @returns {Object} x,y coordinates to shift `toAlign` by, from the top/left
   */
  function alignRectangles(rect, toAlign){
    // get the top left coordinate of the `toAlign` rect
    return {x: rect.x + (rect.width/2) - (toAlign.width/2),
	    y: rect.y + (rect.height/2) - (toAlign.height/2)};
  }
  
  function GraphingCalculator(K, Utils){
   
    function stops(a,b){
      var length = b - a,
	  // get even stops; look for one lower order of magnitude
	  // than our length - TODO: too many stops for lower ranges
	  interval = Utils.orderOfMagnitude(length, -1),
	  // want to capture some margin
	  start = Math.max(a - interval,0),
	  end = b+interval;
      // want to be start/end INCLUSIVE so pad the end
      return {
	interval: interval,
	stops: _.range(start, end+interval, interval)
      };
    }

    function xaxis(opts){
      var m = chartMetrics(opts),
	  s = stops(opts.range[0], opts.range[1]),
	  length = opts.range[1] - opts.range[0],
	  interval = s.interval,
	  ss = s.stops,
	  // want to capture some margin
	  start = ss[0],
	  end = _.last(ss),
	  xAxisTicks = ss.length,
	  ticks = ss.map(makeTicks),
	  scaleX = m.graph.w / length,
	  // blow the text back up, we don't want it scaled
	  labelScaleX = 1/scaleX,
	  labels = ss.map(makeLabels)
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

      function makeTicks(lbl, n){
	var x = parseInt(start + n*interval);
	return {
	  x: x, y: 0,
	  x2: x, y2: opts.profileTicks
	};
      }

      function makeLabels(lbl, n){
	return {
	  x: ticks[n].x, y: opts.profileTicks,
	  coord: ticks[n].x,
	  text: lbl,
	  scaleX: labelScaleX
	};
      }
    }

    function zoom(obj){
      // zoom while retaining our position, so if you zoom in on gene
      // 2000, after the zoom the mouse cursor is still on 2000, but
      // the rest of the graph has shifted. Treats the genes as a
      // number line.

      // s = initial start of the gene space we're showing
      // L = the low end of the visible graph
      // H = the high end of the visible graph
      // Z = where the user is zooming into
      // o = offset (not shown in the diagrams), how much the user has
      //     scrolled. Represented as the space between `s` and `L`

      // -------------------------------------
      //               BEFORE
      // -------------------------------------
      //       s        Z                H # no pan (zero offset)
      //    s  L        Z                H # panned left (positive offset)
      //       L   s    Z                H # panned right (negative offset)
      // -------------------------------------

      // after the zoom, we want Z still in the same spot on the
      // screen, with everything else shifting. Offset increases

      // -------------------------------------
      //               AFTER
      // -------------------------------------
      //    s  L        Z                H # no pan
      //  s    L        Z                H # panned left
      //       L s      Z                H # panned right
      // -------------------------------------

      // some algebra (LZ means length between L and Z):
      // o = L - s
      // LZ = LH * LZ_pct
      // ZH = LH * (1-LZ_pct)
      // L = LH - (1-LZ_pct)*Z
      // Z = LG * LZ_pct

      // we deal with 3 coordinate systems:
      //  * gene space (gn suffix)
      //  * before-zoom graph space (bpx suffix)
      //  * after-zoom graph space (apx suffix)

      // some aliases
      var zoomLevel = obj.zoom,
	  baseScaleX = obj.baseScaleX,
	  s_gn = obj.start_gn,
	  // -------------------------------------
	  // calculate some pre-zoom values
	  oldScaleX = baseScaleX*obj.oldZoom,
	  oldLength_gn = obj.length_gn / obj.oldZoom,
	  oldLength_bpx = oldLength_gn * oldScaleX,

      	  // new scale factor
	  scaleX = baseScaleX * zoomLevel,
	  // where theaccount for offset
	  Z_bpx = obj.centerOn_px + obj.offsetX,
	  // convert graph space to gene space
	  Z_gn = (Z_bpx / oldScaleX) + s_gn,
	  // length of new visible area
	  LH_gn = obj.length_gn / zoomLevel,

	  // -------------------------------------
	  // find the pre-zoom length of LZ, as a percentage of LH,
	  // using graph coordinates
	  LZ_bpx = obj.centerOn_px,
	  LH_bpx = oldLength_bpx,
	  LZ_pct = LZ_bpx / LH_bpx,

	  // -------------------------------------
	  // based on Z_gn, LH_gn and LZ_pct, find the new offset
	  LZ_gn = LH_gn * LZ_pct,
	  L_gn = Z_gn - LZ_gn,
	  offset_gn = L_gn - s_gn,
	  // convert gene space to graph space
	  offset_apx = offset_gn * scaleX
	  ;

      return {
	scaleX : scaleX,
	offsetX: offset_apx,
	textScaleX: 1/scaleX
      };
    };

    return {
      chart:chartMetrics,
      xaxis:xaxis,
      alignRectangles:alignRectangles,
      stops:stops,
      zoom:zoom
    };
  }


  angular.module('npact')
    .factory('GraphingCalculator', GraphingCalculator)
  ;

}());

