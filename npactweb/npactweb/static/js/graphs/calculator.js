angular.module('npact')
  .service('GraphingCalculator', function GraphingCalculator(K, Utils) {
    'use strict';
    var self = this;

    self.stops = function(a,b){
      var length = b - a,
          // get even stops; look for one lower order of magnitude
          // than our length - TODO: too many stops for lower ranges
          interval = Utils.orderOfMagnitude(length, -1),
          // want to get ticks aligned on `interval`. a and b might
          // have weird offsets
          start = Math.floor(a/interval)*interval,
          end = Math.floor(b/interval)*interval;

      return {
        interval: interval,
        // capture some margin
        stops: _.range(
          Math.max(0, start - interval),
          // want to be start/end INCLUSIVE so pad the end
          end + 2*interval, interval)
      };
    };

    self.xaxis = function(opts){

      var makeTicks = function(lbl, n){
        var x = parseInt(start + n*interval);
        return {
          x: x, y: 0,
          x2: x, y2: opts.profileTicks
        };
      };

      var makeLabels = function(lbl, n){
        return {
          x: ticks[n].x, y: opts.profileTicks,
          coord: ticks[n].x,
          text: lbl,
          scaleX: labelScaleX
        };
      };

      var m = self.chart(opts),
          s = self.stops(opts.startBase, opts.endBase),
          length = opts.endBase - opts.startBase,
          interval = s.interval,
          ss = s.stops,
          // want to capture some margin
          start = ss[0],
          end = _.last(ss),
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
    };


    /**
     * determine new offset and basesPerGraph for a zoom, retaining
     * relative position
     *
     * so if you zoom in on gene 2000, after the zoom the mouse cursor
     * is still on 2000, but the rest of the graph has shifted. Treats
     * the genes as a number line.
     *
     * @param {Number} startBase - low end of the graph that we're zooming on
     * (gene space)
     * @param {Number} zoomOnPct - where we're zooming to, as a percent of
     * the graph
     * @param {Number} basesPerGraph - how many bases per graph,
     * indicates current zoom level
     * @param {Number} offset - the current offset for the graphs
     * (gene space)
     * @param {Boolean} zoomingOut - true for zooming out, false for
     * zooming in
     * @returns {Object} new offset and basesPerGraph such that we remain focused on the selected gene
     */
    self.zoom = function(opts){
      //
      // O = offset from L such that Z is positionally at the same
      // point on the screen
      // LH = distance between L and H, a measure of how zoomed in we are
      //
      // all units are in gene space

      var zoomToBase = opts.startBase + Math.floor(opts.zoomOnPct*opts.basesPerGraph),
          row = Math.floor((opts.startBase - opts.offset)/opts.basesPerGraph),
          zoomFactor = opts.zoomingOut ? 2 : 0.5,
          newbasesPerGraph = opts.basesPerGraph * zoomFactor,
          newstartBase = zoomToBase - Math.floor(newbasesPerGraph * opts.zoomOnPct),
          newoffset = newstartBase - (row*newbasesPerGraph);

      return {offset: newoffset, basesPerGraph: newbasesPerGraph};
    };

    /**
     * calculate measurements about the chart
     */
    self.chart = function(opts){
      var yStops = [100, 80, 60, 40, 20, 0],
          yAxisTicks = yStops.length - 1,
          profileHeight = opts.height - opts.headerY - opts.axisLabelFontsize -
            3*opts.profileTicks,

          // the line graph, excluding tick marks and axis labels
          g = {
            x: opts.leftPadding,
            y: opts.height - profileHeight -
              opts.axisLabelFontsize - // x-axis labels
              2*opts.profileTicks, // tick marks
            h: profileHeight,
            w: opts.width - opts.leftPadding - opts.rightPadding
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
    };

    /**
     * find a position that will align the two rectangles on their center
     *
     * @param {Object} rect - the rectangle to align to {width,height}
     * @param {Object} toAlign - the rectangle you want to move {width,height}
     * @returns {Object} x,y coordinates to shift `toAlign` by, from the top/left
     */
    self.alignRectangles = function(rect, toAlign){
      // get the top left coordinate of the `toAlign` rect
      return {x: rect.x + (rect.width/2) - (toAlign.width/2),
              y: rect.y + (rect.height/2) - (toAlign.height/2)};
    };

    /**
     * determine background shading
     *
     * @param {Object} opts - {startBase, endBase, interval}, where
     * interval is the space between axis labels
     */
    self.shades = function(opts) {
      var width = opts.interval / 2,
          // want the X to start on 0 or multiples of `interval`
          nearestEvenStartBase = Math.floor(opts.startBase / opts.interval) *
            opts.interval,
          xes = _.range(nearestEvenStartBase, opts.endBase, opts.interval)
      ;
      return _.map(xes, function(x) {
        var w = width;
        // handle shades falling off the left side
        if (x < opts.startBase){
          w -= opts.startBase - x;
          x = opts.startBase;
        }
        // handle shades falling off the right side
        if(x + width > opts.endBase){
          w = opts.endBase - x;
        }
        return { width: w, x: x };
      });
    };
  })
;
