angular.module('npact')
  .service('GraphingCalculator', function GraphingCalculator(K, Utils, npactConstants) {
    'use strict';
    var self = this;

    this.partition = function(opts) {
      var offset = opts.offset || 0,
          startBase = opts.startBase + offset,
          endBase = opts.endBase,
          basesPerGraph = opts.basesPerGraph;
      if(basesPerGraph)
        return _.range(startBase, endBase, basesPerGraph);
      else
        return [];
    };


    self.stops = function(a, b){
      var length = b - a,
          // get even stops; look for one lower order of magnitude
          // than our length
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


    self.trackSizeCalc = function(activeTracks) {
      var offset = 0;
      var tracks = _.map(activeTracks, function(cfg) {
        var h = npactConstants.trackHeights[cfg.lineType],
            y = offset;
        offset += h;
        return {
          text: cfg.text,
          lineType: cfg.lineType,
          y: y,
          height: h
        };
      });
      return {
        tracks: tracks,
        totalTrackHeight: offset
      };
    };

    /**
     * calculate measurements about the chart
     */
    self.chart = function(opts) {
      var style = npactConstants.graphStyle;

      var pstyle = style.profile;
      var graphTop = opts.totalTrackHeight + style.paddingUnit,
          xAxisTop = graphTop + pstyle.height,
          xAxisHeight = pstyle.axis.text.fontSize + pstyle.tickLength,
          totalHeight = xAxisTop + xAxisHeight;

      return {
        height: totalHeight,
        graph: {
          x: style.leftPadding,
          y: graphTop,
          h: style.profile.height,
          w: opts.width - style.leftPadding
        },
        xaxis: {
          height: xAxisHeight,
          y: xAxisTop
        },
        yaxis: {
          height: style.profile.height,
          y: graphTop,
          x: style.leftPadding
        }
      };
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
  })
;
