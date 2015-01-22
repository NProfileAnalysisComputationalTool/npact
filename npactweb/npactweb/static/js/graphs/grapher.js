angular.module('npact')

  .service('Tooltip', function($log, $rootScope, $compile) {
    'use strict';
    this.show = function ($el, extract, pageX, pageY) {
      this.clearAll();
      var scope = $rootScope.$new(),
          tpl = '<div npact-extract="extract"></div>';
      scope.extract = extract;
      $el.qtip({
        content: {text: $compile(tpl)(scope)},
        position: {
          my: scope.extract.complement === 0 ? 'top center' : 'bottom center',
          target: [pageX, pageY]
        },
        show: {event: 'tooltipShow.npact'},
        hide: {event: 'tooltipHide.npact'}
      });
      $el.trigger('tooltipShow.npact');
    };
    this.clearAll = function() {
      jQuery('.qtip').qtip('destroy');
    };
  })

  .factory('Grapher', function(K, $log, GraphingCalculator, Tooltip, NProfiler, TrackReader, $q, Evt, Utils) {
    'use strict';

    function addMany(container, children) {
      if(children && children.length) {
        container.add.apply(container, children);
      }
    }

    function Grapher(element, opts) {
      this.$element = jQuery(element);
      angular.extend(this, opts);
      // invariants: startBase, endBase
      this.margin = Utils.orderOfMagnitude(this.endBase - this.startBase, -1);
      this.startBaseM = Math.max(this.startBase - this.margin, 0);
      this.endBaseM = this.endBase + this.margin;

      this.stage = new K.Stage({
        container: element,
        height: this.height,
        width: this.width
      });
    }
    var GP = Grapher.prototype;

    GP.destroy = function() {
      if(this.stage) { this.stage.destroy(); }
    };

    GP._onProfilePoints = null;
    GP.getProfilePoints = function() {
      if(!this._onProfilePoints) {
        // start slicing the NProfile into Kinetic-compatible [x1, y1,
        // x2, y2, ...] lists, make a promise for the completed group of
        // points
        var r= [], g= [], b= [];
        this._onProfilePoints = NProfiler
          .slice({ startBase: this.startBaseM, endBase: this.endBaseM,
                   onPoint: function(coord, rv, gv, bv) {
                     //invert because drawing is from the top so 100% is 0 pix
                     r.push(coord); r.push(100.0 - rv);
                     g.push(coord); g.push(100.0 - gv);
                     b.push(coord); b.push(100.0 - bv);
                   }})
          .then(function(opts) { return {r: r, g: g, b: b}; });
      }
      return this._onProfilePoints;
    };

    GP.drawAxisTicks = function(ticks) {
      var tickOpts = {x: 0, y: 0,
                      stroke: this.borderColor,
                      strokeScaleEnabled: false};
      return _.map(ticks, function(t) {
        tickOpts.points = [t.x, t.y, t.x2, t.y2];
        return new K.Line(tickOpts);
      });
    };

    GP.drawAxisLabels = function(labels, textOpts) {
      // draw labels at the right spacing
      var defaultTextOpts = {
        fontSize: this.axisLabelFontsize,
        fill: this.axisFontcolor
      };

      return _.map(labels, function(lbl) {
        var txtOpts = angular.extend({}, lbl, defaultTextOpts, textOpts);
        return new K.Text(txtOpts);
      });
    };

    GP.headerGroup = function(g) {
      // TODO: derive opts.leftPadding
      // * add all labels
      // * loop over K.Text objects, find max width
      // * leftPadding = maxWidth + headerLabelPadding
      // * loop over K.Text objects, set width to leftPadding

      var defaultTextOpts = {
        align: 'right', x: 0,
        fontSize: this.headerLabelFontsize,
        width: this.leftPadding - this.headerLabelPadding,
        fill: this.headerLabelFontcolor
      };

      addMany(g, _.map(this.headers, function(header) {
        return new K.Text(angular.extend(defaultTextOpts, header));
      }));
      return g;
    };

    GP.yAxisGroup = function(g) {
      var m = this.m;

      addMany(g, this.drawAxisTicks(m.yaxis.ticks));
      // draw labels at the right aligned
      addMany(g, this.drawAxisLabels(m.yaxis.labels, {align: 'right'}));

      // the title
      var title = new K.Text({
        y: 0, x: 0, // reposition this below
        rotation: -90,
        fill: this.axisFontcolor,
        fontSize: this.axisTitleFontsize,
        text: this.axisTitle
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

    var leftLayerCache = {};  // Shared amongst all Graphers
    GP.leftLayer = function(stage) {
      var key = angular.toJson(this.headers) + this.axisTitle,
          opts = {
            width: this.leftPadding,
            x: 0, y: 0,
            height: this.height
          },
          layer = new K.FastLayer(opts);
      stage.add(layer);

      if(leftLayerCache[key]) {
        return leftLayerCache[key].then(function(image) {
          layer.add(new K.Image(angular.extend(opts, {image: image})));
          layer.draw();
        });
      }
      else {
        this.yAxisGroup(layer);
        this.headerGroup(layer);
        var d = $q.defer();
        layer.toImage(angular.extend(opts, {callback: function(v) {
          d.resolve(v);
          layer.draw();
        }}));
        return (leftLayerCache[key] = d.promise);
      }
    };

    GP.genomeLayer = function(stage) {
      //Everything in this layer is in the coordinate system of the
      //genome and will be scaled to pixels by Kinetic
      var l = new K.Layer({
        x: this.m.graph.x,
        clip: {
          x: 0, y: 0,
          width: this.m.graph.w,
          height: 1000
        }
      });
      stage.add(l);
      var dg = new K.Group({
        scaleX: this.xaxis.scaleX,
        x: 0,
        draggable: true,
        dragBoundFunc: function(pos) {
          pos.y = 0; //Disallow vertical movement
          return pos;
        }
      });
      l.add(dg);
      dg.on('dragstart dragend', Tooltip.clearAll);
      dg.on('dragend', _.bind(this.onDragEnd, this));
      dg.on('mouseover', function() { document.body.style.cursor = 'pointer'; });
      dg.on('mouseout', function() { document.body.style.cursor = 'default'; });
      dg.on('dblclick', _.bind(this.onDblClick, this));

      // need a shape that can be clicked on to allow dragging the
      // entire canvas
      dg.add(new K.Rect({x: 0, y: 0, //fill: '#BDD',
                          width: this.xaxis.length,
                          height: this.height}));

      //Everything in this layer needs to be part of that draggable.
      dg.add(this.xAxisGroup());

      var p1 = this.profileGroup().then(function(pg) { dg.add(pg); });
      var p2 = $q.all(_.map(this.headers, this.makeHeader, this))
            .then(function(list) { addMany(dg, list); });
      return $q.all([p1, p2]).then(function() {
        l.draw();
        return l;
      });
    };

    GP.onDragEnd = function(evt) {
      $log.log('dragEnd', evt.target.x(), evt.target.getScaleX());
      this.onPan(-evt.target.x() / evt.target.getScaleX());
    };

    GP.onDblClick = function(evt) {
      Tooltip.clearAll();

      var zoomOnPx = evt.evt.layerX - this.m.graph.x,
          zoomOnPct = zoomOnPx / this.m.graph.w;
      // tell the world
      this.onZoom({
        evt: evt,
        // these keys must match what's expected by `GraphingCalculator.zoom`
        startBase: this.startBase,
        zoomOnPct: zoomOnPct,
        zoomingOut: evt.evt.shiftKey
      });
    };

    GP.frameLayer = function(stage) {
      var m = this.m,
          l = new K.FastLayer(),
          // frame around the graph
          border = new K.Rect({
            x: m.graph.x,
            y: m.graph.y,
            width: m.graph.w,
            height: m.graph.h,
            stroke: this.borderColor
          });
      stage.add(l);

      l.add(border);
      l.draw();
      return $q.when(l);
    };

    function centerXLabel(txt, scaleX) {
      // center on the tick marks, at this point we know how wide
      // this text is
      var w = txt.getWidth() / scaleX,
          x = txt.getAttr('coord'),
          newX = x - w/2;
      txt.x(newX);
    }

    GP.xAxisGroup = function() {
      var m = this.m,
          xaxis = this.xaxis,
          g = new K.Group({
            x: 0, y: m.graph.h + m.graph.y,
            width: xaxis.length,
            offsetX: this.startBase
          });

      addMany(g, this.drawAxisTicks(xaxis.ticks));
      _.forEach(this.drawAxisLabels(xaxis.labels), function(lbl) {
        g.add(lbl);
        centerXLabel(lbl, xaxis.scaleX);
      });
      return g;
    };

    GP.profileGroup = function() {
      var m = this.m, xaxis = this.xaxis,
          colors = this.colors,
          g = new K.Group({
            x: 0, y: m.graph.y,
            height: m.graph.h, width: xaxis.length,
            // convert % to px
            scaleY: m.graph.h / 100,
            offsetX: this.startBase
          }),
          buildLine = function(points, color) {
            return new K.Line({
              points: points,
              stroke: color,
              strokeWidth: 1,
              strokeScaleEnabled: false
            });
          },
          shades = GraphingCalculator.shades({
            startBase: this.startBaseM,
            endBase: this.endBaseM,
            interval: xaxis.interval
          }),
          shadeOpts = {
            y: 2, height: 96,
            fill: this.profileShadeColor
          };

      _.each(shades, function(shd) {
        shadeOpts.width = shd.width;
        shadeOpts.x = shd.x;
        g.add(new K.Rect(shadeOpts));
      });

      return this.getProfilePoints()
        .then(function(points) {
          _.forEach(points, function(v, k) {
            g.add(buildLine(v, colors[k]));
          });
          return g;
        });
    };

    /**
     * get the relevant slice of track data for the given header
     *
     * Keeps a cache on this to avoid duplicate computation
     *
     * @return {Promise}
     */
    GP.trackSlice = function(name) {
      if(!this.trackSliceCache) { this.trackSliceCache = {}; }

      return this.trackSliceCache[name] ||
        (this.trackSliceCache[name] = TrackReader.slice({
          name: name,
          startBase: this.startBaseM,
          endBase: this.endBaseM
        }));
    };

    GP.redraw = function(newOpts) {
      var t1 = new Date();
      angular.extend(this, newOpts);
      this.stage.destroyChildren();
      this.stage.setWidth(this.width);
      this.stage.setHeight(this.m.height);
      var llp = this.leftLayer(this.stage);
      var flp = this.frameLayer(this.stage);
      var glp = this.genomeLayer(this.stage);
      return $q.all([llp, flp, glp])
        .then(function() {
          $log.log("Finished draw at", newOpts.startBase, "in", new Date() - t1);
        });
    };

    function centerExtractLabel(txt, scaleX) {
      // now that lbl is on the canvas, we can see what it's
      // height/width is
      var arrowBounds = txt.getAttr('arrowBounds'),
          pos = GraphingCalculator.alignRectangles(
            arrowBounds,
            {
              // convert to gene space
              width: txt.getWidth() / scaleX,
              height: txt.getHeight()
            });
      pos.x = Math.max(pos.x, arrowBounds.x+1);
      txt.position(pos);
    }

    GP.makeHeader = function(hdr) {
      return this.trackSlice(hdr.text)
        .then(_.bind(function(data) {
          switch(hdr.lineType) {
          case 'extracts':
            return this.cdsGroup(hdr, data);
          case 'hits':
            return this.drawHit(hdr, data);
          default:
            throw new Error("don't know how to draw " + hdr);
          }
        }, this));
    };
    GP.cdsGroup = function(header, cds) {
      var xaxis = this.xaxis, $el = this.$element,
          colors = this.colors,
          g = new K.Group({
            x: 0, y: 0,
            offsetX: this.startBase
          }),
          colorNames = 'rgb',
          y = header.y,
          ahHalfHeight = this.headerArrowHeight/2,
          ahw = this.headerArrowWidth/xaxis.scaleX,
          textOpts = {
            fontSize: this.headerArrowFontsize,
            fill: this.axisFontcolor,
            scaleX: 1/xaxis.scaleX,
            strokeScaleEnabled: false
          },
          arrowOpts = {
            x: 0, y: 0,
            closed: true,
            strokeWidth: 1,
            strokeScaleEnabled: false
          }
      ;

      _.forEach(cds, function(x) {
          var isComplement = x.complement === 1,
              c = colors[colorNames[x.phase]],
              baseY = isComplement ? y + this.headerArrowHeight : y,
              arrowPointY = baseY + ahHalfHeight,
              arrowMaxY = baseY + this.headerArrowHeight,
              shape = isComplement ?
                [
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
              arrowBounds = isComplement ?
                {
                  x: x.start+ahw, y: baseY,
                  width: x.end-x.start-ahw, height: this.headerArrowHeight
                } : {
                  x: x.start, y: baseY,
                  width: x.end-x.start-ahw, height: this.headerArrowHeight
                },
              line = new K.Line(angular.extend({
                extract: x,
                points: shape,
                stroke: c
              }, arrowOpts)),
              // render the name, too
              lbl = new K.Text(angular.extend({
                extract: x,
                arrowBounds: arrowBounds,
                text: x.name
              }, textOpts)),
              // need a dummy group for clipping, `Text` doesn't
              // support clip directly
              lblGroup = new K.Group({clip: arrowBounds});

          lblGroup.add(lbl);
          g.add(line, lblGroup);
          // now that lbl is on the canvas, we can see what it's
          // height/width is
          centerExtractLabel(lbl, xaxis.scaleX);
      }, this);

      g.on('click', function(evt) {
        Tooltip.show($el, evt.target.getAttrs().extract, evt.evt.pageX, evt.evt.pageY);
      });
      return g;
    };

    GP.drawHit = function(header, hits) {
      var xaxis = this.xaxis,
          startBase = this.startBase,
          endBase = this.endBase,
          colors = this.colors,
          colorNames = 'rgb',
          offset = header.height / 4,
          hitStrokeWidth = offset / 2,
          guideYOffset = 2,
          // arrow sticks out ~1%
          guideArrowXOffset = Math.floor(0.01 * (endBase - startBase)),
          baseY = header.y + (header.height / 2),
          g = new K.Group({
            x: 0, y: 0,
            offsetX: startBase
          }),
          guideLineOpts = {
            stroke: '#ddd'
          }
      ;
      // draw the guide lines
      g.add(new K.Line(angular.extend({
        points: [startBase, baseY + guideYOffset,
                 endBase, baseY + guideYOffset,
                 endBase - guideArrowXOffset, header.y + header.height]
      }, guideLineOpts)));
      g.add(new K.Line(angular.extend({
        points: [startBase + guideArrowXOffset, header.y,
                 startBase, baseY - guideYOffset,
                 endBase, baseY - guideYOffset]
      }, guideLineOpts)));

      // draw each hit
      _.forEach(hits, function(hit) {
        var y = baseY + (hit.complement ? offset : -offset);

        g.add(
          new K.Line({
            points: [hit.start, y, hit.end, y],
            stroke: colors[colorNames[hit.phase]],
            strokeWidth: hitStrokeWidth
          }));
      });

      return g;
    };
    return Grapher;
  });
