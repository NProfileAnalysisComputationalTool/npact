angular.module('npact')
  .factory('Grapher', function(K, $log, GraphingCalculator, $rootScope, $compile, NProfiler, TrackReader, $q, Evt) {
    'use strict';
    function addMany(container, children) {
      if(children && children.length) {
        container.add.apply(container, children);
      }
    }

    function Grapher(element, opts) {
      this.$element = jQuery(element);
      angular.extend(this, opts);
      // invariants: startBase, endBase, margin
      this.startBaseM = Math.max(this.startBase - this.margin, 0);
      this.endBaseM = this.endBase + this.margin;

      this.stage = new K.Stage({
        container: element,
        height: this.height,
        width: this.width
      });

      this.stage.add(this.leftLayer(), this.chartLayer());
    }
    var GP = Grapher.prototype;

    //Tie margin to 0 until we figure out clipping issues (not sure margin helped anyways.)
    GP.margin = 0;

    GP.destroy = function() {
      if(this.stage) {
        this.stage.destroy();
      }
    };

    GP._onProfilePoints = null;
    GP.getProfilePoints = function() {
      if(!this._onProfilePoints) {
        // start slicing the NProfile into Kinetic-compatible [x1, y1,
        // x2, y2, ...] lists, make a promise for the completed group of
        // points
        var r= [], g= [], b= [];
        this.onProfilePoints = NProfiler
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
      },
          makeLabel = function(header) {
            var txtOpts = angular.extend({}, defaultTextOpts, header);
            return new K.Text(txtOpts);
          },
          lbls =  _.map(this.headers, makeLabel);

      addMany(g, lbls);
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

    var leftLayerCache = {};
    GP.leftLayer = function() {
      var key = angular.toJson(this.headers) + this.axisTitle,
          opts = {
            width: this.leftPadding,
            x: 0, y: 0,
            height: this.height
          },
          layer = new K.FastLayer(opts);

      if(leftLayerCache[key]) {
        leftLayerCache[key].then(function(image) {
          layer.add(new K.Image(angular.extend(opts, {image: image})));
          layer.draw();
        });
      }else{
        this.yAxisGroup(layer);
        this.headerGroup(layer);
        leftLayerCache[key] = $q(function(resolve) {
          $log.log('saving leftLayer image');
          layer.toImage(angular.extend(opts, {callback: resolve}));
        });
      }
      return layer;
    };

    GP.genomeGroup = function() {
      var self = this,
          gx = this.m.graph.x,
          dragRes = {x: gx, y: 0},
          g = new K.Group({
            x: gx,
            draggable: true,
            dragBoundFunc: function(pos) {
              // pos is the distance dragged, except `pos.x` always
              // starts at `gx`. Something due to container coordinate
              // system vs stage coordinate system. `pos.y` is not
              // similarly affected.
              var nextOffset = - (pos.x - gx);
              // change position via offsetX
              this.offsetX(nextOffset);
              // never change position via this return value
              return dragRes;
            }
          });
      // remember where we are at the start/end of a drag, so dragging can
      // match our "scroll"
      g.on('dragstart dragend', function() {
        // kill ALL tooltips, not just the ones on this graph
        jQuery('.qtip').qtip('destroy');
      });

      g.on('dragend', function(evt) {
        var oldStartBase = self.startBase,
            newStartBase = (this.offsetX() / self.xaxis.scaleX) + oldStartBase;
        // tell the world
        self.onPan(
          {oldStartBase: oldStartBase, newStartBase: newStartBase, evt: evt});
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
      g.add(new K.Rect({x: 0, y: this.m.graph.y,
                        width: this.m.graph.w,
                        height: this.m.graph.h}));
      return (this._genomeGroup = g);
    };

    GP.onDblClick = function(evt) {
      jQuery('.qtip').qtip('destroy');

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

    GP.chartLayer = function() {
      var m = this.m,
          l = new K.Layer({
            clip: {
              x: m.graph.x-1, y: 0,
              width: m.graph.w+2,
              height: 1000
            }
          }),
          // frame around the graph
          border = new K.Rect({
            x: m.graph.x,
            y: m.graph.y,
            width: m.graph.w,
            height: m.graph.h,
            stroke: this.borderColor
          });

      l.add(border, this.genomeGroup());
      return l;
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
            scaleX: xaxis.scaleX,
            offsetX: this.startBase
          });

      addMany(g, this.drawAxisTicks(xaxis.ticks));
      var labels = this.drawAxisLabels(xaxis.labels);
      addMany(g, labels);

      // call `centerXLabel(txt, xaxis.scaleX)`
      _.map(labels, function(txt) {
        centerXLabel(txt, xaxis.scaleX);
      });

      return (this._xAxisGroup = g);
    };

    GP.profileGroup = function() {
      var m = this.m, xaxis = this.xaxis,
          colors = this.colors,
          g = new K.Group({
            x: 0, y: m.graph.y,
            height: m.graph.h, width: xaxis.length,
            scaleX: xaxis.scaleX,
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
          }
      ;
      _.each(shades, function(x) {
        g.add(new K.Rect(angular.extend(x, shadeOpts)));
      });


      this.getProfilePoints()
        .then(function(points) {
          if(!g.getLayer()) {return;}
          _.forEach(points, function(v, k) {
            g.add(buildLine(v, colors[k]));
          });
          g.draw();
        });

      return (this._profileGroup = g);
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
      this.stage.add(this.chartLayer(), this.leftLayer());

      var gg = this._genomeGroup;
      gg.add(this.xAxisGroup(), this.profileGroup());

      var drawings = _.map(this.headers, function(hdr) {
        return this.trackSlice(hdr.text)
          .then(function(data) {
            switch(hdr.lineType) {
            case 'extracts':
              return this.cdsGroup(hdr, data);
            case 'hits':
              return this.drawHit(hdr, data);
            default:
              throw new Error("don't know how to draw " + hdr);
            }
          }.bind(this)).then(function(img) { gg.add(img); });
      }, this);

      return $q.all(drawings).then(function() {
        this.stage.draw();
        $log.log('drew npactGraph', this.startBase,
                 'took', new Date() - t1, 'ms');
      }.bind(this));
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

    GP.cdsGroup = function(header, cds) {
      var xaxis = this.xaxis, $el = this.$element,
          colors = this.colors,
          g = new K.Group({
            x: 0, y: 0,
            scaleX: xaxis.scaleX,
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


      // TODO: move tooltip control to the `npactGraph`, just throw an
      // event or callback here with the extract, let `npactGraph`
      // handle this stuff
      g.on('click', function(evt) {
        var scope = $rootScope.$new(),
            tpl = '<div npact-extract="extract"></div>';
        scope.extract = evt.target.getAttrs().extract;
        $el.qtip({
          content: {text: $compile(tpl)(scope)},
          position: {
            my: scope.extract.complement === 0 ? 'top center' : 'bottom center',
            target: [evt.evt.pageX, evt.evt.pageY]
          },
          show: {event: 'tooltipShow.npact'},
          hide: {event: 'tooltipHide.npact'}
        });
        $el.trigger('tooltipShow.npact');
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
            scaleX: xaxis.scaleX,
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
        var y = baseY + (hit.complement ? -offset : offset);

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
