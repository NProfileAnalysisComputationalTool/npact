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

  .factory('Grapher', function(K, $log, GraphingCalculator, Tooltip, NProfiler, TrackReader, $q, Evt, Utils, npactConstants) {
    'use strict';

    var style = npactConstants.graphStyle;

    function addMany(container, children) {
      if(children && children.length) {
        container.add.apply(container, children);
      }
    }

    function boundingBox(container) {
      var children = [];
      if (container.hasChildren && container.hasChildren()) {
        children = container.getChildren();
      }
      else if (_.isArray(container)) {
        children = container;
      }

      var l,r,t,b;
      _.forEach(children, function (child, idx) {
        var cr,cb;
        var cl = child.x() - child.offsetX();
        var ct = child.y() - child.offsetY();

        if(child.hasChildren()) {
          child = new Kinetic.Shape(boundingBox(child));
          //already applied coords of group inside our current container,
          //but the items in there might have additional offset (i.e. not
          //be at 0,0)
          cl += child.x();
          ct += child.y();
        }
        cr = cl + child.width();
        cb = ct + child.height();

        if (idx === 0) {
          l = cl;  r = cr;
          t = ct;  b = cb;
        } else {
          l = Math.min(l, cl);
          r = Math.max(r, cr);
          t = Math.min(t, ct);
          b = Math.max(b, cb);
        }
      });
      return  {
        x: l,
        y: t,
        width: r - l,
        height: b - t
      };
    }

    function digitCount(x) {
      return (1 + Math.floor(Math.log(x) / Math.LN10));
    }

    function Grapher(element, opts) {
      this.$element = jQuery(element);
      angular.extend(this, opts);
      // invariants: startBase, endBase
      var length = this.endBase - this.startBase;
      this.xaxis = {
        length: length,
        scaleX: this.m.graph.w / length
      };
      this.margin = Utils.orderOfMagnitude(length, -1);
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


    GP.headerGroup = function(g) {
      // TODO: derive opts.leftPadding
      // * add all labels
      // * loop over K.Text objects, find max width
      // * leftPadding = maxWidth + headerLabelPadding
      // * loop over K.Text objects, set width to leftPadding

      var defaultTextOpts = _.defaults({
        align: 'right', x: 0,
        width: style.leftPadding - 2 * style.paddingUnit
      }, style.tracks.text);

      addMany(g, _.map(this.tracks, function(track) {
        defaultTextOpts.text = track.text;
        defaultTextOpts.y = track.y;
        defaultTextOpts.height = track.height;
        return new K.Text(defaultTextOpts);
      }));
      return g;
    };


    GP.yAxisTicks = function(height) {
      // Draw marks for percentage on domain of [0,100] then scale
      // that to the height we need
      var axisg = new K.Group({scaleY: height / 100, strokeScaleEnabled: false});
      var profileSpec = style.profile,
          axisSpec = style.profile.axis;

      _.forEach(profileSpec.yStops, function(ystop) {
        // The tick itself
        axisg.add(new K.Line({
          points: [-axisSpec.tickLength, ystop, 0, ystop],
          stroke: profileSpec.borderColor, strokeScaleEnabled: false}));

        //label for the tick
        var width  = 2 * axisSpec.text.fontSize;
        axisg.add(new K.Text(_.defaults({
          text: ystop,
          x: - width - style.paddingUnit * 2,
          width: width,
          y: 100 - ystop, // draw from the top
          offsetY: axisSpec.text.fontSize / 2, // center the text at that point
          scaleY: 100 / height, //text itself needs to be unscaled
          align: 'right'
        }, axisSpec.text)));
      }, this);
      return axisg;
    };
    GP.yAxisTitle = function(height) {
      var title = new K.Text(_.defaults({
        rotation: -90,
        fontSize: style.profile.titleFontSize,
        text: this.axisTitle
      }, style.profile.axis.text));
      //All the width height x/y calcuations go on pre-rotation, so
      // swap, align vertically and then stick in a group sized right
      // so elsewhere we don't have to think about htat.
      var w = title.getHeight(), h = title.getWidth();
      title.y((height + h) / 2); // Align vertically in our space
      title.x(-w);               // we're drawing leftwards
      var g = new K.Group({width: w, height: h});
      g.add(title);
      return g;
    };

    GP.yAxisGroup = function() {
      //This group is drawn with m.graph.x as x=0 and we draw left from there.
      var m = this.m;
      var g = new K.Group({x: m.yaxis.x, y: m.yaxis.y});
      var ticks = this.yAxisTicks(m.yaxis.height);
      var title = this.yAxisTitle(m.yaxis.height);
      title.offsetX(boundingBox(ticks).width + 2 * style.paddingUnit);
      g.add(ticks); g.add(title);
      return g;
    };

    GP._leftLayerImage = function() {
      return  $q(_.bind(function(resolve) {
        var opts = {
          width: style.leftPadding,
          x: 0, y: 0,
          height: this.m.height
        };
        var layer = new K.Layer(opts);
        layer.add(this.yAxisGroup());
        this.headerGroup(layer);
        layer.toImage(angular.extend(opts, {callback: resolve}));
      }, this));
    };

    var leftLayerCache = {};  // Shared amongst all Graphers
    GP.leftLayer = function(stage) {
      var key = angular.toJson(this.tracks) + this.axisTitle;
      var layer = new K.FastLayer();
      stage.add(layer);
      return (leftLayerCache[key] || (leftLayerCache[key] = this._leftLayerImage()))
        .then(function(image) {
          layer.add(new K.Image( {image: image}));
          layer.draw();
        });
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
        offsetX: this.startBase,
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

      var dgAdd = _.bind(dg.add, dg);
      // need a shape that can be clicked on to allow dragging the
      // entire canvas
      dgAdd(new K.Rect({x: this.startBase, y: 0,// fill: '#BDD',
                         width: this.xaxis.length,
                         height: this.m.height}));

      var p3 = this.xAxisGroup().then(dgAdd);
      var p2 = this.tracksGroup().then(dgAdd);
      var p1 = this.profileGroup().then(dgAdd);

      return $q.all([p1, p2, p3]).then(function() {
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
            x: m.graph.x, y: m.graph.y,
            width: m.graph.w, height: m.graph.h,
            stroke: style.profile.borderColor
          });
      stage.add(l);
      l.add(border);
      l.draw();
      return $q.when(l);
    };

    GP.xAxisGroup = function() {
      var xaxis = this.xaxis,
          stops = GraphingCalculator.stops(this.startBaseM, this.endBaseM),
          g = new K.Group({
            width: xaxis.length,
            x: 0, y: this.m.xaxis.y
          }),
          tickOpts = {
            stroke: style.profile.borderColor, strokeScaleEnabled: false
          },
          labelOpts = _.assign({
            scaleX: 1/xaxis.scaleX,
            y: style.profile.tickLength
          }, style.profile.axis.text),
          shadeOpts = {
            height: this.m.graph.h,
            offsetY: this.m.graph.h,
            fill: style.profile.shadeColor,
            width: stops.interval / 2
          };

      _.forEach(stops.stops, function(stop) {
        tickOpts.points = [stop, 0, stop, style.profile.tickLength];
        labelOpts.text = stop;
        labelOpts.x = stop;
        labelOpts.offsetX =
          (digitCount(stop) * style.profile.axis.text.fontSize / 4);
        shadeOpts.x = stop;
        g.add(new K.Line(tickOpts));
        g.add(new K.Text(labelOpts));
        g.add(new K.Rect(shadeOpts));
      });
      return $q.when(g);
    };

    GP.profileGroup = function() {
      var m = this.m, xaxis = this.xaxis,
          colors = this.colors,
          g = new K.Group({
            x: 0, y: m.graph.y,
            height: m.graph.h, width: xaxis.length,
            // convert % to px
            scaleY: m.graph.h / 100
          }),
          buildLine = function(points, color) {
            return new K.Line({
              points: points,
              stroke: color,
              strokeWidth: 1,
              strokeScaleEnabled: false
            });
          };

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

      var glp = this.genomeLayer(this.stage);
      var llp = this.leftLayer(this.stage);
      var flp = this.frameLayer(this.stage);
      return $q.all([llp, flp, glp])
        .then(function() {
          $log.log("Finished draw at", newOpts.startBase, "in", new Date() - t1);
        });
    };

    /**
     * Paint a green border around the graph to help verify the height and width
     */
    GP.testFrame = function() {
      var testLayer = new K.Layer(
        {x: 0, y: 0, width: this.width, height: this.m.height});
      testLayer.add(
        new K.Rect({x: 0, y: 0, width: this.width, height: this.m.height,
                    stroke: 'green', strokeWidth:1
                   }));
      this.stage.add(testLayer);
      testLayer.draw();
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

    GP.tracksGroup = function() {
      var g = new K.Group();
      return $q.all(_.map(this.tracks, function(track) {
        return this.trackSlice(track.text)
          .then(_.bind(function(data) {
            switch(track.lineType) {
            case 'extracts':
              return this.cdsGroup(track, data);
            case 'hits':
              return this.drawHit(track, data);
            default:
              throw new Error("don't know how to draw " + track);
            }
          }, this));
      }, this))
        .then(function(list) {
          addMany(g, list);
          return g;
        });
    };
    GP.cdsGroup = function(track, cds) {
      var xaxis = this.xaxis, $el = this.$element,
          colors = this.colors,
          g = new K.Group({
            x: 0, y: 0
          }),
          colorNames = 'rgb',
          y = track.y,
          arrowHeight = style.tracks.arrow.height,
          ahw = style.tracks.arrow.width / xaxis.scaleX,
          textOpts = _.assign({
            scaleX: 1/xaxis.scaleX,
            strokeScaleEnabled: false
          }, style.tracks.text),
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
              baseY = isComplement ? y + arrowHeight : y,
              arrowPointY = baseY + arrowHeight / 2,
              arrowMaxY = baseY + arrowHeight,
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

    GP.drawHit = function(track, hits) {
      var startBase = this.startBase,
          endBase = this.endBase,
          colors = this.colors,
          colorNames = 'rgb',
          offset = track.height / 4,
          hitStrokeWidth = offset / 2,
          guideYOffset = 2,
          // arrow sticks out ~1%
          guideArrowXOffset = Math.floor(0.01 * (endBase - startBase)),
          baseY = track.y + (track.height / 2),
          g = new K.Group({ x: 0, y: 0 }),
          guideLineOpts = { stroke: '#ddd' }
      ;
      // draw the guide lines
      g.add(new K.Line(angular.extend({
        points: [startBase, baseY + guideYOffset,
                 endBase, baseY + guideYOffset,
                 endBase - guideArrowXOffset, track.y + track.height]
      }, guideLineOpts)));
      g.add(new K.Line(angular.extend({
        points: [startBase + guideArrowXOffset, track.y,
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
