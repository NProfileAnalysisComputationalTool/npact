angular.module('npact')

  .factory('Grapher', function($log, $q, $timeout, GraphConfig, GraphingCalculator,
                        Tooltip, NProfiler, Utils, npactConstants) {
    'use strict';
    var K = Konva;
    var style = npactConstants.graphStyle;

    function Grapher($element, $scope, opts) {
      this.$element = $element;
      this.$scope = $scope;
      _.assign(this, opts);
      //$element.find('.coords').detach(); $element.prepend("<div class='coords'> start:"+opts.startBase+" end"+opts.endBase+"</div>");
      // invariants: startBase, endBase

      if(this.endBase === undefined) {
        throw new Exception("Undefined endBase");
      }
      var length = this.endBase - this.startBase;

      this.margin = Utils.orderOfMagnitude(length, -1);
      this.startBaseM = Math.max(this.startBase - this.margin, 0);
      this.endBaseM = Math.min(this.endBase + this.margin, GraphConfig.endBase);
    }
    var GP = Grapher.prototype;

    GP.destroy = function() { if(this.stage) { this.stage.destroy(); } };

    GP.draw = function(newOpts) {
      var t1 = new Date();
      _.assign(this, newOpts);
      this.m.xaxis.length = this.endBase - this.startBase;
      this.m.xaxis.scaleX = this.m.graph.w / this.m.xaxis.length;

      this.colors = GraphConfig.colorBlindFriendly ?
        npactConstants.colorBlindLineColors :
        npactConstants.lineColors;

      if(!this.stage) {
        this.stage = new K.Stage({
          container: this.$element[0],
          height: this.m.height,
          width: this.width
        });
      }
      else {
        this.stage.destroyChildren();
        this.stage.setWidth(this.width);
        this.stage.setHeight(this.m.height);
      }

      var llp = this.leftLayer(this.stage);
      var flp = this.frameLayer(this.stage);
      var glp = this.genomeLayer(this.stage);
      return $q.all([llp, flp, glp]);
        /*.then(_.bind(function() {
          $log.log("Finished draw at", this.startBase, "in", new Date() - t1);
        }, this));*/
    };

    GP._onProfilePoints = null;
    GP.clearProfilePoints = function() { this._onProfilePoints = null; };
    GP.getProfilePoints = function() {
      if(!this._onProfilePoints) {
        // start slicing the NProfile into Kinetic-compatible [x1, y1,
        // x2, y2, ...] lists, make a promise for the completed group of
        // points
        var p0 = [], p1 = [], p2 = [];
        this._onProfilePoints = NProfiler.slice({
          startBase: this.startBaseM, endBase: this.endBaseM,
          onPoint: function(coord, p0v, p1v, p2v) {
            //invert because drawing is from the top so 100% is 0 pix
            p0.push(coord); p0.push(100.0 - p0v);
            p1.push(coord); p1.push(100.0 - p1v);
            p2.push(coord); p2.push(100.0 - p2v);
          }})
          .then(function(opts) { return [p0, p1, p2]; });
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
        defaultTextOpts.text = track.name;
        defaultTextOpts.y = track.y + track.style.height/4;
        defaultTextOpts.height = track.style.height;
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
        text: GraphConfig.profileTitle()
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
      $log.debug("Building leftLayer image");
      var t1 = new Date();
      return  $q(_.bind(function(resolve) {
        var opts = {
          x: 0, y: 0,
          height: this.m.height, width: style.leftPadding,
          listening: false
        };
        var layer = new K.Layer();
        layer.hitGraphEnabled(false);
        layer.add(this.yAxisGroup());
        this.headerGroup(layer);
        layer.toImage(angular.extend(opts, {callback: resolve}));
      }, this))
        .then(function (image) {
          $log.log("Finished Building leftLayer image", new Date() - t1);
          return image;
        });
    };

    var leftLayerCache = null;  // Shared amongst all Graphers
    GP.leftLayer = function(stage) {
      var layer = new K.FastLayer();
      stage.add(layer);
      return (leftLayerCache || (leftLayerCache = this._leftLayerImage()))
        .then(function(image) {
          layer.add(new K.Image( {image: image}));
          layer.draw();
        });
    };
    GP.clearLeftLayerCache = function () { leftLayerCache = null; };

    GP.genomeLayer = function(stage) {
      //Everything in this layer is in the coordinate system of the
      //genome and will be scaled to pixels by Kinetic. Quirk: apply
      //the scale to the first group inside the layer so that the
      //clipping on the layer itself (which makes the viewport) isn't
      //scaled
      var metrics = this.m,
          scaleX = metrics.xaxis.scaleX,
          l = new K.Layer({
            x: metrics.graph.x,
            scaleX: scaleX,
            clip: {
              x: 0, width: metrics.graph.w / scaleX,
              y: 0, height: metrics.height
            }
          });
      stage.add(l);
      var dg = new K.Group({
        offsetX: this.startBase,
        draggable: true,
        dragBoundFunc: function(pos) {
          pos.y = 0; //Disallow vertical movement
          return pos;
        }
      });
      l.add(dg);
      this.dragEventsHandlers(dg, metrics);
      var dgAdd = _.bind(dg.add, dg);
      var p3 = this.xAxisGroup().then(dgAdd);
      var p2 = this.tracksGroup().then(dgAdd);
      var p1 = this.profileGroup().then(dgAdd);
      var maybeGoto = _.bind(function() {
        if(GraphConfig.gotoBase &&
           GraphConfig.gotoBase > this.startBaseM && GraphConfig.gotoBase < this.endBaseM) {
          dgAdd(new K.Line({stroke: style.profile.axis.text.fill, strokeWidth: 1,
                            strokeScaleEnabled: false,
                            points: [
                              GraphConfig.gotoBase, 0, GraphConfig.gotoBase, metrics.height]}));
        }

        return l;
      }, this);

      return $q.all([p1, p2, p3]).then(function() {
        maybeGoto();
        l.draw();
        return l;
      });
    };

    GP.dragEventsHandlers = function (dg, metrics, $scope) {
      var dragging = false,    //indicates whether *this* row is being dragged
          selectionRect = null,
          onDragMove = this.onDragMove,
          onDragEnd = this.onDragEnd,
          onRegionSelected = this.onRegionSelected,
          scaleX = metrics.xaxis.scaleX,
          toScaledDgX = function (layerX) {
            return dg.offsetX() + (layerX - metrics.graph.x) / scaleX;
          };

      this.offset = function(dx) {
        if(dragging === false) {  //don't double move the one the mouse is on.
          dg.offsetX(dg.offsetX() + dx);
          dg.getLayer().batchDraw();
        }
      };

      // need a shape that can be clicked on to allow dragging the
      // entire layer
      dg.add(new K.Rect({x: this.startBaseM, y: 0, //fill: '#BDD',
                        width: metrics.xaxis.length * 2,
                        height: metrics.height}));

      if(onDragMove) {
        dg.on('dragstart', function (evt) { dragging = 0; });
        dg.on('dragmove', function(evt) {
          var delta = evt.target.x() - dragging;
          onDragMove(-delta);
          dragging = evt.target.x();
        });
      }
      if(onDragEnd) {
        dg.on('dragend', function(evt) {
          dragging = false;
          onDragEnd(-evt.target.x());
        });
      }

      if(onRegionSelected) {
        var expandSelection = function (evt) {
          var w = toScaledDgX(evt.evt.layerX) - selectionRect.x();
          selectionRect.setWidth(w);
          dg.getLayer().batchDraw();
        };
        var finishSelection = function (evt) {
          dg.setDraggable(true);
          dg.off('mousemove', expandSelection);
          dg.off('mouseup', finishSelection);
          if(!selectionRect) return; //shouldn't happen
          var start = Math.floor(selectionRect.x()),
              end = Math.ceil(selectionRect.x() + selectionRect.width());
          onRegionSelected({
            type: 'region',
            start: Math.min(start, end),
            end: Math.max(start, end)
          });
          selectionRect.remove();
          selectionRect = null;
        };
        dg.on('mousedown', function (evt) {
          if(evt.evt.shiftKey){
            dg.setDraggable(false);
            dg.add(selectionRect = new K.Rect({
              x: toScaledDgX(evt.evt.layerX), y: 0,
              width:1, height: metrics.height,
              fill: '#bbb', stroke: '#999', opacity: 0.2
            }));
            dg.on('mousemove', expandSelection);
            dg.on('mouseup', finishSelection);
          }
        });
      }
      dg.on('mouseover', function() { document.body.style.cursor = 'pointer'; });
      dg.on('mouseout', function() { document.body.style.cursor = 'default'; });
      dg.on('dblclick', _.bind(this.onZoom, this));
    };

    GP.onZoom = function(evt) {
      var zoomOnPx = evt.evt.layerX - this.m.graph.x,
          zoomOnPct = zoomOnPx / this.m.graph.w,
          opts = {
            // these keys must match what's expected by `GraphingCalculator.zoom`
            startBase: this.startBase,
            zoomOnPct: zoomOnPct,
            basesPerGraph: GraphConfig.basesPerGraph,
            offset: GraphConfig.offset,
            zoomingOut: evt.evt.shiftKey
          };
      $log.log('Zoom event:', opts);
      this.$scope.$evalAsync(function() {
        //updates `offset`, and `basesPerGraph`
        angular.extend(GraphConfig, GraphingCalculator.zoom(opts));
      });
    };

    GP.frameLayer = function(stage) {
      // This holds the frame around the profile and the guide lines
      // for any hits tracks; it is fixed, it does not scroll left and
      // right.
      var l = new K.FastLayer();
      stage.add(l);

      var mg = this.m.graph,
          // frame around the graph
          border = new K.Rect({
            x: mg.x-1, y: mg.y,
            width: mg.w + 2, height: mg.h,
            stroke: style.profile.borderColor, strokeWidth: 1
          });
      l.add(border);
      var guides = _(this.tracks)
            .filter({type: 'hits'})
            .map(function(track) {
              return this.hitsTrackGuideLinesImage(track.style.height, mg.w)
                .then(function(img) {
                  l.add(new K.Image({
                    image: img,
                    x: mg.x, y: track.y,
                    height: track.style.height, width: mg.w
                  }));
                });
            }, this).value();
      return $q.all(guides).then(function() { l.draw(); });
    };

    var hitsTrackGuideLinesCache = {};
    GP.hitsTrackGuideLinesImage = function(height, width) {
      var key = _([height, width]).toString();
      return hitsTrackGuideLinesCache[key] ||
        (hitsTrackGuideLinesCache[key] = $q(function(resolve) {
          var layer = new K.FastLayer(),
              midY = (height / 2),
              offset = 2,  //how far off midline
              guideArrowXOffset = 8,
              guideLineOpts = style.tracks.guidelines;
          layer.add(new K.Line(angular.extend({
            points: [0, midY - offset,
                     width, midY - offset,
                     width - guideArrowXOffset, 0]

          }, guideLineOpts)));
          layer.add(new K.Line(angular.extend({
            listening: false,
            points: [0 + guideArrowXOffset, height,
                     0, midY + offset,
                     width, midY + offset]
          }, guideLineOpts)));

          layer.toImage({height: height, width: width, callback: resolve});
        }));
    };


    GP.xAxisGroup = function() {
      var xaxis = this.m.xaxis,
          stops = GraphingCalculator.stops(this.startBaseM, this.endBaseM, xaxis.length),
          g = new K.Group({
            width: xaxis.length,
            x: 0, y: this.m.xaxis.y,
            listening: false
          }),
          tickOpts = {
            stroke: style.profile.borderColor, strokeScaleEnabled: false,
            listening: false
          },
          labelOpts = _.assign({
            scaleX: 1/xaxis.scaleX,
            y: style.profile.tickLength,
            listening: false
          }, style.profile.axis.text),
          shadeOpts = {
            // a touch smaller to sit inside the border
            height: this.m.graph.h - 2,
            offsetY: this.m.graph.h - 1,
            fill: style.profile.shadeColor,
            width: stops.interval / 2,
            listening: false
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
      var m = this.m, xaxis = this.m.xaxis,
          colors = this.colors,
          g = new K.Group({
            x: 0, y: m.graph.y,
            height: m.graph.h, width: xaxis.length,
            // convert % to px
            scaleY: m.graph.h / 100,
            listening: false
          }),
          buildLine = function(points, color) {
            return new K.Line({
              points: points,
              stroke: color,
              strokeWidth: 1,
              strokeScaleEnabled: false,
              listening: false
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

    GP.tracksGroup = function() {
      var g = new K.Group();
      return $q.all(_.map(this.tracks, this.oneTrack, this))
        .then(function(list) {
          addMany(g, list);
          return g;
        });
    };

    GP.oneTrack = function(track) {
      return track.slice({
          startBase: this.startBaseM,
          endBase: this.endBaseM
        })
        .then(_.bind(function(data) {
          switch(track.type) {
          case 'neworfs':
          case 'modified':
          case 'extracts':
            return this.drawORFsTrack(track, data);
          case 'hits':
            return this.drawHitsTrack(track, data);
          default:
            throw new Error("don't know how to draw " + track);
          }
        }, this));
    };

    var arrowHeight = style.tracks.arrow.height,
        arrowHalfHeight = arrowHeight / 2,
        leftArrow = function(width, headWidth, tailWidth) {
          //Starting at the tip, clockwise
          return [ 0, arrowHalfHeight,
                   headWidth, 0,
                   width, 0,
                   width, arrowHeight,
                   headWidth, arrowHeight ];
        },
        rightArrow = function(width, headWidth, tailWidth) {
          //Starting at the tip, clockwise
          return [ width, arrowHalfHeight,
                   tailWidth, arrowHeight,
                   0, arrowHeight,
                   0, 0,
                   tailWidth, 0 ];
        };
    GP.drawORFsTrack = function(track, orfs) {
      var xaxis = this.m.xaxis,
          colors = this.colors,
          trackGroup = new K.Group({
            y: track.y,
            track: track,
            name: 'orf-track'
          }),
          arrowHeadWidth = style.tracks.arrow.width / xaxis.scaleX;

      trackGroup.add(new K.Rect({
        name: 'orf-track-hit-rectangle',
        x: 0, width: this.endBaseM,
        y: 0, height: track.style.height
      }));

      if(orfs && orfs.length > 0) {
        this.orfTrackEvents(trackGroup, track);
      }
      // Go through the list of genes in the track.
      _.forEach(orfs, function(x) {
        //x = _.clone(x);
        x.track = track.filename;
        var width = x.end - x.start,
            baseY = 0, shape,
            headWidth = Math.min(width, arrowHeadWidth),
            tailWidth = Math.max(width - headWidth, 0),
            g = new K.Group({
              draggable: true,
              orf: x,
              name: 'orf'
            });
        if(x.complement) {
          shape = leftArrow(width, headWidth, tailWidth);
          baseY = arrowHalfHeight;       // complement gets drawn lower
        }
        else {
          shape = rightArrow(width, headWidth, tailWidth);
        }
        g.add(new K.Line({
          x: x.start, y: baseY,
          points: shape, closed: true,
          stroke: track.style.light ? shadeBlend(0.7, colors[x.phase]) : colors[x.phase],
          fill: x.selected ? shadeBlend(0.75, colors[x.phase]) : null,
          strokeWidth: track.style.strokeWidth,
          strokeScaleEnabled: false
        }));
        // Only put a label in it if there is any room
        if(tailWidth && x.name) {
          // need a dummy group for clipping, `Text` doesn't
          // support clip directly
          var textBounds =
              {
                x: x.complement === 1 ? headWidth : 0, y:0,
                width: tailWidth, height: arrowHeight
              },
              lblGroup = new K.Group({
                x: x.start, y: baseY,
                clip: textBounds,
                listening: false
              }),
              lbl = new K.Text(_.assign({
                x: textBounds.x + 1,
                text: x.name,
                listening: false,
                scaleX: 1/xaxis.scaleX //undo parent scaling for readable txt
              }, style.tracks.text)),
              lblOffsetX = (lbl.getWidth() - (tailWidth * xaxis.scaleX)) / 2,
              lblOffsetY =  (lbl.getHeight() - arrowHeight) / 2;
          lblGroup.add(lbl);
          g.add(lblGroup);
          lbl.setOffset({ x: Math.min(lblOffsetX, 0), y: lblOffsetY });
        }
        trackGroup.add(g);
      }, this);

      return trackGroup;
    };

    GP.orfTrackEvents = function (trackGroup, track) {
      var self = this;
      var orfClick = function(evt) {
        var orf = getScopedAttr(evt.target, 'orf');
        //$log.log(orf);
        self.onOrfSelected({
          type: 'ORF',
          track: track,
          item: orf
        });
      };
      trackGroup.on('click', orfClick);
      var startPos =null;
      var trackLayer = null, dragLayer = null;
      var isShortDrag = function(start, end){
        var dy = Math.abs(end.y-start.y), dx = Math.abs(end.x-start.x);
        return dy<8 && dx < 25;
      };
      var dragend = function (e) {
        var pos = self.stage.getPointerPosition(),
            dropshape = trackLayer.getIntersection(pos),
            targetTrack = getScopedAttr(dropshape, 'track'),
            orf = getScopedAttr(e.target, 'orf'),
            parentTrack = _.find(GraphConfig.tracks, {filename: orf.track});
        if(isShortDrag(startPos, pos)){
          console.log('short drag, click instead');
          orfClick(e);
        }
        else if(targetTrack && track !== targetTrack
                && track.type === targetTrack.type) {
          console.log('copying', orf, e, parentTrack);
          parentTrack.remove(orf);
          parentTrack.save();
          targetTrack.add(orf);
          targetTrack.save();
        }
        $timeout(_.bind(self.draw, self), 0, false);
        return true;
      };
      var dragstart = function (e) {
        startPos = self.stage.getPointerPosition();
        trackLayer = e.target.getLayer();
        dragLayer = new Konva.Layer({
          name: 'dragLayer',
          scaleX: trackLayer.scaleX()
        });
        dragLayer.on('dragend', dragend);
        e.target.getStage().add(dragLayer);
        e.target.moveTo(dragLayer);
        trackLayer.batchDraw();
        dragLayer.batchDraw();
      };
      trackGroup.on('dragstart', dragstart);
    };

    GP.drawHitsTrack = function(track, hits) {
      var midY = (track.style.height / 2),
          offset = 2,  //how far off midline
          g = new K.Group({
            name: 'hits-track',
            x: 0, y: track.y,
            track: track}),
          colors = {
            'H': this.colors,
            //G type hits should be lighter
            'G': _.map(this.colors, function(c) { return shadeBlend(0.5, c); })}
      ;
      // draw each hit
      _.forEach(hits, function(hit, i) {
        hit.track = track.filename;
        var type = hit.name[0],  // {G,H}
            confidence = _.parseInt(hit.name[1]); // {0,1,2,3}
        var height =  16,
            c =colors[type][hit.phase],
            cfill = c,
            cstroke = hit.selected ? "#000000" : shadeBlend(0.1, c);

        if(confidence == 1){
          height = 12;
          cfill = shadeBlend(0.1, c);
        }
        if(confidence == 0){
          height = 8;
          cfill = shadeBlend(0.5, c);
        }

        //if( i<20 ) console.log("c:", c, "cfill:", cfill, 'stroke:', cstroke, 'height', height, 'c',  confidence);
        g.add(new K.Rect({
          hit: hit,
          x: hit.start,
          y:  midY,
          width:hit.end-hit.start,
          height: height,
          // set it outside the guidelines instead of ontop:
          offsetY: (hit.complement ? 0 : +height),
          fill: cfill,
          stroke: cstroke,
          strokeWidth: hit.selected ? 2 : 1
        }));
      });
      g.on('click', _.bind(function(evt) {
        var hit = evt.target.getAttrs().hit;
        GraphConfig.clearORFSelection();
        hit.selected = true;
        this.onHitSelected({
          type: 'hit',
          track: _.find(GraphConfig.tracks, {filename: hit.track}),
          item: hit
        });
      }, this));

      return g;
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

    /**
     * Capture the fully rendered graph row as a static image and
     * replace the canvas element with it.
     */
    GP.replaceWithImage = function() {
      var self = this;
      var $el = this.$element;
      return $q(function(resolve) {
        if(self.stage === null) {
          //already converted to an image
          resolve();
          return;
        }
        self.stage.toImage({
          mimeType: 'image/png',
          callback: function(image) {
            self.stage.destroy();
            self.stage = null;
            $el.append(image);
            resolve();
          }});
      });
    };

    function addMany(container, children) {
      if(children && children.length) {
        container.add.apply(container, children);
      }
    }
    function getScopedAttr(node, attr) {
      if(!node) return null;
      var val = node.getAttr(attr);
      return val !== undefined ? val : getScopedAttr(node.getParent(), attr);
    }

    function getContainingGroup(node) {
      return node.nodeType === "Group" ? node : containingGroup(node.getParent());
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
      //avoid log of 0
      return x ? 1 + Math.floor(Math.log(x) / Math.LN10) : 0;
    }

    // From http://stackoverflow.com/questions/5560248/programmatically-lighten-or-darken-a-hex-color-or-rgb-and-blend-colors
    function shadeBlend(p,c0,c1) {
      var n=p<0?p*-1:p,u=Math.round,w=parseInt;
      if(c0.length>7){
        var f=c0.split(","),t=(c1?c1:p<0?"rgb(0,0,0)":"rgb(255,255,255)").split(","),R=w(f[0].slice(4)),G=w(f[1]),B=w(f[2]);
        return "rgb("+(u((w(t[0].slice(4))-R)*n)+R)+","+(u((w(t[1])-G)*n)+G)+","+(u((w(t[2])-B)*n)+B)+")";
      }else{
        var f=w(c0.slice(1),16),t=w((c1?c1:p<0?"#000000":"#FFFFFF").slice(1),16),R1=f>>16,G1=f>>8&0x00FF,B1=f&0x0000FF;
        return "#"+(0x1000000+(u(((t>>16)-R1)*n)+R1)*0x10000+(u(((t>>8&0x00FF)-G1)*n)+G1)*0x100+(u(((t&0x0000FF)-B1)*n)+B1)).toString(16).slice(1);
      }
    }

    return Grapher;
  })

  .controller('ExtractTooltipCtrl', function($scope, extract, $log) {
    'use strict';
    $log.log('extract', extract);
    $scope.extract = extract;
  })

  .service('Tooltip', function($uibModal) {
    'use strict';
    this.show = function (extract, pageX, pageY) {
      var modalInstance = $uibModal.open({
        template: '<div npact-extract="extract"></div>',
        controller: 'ExtractTooltipCtrl',
        size: 'lg',
        resolve: {
          extract: function() {return extract;}
        }
      });
      return modalInstance.result;
    };
  });
