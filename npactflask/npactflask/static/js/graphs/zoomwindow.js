angular.module('npact')

// {type: 'region', start: 0, end: 10}
// {type: 'hit}

  .service('FocusData', function ($log, $location, GraphConfig) {
    /**
     * focusData interface:
     *
     *  type       : "(hit|orf|region)"
     *  track      : Reference to a Track object
     *  item       : Data Object inside a track
     *  start      : start base of zoom region
     *  end        : end base of zoom region
     *  complement : (0|1), indicates strand
     *  phase      : (0|1|2)
     *  name       : name of item
     *
     */
    this.deserialize = function () {
      if(!(GraphConfig.zoomTrack && GraphConfig.zoomIdx >= 0))
        return null;

      var track = GraphConfig.findTrack(GraphConfig.zoomTrack);
      var cds = track.data[GraphConfig.zoomIdx];
      //$log.log('loading', GraphConfig.zoomTrack, GraphConfig.zoomIdx,track, cds);
      var focusData = {
        type: 'orf',
        track: track,
        start: cds.start,
        end: cds.end,
        complement: cds.complement,
        phase: cds.phase,
        name: cds.name,
        item: cds
      };
      return focusData;
    };
    this.serialize = function (focusData) {
      $log.log('Serializing: ', focusData, focusData.track.filename, focusData.item.cdsidx);
      if(!focusData) return $location.absUrl();
      if(focusData.track)
        GraphConfig.zoomTrack = focusData.track.filename;
      if(focusData.item && focusData.item.cdsidx>=0)
        GraphConfig.zoomIdx = focusData.item.cdsidx;
      return $location.absUrl();
    };
    this.clearQuery = function () {
      //$log.log("Clearing FocusData from querystring");
      GraphConfig.zoomTrack=null;
      GraphConfig.zoomIdx = null;
    };
  })
  .service('ZoomWindowHandler', function ($log, $uibModal, FocusData, STATIC_BASE_URL, Evt) {

    this.register = function ($scope) {
      $scope.$on('ORF-selected',
                 _.bind(function (evt, data) { this.popup(data,null, $scope); }, this));
      $scope.$on('hit-selected',
                 _.bind(function (evt, data) { this.popup(data,null, $scope); }, this));
      $scope.$on('region-selected',
                 _.bind(function (evt, data) { this.popup(data,null, $scope); }, this));
    };

    this.maybePopup = function ($scope) {
      var fd = FocusData.deserialize();
      //console.log('maybe popup', fd);
      if(fd) return this.popup(fd,null, $scope);
      return null;
    };
    this.popup = function (focusData, modalopts, $scope) {
      $log.log("popping!", focusData, modalopts);
      FocusData.serialize(focusData);
      if(GraphConfig.zoomwindow){
        GraphConfig.zoomwindow.$scope.data = focusData;
        $scope.$broadcast(Evt.REDRAW);
        return GraphConfig.zoomwindow.result;
      }

      var child = $scope.$new();
      child.data = focusData;
      var modalDefaults = {
        templateUrl: STATIC_BASE_URL + 'js/graphs/zoomwindow.html',
        controller: 'ZoomWindowCtrl',
        controllerAs: 'zw',
        bindToController: true,
        scope: child,
        animation: true,
        size: 'lg'
      };
      GraphConfig.zoomwindow = $uibModal.open(_.assign(modalDefaults, modalopts));
      GraphConfig.zoomwindow.result.catch(function() {
        // when exiting the window, data may have changed that requires a redraw
        $scope.$broadcast(Evt.REDRAW);
      }).finally( function() {
        FocusData.clearQuery();
        GraphConfig.zoomwindow = null;
      });
      return GraphConfig.zoomwindow.result;
    };
  })
  .controller('ZoomWindowCtrl', function ($scope, $log, FocusData, getPhase,
                                   GraphConfig, Utils, $location, $uibModalInstance,
                                   $window ) {
    //$log.log('ZoomWindowCtrl', focusData);
    GraphConfig.zoomwindow.$scope = $scope;
    var type = $scope.data.type;
    var self = this;
    this.cancel = function() {
      console.log('Cancelling track edit');
      $window.location.reload();
    };
    this.save = function(track) {
      console.log('saving track', track);
      $uibModalInstance.dismiss();
      return track.save();
    };
    this.delete = function(orf, track) {
      console.log('deleting orf', orf, track);
      self.track.remove(orf);
      return self.save(self.track);
    };
    if(type === 'region') {
      $scope.data = {
        name: "New",
        start: Utils.floor3($scope.data.start),
        end: Utils.ceil3($scope.data.end) -1, //the end is an inclusive range
        complement: 0
      };
    }
    else {
      this.item = $scope.data.item;
    }
    this.track = $scope.data.track || _.first(GraphConfig.activeTracks);

    $scope._setGraphBounds = function() {
      var margin = Utils.ceil3(Utils.orderOfMagnitude(
        $scope.data.item.end - $scope.data.item.start, -1));
      $scope.startBase = Math.max($scope.data.item.start - margin, 0);
      $scope.endBase = Math.min($scope.data.item.end + margin, GraphConfig.endBase);
    };
    $scope._setGraphBounds();

    $scope.$watch(_.partial(getPhase, $scope.data.item),
                  _.bind(function (phase) { $scope.data.item.phase = phase; }, this));

    //$log.log("Focusing on", focusData);
    $scope.$on('offset', _.bind(function ($evt, dx) {
      $evt.stopPropagation();
      $log.log("Updating offsets", dx);
      dx = Utils.round3(dx);
      $scope.startBase = $scope.startBase + dx;
      $scope.endBase = $scope.endBase + dx;
      $scope.$apply();
    }, this));
    var redraw = function () {
      $scope._setGraphBounds();
      $scope.$broadcast('redraw');
    };
    $scope.$watchCollection("data.item", redraw);
    $scope.$watch('zw.track', redraw);
    //Keep phase up to date with the end

    this.buildPermalink = function () {
      //TODO: Params to this function
      return FocusData.serialize(focusData);
    };

  })

  .directive('npactSingleGraph', function($log, $timeout,
                                   GraphConfig, GraphingCalculator, Grapher, Evt) {
    'use strict';
    return {
      restrict: 'A',
      scope: {
        startBase: '=',
        endBase: '='
      },
      link: function($scope, $element, $attrs) {
        var g = null;
        var draw = function() {
          $log.debug("Redrawing single graph");
          var opts = {
            width: $element.width(),
            m: null,
            tracks: GraphConfig.activeTracks(),
            startBase: $scope.startBase,
            endBase: $scope.endBase,
            onDragEnd: function (dx) { $scope.$emit('offset', dx); },
            onRegionSelected: function (data) { $scope.$emit('region-selected', data); },
            onOrfSelected: function (data) { $scope.$emit('ORF-selected', data); },
            onHitSelected: function (data) { $scope.$emit('hit-selected', data); }
          };
          opts.m = GraphingCalculator.chart(opts);
          $scope.graphHeight = opts.m.height;
          $log.debug("Redrawing with params:", opts);
          if(g) { g.destroy(); }
          g = new Grapher($element, $scope, opts);
          return g.draw();
        };
        var schedule = _.debounce(function () {$timeout(draw);},400);
        $scope.$on('redraw', schedule);
        $scope.$watchGroup(['startBase', 'endBase'], schedule);
        $scope.$watch(function() {return $element.width();}, schedule);
      }
    };
  })

  .directive('npactProteinTranslation', function ($log, TranslatePath) {
    'use strict';
    return {
      restrict: 'E',
      scope: {
        item: '=',
        start: '=',
        end: '=',
        complement: '='
      },
      link: function ($scope, $element) {
        $scope.$watchGroup(
          ['start', 'end', 'complement'],
          function () {
            $element.text('');
            TranslatePath($scope.start, $scope.end, $scope.complement)
              .then(function (data) {
                if($scope.item){
                  $scope.item.ddnaP = data.trans;
                  $scope.item.ddna = data.seq;
                }
                $element.text(data.trans);
              });
          }
        );
      }
    };
  })
;
