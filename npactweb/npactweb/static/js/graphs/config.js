angular.module('npact')
  .service('GraphConfig', function(Err, npactConstants, Evt, $log) {
    var self = this;
    self.tracks = [];
    self.colorBlindFriendly = false;
    self.basesPerGraph = 10000;
    self.nucleotides = ['C', 'G'];
    self.offset = 0; // how much the graph is panned left/right

    /**
     * what's the right title for the current nucleotides?
     */
    self.profileTitle = function() {
      return '% ' + self.nucleotides.join('');
    };

    self.activeTracks = function(){
      return _.filter(self.tracks, 'active');
    };

    /**
     * do we have a track with a given name?
     */
    this.hasTrack = function(name){
      return _.some(self.tracks, {text: name});
    };

    /**
     * register a track to be displayed on the graph
     */
    this.loadTrack = function(name, type) {
      //TODO: This needs to replace existing track on name match
      if(self.hasTrack(name)){ throw new Err.TrackAlreadyDefined(); }
      $log.log('loading track', name, type);
      var newTrack = { text: name, lineType: type, active: true },
          isHit = function(t) { return t.lineType === 'hits'; },
          hitIdx = self.tracks.indexOf(_.find(self.tracks, isHit))
      ;

      // ensure any hits are last, which will display them close to
      // the graph
      if(!isHit(newTrack) && hitIdx !== -1){
        self.tracks.splice(hitIdx, 0, newTrack);
      }else{
        self.tracks.push(newTrack);
      }
    };

    /**
     * calcuate the header information
     */
    this.partition = function() {
      var idx = 0,
          offset = self.offset || 0,
          startBase = Math.max(self.startBase + offset, 0),
          endBase = self.endBase,
          bpg = self.basesPerGraph;
      if(!(startBase >= 0 && endBase >= 0 && bpg >=0)) { return []; }
      var g = new Array(Math.ceil((endBase - startBase) / bpg));
      while(startBase < endBase) {
        g[idx] = {
          startBase: startBase,
          endBase: startBase + bpg - 1
        };
        idx++;
        startBase = startBase + bpg;
      }
      return g;
    };
  })
  .directive('npactGraphConfig', function npactGraphConfig(STATIC_BASE_URL) {
    'use strict';
    return {
      restrict: 'A',
      templateUrl: STATIC_BASE_URL + 'js/graphs/config.html',
      controller: 'npactGraphConfigCtrl as gcctrl'
    };
  })
  .controller('npactGraphConfigCtrl', function($scope, GraphConfig, PredictionManager) {
    'use strict';
    $scope.gc = GraphConfig;
    $scope.$watch('gc.significance', PredictionManager.onSignificanceChange);
  })

  .factory('headerSpecCalc', function(npactConstants) {
    'use strict';
    return function(activeTracks) {
      var offset = npactConstants.graphSpecDefaults.headerY;
      var headers = _.map(activeTracks, function(cfg) {
        var h = npactConstants.headerSizes[cfg.lineType],
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
        headers: headers,
        headerY: offset
      };
    };
  })
;
