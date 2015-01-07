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

    var activeTracks = function(){
      return _.filter(self.tracks, {active: true}, 'active');
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
      self.tracks.push({
        text: name,
        lineType: type,
        active: true
      });
    };

    /**
     * calcuate the header information
     */
    this.headerSpec = function() {
      var offset = npactConstants.graphSpecDefaults.headerY;
      var headers = _.map(activeTracks(), function(cfg) {
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
  .directive('npactGraphConfig', function npactGraphConfig(STATIC_BASE_URL, GraphConfig) {
    return {
      restrict: 'A',
      templateUrl:STATIC_BASE_URL+'js/graphs/config.html',
      link: function($scope){
        $scope.gc = GraphConfig;
      }
    };
  });
