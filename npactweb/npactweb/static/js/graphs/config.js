angular.module('npact')
  .service('GraphConfig', function(Err, npactConstants) {
    var self = this;
    self.tracks = [];
    self.colorBlindFriendly = false;
    self.basesPerGraph = 10000;

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
      if(self.hasTrack(name)){ throw Err.TrackAlreadyDefined; }
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

    /**
     * toggle this track from showing/hiding
     */
    this.toggleTrack = function(name){
      if(!self.hasTrack(name)){ throw Err.TrackNotFound; }
      var cfg = _.find(self.tracks, {text:name});
      cfg.active = !cfg.active;
      return cfg;
    };
  })
  .directive('npactGraphConfig', function npactGraphConfig(STATIC_BASE_URL, GraphDealer, $log, GraphConfig, $rootScope, Evt) {
    return {
      restrict: 'A',
      templateUrl:STATIC_BASE_URL+'js/graphs/config.html',
      controller: function($scope){
        $scope.graphConfig = GraphDealer.opts;
        $scope.gc = GraphConfig;
      },
      controllerAs: 'ctrl'
    };
  });
