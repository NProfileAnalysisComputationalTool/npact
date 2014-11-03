angular.module('npact')
  .service('GraphConfig', function(TrackReader, npactConstants) {
    var self = this;
    self.tracks = {};
    self.TrackAlreadyDefined = TrackReader.TrackAlreadyDefined;

    /**
     * register a track to be displayed on the graph
     */
    this.loadTrack = function(name, type) {
      if(_.has(self.tracks, name)){ throw self.TrackAlreadyDefined; }
      self.tracks[name] = type;
    };

    /**
     * calcuate the header information
     */
    this.headerSpec = function() {
      var offset = npactConstants.graphSpecDefaults.headerY;
      var headers = _.map(self.tracks, function(lineType, text) {
        var h = npactConstants.headerSizes[lineType],
            y = offset;
        offset += h;
        return {
          text: text,
          lineType: lineType,
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
  .directive('npactGraphConfig', function npactGraphConfig(STATIC_BASE_URL, GraphDealer, $log) {
    return {
      restrict: 'A',
      templateUrl:STATIC_BASE_URL+'js/graphs/config.html',
      controller: function(){
        this.setZoom = _.debounce(GraphDealer.setZoom, 500);
        this.setColors = _.debounce(GraphDealer.setColors, 250);
      },
      controllerAs: 'ctrl',
      link:function($scope, $element, $attrs){
        $scope.graphConfig = GraphDealer.opts;
      }
    };
  });
