angular.module('npact')
  .constant('PUBLIC_CONFIG_KEYS',
            ['first_page_title', 'following_page_title', 'nucleotides',
             'significance', 'alternate_colors', 'startBase', 'endBase'])


  .service('GraphConfig', function(Err, npactConstants, Evt, PUBLIC_CONFIG_KEYS,
                            $location, $log) {
    var self = this;
    self.tracks = [];
    self.colorBlindFriendly = false;
    self.basesPerGraph = 10000;
    self.nucleotides = ['C', 'G'];
    self.offset = 0; // how much the graph is panned left/right

    //Get values from the querystring
    _.forEach(PUBLIC_CONFIG_KEYS, function(k) {
      if($location.search()[k]) {
        self[k] = $location.search()[k];
      }
    });

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
    this.findTrack = function(name){
      return _.find(self.tracks, {text: name});
    };

    /**
     * register a track to be displayed on the graph
     */
    this.loadTrack = function(name, type, weight) {
      weight = weight || 0;
      var existing = this.findTrack(name);
      if(existing) {
        $log.log("Updating track:", name);
        existing.lineType = type;
        existing.weight = weight;
      }
      else {
        $log.log('loading new track', name, type, weight);
        self.tracks.push({
          text: name, lineType: type,
          active: true, weight: weight
        });
      }
      self.tracks = _.sortBy(self.tracks, 'weight');
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

  /*
   * Run a VERB (e.g. extract, acgt_gamma) on server.
   * Get back a config dictionary (which automatically updates GraphConfig)
   */
  .factory('processOnServer', function(GraphConfig, KICKSTART_BASE_URL,
                                PUBLIC_CONFIG_KEYS,
                                $http, $log, $location) {
    'use strict';
    //The keys the server is going to accept from our GraphConfig
    return function(verb) {
      var url = KICKSTART_BASE_URL;
      var postconfig = angular.extend({verb: verb}, $location.search());
      $log.log('Going to request', verb, url, postconfig);
      PUBLIC_CONFIG_KEYS.forEach(function(k) {
        if(GraphConfig[k])
          postconfig[k] = GraphConfig[k];
      });

      return $http.get(url, { cache: true, params: postconfig })
        .then(function(res) {
          $log.log('Server process successful:', res.data);
          angular.extend(GraphConfig, res.data);
          return res.data;
        });
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
  .controller('npactGraphConfigCtrl', function($scope, $window, $location, $log,
                                        GraphConfig, PredictionManager,
                                        PUBLIC_CONFIG_KEYS) {
    'use strict';
    $scope.gc = GraphConfig;
    $scope.$watch('gc.significance', PredictionManager.onSignificanceChange);

    //  If any of the GraphConfig values change update the querystring
    var gcpubkeys = PUBLIC_CONFIG_KEYS.map(function(k) { return 'gc.' + k; });
    $scope.$watchGroup(gcpubkeys, function(newVals) {
      $location.search(_.object(PUBLIC_CONFIG_KEYS, newVals));
    });
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
