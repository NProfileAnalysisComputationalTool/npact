angular.module('npact')
  .constant('PUBLIC_CONFIG_KEYS',
            ['first_page_title', 'following_page_title', 'nucleotides',
             'significance', 'startBase', 'endBase', 'basesPerGraph','offset'])


  .service('GraphConfig', function(Err, npactConstants, Evt, PUBLIC_CONFIG_KEYS, Track,
                            $location, $log, $rootScope) {
    var self = this;
    self.tracks = [];
    self.colorBlindFriendly = false;
    self.basesPerGraph = 10000;
    self.nucleotides = ['C', 'G'];
    self.offset = 0; // how much the graph is panned left/right

    //Get values from the querystring during intialization
    _.forEach(PUBLIC_CONFIG_KEYS, function(k) {
      var v = $location.search()[k];
      if(v) {
        try {
          self[k] = Number(v);
        }
        catch(e) {
          self[k] = v;
        }
      }
    });

    /**
     * what's the right title for the current nucleotides?
     */
    self.profileTitle = function() {
      return self.nucleotides ? '% ' + self.nucleotides.join('') : null;
    };

    self.activeTracks = function(){
      return _.filter(self.tracks, 'active');
    };

    /**
     * do we have a track with a given name?
     */
    this.findTrack = function(name){
      return _.find(self.tracks, {name: name});
    };

    /**
     * register a track to be displayed on the graph
     */
    this.loadTrack = function(track) {
      $log.log('loading new track', track);
      _.remove(self.tracks, {name: track.name});  //mutates
      self.tracks.push(track);
      self.tracks = _.sortBy(self.tracks, 'weight');
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
      _.forEach(PUBLIC_CONFIG_KEYS, function(k) {
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
    var gcpubkeys = _.map(PUBLIC_CONFIG_KEYS, function(k) { return 'gc.' + k; });
    $scope.$watchGroup(gcpubkeys, function(newVals) {
      $location.search(_.object(PUBLIC_CONFIG_KEYS, newVals));
    });
  })
;
