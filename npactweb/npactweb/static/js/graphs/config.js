angular.module('npact')
  .constant('PUBLIC_CONFIG_KEYS',
            ['first_page_title', 'following_page_title', 'nucleotides',
             'significance', 'startBase', 'endBase', 'basesPerGraph',
             'offset','mycoplasma'])


  .service('GraphConfig', function(Err, npactConstants, Evt, PUBLIC_CONFIG_KEYS, Track,
                            $location, $log, $rootScope, $cookies) {
    var self = this;
    self.cookieBools = ["colorBlindFriendly"];
    self.colorBlindFriendly = false;
    self.cookieInit = function() {
      _.each(self.cookieBools, function(v){
        var cv = $cookies.get(v);
        if( cv == "true" ) self[v] = true;
        else if(cv == "false") self[v] = false;
      });
    };
    self.cookiePersist  = function(){
      _.each(self.cookieBools, function(v){
        if(self[v]) $cookies.put(v, "true");
        else $cookies.put(v, "false");
      });
    };
    self.tracks = [];
    self.cookieInit();

    self.basesPerGraph = 10000;
    self.nucleotides = ['C', 'G'];
    self.offset = 0; // how much the graph is panned left/right

    //Get values from the querystring during intialization
    _.forEach(PUBLIC_CONFIG_KEYS, function(k) {
      var v = $location.search()[k];
      if(v) {
        self[k] = !isNaN(Number(v)) ? Number(v) : v ;
      }
    });

    /**
     * what's the right title for the current nucleotides?
     */
    self.profileTitle = function() {
      return self.nucleotides ? '% ' + self.nucleotides.join('') : null;
    };

    self.activeTracks = function() {
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
                                PUBLIC_CONFIG_KEYS, Utils,
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
      postconfig['x-tics'] = Utils.orderOfMagnitude(GraphConfig.basesPerGraph, -1);

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

  .directive('npactOrfFinder', function(GraphConfig, MessageBus, $q, $log) {
    'use strict';
    return {
      restrict: 'A',
      require: 'ngModel',
      link: function($scope, element, attrs, ngModel) {
        function doSearch(val) {
          $log.log("Searching for:", val);
          $scope.hitsIndex = null;
          if(!val) return $q.when(true);
          function hitInRange(hit) {
            return hit.start >= GraphConfig.startBase &&
              hit.end <= GraphConfig.endBase;
          }
          return $q
            .all(_.map(GraphConfig.activeTracks(), function(t) { return t.findByName(val); }))
            .then(function(values) {
              $scope.hits = _(values).flatten().filter(hitInRange).sortBy('start').value();
              $log.log("Total results", $scope.hits.length);
              if($scope.hits.length > 0) {
                $scope.hitsIndex = 0;
                return true;
              }
              else {
                MessageBus.warning("No matches found", 1000);
                return $q.reject("No matches found for: '" + val + "'");
              }
            });
        }
        $scope.$watchGroup(['gc.startBase', 'gc.endBase'], function() {
          if(GraphConfig.findORF) {
            doSearch(GraphConfig.findORF);
          }
        });
        ngModel.$asyncValidators.matchingOrfs = doSearch;
        $scope.$watch('hitsIndex', function(idx) {
          if(!$scope.hits || $scope.hits.length === 0) return;
          if(idx < 0) { $scope.hitsIndex = 0; }
          if(idx >= $scope.hits.length) { $scope.hitsIndex = $scope.hits.length - 1; }
          GraphConfig.gotoBase = $scope.hits[$scope.hitsIndex].start;
        });

      }
    };
  })
;
