angular.module('npact')
  
  .constant('PUBLIC_CONFIG_KEYS',
            ['first_page_title', 'following_page_title', 'nucleotides',
             'significance', 'startBase', 'endBase', 'basesPerGraph',
             'offset', 'mycoplasma'])


  .service('GraphConfig', function(Err, npactConstants, Evt, PUBLIC_CONFIG_KEYS, Track,
                            $location, $log, $rootScope, $cookies) {
    var self = this;
    self.cookieBools = ["colorBlindFriendly"];
    self.colorBlindFriendly = false;
    self.cookieInit = function() {
      _.each(self.cookieBools, function(v){
        var cv = $cookies.get(v);
        if( cv === "true" ) self[v] = true;
        else if(cv === "false") self[v] = false;
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
        $log.log(k, v);
        if(v === true || v === false) {
          self[k] = v;
        }
        else if(!isNaN(v)) {
          self[k] = Number(v);
        }
        else {
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
  .factory('processOnServer', function(GraphConfig, KICKSTART_URL,
                                PUBLIC_CONFIG_KEYS, Utils,
                                $http, $log, $location) {
    'use strict';
    //The keys the server is going to accept from our GraphConfig
    return function(verb) {
      var url = KICKSTART_URL;
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
                                        $timeout,
                                        GraphConfig, PredictionManager,
                                        PUBLIC_CONFIG_KEYS) {
    'use strict';
    $scope.gc = GraphConfig;
    //  If any of the GraphConfig values change update the querystring
    var gcpubkeys = _.map(PUBLIC_CONFIG_KEYS, function(k) { return 'gc.' + k; });
    $scope.$watchGroup(gcpubkeys, function(newVals) {
      $location.search(_.object(PUBLIC_CONFIG_KEYS, newVals));
    });
    $scope.$watch('gc.mycoplasma', function() {
      $timeout(PredictionManager.start, 50);
    });
  })

  .directive('npactOrfFinder', function(GraphConfig, MessageBus, $q, $log, $timeout) {
    'use strict';
    return {
      restrict: 'A',
      require: 'ngModel',
      link: function($scope, element, attrs, ngModel) {
        function resultInRange(result) {
          return result.start >= GraphConfig.startBase &&
            result.end <= GraphConfig.endBase;
        }

        function doSearch(val) {
          $log.log("Searching for:", val);
          delete $scope.results;
          delete $scope.resultsIndex;
          if(!val) return $q.when(true);
          return $q
            .all(_.map(GraphConfig.activeTracks(), function(t) { return t.findByName(val); }))
            .then(function(values) {
              return $timeout(function() {  // let other code run
                $scope.results = _(values).flatten().filter(resultInRange).sortBy('start').value();
                $log.log("Total results", $scope.results.length);
                if($scope.results.length > 0) {
                  $scope.resultsIndex = 0;
                  GraphConfig.gotoBase = $scope.results[$scope.resultsIndex].start;
                  return true;
                }
                else {
                  MessageBus.warning("No matches found", 1000);
                  return $q.reject("No matches found for: '" + val + "'");
                }
              });
            });
        }
        ngModel.$asyncValidators.matchingOrfs = doSearch;
        $scope.$watchGroup(['gc.startBase', 'gc.endBase'], function() {
          if(GraphConfig.findORF) {
            doSearch(GraphConfig.findORF);
          }
        });
        $scope.$watch('resultsIndex', function(idx) {
          if(!$scope.results || $scope.results.length === 0) return;
          if(idx < 0) { $scope.resultsIndex = 0; }
          if(idx >= $scope.results.length) { $scope.resultsIndex = $scope.results.length - 1; }
          GraphConfig.gotoBase = $scope.results[$scope.resultsIndex].start;
        });
      }
    };
  })

  .controller('TooltipCtrl', function($scope, $sce) {
  })
;
