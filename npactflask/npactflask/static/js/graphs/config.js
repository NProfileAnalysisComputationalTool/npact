angular.module('npact')
  .constant('PUBLIC_CONFIG_KEYS',
            ['first_page_title', 'following_page_title', 'nucleotides',
             'significance', 'startBase', 'endBase', 'basesPerGraph',
             'offset', 'mycoplasma', 'zoomTrack', 'zoomIdx'])

  .service('GraphConfig', function(Err, npactConstants, Evt, PUBLIC_CONFIG_KEYS,
                            $location, $log, $rootScope, $cookies, $window,
                            STATIC_BASE_URL, CodonFinder) {
    'use strict';
    var self = this;
    self.graphMargin = 198;// multiple of 3 ;)
    $window.GraphConfig = this;
    self.baseUrl = STATIC_BASE_URL;
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
    self.tracks = null;
    self.cookieInit();

    self.basesPerGraph = 10000;
    self.nucleotides = ['C', 'G'];
    self.offset = 0; // how much the graph is panned left/right

    var inputConfig = {};
    //Get values from the querystring during intialization
    _.forEach(PUBLIC_CONFIG_KEYS, function(k) {
      var v = $location.search()[k];
      if(v) {
        if(v === true || v === false) {
          inputConfig[k] = v;
        }
        else if (v === "true") {
          inputConfig[k] = true;
        }
        else if (v === "false") {
          inputConfig[k] = false;
        }
        else if(!isNaN(v)) {
          inputConfig[k] = Number(v);
        }
        else {
          inputConfig[k] = v;
        }
      }
      //  Watch for the value changing later
      $rootScope.$watch(function () { return GraphConfig[k]; },
                        function (v) { $location.search(k, v); });
    });
    inputConfig.trackPaths = $location.search().trackPaths;
    //  Watch for the track paths changing changing later
    $rootScope.$watch(
      function () {
        return GraphConfig.tracks ? _.map(GraphConfig.tracks, "filename").join(",")
          : GraphConfig.trackPaths;
      },
      function (v) {
        $location.search('trackPaths', v);
        GraphConfig.trackPaths = v;
      });
    $log.debug("Finished reading config from querystring:", inputConfig);
    _.assign(self, inputConfig);

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
      return _.find(self.tracks,function(tr) {
        return tr.name == name || tr.filename == name;
      });
    };
    self.CodonFinder = CodonFinder;


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
                                        $timeout, GraphConfig, PUBLIC_CONFIG_KEYS) {
    'use strict';
    $scope.gc = GraphConfig;
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
          delete $scope.results;
          delete $scope.resultsIndex;
          if(!val) return $q.when(true);
          $log.log("Searching for:", val);
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
                  MessageBus.warning("No matches found");
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
