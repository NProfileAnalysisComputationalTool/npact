angular.module('npact')
  .constant('PUBLIC_CONFIG_KEYS',
            ['first_page_title', 'following_page_title', 'nucleotides',
             'significance', 'startBase', 'endBase', 'basesPerGraph',
             'offset', 'mycoplasma', 'zoomTrack', 'zoomIdx'])

  .service('GraphConfig', function(Err, npactConstants, Evt, PUBLIC_CONFIG_KEYS,
                            $location, $log, $rootScope, $cookies, $window,
                            STATIC_BASE_URL) {
    'use strict';
    var self = this;
    var GraphConfig = this;

    GraphConfig.graphMargin = 198;// multiple of 3 ;)
    $window.GraphConfig = this;
    GraphConfig.baseUrl = STATIC_BASE_URL;
    GraphConfig.cookieBools = ["colorBlindFriendly"];
    GraphConfig.colorBlindFriendly = false;
    GraphConfig.clearORFSelection =function() {
      _.each(GraphConfig.tracks,function(t) {
        _.each(t.data,function(orf,k) { orf.selected = false; });
      });
    };
    GraphConfig.cookieInit = function() {
      _.each(GraphConfig.cookieBools, function(v){
        var cv = $cookies.get(v);
        if( cv === "true" ) GraphConfig[v] = true;
        else if(cv === "false") GraphConfig[v] = false;
      });
    };
    GraphConfig.cookiePersist  = function(){
      _.each(GraphConfig.cookieBools, function(v){
        if(GraphConfig[v]) $cookies.put(v, "true");
        else $cookies.put(v, "false");
      });
    };
    GraphConfig.tracks = null;
    GraphConfig.cookieInit();

    GraphConfig.basesPerGraph = 10000;
    GraphConfig.nucleotides = ['C', 'G'];
    GraphConfig.offset = 0; // how much the graph is panned left/right

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
    _.assign(GraphConfig, inputConfig);

    /**
     * what's the right title for the current nucleotides?
     */
    GraphConfig.profileTitle = function() {
      return GraphConfig.nucleotides ? '% ' + GraphConfig.nucleotides.join('') : null;
    };

    GraphConfig.activeTracks = function() {
      return _.filter(GraphConfig.tracks, 'active');
    };

    /**
     * do we have a track with a given name?
     */
    GraphConfig.findTrack = function(name){
      return _.find(GraphConfig.tracks,function(tr) {
        return tr.name == name || tr.filename == name;
      });
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
            .all(_.map(GraphConfig.activeTracks(), function(t) {
              var orf = t.findByName(val);

              return orf;
            }))
            .then(function(values) {
              values = _(values).flatten();
              return $timeout(function() {  // let other code run
                //$log.log('Found orfs', values);
                GraphConfig.clearORFSelection();
                _.each(values, function(orf) {
                  orf.selected = true;
                });

                $scope.results = _(values).filter(resultInRange).sortBy('start').value();
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
