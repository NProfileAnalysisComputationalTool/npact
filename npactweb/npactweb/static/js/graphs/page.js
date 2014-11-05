angular.module('npact')
  .directive('npactGraphPage', function(STATIC_BASE_URL) {
    'use strict';
    return {
      restrict: 'A',
      templateUrl: STATIC_BASE_URL + 'js/graphs/page.html',
      controller: 'npactGraphPageCtrl',
      controllerAs: 'ctrl'
    };
  })

  .controller('npactGraphPageCtrl', function($scope, Fetcher, npactConstants, $q, $log, StatusPoller, FETCH_URL, $window, $element, GraphConfig, Evt, TrackReader, GraphingCalculator, ProfileReader, Utils) {
    'use strict';
    var self = this,
        visibleGraphs = 5,
        graphSpecs = [],
        // helper functions
        getWidth = function(){ return $element.width(); },
        getGraphConfig = function() { return GraphConfig; },
        addTrack = function(key, name, data) {
          return TrackReader.load(name, data)
            .then(function() { GraphConfig.loadTrack(name, key); });
        },
        addExtract = function(name, data) {
          return addTrack('extracts', name, data);
        },
        addInputFileCds = function(config) {
          return Fetcher.inputFileCds(config)
            .then(function(data) { return addExtract('Input file CDS', data); });
        },
        addNewCds = function(config) {
          return Fetcher.fetchFile(config['File_of_new_CDSs'])
            .then(function(data) {
              return addExtract('Newly Identified ORFs', data);
            });
        },
        addHits = function(config) {
          return Fetcher.fetchFile(config['File_of_G+C_coding_potential_regions'])
            .then(function(data) { addTrack('hits', 'Hits', data); });
        },
        addProfile = function(config) {
          return Fetcher.nprofile(config)
            .then(ProfileReader.load)
            .then(function(summary) {
              // find a sensible zoom level
              var basesPerGraph = summary.length / visibleGraphs;
              // if we're really short, reset out bases per graph
              if (basesPerGraph < GraphConfig.basesPerGraph) {
                GraphConfig.basesPerGraph = Utils.orderOfMagnitude(basesPerGraph);
              }
              GraphConfig.profileSummary = summary;
            });
        }
    ;

    $scope.miscFiles = [];
    $scope.graphHeight = npactConstants.graphSpecDefaults.height;
    $scope.status = 'Initializing';
    $scope.ready = false;
    $scope.FETCH_URL = FETCH_URL;

    $scope.$watch(getWidth, function(newValue, oldValue){
      if (newValue > 0){ GraphConfig.width = newValue; }
    });

    $scope.$watch(getGraphConfig, function(newValue, oldValue){
      var cmd = newValue.refreshCommand(oldValue);
      $log.log('graph config changed:', cmd);
      switch(cmd){
      case Evt.REBUILD:
        graphSpecs = ProfileReader.partition(GraphConfig);
        $scope.graphSpecs = _.take(graphSpecs, visibleGraphs);
        break;
      case Evt.REDRAW:
        $scope.$broadcast(cmd);
        break;
      }
    }, true); // deep-equality

    $scope.$on(Evt.PAN, function(evt, opts) {
      var offset = Math.floor(opts.newStartBase - opts.oldStartBase);
      GraphConfig.offset += offset;
      // TODO: Event originated from outside ng, but why doesn't
      // `$watch` pick up the `offset` change?
      $scope.$apply();
    });

    $scope.$on(Evt.ZOOM, function(evt, opts) {
      var res = GraphingCalculator.zoom(angular.extend({}, opts, GraphConfig));
      GraphConfig.offset = res.offset;
      GraphConfig.basesPerGraph = res.basesPerGraph;
      // TODO: Event originated from outside ng, but why doesn't
      // `$watch` pick up the `offset` change?
      $scope.$apply();
    });

    // check for scope changes on resize
    angular.element($window).bind('resize', function () { $scope.$apply(); });

    /**
     * add more visible entries to $scope
     */
    self.addMore = function(){
      if($scope.graphSpecs){
        $log.log('scrolling down via infinite scroller');
        Utils.extendByPage(graphSpecs, $scope.graphSpecs, visibleGraphs);
      }
    };

    // start it up
    Fetcher.kickstart()
      .then(function(config) {
        $log.log('Kickstart successful:', config);
        $scope.status = 'Running';
        $scope.config = config;
        // got config, request the first round of results
        $scope.title = config.first_page_title;
        var nprofile = addProfile(config);

        // Non-gbk files don't have CDSs we can extract.
        var inputFileCds = config.isgbk && addInputFileCds(config);

        var extraFileList = Fetcher.acgtGammaFileList(config)
              .then(function(fileList) {
                $scope.miscFiles.push.apply($scope.miscFiles, fileList);
                addNewCds(config);
                addHits(config);
              });
        $q.all([nprofile, inputFileCds, extraFileList])
          .then(function() { delete $scope.status; })
          .catch(function(err) { $scope.status = err; });

        StatusPoller.start(config['pdf_filename'])
          .then(function(pdfFilename) { $scope.miscFiles.push(pdfFilename); });
      });
  })

  .service('Fetcher', function(StatusPoller, $http, FETCH_URL, ACGT_GAMMA_FILE_LIST_URL, KICKSTART_BASE_URL, $window) {
    'use strict';
    var self = this;
    /**
     * download contents from any url
     */
    self.rawFile = function(url) {
      return $http.get(url).then(function(res) { return res.data; });
    };

    /**
     * download contents from a "fetch" path
     */
    self.fetchFile = function(path){ return self.rawFile(FETCH_URL + path); };


    /**
     * poll the server for when `path` is ready, then fetch it
     */
    self.pollThenFetch = function(path) {
      return StatusPoller.start(path).then(self.fetchFile);
    };

    self.kickstart = function(){
      return self.rawFile(KICKSTART_BASE_URL + $window.location.search);
    };

    self.nprofile = function(config) {
      return self.pollThenFetch(config['nprofileData']);
    };

    self.inputFileCds = function(config) {
      return self.pollThenFetch(config['File_of_published_accepted_CDSs']);
    };

    self.acgtGammaFileList = function(config) {
      return StatusPoller.start(config['acgt_gamma_output'])
        .then(function(path) {
          return self.rawFile(ACGT_GAMMA_FILE_LIST_URL + path);
        });
    };
  })


  .service('StatusPoller', function(STATUS_BASE_URL, $q, $http, $timeout, $log) {
    'use strict';
    var POLLTIME = 2500;

    function poller(tid, deferred) {
      // remember our arguments
      var pollAgain = _.partial(poller, tid, deferred);

      $http.get(STATUS_BASE_URL + tid)
        .then(function(res) {
          if(res.data.ready) { deferred.resolve(tid); }
          else { $timeout(pollAgain, POLLTIME); }
        })
        .catch(function(err) {
          $log.error('Error while fetching tid: ', tid, err);
          deferred.reject(err.data.message);
        });

      return deferred.promise;
    }

    this.start = function(tid) {
      if(!tid || tid.length === 0){
        return $q.reject(new Error('Invalid task id'));
      }

      return poller(tid, $q.defer());
    };
  });
