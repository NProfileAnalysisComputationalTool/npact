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
        getWidth = function(){ return $element.width(); },
        getGraphConfig = function() { return GraphConfig; },
        addTrack = function(key, name, data) {
          return TrackReader.load(name, data)
            .then(function() { GraphConfig.loadTrack(name, key); });
        },
        addExtract = function(name, data) {
          return addTrack('extracts', name, data);
        },
        addHits = function(name, data) {
          return addTrack('hits', name, data);
        },
        configureProfile = function(summary) {
          // find a sensible zoom level
          var basesPerGraph = summary.length / visibleGraphs;
          // if we're really short, reset out bases per graph
          if (basesPerGraph < GraphConfig.basesPerGraph) {
            GraphConfig.basesPerGraph = Utils.orderOfMagnitude(basesPerGraph);
          }
          GraphConfig.profileSummary = summary;
        }
    ;

    $scope.miscFiles = [];
    $scope.graphHeight = npactConstants.graphSpecDefaults.height;
    $scope.status = 'Initializing';
    $scope.ready = false;
    $scope.FETCH_URL = FETCH_URL;

    self.addMore = _.debounce(function(){
      if($scope.graphSpecs){
        $log.log('scrolling down via infinite scroller');
        Utils.extendByPage(graphSpecs, $scope.graphSpecs, visibleGraphs);
      }
    }, 250);

    $scope.$watch(getWidth, function(newValue, oldValue){
      if (newValue > 0){
        $log.log('width changed from', oldValue, '->', newValue);
        GraphConfig.width = newValue;
      }
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

    // start it up
    Fetcher.kickstart()
      .then(function(config) {
        $log.log('Kickstart successful:', config);
        $scope.status = 'Running';
        $scope.config = config;
        // got config, request the first round of results
        $scope.title = config.first_page_title;
        var nprofile = Fetcher.nprofile(config)
              .then(ProfileReader.load)
              .then(configureProfile);

        // Non-gbk files don't have CDSs we can extract.
        var inputFileCds = config.isgbk && Fetcher.inputFileCds(config)
              .then(function(data) { return addExtract('Input file CDS', data); });

        var extraFileList = Fetcher.acgtGammaFileList(config)
              .then(function(fileList) {
                $scope.miscFiles.push.apply($scope.miscFiles, fileList);
                Fetcher.fetchFile(config['File_of_new_CDSs'])
                  .then(function(data) {
                    return addExtract('Newly Identified ORFs', data);
                  });
                Fetcher.fetchFile(config['File_of_G+C_coding_potential_regions'])
                  .then(function(data) { addHits('Hits', data); });
              });
        $q.all([nprofile, inputFileCds, extraFileList])
          .then(function() {
            delete $scope.status;
          })
          .catch(function(err) {
            $scope.status = err;
          });

        StatusPoller.start(config['pdf_filename'])
          .then(function(pdfFilename) {
            $scope.miscFiles.push(pdfFilename);
          });
      });
  })

  .service('Fetcher', function(StatusPoller, $http, FETCH_URL, ACGT_GAMMA_FILE_LIST_URL, KICKSTART_BASE_URL, $window) {
    'use strict';
    /**
     * download contents from any url
     */
    function rawFile(url) {
      return $http.get(url).then(function(res) { return res.data; });
    }

    /**
     * download contents from a "fetch" path
     */
    function fetchFile(path){
      return rawFile(FETCH_URL + path);
    }

    this.rawFile = rawFile;
    this.fetchFile = fetchFile;
    this.kickstart = function(){
      return rawFile(KICKSTART_BASE_URL + $window.location.search);
    };
    this.nprofile = function(config) {
      return StatusPoller.start(config['nprofileData'])
        .then(fetchFile);
    };
    this.inputFileCds = function(config) {
      return StatusPoller.start(config['File_of_published_accepted_CDSs'])
        .then(fetchFile);
    };

    this.acgtGammaFileList = function(config) {
      return StatusPoller.start(config['acgt_gamma_output'])
        .then(function(path) {
          return rawFile(ACGT_GAMMA_FILE_LIST_URL + path);
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
