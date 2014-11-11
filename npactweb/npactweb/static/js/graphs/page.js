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

  .controller('npactGraphPageCtrl', function($scope, Fetcher, $q, $log, StatusPoller, FETCH_URL, $window, $element, GraphConfig, TrackReader, ProfileReader, Pynpact) {
    'use strict';
    // helper functions
    var getWidth = function(){ return $element.width(); },
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
          return Fetcher.fetchFile(config[Pynpact.NEW_CDS])
            .then(function(data) {
              return addExtract('Newly Identified ORFs', data);
            });
        },
        addHits = function(config) {
          return Fetcher.fetchFile(config[Pynpact.HITS])
            .then(function(data) { addTrack('hits', 'Hits', data); });
        },
        addProfile = function(config) {
          return Fetcher.nprofile(config)
            .then(ProfileReader.load);
        }
    ;

    $scope.miscFiles = [];
    $scope.status = 'Initializing';
    $scope.ready = false;
    $scope.FETCH_URL = FETCH_URL;

    $scope.$watch(getWidth, function(newValue, oldValue){
      if (newValue > 0){ GraphConfig.width = newValue; }
    });

    // check for scope changes on resize
    angular.element($window).bind('resize', function () { $scope.$apply(); });

    // start it up
    Fetcher.kickstart()
      .then(function(config) {
        $log.log('Kickstart successful:', config);
        $scope.status = 'Running';
        // got config, request the first round of results
        $scope.title = config[Pynpact.TITLE];
        $scope.email = config[Pynpact.EMAIL];
        $scope.configureUrl = config[Pynpact.CONFIGURE_URL];
        var nprofile = addProfile(config);

        // Non-gbk files don't have CDSs we can extract.
        var inputFileCds = config[Pynpact.HAS_CDS] && addInputFileCds(config);

        var extraFileList = Fetcher.acgtGammaFileList(config)
              .then(function(fileList) {
                $scope.miscFiles.push.apply($scope.miscFiles, fileList);
                addNewCds(config);
                addHits(config);
              });
        $q.all([nprofile, inputFileCds, extraFileList])
          .then(function() { delete $scope.status; })
          .catch(function(err) { $scope.status = err; });

        StatusPoller.start(config[Pynpact.PDF])
          .then(function(pdfFilename) { $scope.miscFiles.push(pdfFilename); });
      });
  })

  .service('Fetcher', function(StatusPoller, $http, FETCH_URL, ACGT_GAMMA_FILE_LIST_URL, KICKSTART_BASE_URL, $window, Pynpact) {
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
      return self.pollThenFetch(config[Pynpact.NPROFILE]);
    };

    self.inputFileCds = function(config) {
      return self.pollThenFetch(config[Pynpact.CDS]);
    };

    self.acgtGammaFileList = function(config) {
      return StatusPoller.start(config[Pynpact.ACGT_GAMMA_FILES])
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
