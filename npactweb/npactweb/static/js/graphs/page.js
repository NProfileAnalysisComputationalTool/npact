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

  .controller('npactGraphPageCtrl', function($scope, Fetcher, $q, $log, StatusPoller, FETCH_URL, GraphConfig, Pynpact, KickstartManager, FileManager) {
    'use strict';

    $scope.miscFiles = [];
    $scope.status = 'Initializing';
    $scope.ready = false;
    $scope.FETCH_URL = FETCH_URL;

    KickstartManager.start().then(function(config) {
      $scope.status = 'Running';
      // got config, request the first round of results
      $scope.title = config[Pynpact.TITLE];
      $scope.email = config[Pynpact.EMAIL];
      $scope.configureUrl = config[Pynpact.CONFIGURE_URL];

      // $q.all([nprofile, inputFileCds, extraFileList])
      //     .then(function() { delete $scope.status; })
      //     .catch(function(err) { $scope.status = err; });
    });
    $scope.$watch(FileManager.getFiles, function(val) {
      $scope.miscFiles = val;
    }, true);
  })

  .service('KickstartManager', function(KICKSTART_BASE_URL, ProfileReader, TrackReader, $window, $http, $log, NProfiler, PredictionManager, ExtractManager, FileManager) {
    this.start = function() {
      var url = KICKSTART_BASE_URL + $window.location.search;
      var basePromise = $http.get(url)
        .then(function(res) {
          $log.log('Kickstart successful:', res.data);
          return res.data;
        });
      basePromise.then(NProfiler.start);
      basePromise.then(PredictionManager.start);
      basePromise.then(ExtractManager.start);
      basePromise.then(FileManager.start);
      return basePromise;
    };
  })
  .service('ExtractManager', function(Fetcher, Pynpact, TrackReader, GraphConfig) {
    this.start = function(config) {
      if(config[Pynpact.HAS_CDS]) {
        Fetcher.inputFileCds(config)
          .then(function(data) {
            var name = 'Input file CDS';
            TrackReader.load(name, data)
              .then(function() { GraphConfig.loadTrack(name, 'extracts'); });
          });
      }
    };
  })
  .service('PredictionManager', function(Fetcher, Pynpact, TrackReader, GraphConfig, $log) {
    var self = this;
    self.files = null;
    self.start = function(config) {
      Fetcher.fetchFile(config[Pynpact.NEW_CDS])
        .then(function(data) {
          var name = 'Newly Identified ORFs';
          TrackReader.load(name, data)
            .then(function() { GraphConfig.loadTrack(name, 'extracts'); });
        });
      Fetcher.fetchFile(config[Pynpact.HITS])
        .then(function(data) {
          var type = 'hits';
          var name = 'Hits';
          TrackReader.load(name, data)
            .then(function() { GraphConfig.loadTrack(name, type); });

        });
      Fetcher.acgtGammaFileList(config).then(function(files) {
        $log.log('ACGT Gamma Files ready', files);
        self.files = files;
      });
    };
  })
  .service('NProfiler', function(ProfileReader, Fetcher) {
    this.start = function(config) {
      return Fetcher.nprofile(config).then(ProfileReader.load);
    };
  })
  .service('FileManager', function(PredictionManager, StatusPoller, Pynpact, $log) {
    var pdffile = null;
    this.start = function(config) {
      StatusPoller.start(config[Pynpact.PDF])
        .then(function(pdfFilename) {
          $log.log('PDF ready', pdfFilename);
          pdffile = pdfFilename;
        });
    };
    this.getFiles = function() {
      var list = [];
      if(pdffile) {
        list.push(pdffile);
      }
      if(PredictionManager.files) {
        list.push.apply(list, PredictionManager.files);
      }
      return list;
    };
  })
  ;
