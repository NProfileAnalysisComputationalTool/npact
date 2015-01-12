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

  .controller('npactGraphPageCtrl', function($scope, Fetcher, $q, $log, StatusPoller, FETCH_URL, GraphConfig, Pynpact, FileManager, kickstarter) {
    'use strict';

    $scope.miscFiles = [];
    $scope.status = 'Initializing';
    $scope.ready = false;
    $scope.FETCH_URL = FETCH_URL;

    kickstarter.start().then(function(config) {
      $scope.status = 'Running';
      // got config, request the first round of results
      $scope.title = config[Pynpact.TITLE];
      $scope.email = config[Pynpact.EMAIL];
      $scope.configureUrl = config[Pynpact.CONFIGURE_URL];

      kickstarter.everything
        .catch(function(err) {
          $scope.error = true;
        })
        .then(function() { delete $scope.status; });

    });
    $scope.$watch(FileManager.getFiles, function(val) {
      $scope.miscFiles = val;
    }, true);
  })

  .service('kickstarter', function($q, $log, processOnServer,
                            NProfiler, PredictionManager, ExtractManager, FileManager) {
    'use strict';
    //Kickstart the whole process, start all the main managers
    this.start = function() {
      $log.log('kickstarting');
      this.basePromise = processOnServer( 'extract');

      this.everything = $q.all([
        this.basePromise.then(NProfiler.start),
        this.basePromise.then(PredictionManager.start),
        this.basePromise.then(ExtractManager.start),
        this.basePromise.then(FileManager.start)]);
      return this.basePromise;
    };
  })

  .service('ExtractManager', function(Fetcher, Pynpact, TrackReader, GraphConfig, $log) {
    'use strict';
    this.start = function(config) {
      if(config[Pynpact.CDS]) {
        Fetcher.pollThenFetch(config[Pynpact.CDS])
          .then(function(data) {
            var name = 'Input file CDS';
            TrackReader.load(name, data)
              .then(function() { GraphConfig.loadTrack(name, 'extracts'); });
          });
      }
    };
  })
  .service('PredictionManager', function(Fetcher, StatusPoller, Pynpact, TrackReader, GraphConfig,
                                  processOnServer, $log, ACGT_GAMMA_FILE_LIST_URL) {
    'use strict';
    var self = this;
    self.files = null;
    self.process = function(config) {
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
            .then(function() { GraphConfig.loadTrack(name, type, 100); });

        });
      StatusPoller.start(config[Pynpact.ACGT_GAMMA_FILES])
        .then(function(path) {
          return Fetcher.rawFile(ACGT_GAMMA_FILE_LIST_URL + path);
        })
        .then(function(files) {
          $log.log('ACGT Gamma Files ready', files);
          self.files = files;
        });
    };
    //just match the api of expecting a 'start' method
    self.start = self.process;

    self.onSignificanceChange = function(significance) {
      $log.log("New prediction significance: ", significance);
      self.files = null;
      processOnServer('acgt_gamma').then(self.process);
    };
  })
  .service('FileManager', function(PredictionManager, StatusPoller, Pynpact, $log) {
    'use strict';
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
