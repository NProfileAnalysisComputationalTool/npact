angular.module('npact')
  .directive('npactGraphPage', function(STATIC_BASE_URL) {
    'use strict';
    return {
      restrict: 'A',
      templateUrl: STATIC_BASE_URL + 'js/graphs/page.html',
      controller: 'npactGraphPageCtrl',
      controllerAs: 'pageCtrl'
    };
  })

  .service('MessageBus', function($log) {
    this.error = function(msg) {
      return this.log('error', msg);
    };
    this.info = function(msg) {
      return this.log('info', msg);
    };
    this.log = function(level, msg) {
      $log.log(level, msg);
      /// TODO: Make these user-visible on the screen
      // var msgpane = angular.element('#msgpane');
      // var newmessage = angular.element('<div>' + msg + '</div>')
      //       .addClass('alert').addClass('alert-' + level)
      //       .css({display: 'none'});
      // msgpane.append(newmessage);
      // newmessage.slideDown(1000, function() {
      //   $timeout(function () { newmessage.slideUp(1000); },
      //            MessageDisplayTime);
      // });
    };
  })

  .controller('npactGraphPageCtrl', function($scope, $q, $log, Fetcher, FETCH_URL,
                                      GraphConfig, FileManager, kickstarter) {
    'use strict';

    $scope.miscFiles = [];
    $scope.status = 'Initializing';
    $scope.ready = false;
    $scope.FETCH_URL = FETCH_URL;
    $scope.config = GraphConfig;

    kickstarter.start();
    $scope.$watch(FileManager.getFiles, function(val) {
      $scope.miscFiles = val;
    }, true);

    this.print = function() {
      $log.log("print requested: ", new Date());
      $scope.$broadcast('print');
      $log.log("finished rendering: ", new Date());
    };
  })

  .service('kickstarter', function($q, $log, processOnServer, MessageBus,
                            NProfiler, PredictionManager, ExtractManager, FileManager) {
    'use strict';
    //Kickstart the whole process, start all the main managers
    this.start = function() {
      MessageBus.info('kickstarting');
      this.basePromise = processOnServer('parse');
      this.basePromise.then(NProfiler.start);
      this.basePromise.then(PredictionManager.start);
      this.basePromise.then(ExtractManager.start);
      this.basePromise.then(FileManager.start);
      return this.basePromise;
    };
  })

  .service('ExtractManager', function(Fetcher, Pynpact, TrackReader, GraphConfig,
                               processOnServer, MessageBus, $log) {
    'use strict';
    this.start = function(config) {
      if(config.format != 'genbank') { return; }
      processOnServer('extract').then(function(config) {
        if(config[Pynpact.CDS]) {
          Fetcher.pollThenFetch(config[Pynpact.CDS])
            .then(function(data) {
              var name = 'Input file CDS';
              TrackReader.load(name, data)
                .then(function() { GraphConfig.loadTrack(name, 'extracts'); });
            });
        }
      });
    };
  })
  .service('PredictionManager', function(Fetcher, StatusPoller, Pynpact, TrackReader,
                                  GraphConfig, processOnServer, $log,
                                  ACGT_GAMMA_FILE_LIST_URL) {
    'use strict';
    var self = this;
    self.files = null;
    self.updateFiles = function(path) {
      $log.log('Fetching the ACGT_GAMMA_FILE_LIST from', path);
      Fetcher.rawFile(ACGT_GAMMA_FILE_LIST_URL + path)
        .then(function(files) {
          $log.log('ACGT Gamma Files ready', files);
          self.files = files;
        });
    };
    self.disableTrack = function(nameBase, oldSig) {
      if(oldSig) {
        var oldTrack = GraphConfig.findTrack(nameBase + oldSig);
        if(oldTrack) { oldTrack.active = false; }
      }
    };
    self.newHits = function(waitOn, config, oldSig) {
      var nameBase = 'Hits @';
      self.disableTrack(nameBase, oldSig);
      waitOn.then(function(path) {
        Fetcher.fetchFile(config[Pynpact.HITS])
          .then(function(data) {
            var type = 'hits', name = nameBase + config.significance;
            TrackReader.load(name, data)
              .then(function() { GraphConfig.loadTrack(name, type, 100); });

          });
      });
    };
    self.newCds = function(waitOn, config, oldSig) {
      var nameBase = 'New ORFs @';
      self.disableTrack(nameBase, oldSig);
      waitOn.then(function(path) {
        Fetcher.fetchFile(config[Pynpact.NEW_CDS])
          .then(function(data) {
            var name = nameBase + config.significance;
            TrackReader.load(name, data)
              .then(function() { GraphConfig.loadTrack(name, 'extracts', 15); });
          });
      });
    };

    self.onSignificanceChange = function(significance, oldSig) {
      if(significance) {
        $log.log("New prediction significance: ", significance);
        self.files = null;
        processOnServer('acgt_gamma').then(function(config) {
          var waitOn = StatusPoller.start(config[Pynpact.ACGT_GAMMA_FILES]);
          self.newHits(waitOn, config, oldSig);
          self.newCds(waitOn, config, oldSig);
          waitOn.then(self.updateFiles);
          return waitOn;
        });
      }
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
