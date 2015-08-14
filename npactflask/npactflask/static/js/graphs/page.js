angular.module('npact')
  .controller('Results', function ($scope) {
    $scope.status = {
      isFirstOpen: true,
      isFirstDisabled: false
    };
  })

  .directive('npactGraphPage', function(STATIC_BASE_URL) {
    'use strict';
    return {
      restrict: 'A',
      templateUrl: STATIC_BASE_URL + 'js/graphs/page.html',
      controller: 'npactGraphPageCtrl',
      controllerAs: 'pageCtrl'
    };
  })

  .controller('npactGraphPageCtrl', function($scope,$q, $window, $log, PrintModal,
                                      Fetcher, BASE_URL, PATH,
                                      FETCH_BASE_URL, EmailBuilder,
                                      STATIC_BASE_URL, GraphConfig,
                                      kickstarter, processOnServer, $location) {
    'use strict';

    $scope.FETCH_BASE_URL = FETCH_BASE_URL;
    $scope.config = GraphConfig;
    $scope.email = EmailBuilder.send;
    $scope.BASE_URL = BASE_URL;
    $scope.PATH = PATH;
    
    kickstarter.start();

    var _doPrint = function() {
      var t1 = new Date();
      $log.log("print requested: ", t1);
      $scope.$broadcast('printresize', true);
      $scope.printCounter = 0;
      var waitList = [];
      var evt = $scope.$broadcast('print', function(promise) {
        $scope.printCounter++;
        waitList.push(promise);
      });
      $q.all(waitList).then(function() {
        $scope.printCounter = 0;
        $log.log('Finished rendering', new Date() - t1);
        $window.print();
        $scope.$broadcast('printresize', false);
      });
    };

    this.print = function() {
      PrintModal.show().then(_doPrint);
    };

    this.requestPDF = function(){
      $scope.PARAMS = $.param($location.search());
      $window.open(Fetcher.buildUrl('getpdf'), '_blank');
    };
  })

  .controller('DownloadsCtrl', function($scope, $log, PredictionManager, MessageBus, Pynpact, StatusPoller, GraphConfig, dialogService) {
    'use strict';
    $scope.$watch( function() { return PredictionManager.files; },
                   function(val) { $scope.predictionFiles = val; });
  })

  .service('PrintModal', function(STATIC_BASE_URL, $modal) {
    var printTemplate = STATIC_BASE_URL + 'js/graphs/printConfirm.html';
    this.show = function() {
      var modalInstance = $modal.open({
        animation: true,
        templateUrl: printTemplate,
        controller: 'ModalInstanceCtrl'
      });
      return modalInstance.result;
    };
  })

.controller('ModalInstanceCtrl', function($scope, $modalInstance) {
    $scope.proceed = function() {
      $modalInstance.close();
    };

    $scope.cancel = function() {
      $modalInstance.dismiss();
    };
  })

  .service('kickstarter', function($q, $log, processOnServer, MessageBus,
                            NProfiler, PredictionManager, ExtractManager) {
    'use strict';
    //Kickstart the whole process, start all the main managers
    this.start = function() {
      this.basePromise = processOnServer('parse');
      var safePromise = this.basePromise.catch(function(e) {
        $log.error(e);
        MessageBus.danger('There\'s been an error starting up; please try starting over.');
      });
      MessageBus.info('kickstarting', safePromise);
      this.basePromise.then(NProfiler.start);
      this.basePromise.then(PredictionManager.start);
      this.basePromise.then(ExtractManager.start);
      return this.basePromise;
    };
  })

  .service('ExtractManager', function(Fetcher, Pynpact, Track, GraphConfig,
                               processOnServer, MessageBus, $log) {
    'use strict';
    this.start = function(config) {
      if(config.format != 'genbank') { return; }
      var p = processOnServer('extract').then(function(config) {
        if(config[Pynpact.CDS]) {
          Fetcher.pollThenFetch(config[Pynpact.CDS])
            .then(function(data) {
              return GraphConfig.loadTrack(new Track('Input file CDS', data, 'extracts'));
            });
        }
      });
      MessageBus.info(
        "Fetching extract data from server",
        p.catch("Failure while extracting known genes."));
    };
  })
  .service('PredictionManager', function(Fetcher, StatusPoller, Pynpact, Track,
                                  GraphConfig, processOnServer, $log,
                                  MessageBus, PATH, $location,
                                  ACGT_GAMMA_FILE_LIST_BASE_URL) {
    'use strict';
    var self = this;
    self.files = null;
    var results = {}; //hash keyed on significance of already requested results.

    self.start = function(_config) {
      var url = Fetcher.buildUrl('acgt_gamma');
      var acgt_gamma_promise = Fetcher.rawFile(url)
          .then(function(response) {
            self.files = response.files;
            var config = response.config;
            angular.extend(GraphConfig, config);
            return GraphConfig;
          });
      acgt_gamma_promise.then(self.fetchHits);
      acgt_gamma_promise.then(self.fetchModifiedOrfs);
      acgt_gamma_promise.then(self.fetchNewOrfs);

      MessageBus.info(
        'Identifying significant 3-base periodicities',
        acgt_gamma_promise.catch(function(e) {
          MessageBus.danger('Failure while identifying significant 3-base periodicities');
        }));
    };

    self.fetchHits = function(config) {
      config.hits = Fetcher.fetchFile(config[Pynpact.HITS])
          .then(function(data) {
            var name = 'Hits',
                track = new Track(name, data, 'hits', 100);
            GraphConfig.loadTrack(track);
            return track;
          });
    };

    self.fetchNewOrfs = function(config) {
      config.newORFs = Fetcher.fetchFile(config[Pynpact.NEW_ORFS])
        .then(function(data) {
          var name = 'New ORFs',
              track = new Track(name, data, 'neworfs', 15);
          GraphConfig.loadTrack(track);
          return track;
        });
    };

    self.fetchModifiedOrfs = function(config) {
      config.newORFs = Fetcher.fetchFile(config[Pynpact.MODIFIED])
        .then(function(data) {
          var name = 'Modified ORFs',
              track = new Track(name, data, 'modified', 10);
          GraphConfig.loadTrack(track);
          return track;
        });
    };
  });
