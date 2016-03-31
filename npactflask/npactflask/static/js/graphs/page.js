angular.module('npact')
  .controller('ResultsCtrl', function ($scope, $location, $anchorScroll, $log) {
    'use strict';
    var self = this;
    self.downloadsopen = false;
    self.graphsopen = true;
    $scope.$watch(
      function () {return $location.hash();     },
      function (hash) {
        $log.log("Hash changed to", hash);
        if(_.startsWith(hash, "downloads")) {
          self.downloadsopen = true;
          $anchorScroll(hash);
        }
        else if(_.startsWith(hash, "graphs")) {
          self.graphsopen = true;
          $anchorScroll(hash);
        }
      });
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

  .controller('npactGraphPageCtrl', function($scope,$q, $window, $log, $location, $timeout,
                                      Fetcher, BASE_URL, PATH, Evt,
                                      FETCH_BASE_URL, EmailBuilder,
                                      STATIC_BASE_URL, GraphConfig, GraphingCalculator,
                                      kickstarter, processOnServer,
                                      PrintModal, ZoomWindowHandler) {
    'use strict';

    $scope.FETCH_BASE_URL = FETCH_BASE_URL;
    $scope.config = GraphConfig;
    $scope.email = EmailBuilder.send;
    $scope.BASE_URL = BASE_URL;
    $scope.PATH = PATH;

    var bootp = kickstarter.start();

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

    this.print = function() { PrintModal.show().then(_doPrint); };
    this.requestPDF = function(){ $window.open(Fetcher.buildUrl('getpdf'), '_blank'); };

    ZoomWindowHandler.register($scope);
    bootp
      .then(function () {
        var p = ZoomWindowHandler.maybePopup($scope);
        $log.log("Everything finished initial boot");
      });

    /*** Watch the config for changes we care about ***/
    $scope.$watchCollection(
      function() { return [GraphConfig.basesPerGraph,
                    GraphConfig.offset, GraphConfig.startBase, GraphConfig.endBase]; },
      function() {
        // basic row geometry changed, repartition and rebuild
        if(isNaN(GraphConfig.startBase) ||
           isNaN(GraphConfig.endBase) ||
           isNaN(GraphConfig.basesPerGraph)) { return; }

        $scope.graphSpecs = GraphingCalculator.partition(GraphConfig);
        $log.log('Partitioned into', $scope.graphSpecs.length, 'rows.');
        //$timeout(function () { $log.debug("Rebuilding"); $scope.$broadcast(Evt.REBUILD); });
      });
  })

  .controller('DownloadsCtrl', function($scope, $log, PredictionManager, MessageBus, Pynpact, StatusPoller, GraphConfig, FETCH_BASE_URL, Fetcher) {
    'use strict';
    $scope.FETCH_BASE_URL = FETCH_BASE_URL;
    $scope.$watch( function() { return PredictionManager.files; },
                   function(val) { $scope.predictionFiles = val; });
    $scope.GraphConfig = GraphConfig;
    this.buildGBKDownload = function () {
      var trackPaths = _(GraphConfig.tracks)
          .filter({active: true, type: 'extracts'})
          .pluck('filename').join(',');
      return Fetcher.buildUrl('build_gbk', { trackPaths: trackPaths });
    };
  })
  .service('PrintModal', function(STATIC_BASE_URL, $uibModal) {
    'use strict';
    var printTemplate = STATIC_BASE_URL + 'js/graphs/printConfirm.html';
    this.show = function() {
      var modalInstance = $uibModal.open({
        animation: true,
        templateUrl: printTemplate
      });
      return modalInstance.result;
    };
  })

  .service('kickstarter', function($q, $log, processOnServer, MessageBus, $http,
                            NProfiler, PredictionManager, Fetcher, Track,
                            GraphConfig, CodonFinder) {
    'use strict';
    //Kickstart the whole process, start all the main managers
    this.start = function() {
      this.basePromise = processOnServer('parse');
      var safePromise = this.basePromise.catch(function(e) {
        $log.error(e);
        MessageBus.danger('There\'s been an error starting up; please try starting over.');
      });
      MessageBus.info('kickstarting', safePromise);
      return $q.all([
        this.basePromise.then(NProfiler.start),
        this.basePromise.then(function() {
          if(!GraphConfig.trackPaths) {
            console.log('Starting PredictionManager, because: ', GraphConfig.trackPaths)
            PredictionManager.start();
          }
          else{
            var pths = GraphConfig.trackPaths;
            if(typeof(pths) === 'string') pths=pths.split(',');
            Track.fetchAllTracks(pths)
              .then(function(tracks) {
                GraphConfig.tracks = tracks;
                GraphConfig.clearORFSelection();
              });
          }
          return null;
        }),
        this.basePromise.then(function() {
          //console.log('Getting DDNA string');
          Fetcher.fetchFile(GraphConfig.ddna).then(function(res){
            GraphConfig.ddnaString = res;
            CodonFinder.reindex();
          });
        })
      ]);
    };
  })

  .service('PredictionManager', function(Fetcher, StatusPoller, Pynpact, Track,
                                  $rootScope, $timeout, $log, $q,
                                  GraphConfig, processOnServer, MessageBus) {
    'use strict';
    var self = this;
    self.files = null;
    self.predictionTracks = null;
    self.doPrediction = function() {
      var url = Fetcher.buildUrl('acgt_gamma');
      $log.log("Querying acgt_gamma at", url);
      var acgt_gamma_promise = Fetcher.rawFile(url)
          .then(function(response) {
            self.files = response.files;
            Track.fetchAllTracks(response.trackPaths)
              .then(function(tracks) {
                $log.log('Fetched tracks', tracks);
                var tracksToAdd = _.difference(tracks, self.predictionTracks);
                var tracksToRem = _.difference(self.predictionTracks, tracks);
                GraphConfig.tracks = _(GraphConfig.tracks)
                  .difference(tracksToRem)
                  .union(tracksToAdd).value();

                self.predictionTracks = tracks;

                return;
              });
          });

        MessageBus.info(
          'Identifying significant 3-base periodicities',
          acgt_gamma_promise.catch(function(e) {
            MessageBus.danger('Failure while identifying significant 3-base periodicities');
          }));
      return acgt_gamma_promise;
    };
    self.start = function(_config) {
      $log.log("PredictionManager.start");
      if(GraphConfig.trackPaths){
        $log.log('skipping because we have explicit graphs, how are we even getting here');
        return $q.when(true);
      }
      else{
        var acgt_gamma_promise = self.doPrediction();
        return acgt_gamma_promise;
      }
    };

    $rootScope.$watch(
      function () { return GraphConfig.mycoplasma; },
      function (val, old) {
        if(old !== undefined && old != val) {
          // Wait so that the above update to $location has a chance to complete
          $timeout(function() {
            $log.debug("mycoplasma changed from", old, "to", val);
            self.doPrediction();
          }, 50);
        }
      });
  });
