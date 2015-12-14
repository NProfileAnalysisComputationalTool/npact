angular.module('npact')
  .controller('Results', function ($scope) {
    'use strict';
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
        ZoomWindowHandler.maybePopup();
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
        $timeout(function () { $log.debug("Rebuilding"); $scope.$broadcast(Evt.REBUILD); });
      });
  })

  .controller('DownloadsCtrl', function($scope, $log, PredictionManager, MessageBus, Pynpact, StatusPoller, GraphConfig, FETCH_BASE_URL, Fetcher) {
    'use strict';
    $scope.FETCH_BASE_URL = FETCH_BASE_URL;
    $scope.$watch( function() { return PredictionManager.files; },
                   function(val) { $scope.predictionFiles = val; });
    $scope.GraphConfig = GraphConfig;
    this.buildGBKDownload = function () {
      return Fetcher.buildUrl('BuildGBK', {trackPaths: _.map(GraphConfig.activeTracks(), 'filename').join(',')});
    };
  })
  .directive('lazyHref', function ($log) {
    return {
      scope: {
        lazyHref: '&'
      },
      link: function($scope, $element, attrs) {
        $element.on('mouseenter focus', function () {
          var newUrl = $scope.lazyHref();
          $log.debug("Updated lazyHref to", newUrl);
          $element.attr('href', newUrl);
        });
      }
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

  .service('kickstarter', function($q, $log, processOnServer, MessageBus,
                            NProfiler, PredictionManager, Fetcher, Track, GraphConfig) {
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
            PredictionManager.start();
          }
          else{
            var pths = GraphConfig.trackPaths;
            if(typeof(pths) === 'string') pths=pths.split(',');
            Track.fetchAllTracks(pths)
              .then(function(tracks) { GraphConfig.tracks = tracks; });
          }
          return null;
        })
      ]);
    };
  })

  .service('PredictionManager', function(Fetcher, StatusPoller, Pynpact, Track,
                                  $rootScope, $timeout, $log,
                                  GraphConfig, processOnServer, MessageBus) {
    'use strict';
    var self = this;
    self.files = null;
    self.predictionTracks = null;

    self.start = function(_config) {
      var url = Fetcher.buildUrl('acgt_gamma');
      $log.log("Querying acgt_gamma at", url);
      var acgt_gamma_promise = Fetcher.rawFile(url)
          .then(function(response) {
            self.files = response.files;
            Track.fetchAllTracks(response.trackPaths)
              .then(function(tracks) {
                var tracksToAdd = _.difference(tracks, self.predictionTracks);
                var tracksToRem = _.difference(self.predictionTracks, tracks);
                GraphConfig.tracks = _(GraphConfig.tracks)
                  .difference(tracksToRem)
                  .union(tracksToAdd).value();

                self.predictionTracks = tracks;

                return;
                //////////
                _.each(tracks, function(track) {
                  track.active = true;
                  if(!_.find(GraphConfig.tracks, 'filename', track.filename)) {
                    GraphConfig.tracks.push(track);
                  }
                });
              });
          });

      MessageBus.info(
        'Identifying significant 3-base periodicities',
        acgt_gamma_promise.catch(function(e) {
          MessageBus.danger('Failure while identifying significant 3-base periodicities');
        }));
    };

    $rootScope.$watch(
      function () { return GraphConfig.mycoplasma; },
      function (val, old) {
        if(old !== undefined) {
          $log.debug("mycoplasma changed from", old, "to", val);
          // Wait so that the above update to $location has a chance to complete
          $timeout(self.start, 50);
        }
      });
  });
