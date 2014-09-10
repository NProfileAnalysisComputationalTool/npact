angular.module('npact')
  .directive('npactGraphPage', function(STATIC_BASE_URL, Utils, GraphDealer) {
    return {
      restrict: 'A',
      templateUrl: STATIC_BASE_URL + 'js/graphs/page.html',
      controller: 'npactGraphPageCtrl',
      link: function($scope, $element, $attrs) {
        // TODO: watch for changes in width
        // http://stackoverflow.com/questions/23044338/window-resize-directive
        Utils.widthAvailable($element).then(GraphDealer.setWidth);
      }
    };
  })

  .controller('npactGraphPageCtrl', function($scope, GraphDealer, Fetcher, npactConstants, $q, $log) {

    $scope.miscFiles = [];
    $scope.graphHeight = npactConstants.graphSpecDefaults.height;
    $scope.status = 'Initializing';
    // start it up
    Fetcher.kickstart()
      .then(function(config) {
        $log.log('Kickstart successful:', config);
        $scope.status = 'Running';
        $scope.config = config;
        // got config, request the first round of results
        $scope.title = config.first_page_title;
        var nprofile = Fetcher.nprofile(config).then(GraphDealer.setProfile);
        var inputFileCds = Fetcher.inputFileCds(config)
          .then(function(data) {
            GraphDealer.addExtract({name: 'Input file CDS', data: data});
          });

        var extraFileList = Fetcher.acgtGammaFileList(config)
              .then(function(fileList) {
                $scope.miscFiles.push.apply($scope.miscFiles, fileList);
                Fetcher.fetchFile(config['File_of_new_CDSs'])
                  .then(function(data) {
                    GraphDealer.addExtract({name: 'Newly Identified ORFs', data: data});
                  });

                //TODO: Hits line: config['File_of_G+C_coding_potential_regions']
                // Fetcher.fetchFile(config['File_of_G+C_coding_potential_regions'])
                //   .then(hitsParser)
                //   .then(function(data) {
                //     //TODO: GraphDealer.addHits
                //     //GraphDealer.addHits({name: 'Hits', data: data});
                //   });
              });
        $q.all([nprofile, inputFileCds, extraFileList]).then(function() {
          $scope.status = 'Finished';
        });
      });
  })

  .service('Fetcher', function(StatusPoller, $http, FETCH_URL, ACGT_GAMMA_FILE_LIST_URL, KICKSTART_BASE_URL, $window) {

    /**
     * download contents from any url
     */
    function rawFile(url) {
      return $http.get(url).then(function(res) { return res.data; });
    };

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


  .service('StatusPoller', function(STATUS_BASE_URL, $q, $http, $timeout) {
    var POLLTIME = 2500;

    function poller(tid, deferred) {
      // remember our arguments
      var pollAgain = _.partial(poller, tid, deferred);

      $http.get(STATUS_BASE_URL + tid)
        .then(function(res) {
          if(res.data.ready) { deferred.resolve(tid); }
          else { $timeout(pollAgain, POLLTIME); }
        })
        .catch(deferred.reject);

      return deferred.promise;
    }

    this.start = function(tid) {
      if(!tid || tid.length == 0)
        return $q.reject(new Error("Invalid task id"));

      return poller(tid, $q.defer());
    };
  });
