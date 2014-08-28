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

  .controller('npactGraphPageCtrl', function($scope, $http, KICKSTART_BASE_URL, GraphDealer, $window, Fetcher, ExtractParser, npactConstants) {

    $scope.miscFiles = [];
    $scope.graphHeight = npactConstants.graphSpecDefaults.height;

    $http.get(KICKSTART_BASE_URL + $window.location.search)
      .then(function(config) {
        $scope.title = config.first_page_title;
        Fetcher.nprofile(config).then(GraphDealer.setProfile);
        Fetcher.inputFileCds(config)
          .then(ExtractParser.parse)
          .then(function(data) {
            GraphDealer.addExtract({name: 'Input file CDS', data: data});
          });

        Fetcher.acgtGammaFileList(config).then(function(fileList) {
          $scope.miscFiles.push.apply($scope.miscFiles, fileList);
          Fetcher.rawFile(config['File_of_new_CDSs'])
            .then(ExtractParser.parse)
            .then(function(data) {
              GraphDealer.addExtract({name: 'Newly Identified ORFs', data: data});
            });

          //TODO: Hits line: config['File_of_G+C_coding_potential_regions']
          // Fetcher.rawFile(config['File_of_G+C_coding_potential_regions'])
          //   .then(hitsParser)
          //   .then(function(data) {
          //     //TODO: GraphDealer.addHits
          //     //GraphDealer.addHits({name: 'Hits', data: data});
          //   });
        });
      });
  })

  .service('Fetcher', function(StatusPoller, $http, GraphDealer, FETCH_URL, ACGT_GAMMA_FILE_LIST_URL) {
    function rawFile(path) {
      var url = FETCH_URL + path;
      return $http.get(url)
        .then(function(res) {
          return res.data;
        });
    };
    this.rawFile = rawFile;

    this.nprofile = function(config) {
      return StatusPoller.start(config['nprofileData'])
        .then(rawFile);
    };
    this.inputFileCds = function(config) {
      return StatusPoller.start(config['Input file CDS'])
        .then(rawFile);
    };

    this.acgtGammaFileList = function(config) {
      return StatusPoller.start(config['acgt_gamma_output'])
        .then(function(x) {
          var url = ACGT_GAMMA_FILE_LIST_URL + config['acgt_gamma_output'];
          return $http.get(url)
            .then(function(res) {
              return res.data;
            });
        });
    };
  })


  .service('StatusPoller', function(STATUS_BASE_URL, $q, $http, $timeout) {
    var POLLTIME = 2500;

    function poller(tid, deferred) {
      $http.get(STATUS_BASE_URL + tid)
        .then(function(res) {
          if(res.ready)
            deferred.resolve(tid);
          else
            $timeout(poller, POLLTIME);
        })
        .catch(function(err) {
          deferred.reject(err);
        });
    }

    this.start = function(tid) {
      var deferred = $q.defer();
      if(!tid || !tid.length == 0)
        deferred.reject(new Error("Invalid task id"));
      poller(tid, deferred);
      return deferred.promise;
    };
  });
