angular.module('npact')
  .service('Utils', function($q, $timeout, $log) {
    'use strict';
    // sigil value for quitting forEachAsync
    var STOP_ITERATING = {};
    this.STOP_ITERATING = STOP_ITERATING;

    /**
     * find the nearest order of magnitude to `x`
     *
     * orderOfMagnitude(1012) => 1000
     * orderOfMagnitude(1012, -1) => 100
     *
     * @param {Number} x
     * @param {Number} [exponent] shift the result by this many digits
     */
    this.orderOfMagnitude = function(x, exponent){
      var k = Math.floor(Math.log(x) / Math.LN10),
          t = Math.floor(0.5 + x / Math.pow(10, k));
      return t * Math.pow(10, k + (exponent || 0));
    };

    /**
     * Loop over the list from back to front, asyncronously
     *
     * stop iteration early with `throw Utils.STOP_ITERATING`
     *
     * @returns {Promise} when the iteration is complete
     */
    this.forEachAsync = function(list, fn){
      var d = $q.defer(),
          p = d.promise;

      if(list === null){
        d.reject();
        return p;
      }

      var idx = 0,
          opts = {batchSize:512, delay: 0},
          batches = 0,
          len = list.length
      ;

      function iterate(){
        var t1 = new Date();
        try{
          for(var c = opts.batchSize; idx < len && c >= 0; idx++, c--){
            fn(list[idx], idx);
          }
        }catch(e){
          if(e === STOP_ITERATING){
            idx = len;
          }else{
            d.reject({error:e, line: list[idx]});
          }
        }
        var elapsed = new Date() - t1;
        batches++;
        if(idx < len){
          d.notify({batches:batches, idx:idx, len:len,
                    elapsed:elapsed,
                    batchSize:opts.batchSize});
          // if the loop was very fast, double the batch size, if slow, halve it
          var goalms = 10,
              rowsPerms = opts.batchSize / elapsed;

          opts.batchSize = rowsPerms * goalms;
          // queue the next batch
          $timeout(iterate, opts.delay, false);
        }else{
          d.resolve(batches);
        }
      }
      // start it off soon
      $timeout(iterate, opts.delay, false);
      return p;
    };

    /**
     * adds a "page" of data from src to dst
     */
    this.extendByPage = function(src, dst, pageSize) {
      // build up a list of splice args
      var args = src.slice(dst.length, dst.length + pageSize);
      args.unshift(0);
      args.unshift(dst.length);
      dst.splice.apply(dst, args);
    };
  })
  .service('LineParser', function(Utils, $q) {
    'use strict';
    /**
     * parses many lines
     *
     * @param {string} text - multi-line text
     * @param {function} parseLine - callback to parse a line
     * @returns {array} extract objects
     */
    this.parse = function(text, parseLine){
      if (!_.isString(text)) { return text; }
      var lines = text.split('\n');
      return _.map(lines, parseLine);
    };
    /**
     * parses many lines asynchronously
     *
     * @param {string} text - multi-line text
     * @param {function} parseLine - callback to parse a line
     * @returns {Promise} for array of parsed objects
     */
    this.parseAsync = function(text, parseLine){
      if (!_.isString(text)) {
        throw new Error("Can't parse non-text " + typeof text);
      }

      var results = [];
      // async loop, gathering results
      return Utils.forEachAsync(
        text.split('\n'),
        function(line, idx){
          if(line && line.length > 0) {
            results[idx] = parseLine(line);
          }
        })
        .then(function(){ return results; });
    };
  })
// creates parsers, suitable for returning from factory declarations
  .service('ParserFactory', function(LineParser){
    'use strict';
    this.create = function(parseLineFn) {
      return {
        parseLine: parseLineFn,
        parse: _.partialRight(LineParser.parse, parseLineFn),
        parseAsync:  _.partialRight(LineParser.parseAsync, parseLineFn)
      };
    };
  })
  .factory('ExtractParser', function(ParserFactory){
    'use strict';
    /**
     * Extract format
     *
     * The grammar is roughly:
     *   <extract>    = <name> <address>
     *   <address>    = <complement> | <range>
     *   <complement> = complement(<range>)
     *   <range>      = <base>..<base>
     *   <base>       = any integer
     */
    var EXTRACT_REGEX = /([^ ]+) (complement\()?(<)?(\d+)\.\.(>)?(\d+)/,
        // some indexes for where our regex groups will show
        NAME = 1, COMP = 2, START_APPROX = 3, START = 4, END_APPROX = 5, END = 6;

    /**
     * parses one line as an extract
     *
     * COLORING INDICATES THE COLOR OF THE GENOME PHASE (0,1, OR 2) OF
     * THE THIRD CODON POSITIONS OF GENES, IN INTERNAL COORDINATES. IF A
     * GENE IS ENCODED IN THE DIRECT STRAND, ITS THIRD CODON POSITIONS
     * ARE IN THE SAME PHASE AS THE SECOND COORDINATE. E.G., dnaN
     * 1460..2668, WILL HAVE COLOR CORRESPONDING TO (2668 – 1) % 3 = 0
     * (RED)). INSTEAD, IF IT IS ENCODED ON THE COMPLEMENTARY STRAND, IT
     * WILL HAVE COLOR CORRESPONDING TO THE POSITION INDICATED BY THE
     * FIRST COORDINATE. E.G., radA complement(23561..24766) HAS COLOR
     * (23561 – 1) % 3 = 1 (GREEN). ONE IS SUBTRACTED FROM POSITION TO
     * TRANSLATE OUTPUT/INPUT COORDINATES INTO INTERNAL COORDINATES.
     *
     * @param {string} line - single extract line
     * @returns {Object} extract
     */
    function parseExtract(line){
      var parts = EXTRACT_REGEX.exec(line),
          res = {
            start: parseInt(parts[START]),
            end: parseInt(parts[END]),
            complement: parts[COMP] ? 1 : 0,
            name: parts[NAME],
            approximate: (parts[START_APPROX] || parts[END_APPROX]) ? true : false
          },
          phaseCoordinate = res.complement ? res.start : res.end;
      res.phase = (phaseCoordinate - 1) % 3;
      return res;
    }

    return ParserFactory.create(parseExtract);
  })

  .factory('Translater', function($http, $log, TRANSLATE_URL) {
    'use strict';
    return function (dna, mycoplasma, complement) {
      return $http({
        method: 'POST',
        url: TRANSLATE_URL,
        data:{
          seq: dna,
          complement: complement,
          mycoplasma: mycoplasma},
        headers: {'Content-Type': 'application/json'}
      });
    };
  })

  .service('Fetcher', function(StatusPoller, $http, FETCH_BASE_URL, BASE_URL, PATH, $location) {
    'use strict';
    var self = this;
    self.buildUrl = function(verb) {
      return BASE_URL + '/' + verb + '/'
        + PATH + '?' + $.param($location.search());
    };

    /**
     * download contents from any url
     */
    self.rawFile = function(url) {
      return $http.get(url).then(function(res) { return res.data; });
    };



    /**
     * download contents from a "fetch" path
     */
    self.fetchFile = function(path){
      if(!path) {
        throw new Error("Path is undefined");
      }
      return self.rawFile(FETCH_BASE_URL + path);
    };


    /**
     * poll the server for when `path` is ready, then fetch it
     */
    self.pollThenFetch = function(path) {
      if(!path) { throw new Error("Path is undefined"); }
      var blockonurl = BASE_URL + '/blockon/' + path;
      return $http.get(blockonurl).then(function(res) {
        return res.data;
      });
    };
  })


  .service('StatusPoller', function(STATUS_BASE_URL, FETCH_BASE_URL,
                                    $q, $http, $timeout, $log) {
    'use strict';
    var initialPollTime = 200;

    function poller(tid, deferred, time) {
      // remember our arguments
      var pollAgain = _.partial(poller, tid, deferred, time * 1.5);

      $http.get(STATUS_BASE_URL + tid)
        .then(function(res) {
          if(res.data.ready) { deferred.resolve(tid); }
          else { $timeout(pollAgain, time); }
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
      return poller(tid, $q.defer(), initialPollTime);
    };
  })
;
