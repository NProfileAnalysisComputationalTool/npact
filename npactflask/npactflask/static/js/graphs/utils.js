angular.module('npact')

  .directive('npConfirmClick', function(){
    //http://stackoverflow.com/questions/18313576/confirmation-dialog-on-ng-click-angularjs
    return {
      link: function (scope, element, attr) {
        var msg = attr.npConfirmClick || "Are you sure?";
        var clickAction = attr.confirmedClick;
        element.bind('click',function (event) {
          if ( window.confirm(msg) ) {
            //console.log('applying confirmedclick', clickAction);
            scope.$eval(clickAction);
          }
        });
      }
    };
  })
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
    this.orderOfMagnitude = function(x, exponent) {
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
    this.forEachAsync = function(list, fn) {
      var d = $q.defer(),
          p = d.promise;

      if(list === null) {
        d.reject();
        return p;
      }

      var idx = 0,
          opts = {batchSize: 512, delay: 0},
          batches = 0,
          len = list.length;

      function iterate() {
        var t1 = new Date();
        try {
          for(var c = opts.batchSize; idx < len && c >= 0; idx++, c--) {
            fn(list[idx], idx);
          }
        }catch(e) {
          if(e === STOP_ITERATING) {
            idx = len;
          } else{
            d.reject({error:e, line: list[idx]});
          }
        }
        var elapsed = new Date() - t1;
        batches++;
        if(idx < len) {
          d.notify({batches: batches, idx: idx, len: len,
                    elapsed: elapsed,
                    batchSize: opts.batchSize});
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


    // a lot of our math is based on codons which are three in length
    // so we need to find the nearest multiple of 3 based on some
    // start
    this.round3 = function(x) {
      x = Number(x);
      return x + 1.5 - (x + 1.5) % 3; };
    this.ceil3 = function(x) {
      x = Number(x);
      var r = x%3; return r ? x + 3 - r : x; };
    this.floor3 = function (x) {
      if(!x || isNaN(x)) return null;
      x = Number(x);
      return x - x%3; };

    /*
    this.rangeModulo3 = function(startIdx, stopIdx) {
      // start is 161
      // stop codon index is 236 for indices 236,237,238
      // 238 is part of the codon, so we need to include it
      return (((stopIdx+1) - startIdx) % 3);
    };
    */
    this.isValidRange = function(start, stop) {
      return 0 === this.rangeModulo3(start, stop);
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
    this.parse = function(text, parseLine) {
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
    this.parseAsync = function(text, parseLine) {
      if (!_.isString(text)) {
        throw new Error("Can't parse non-text " + typeof text);
      }

      var results = [];
      // async loop, gathering results
      return Utils.forEachAsync(
        text.split('\n'),
        function(line, idx) {
          if(line && line.length > 0) {
            results[idx] = parseLine(line);
          }
        })
        .then(function() { return results; });
    };
  })

   /**
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
    */
  .factory('getPhase', function () {
    return function (orf) {
      var phaseCoord = orf.complement ? orf.start : orf.end;
      return (phaseCoord - 1) % 3;
    };
  })

// creates parsers, suitable for returning from factory declarations
  .service('ParserFactory', function(LineParser) {
    'use strict';
    this.create = function(parseLineFn) {
      return {
        parseLine: parseLineFn,
        parse: _.partialRight(LineParser.parse, parseLineFn),
        parseAsync:  _.partialRight(LineParser.parseAsync, parseLineFn)
      };
    };
  })
  .factory('ExtractParser', function(ParserFactory, getPhase) {
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
     *
     * @param {string} line - single extract line
     * @returns {Object} extract
     */
    function parseExtract(line) {
      if(line && line[0] === '#'){
        pair = line.split(':');
        return {
          key: pair[0],
          value: pair[1],
          type: 'META'
        };
      }
      var parts = EXTRACT_REGEX.exec(line),
          res = {
            start: parseInt(parts[START]),
            end: parseInt(parts[END]),
            complement: parts[COMP] ? 1 : 0,
            name: parts[NAME],
            approximate: (parts[START_APPROX] || parts[END_APPROX]) ? true : false,
            type:'CDS'
          };
      res.phase = getPhase(res);
      return res;
    }
    var parser = ParserFactory.create(parseExtract);
    return parser;
  })

  .factory('Translater', function($http, $log, Fetcher) {
    'use strict';
    return function (dna, mycoplasma, complement) {
      return $http({
        method: 'POST',
        url: Fetcher.buildUrl('translate'),
        data: {
          seq: dna,
          complement: complement,
          mycoplasma: mycoplasma},
        headers: {'Content-Type': 'application/json'}
      })
        .then(function (res) {
          return res.data;
        });
    };
  })

  .factory('TranslatePath', function ($http, $log, Fetcher) {
    'use strict';
    return function (start, end, complement) {
      var url = Fetcher.buildUrl('translate', {
        startBase: start, endBase: end, rc: complement ? 'true' : 'false'});
      // $log.log("Translating @ ", url);
      return $http.get(url, {responseType: 'json'})
        .then(function (res) {
          return res.data;
        });
    };
  })
  .service('CodonFinder', function(GraphConfig,
                           $log) {
    /* The DNA string is 0 indexed, the dna positions are 1 indexed)
     * CodonFinder deals in string indexes (0 based)
     */
    var self = this;
    GraphConfig.CodonFinder = self;
    self.STOP="stop";
    self.START="start";
    self.ALTSTART="alt";
    self.index = [];
    self.indexes={start:[], alt:[], stop:[]};
    self.comIndexes={start:[], alt:[], stop:[]};

    self.codonRegexes = {
      start: /ATG|GTG/gim,
      altStart: /TTG|CTG|ATT/gim,
      end: /TAG|TAA/gim,
      nonMycoEnd: /TAG|TAA|TGA/gim
    };
    self.logging=false;
    self.classifyCodon = function(codon){
      var stopRe = GraphConfig.mycoplasma ?
          self.codonRegexes.end : self.codonRegexes.nonMycoEnd,
        startRe = self.codonRegexes.start,
        altstartRe = self.codonRegexes.altStart;
      //stops
      if(codon.search(stopRe)===0){
        return self.STOP;
      }
      // start
      else if(codon.search(startRe)===0){
        return self.START;
      }
      // nonstandard start
      else if(codon.search(altstartRe)===0){
        return self.ALTSTART;
      }
      return false;
    };
    self.classifyIdx = function(i, complement, codon) {
      if(i >= 9250 && i <= 9260) console.log('Classifying', arguments);
      if(arguments[2] === undefined){
        complment = GraphConfig.complement;
      }
      if(!codon){
        if(complement){
          codon = self.strRevCom( GraphConfig.ddnaString.substr(i-2, 3) );
        }else {
          codon = GraphConfig.ddnaString.substr(i, 3);
        }
      }
      var idxs = complement ? self.comIndexes : self.indexes ;
      var classification = self.classifyCodon(codon);
      if(classification){
        idxs[classification][i] = codon;
        if(GraphConfig.logging)$log.log('Classified Codon', codon, classification);
      }
      return classification;
    };
    self.strRevCom = function(s) {
      for (var i = s.length - 1, o = ''; i >= 0; i--) {
        var c = s[i];
        switch(c){
          case "A": o+= 'T'; break;
          case "T": o+= 'A'; break;
          case "G": o+= 'C'; break;
          case "C": o+= 'G'; break;
          default: o+= c;    break;
        }
      }
      return o;
    };
    self.reindex = _.debounce(function() {
      self._indexing = true;
      self.indexes={start:[], alt:[], stop:[]};
      self.comIndexes={start:[], alt:[], stop:[]};

      // I tried different methods including a hand coded matcher
      // and reversing the DNA string to match complements
      // this was the best I could come up with at about 1.2sec
      // to index
      var indexByRe = function() {
        var match,
            res = 'TAG|TAA|ATG|GTG|TTG|CTG|ATT|TGA',
            re      = new RegExp(res, 'img'),
            //rev   = /GAT|AAT|GTA|GTG|GTT|GTC|TTA|AGT/img;
            // reverse and complment the codes
            comre   = new RegExp(self.strRevCom(res),'img');
        while((match = re.exec(GraphConfig.ddnaString))){
          re.lastIndex = match.index+1;
          self.classifyIdx(match.index, false, match[0]);
        }
        while((match = comre.exec(GraphConfig.ddnaString))){
          comre.lastIndex = match.index+1;
          self.classifyIdx(match.index+2, true, self.strRevCom(match[0]));
        }
      };
      if(GraphConfig.logging) console.time('reindex re');
      indexByRe();
      if(GraphConfig.logging) console.timeEnd('reindex re');
      self._indexing = false;
    }, 50);
    self.isStart = function(cidx, complement){
      var idxs = complement ? self.comIndexes : self.indexes ;
      return idxs.start[cidx];
    };
    self.isAltStart = function(cidx, complement){
      var idxs = complement ? self.comIndexes : self.indexes ;
      return idxs.alt[cidx];
    };
    self.isStop = function(cidx, complement){
      var idxs = complement ? self.comIndexes : self.indexes ;
      var r = idxs.stop[cidx];
      if(GraphConfig.logging)
        $log.log('isStop Codon', cidx, complement, r);
      return r;
    };
    self.findNextStopCodon = function(cidx, complement){
      var idxs = complement ? self.comIndexes : self.indexes ;
      for(var i = cidx, rightbound=GraphConfig.ddnaString.length;
          i < rightbound ; i+=3){
        if(idxs.stop[i]) return i;
      }
      return null;
    };
    self.findPrevStopCodon = function(cidx, complement){
      var idxs = complement ? self.comIndexes : self.indexes ;
      for(var i = cidx-3; i >= 0 ; i-=3){
        if(idxs.stop[i]) return i;
      }
      return null;
    };
    self.findPrevStartCodon = function(cidx, complement, prevstopidx){
      var idxs = complement ? self.comIndexes : self.indexes ;
      var leftbound = prevstopidx||0;
      for(var i = cidx-3; i >= leftbound ; i-=3){
        if(idxs.start[i]) return i;
      }
      return null;
    };
    self.startOfNextCodon = function(idx){
      if(!idx || isNaN(idx)) return null;
      idx = Number(idx);
      return idx+3;
    };
    self.endOfThisCodon = function(idx){
      if(!idx || isNaN(idx)) return null;
      idx = Number(idx);
      return idx+2;
    };
  })
  .service('Fetcher', function($location, $http, GraphConfig,
                        StatusPoller, FETCH_BASE_URL, BASE_URL, PATH, PUBLIC_CONFIG_KEYS) {
    'use strict';
    var self = this;
    self.buildUrl = function(verb, params) {
      params = _.assign({},
                        $location.search(),
                        _.pick(GraphConfig, PUBLIC_CONFIG_KEYS),
                        params);
      return BASE_URL + '/' + verb + '/' + PATH + '?' + $.param(params);
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
    self.fetchFile = function(path) {
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
      if(!tid || tid.length === 0) {
        return $q.reject(new Error('Invalid task id'));
      }
      return poller(tid, $q.defer(), initialPollTime);
    };
  })
;
