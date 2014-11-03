angular.module('npact')
  .service('Utils', function($q, $timeout, $log) {
    'use strict';
    // sigil value for quitting forEachAsync
    var STOP_ITERATING = {},
        log10 = Math.log(10);
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
      return Math.pow(10, Math.round(Math.log(x) / log10) + (exponent || 0));
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
          len = list.length,
          t1 = new Date();

      p.then(function(){
        $log.log('iterated over', len,
                 'items in', new Date() - t1, 'ms, batchcount: ', batches);
      });

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
            d.reject(e);
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
      return lines.map(parseLine);
    };
    /**
     * parses many lines asynchronously
     *
     * @param {string} text - multi-line text
     * @param {function} parseLine - callback to parse a line
     * @returns {Promise} for array of parsed objects
     */
    this.parseAsync = function(text, parseLine){
      if (!_.isString(text)) { return $q.when(text); }

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
    var EXTRACT_REGEX = /([^ ]+) (complement\()?(\d+)\.\.(\d+)/,
        // some indexes for where our regex groups will show
        NAME = 1, COMP = 2, START = 3, END = 4;

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
            name: parts[NAME]
          },
          phaseCoordinate = res.complement ? res.start : res.end;
      res.phase = (phaseCoordinate - 1) % 3;
      return res;
    }

    return ParserFactory.create(parseExtract);
  })
;
