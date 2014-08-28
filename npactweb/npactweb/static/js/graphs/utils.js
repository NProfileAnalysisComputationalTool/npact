(function(){

  var log10 = Math.log(10);

  /**
   * find the nearest order of magnitude to `x`
   *
   * orderOfMagnitude(1012) => 1000
   * orderOfMagnitude(1012, -1) => 100
   *
   * @param {Number} x
   * @param {Number} [exponent] shift the result by this many digits
   */
  function orderOfMagnitude(x, exponent){
    return Math.pow(10, Math.round(Math.log(x) / log10) + (exponent || 0));
  }

  // sigil value for quitting forEachAsync
  var STOP_ITERATING = {};


  function Utils($q, $timeout, $log, $interval){
    // public interface
    return {
      STOP_ITERATING: STOP_ITERATING,
      orderOfMagnitude: orderOfMagnitude,
      forEachAsync: forEachAsync,
      widthAvailable: widthAvailable
    };

    /**
     * polls for when the element has a width
     *
     * @param {Element} element - what element to scan
     * @returns {Promise} - the width of the element
     */
    function widthAvailable(element){
      var d = $q.defer(),
	  p = d.promise,
	  t1 = new Date(),
	  task = $interval(checkForWidth, 100);
      // stop polling once we've found a value
      p.then(function(){$interval.cancel(task);});
      return p;

      function checkForWidth(){
	var w = element.width();
	// indicate progress via the promise
	d.notify(new Date() - t1);
	if(w > 0){
	  d.resolve(w);
	}
      }
    };

    /**
     * Loop over the list from back to front, asyncronously
     *
     * stop iteration early with `throw Utils.STOP_ITERATING`
     *
     * @returns {Promise} when the iteration is complete
     */
    function forEachAsync(list, fn){
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
	  $log.log('iterated over', len, 'items in', new Date() - t1, 'ms');
	});

	function iterate(){
	  var t1 = new Date();
	  try{
	    for(var c = opts.batchSize; idx < len && c >= 0; idx++, c--){
	      fn(list[idx]);
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
	  if(idx === len){
	    d.notify({batches:batches, idx:idx,
		      elapsed:elapsed,
		      batchSize:opts.batchSize});
	    // if the loop was very fast, double the batch size, if slow, halve it
	    opts.batchSize = (elapsed < 10)
	      ? opts.batchSize << 1
	      : opts.batchSize >> 1;
	    // queue the next batch
	    $timeout(iterate, opts.delay, false);
	  }else{
	    d.resolve(batches);
	  }
	}
	// start it off soon
	$timeout(iterate, opts.delay, false);
	return p;
    }
  }

  /**
   * parses the Extract format
   *
   * The grammar is roughly:
   *   <extract>    = <name> <address>
   *   <address>    = <complement> | <range>
   *   <complement> = complement(<range>)
   *   <range>      = <base>..<base>
   *   <base>       = any integer
   */
  function ExtractParser(){
    var re = /([^ ]+) (complement\()?(\d+)\.\.(\d+)/,
	// some indexes for where our regex groups will show
	NAME = 1, COMP = 2, START = 3, END = 4;

    return {
      parse:parse
    };

    /**
     * parses many lines
     *
     * @param {string} text - multi-line text
     * @returns {array} extract objects
     */
    function parse(text){
      if (_.isObject(text)) return text;
      var lines = text.split('\n');
      return lines.map(parseLine);
    }

    /**
     * parses one line
     *
     * @param {string} line - single extract like
     * @returns {array} extract objects
     */
    function parseLine(line){
      var parts = re.exec(line);
      return {
	start: parseInt(parts[START]),
	end: parseInt(parts[END]),
	complement: parts[COMP] ? 1 : 0,
	name: parts[NAME]
      };
    }
  }

  angular.module('npact')
    .factory('Utils', Utils)
    .factory('ExtractParser', ExtractParser)
  ;

}())
