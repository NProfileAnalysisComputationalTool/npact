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

  
  function Utils($q, $timeout, $log){

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
    
    return {
      STOP_ITERATING: STOP_ITERATING,
      orderOfMagnitude: orderOfMagnitude,
      forEachAsync: forEachAsync
    };
  }

  angular.module('npact')
    .factory('Utils', Utils)
  ;
  
}())
