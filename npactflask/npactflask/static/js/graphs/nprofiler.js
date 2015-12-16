angular.module('npact')
  .service('DDNA', function(Fetcher, Pynpact, $log, $q, $timeout) {
    'use strict';
    var self = this;
    var defered = $q.defer();
    this.fetching = defered.promise;

    self.start = function(config) {
      $log.log("Starting DDNA fetch");
      self.config = config;
      defered.resolve(Fetcher.fetchFile(config[Pynpact.DDNA_FILE]).then(function(ddna) {
        $log.debug('Got back a ddna of length: ', ddna.length);
        return (self.ddna = ddna);
      }));
      return self.fetching;
    };

    self.sliceForExtract = function(extract) {
      return self.fetching.then(function(ddna) {
        // NProfiler.ddna is 0 indexed; the dna by convention (e.g. from the C or NCBI)
        // C/NCBI, is 1 indexed.
        var sliceStart = extract.start - 1;
        // The end index here is inclusive but array.slice isn't so we
        // don't need to subtract 1
        var sliceEnd = extract.end;
        return ddna.slice(sliceStart, sliceEnd);
      });
    };
  })

  .service('NProfiler', function(DDNA, GraphConfig, Utils, $log, $q, $timeout) {
    'use strict';
    var self = this;
    var defered = $q.defer();
    this.fetching = defered.promise;
    self.start = function(config) {
      self.config = config;
      defered.resolve(DDNA.start(config));
      return self.fetching;
    };



    self.defaultStepSize = function(len) {
      //Targetting about 200 datapoints for one graph line
      return Utils.floor3(len / 200);
    };

    self.defaultWindowSize = function(step) {
      //
      // This will be roughly this graph: http://www.wolframalpha.com/input/?i=plot+50*ln%28x%2F200%29+from+1000+to+100000
      // At  2000 it is 111,
      // At 10000 it is 198
      // At 50000 it is 276
      return Utils.floor3(50 * Math.log(step));
    };

    self.slice = function(opts) {
      var d = $q.defer();
      if(opts === undefined ||
         opts.startBase === undefined || opts.endBase === undefined ||
         opts.startBase < 0 || opts.startBase > opts.endBase) {
        d.reject("Incomplete/invalid startBase,endBase config for NProfiler.slice");
        return d.promise;
      }
      // The period is hardcoded to 3 (r,g,b) due to assumptions
      // throughout this graph system
      opts.period = 3;
      var len = GraphConfig.basesPerGraph;
      if(!opts.step)        { opts.step = self.defaultStepSize(len); }
      if(!opts.window)      { opts.window = self.defaultWindowSize(opts.step); }
      if(!opts.nucleotides) { opts.nucleotides = GraphConfig.nucleotides; }

      if(opts.window % opts.period) {
        d.reject("Window size (" + opts.window +
                 ") must be divisable by the period of frames (" +
                 opts.period + ")");
        return d.promise;
      }
      // slice asynchronously once we have the data
      return self.fetching.then(function(ddna) {
        return $timeout(function() {
          self._slice(ddna, opts);
          return opts;
        });
      });
    };

    this._slice = function(ddna, opts) {
      //This function assumes all the error checking and prep has been
      //done by the `slice` function

      // get local references to avoid dictionary lookups in the loop
      var end = Math.min(opts.endBase, ddna.length - 1),
          nucl = opts.nucleotides,
          period = opts.period,
          step = opts.step,
          onPoint = opts.onPoint,
          window = opts.window;

      var n = 0, // number of bases examined
          idx = opts.startBase, // current genome idx (0 based)
          box = new Array(window),
          profile = [0, 0, 0], // the current profile values
          normfactor = 100.0 / (window / period),
          win2 = Math.floor(window / 2);

      for(idx=opts.startBase; idx < end; idx++) {
        //is this base one we're searching for?
        if(nucl.indexOf(ddna[idx]) >= 0) {
          if(! box[n % window]) {
            //the flag wasn't set but we have a hit
            ++profile[idx % period];
            box[n % window] = true;
          }
        }
        else {
          if(box[n % window]) {
            --profile[idx % period];
            box[n % window] = false;
          }
        }
        n++;
        if(n >= window && ((n-window) % step) === 0) {
          //idx+1: we 1 index our bases traditionally
          //win2: we want report the coordinate at the center of the window
          onPoint((idx+1) - win2,
                  normfactor * profile[0],
                  normfactor * profile[1],
                  normfactor * profile[2]);
        }
      }
    };
  })
;
