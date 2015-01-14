angular.module('npact')
  .service('NProfiler', function(Fetcher, Pynpact, $log, GraphConfig, $q, $timeout) {
    'use strict';
    var self = this;
    self.start = function(config) {
      $log.log("Starting nprofiler", self);
      self.config = config;

      self.fetching = Fetcher.fetchFile(config[Pynpact.DDNA_FILE]).then(function(ddna) {
        $log.debug('Got back a ddna of length: ', ddna.length);
        self.ddna = ddna;
      });
      return self.fetching;
    };

    var round3 = function(num) {
      //round to multiple of 3
      return num + 1.5 - (num + 1.5) % 3;
    };

    self.defaultStepSize = function(len) {
      //Targetting about 200 datapoints for one graph line
      return round3(len / 200);
    };

    self.defaultWindowSize = function(step) {
      //
      // This will be roughly this graph: http://www.wolframalpha.com/input/?i=plot+50*ln%28x%2F200%29+from+1000+to+100000
      // At  2000 it is 111,
      // At 10000 it is 198
      // At 50000 it is 276
      return round3(50 * Math.log(step));
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
      var len = opts.endBase - opts.startBase;
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
      return self.fetching.then(function() {
        return $timeout(function() {
          self._slice(opts);
          return opts;
        });
      });
    };

    self._slice = function(opts) {
      //This function assumes all the error checking and prep has been
      //done by the `slice` function

      // get local references to avoid dictionary lookups in the loop
      var t1 = new Date(),
          ddna = self.ddna,
          end = Math.min(opts.endBase, ddna.length - 1),
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
      $log.debug('Slice @', opts.startBase, 'step,win', step,window,
                 'took', new Date() - t1);
    };
  })
;
