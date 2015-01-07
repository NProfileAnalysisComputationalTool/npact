angular.module('npact')
  .service('NProfiler', function(Fetcher, Pynpact, $log, GraphConfig, $q, $timeout) {
    'use strict';
    var self = this;
    self.start = function(config) {
      $log.log("Starting nprofiler", self);
      self.config = config;

      return self.fetching = Fetcher.fetchFile(config[Pynpact.DDNA_FILE]).then(function(ddna) {
        $log.debug('Got back a ddna of length: ', ddna.length);
        self.ddna = ddna;
      });
    };

    self.avgParams = function(len) {
      return {
        window: 201,
        step: 51,
        nucleotides: GraphConfig.nucleotides
      };
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
      opts = angular.extend(self.avgParams(), opts);

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
          onPoint((idx+1) - win2,
                  normfactor * profile[0],
                  normfactor * profile[1],
                  normfactor * profile[2]);
        }
      }
      $log.debug('Slice took ', new Date() - t1);
    };
  })
;
