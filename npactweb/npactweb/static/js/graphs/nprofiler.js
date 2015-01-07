angular.module('npact')
  .service('NProfiler', function(Fetcher, Pynpact, $log) {
    var self = this;
    self.start = function(config) {
      $log.log("Starting nprofiler");
      self.config = config;
      return Fetcher.fetchFile(config[Pynpact.DDNA_FILE]).then(function(ddna) {
        $log.debug('Got back a ddna of length: ', ddna.length);
        self.ddna = ddna;
      });
    };

    self.slice = function(opts) {
      if(opts.startBase === undefined || opts.endBase === undefined) {
        throw new Error('Incomplete options to NProfiler.startBase');
      }

    };
  })
;
