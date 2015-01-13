describe('NProfiler', function(){
  'use strict';
  beforeEach(module('assets'));
  beforeEach(module('npact', function($provide) {
    $provide.value('$log', console);
    $provide.service('Fetcher', function($q, $log, $templateCache) {
      this.fetchFile = function(path) {
        $log.debug('Fetching ' + path);
        if(sampleConfig.ddna !== path) {
          throw new Error("fetchFile called with wrong path");
        }
        return $q.when(sampleDdnaFile);
      };
      this.nprofile = function() {
        fail("Shouldn't call nprofile");
      };
    });
  }));




  var NP, Pynpact, $digest, $timeout, sampleDdnaFile,
      sampleConfig = {
        length: 11141,
        ddna: 'sampleDdnaFile.ddna'
      };

  beforeEach(inject(function(NProfiler, _Pynpact_, $rootScope, _$timeout_, $templateCache) {
    NP = NProfiler;
    Pynpact = _Pynpact_;
    $digest = $rootScope.$digest;
    $timeout = _$timeout_;
    sampleDdnaFile = $templateCache.get('/js/test-data/' + sampleConfig.ddna);
  }));

  describe('start', function() {
    it('should fetch the ddna file', function(done) {
      NP.start(sampleConfig).then(function() {
        expect(NP.ddna).toBe(sampleDdnaFile);
        done();
      });
      $digest();
    });
  });

  describe('.stepSize', function() {
    it('should work', function() {
      var ss = NP.stepSize({startBase: 0, endBase: 10000});
      expect(ss).toBe(51);
    });
    it('should always be multiple of 3',function() {
      expect(NP.stepSize({startBase: 0, endBase: 9040}) % 3).toBe(0);
    });
  });

  describe('.slice', function() {
    beforeEach(inject(function($q) {
      NP.ddna = sampleDdnaFile;
      NP.fetching = $q.when(true);
    }));
    it('should be defined', function() {
      expect(NP.slice).toBeDefined();
    });
    it('should calculate the nprofile of a section', function(done) {
      // This is using the same parameters as the existing nprofile.c
      // and its output to compare
      var expected = [
        [101,38.8,41.8,28.4],
        [152,41.8,40.3,28.4],
        [203,46.3,40.3,29.9],
        [254,50.7,38.8,31.3],
        [305,50.7,40.3,32.8],
        [356,50.7,38.8,38.8],
        [407,47.8,38.8,40.3],
        [458,41.8,38.8,46.3],
        [509,43.3,31.3,50.7],
        [560,41.8,28.4,49.3]
      ];
      NP.slice({
        startBase: 0, endBase: 700,
        window: 201, step: 51,
        bases: ['C', 'G'],
        onPoint: function(coord, r, g, b) {
          var exp = expected.shift();
          expect(coord).toBe(exp[0]);
          expect(r).toBeCloseTo(exp[1], 1);
          expect(g).toBeCloseTo(exp[2], 1);
          expect(b).toBeCloseTo(exp[3], 1);
        }
      })
        .then(function() { expect(expected.length).toBe(0); })
        .finally(done);
      $timeout.flush();
    });
    it('Shouldn\'t fail with an end pass the length', function(done) {
      var count = 0;
      NP.slice({startBase: 0, endBase: 50000,
                step: 51, window: 201,
                onPoint: function() { count++; }})
        .then(function() {
          expect(count).toBe(215);
        }).finally(done);
      $timeout.flush();
    });
    it('returns a rejected promise when there\'s no ddna', function(done) {
      NP.ddna = null;
      var flag = false;
      NP.slice()
        .catch(function() { flag = true; })
        .finally(function() { expect(flag).toBe(true); done(); });
      $timeout.flush();
    });
    it('returns a rejected promise when there\'s no start/end base', function(done) {
      var flag = false;
      NP.slice()
        .catch(function() { flag = true; })
        .finally(function() { expect(flag).toBe(true); done(); });
      $timeout.flush();
    });
    it('returns a rejected promise when the window size is invalid', function(done) {
      var flag = false;
      NP.slice({startBase:1, endBase:50, window: 50})
        .catch(function() { flag = true; })
        .finally(function() { expect(flag).toBe(true); done(); });
      $timeout.flush();
    });

    var endBase = Math.max(2000,
                           Math.round(Math.random() * sampleConfig.length));
    it('has about 200 data points between 0 and ' + endBase, function(done) {
      var count = 0;
      NP.slice({startBase: 0, endBase: endBase, onPoint: function() { count++; }})
        .then(function() {
          expect(count).toBeGreaterThan(170);
          expect(count).toBeLessThan(230);
          done();
      });
      $timeout.flush();
    });
  });
});
