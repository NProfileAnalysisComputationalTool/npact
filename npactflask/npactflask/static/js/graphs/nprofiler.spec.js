describe('NProfiler', function(){
  'use strict';

  beforeEach(module('assets'));
  beforeEach(module('npact', function($provide) {
//    $provide.value('$log', console);
    $provide.service('Fetcher', function($q, $log) {
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


  var DDNA, NP, Pynpact, $digest, $timeout, sampleDdnaFile,
      sampleConfig = {
        length: 11141,
        ddna: 'sampleDdnaFile.ddna'
      };

  beforeEach(inject(function(_DDNA_, NProfiler, _Pynpact_, $rootScope, _$timeout_, $templateCache) {
    NP = NProfiler;
    DDNA = _DDNA_;
    Pynpact = _Pynpact_;
    $digest = $rootScope.$digest;
    $timeout = _$timeout_;
    sampleDdnaFile = $templateCache.get('/js/test-data/' + sampleConfig.ddna);
  }));

  describe('DDNA start', function() {
    it('should fetch the ddna file', function(done) {
      DDNA.start(sampleConfig).then(function(ddna) {
        expect(ddna).toBe(sampleDdnaFile);
        expect(DDNA.ddna).toBe(sampleDdnaFile);
        expect(DDNA.fetching).toBeTruthy();
        done();
      });
      $digest();
    });
  });
  describe('NProfiler start', function() {
    it('should call DDNA.start preserving the old interface', function(done) {
      NP.start(sampleConfig).then(function(ddna) {
        expect(ddna).toBe(sampleDdnaFile);
        expect(NP.fetching).toBeTruthy();
        done();
      });
      $digest();
    });
  });

  describe('.defaultStepSize', function() {
    it('should work', function() {
      var ss = NP.defaultStepSize(10000);
      expect(ss).toBe(51);
    });
    it('should always be multiple of 3',function() {
      expect(NP.defaultStepSize(9040) % 3).toBe(0);
      expect(NP.defaultStepSize(10000) % 3).toBe(0);
      expect(NP.defaultStepSize(20000) % 3).toBe(0);
      expect(NP.defaultStepSize(50000) % 3).toBe(0);
    });
  });

  describe('DDNA.slice', function() {
    beforeEach(inject(function($q) {
      DDNA.fetching = NP.fetching = $q.when(sampleDdnaFile);
    }));
    it('should slice length', function(done) {
      DDNA.sliceForExtract({start:1, end: 9}).then(function(slice) {
        expect(slice.length).toEqual(9);
        expect(slice[0]).toEqual(sampleDdnaFile[0]);
        done();
      });
      $digest();
    });
  });
  describe('NP.slice', function() {
    beforeEach(inject(function($q) {
      DDNA.fetching = NP.fetching = $q.when(sampleDdnaFile);
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
  });
});
