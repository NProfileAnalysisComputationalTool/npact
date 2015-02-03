describe('ITrackReader', function() {
  'use strict';
  beforeEach(module('npact', function($provide) {
    $provide.value('$log', console);
  }));
  beforeEach(module('assets'));

  var $timeout, Err;
  beforeEach(inject(function(_$timeout_, _Err_) {
    $timeout = _$timeout_;
    Err = _Err_;
  }));
  var extract = ['H-51*G complement(57104..57904)',
                 'H-53*A complement(58013..59380)',
                 'H-64-C 71945..72100',
                 'G-125-G 88111..88275',
                 'G-124*t 88544..88906',
                 'H-102-C 112734..112931',
                 'H-103*a 112864..113175',
                 'H-115-a complement(129294..129599)',
                 'H-116-G complement(131053..131163)'].join('\n');

  var parsedData;
  beforeEach(inject(function(ExtractParser) {
    parsedData = ExtractParser.parse(extract);
  }));

  _.map(
    ['TrackRTreeIndex', 'TrackNaiveIndex', 'TrackBinSearchIndex'],
    function (TrackIndexName) {
      describe(TrackIndexName, function() {
        var TrackIndex, tIndex;
        beforeEach(inject([TrackIndexName, function(tr) {
          TrackIndex = tr;
          tIndex = new TrackIndex(parsedData);
        }]));
        it('instantiates', function() {
          expect(tIndex).toEqual(jasmine.any(Function));
        });

        describe('slices', function() {
          it('slices', function() {
            tIndex({name:'test', startBase:72000, endBase:88600})
              .then(function(slice) {
                expect(slice)
                  .toContain({start: 71945, end: 72100, complement: 0, name: 'H-64-C', phase: 0, approximate: false});
                expect(slice)
                  .toContain({start: 88544, end: 88906, complement: 0, name: 'G-124*t', phase: 0, approximate: false});
                expect(slice)
                  .toContain({start: 88111, end: 88275, complement: 0, name: 'G-125-G', phase: 2, approximate: false});
              });
            $timeout.flush();
          });
          it('returns empty array on no matches', function() {
            tIndex({name:'test', startBase:0, endBase:100})
              .then(function(slice) {
                expect(slice).toEqual([]);
              });
            $timeout.flush();
          });
        });
      });
    });

  describe('Track', function() {
    var track;
    beforeEach(inject(function(Track) {
      track = new Track('test', extract, 'extracts', 0);
    }));
    it('should load and parse', function() {
      expect(track).toBeDefined();
      expect(track.name).toEqual('test');
      expect(track.weight).toEqual(0);
      expect(track.type).toEqual('extracts');
      expect(track.height).toBeTruthy();
    });

    it('should throw on unknown track type', function() {
      expect(function() { new Track('test', null, 'unknown'); }).toThrow();
    });

    it('should be active by default', function() {
      expect(track.active).toBe(true);
    });

    it('should be sliceable', function(done) {
      track.slice({name:'test', startBase:72000, endBase:88600}).then(function(results) {
        expect(results.length).toEqual(3);
      }).then(function() {
        // Ensure it works more than once
        track.slice({name:'test', startBase:72000, endBase:88600}).then(function(results) {
          expect(results.length).toEqual(3);
          done();
        });
      });
      $timeout.flush();
    });
  });
});
