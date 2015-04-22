describe('ITrackReader', function() {
  'use strict';
  beforeEach(module('npact', function($provide) {
  //  $provide.value('$log', console);
  }));
  beforeEach(module('assets'));

  var $timeout, Err;
  beforeEach(inject(function(_$timeout_, _Err_) {
    $timeout = _$timeout_;
    Err = _Err_;
  }));
  var extract = ['Adeh_0001 22..1395',
                 'H-51*G complement(57104..57904)',
                 'H-53*A complement(58013..59380)',
                 'H-64-C 71945..72100',
                 'H complement(78698..81685)',
                 'H 81578..81631',
                 'H 81894..83693',
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

  //Sanity check test
  describe('parsedData', function() {
    it('should have the right number of elements', function() {
      expect(parsedData.length).toEqual(13);
    });
  });

  _.map(
    ['TrackRTreeIndex', 'TrackNaiveIndex', 'TrackBinSearchIndex'],
    function (TrackIndexName) {
      describe(TrackIndexName, function() {
        var TrackIndex, tIndex;
        beforeEach(inject([TrackIndexName, function(tr) {
          TrackIndex = tr;
          tIndex = new TrackIndex(parsedData);
        }]));
        it('instantiates', function() { expect(tIndex).toEqual(jasmine.any(Function)); });

        describe('slices', function() {
          it('slices', function() {
            tIndex({startBase: 0, endBase: 1200}).then(function(slice) {
              expect(slice).toContain({start: 22, end: 1395, complement: 0, name: 'Adeh_0001', phase: 2, approximate: false});
            });
            tIndex({startBase:72000, endBase:88600})
              .then(function(slice) {
                expect(slice)
                  .toContain({start: 71945, end: 72100, complement: 0, name: 'H-64-C', phase: 0, approximate: false});
                expect(slice)
                  .toContain({start: 88544, end: 88906, complement: 0, name: 'G-124*t', phase: 0, approximate: false});
                expect(slice)
                  .toContain({start: 88111, end: 88275, complement: 0, name: 'G-125-G', phase: 2, approximate: false});
              });
          });
          it('returns empty array on no matches', function() {
            tIndex({startBase:0, endBase:10})
              .then(function(slice) { expect(slice).toEqual([]); });
            tIndex({startBase: 132000, endBase: 135000 })
              .then(function(slice) { expect(slice).toEqual([]); });
          });
          it('doesn\'t match when abutting a gene border from the outside', function() {
            //We have gene at 22..1395
            tIndex({startBase:0, endBase: 22})
              .then(function(slice) { expect(slice).toEqual([]); });
            tIndex({startBase:1395, endBase: 1400})
              .then(function(slice) { expect(slice).toEqual([]); });
          });
          it('does match when abutting a gene border from the inside', function() {
            //We have gene at 22..1395
            tIndex({startBase:22, endBase: 50}).then(function(slice) {
              expect(slice).toContain({start: 22, end: 1395, complement: 0, name: 'Adeh_0001', phase: 2, approximate: false});
              tIndex({startBase:16, endBase: 1395}).then(function(slice) {
                expect(slice).toContain({start: 22, end: 1395, complement: 0, name: 'Adeh_0001', phase: 2, approximate: false});
              });
            });
          });

          it('should match this one complement that it is not right now', function() {
            // `H complement(78698..81685)` is hidden sometimes
            // because it contains `H 81578..81631`
            tIndex({startBase: 73000, endBase: 91000}).then(function(slice) {
              expect(slice).toContain(
                { start: 78698, end: 81685, complement: 1, name: 'H', approximate: false, phase: 1 }
              );
            });
          });
          afterEach(function() { $timeout.flush(); });
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
        expect(results.length).toEqual(6);
      }).then(function() {
        // Ensure it works more than once
        track.slice({name:'test', startBase:72000, endBase:88600}).then(function(results) {
          expect(results.length).toEqual(6);
          done();
        });
      });
      $timeout.flush();
    });
  });
});
