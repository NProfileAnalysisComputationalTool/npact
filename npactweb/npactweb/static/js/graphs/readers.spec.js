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

  _.map(
    ['TrackReaderNaive', 'TrackReaderRTreeIndex', 'TrackBinarySearch'],
    function testATrackReader(TrackReaderName) {
      describe(TrackReaderName, function() {
        var TrackReader;
        beforeEach(inject([TrackReaderName, function(tr) { TrackReader = tr; }]));

        describe('.load',function() {
          it('loads', function(done) {
            TrackReader.load('test', extract).then(done);
            expect(TrackReader.tracks.test).toBeDefined();
            $timeout.flush();
          });
          it('throws on bad parse', function() {
            expect(function() {
              TrackReader.load('test', null);
            }).toThrow();
          });
        });
        describe('.slice',function() {
          beforeEach(function() {
            TrackReader.load('test', extract);
            $timeout.flush();
          });

          it('slices', function() {
            TrackReader.slice({name:'test', startBase:72000, endBase:88600})
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
          it('throws on bad names', function() {
            expect(function() {
              TrackReader.slice({name:'test2', startBase:0, endBase:0});
            }).toThrow(new Err.TrackNotFound());
          });
          it('returns empty array on no matches', function() {
            TrackReader.slice({name:'test', startBase:0, endBase:100})
              .then(function(slice) {
                expect(slice).toEqual([]);
              });
            $timeout.flush();
          });
        });
      });
    });
});
