describe('Graphs', function(){
  'use strict';

  beforeEach(module('npact', function($provide) {
    $provide.value('$log', console);
  }));
  beforeEach(module('assets'));

  var ng = {}, $scope, $compile, $q, $httpBackend, Err;
  beforeEach(inject(function (_$rootScope_, _$compile_, _$q_, _$httpBackend_, _$timeout_, _Err_) {
    $scope = _$rootScope_;
    $compile = _$compile_;
    $q = _$q_;
    $httpBackend = _$httpBackend_;
    ng.$timeout = _$timeout_;
    Err = _Err_;
  }));

  describe("Enums", function() {
    // Setup the custom matcher
    beforeEach(function() {
      jasmine.addMatchers({
        toHaveUniqueValues: function(util, customEqualityTesters) {
          return {
            compare: function(actual, expected) {
              var vals = _.values(actual);
              return {
                pass: vals.length === _.uniq(vals).length,
                message: "Expected the values of " + actual + " to be unique."
              };
            }
          };
        }
      });
    });
    it('Evt has unique values', inject(function(Evt) {
      expect(Evt).toHaveUniqueValues();
    }));
    it('Pynpact has unique values', inject(function(Pynpact) {
      expect(Pynpact).toHaveUniqueValues();
    }));
  });


  function make(html){
    var el = $compile(html)($scope);
    $scope.$digest();
    return el;
  }

  describe('Utils', function(){
    var U;

    beforeEach(inject(function(Utils){
      U = Utils;
    }));

    describe('.orderOfMagnitude', function() {
      it('handles no offset', function(){
        expect(U.orderOfMagnitude(1012)).toBe(1000);
        expect(U.orderOfMagnitude(1012, -1)).toBe(100);
      });
      it('handles negative offset', function(){
        expect(U.orderOfMagnitude(1012, -1)).toBe(100);
      });
      it('handles positive offset', function(){
        expect(U.orderOfMagnitude(1012, 1)).toBe(10000);
      });
      it('shouldn\'t jerk badly', function() {
        expect(U.orderOfMagnitude(31012, -1)).toBe(3000);
        expect(U.orderOfMagnitude(30012, -1)).toBe(3000);
        expect(U.orderOfMagnitude(33012, -1)).toBe(3000);
        expect(U.orderOfMagnitude(39012, -1)).toBe(4000);
      });
    });

    describe('.extendByPage', function() {
      it('works', function() {
        var src = _.range(100),
            dst = _.take(src, 5);
        U.extendByPage(src, dst, 5);
        expect(dst).toEqual(_.range(10));
        U.extendByPage(src, dst, 5);
        expect(dst).toEqual(_.range(15));
      });
    });

  });

  describe('GraphingCalculator', function() {
    var GC;
    beforeEach(inject(function(GraphingCalculator){
      GC = GraphingCalculator;
    }));

    describe('.partition', function() {
      it('partitions', function(){
        var p = GC.partition({startBase: 0, endBase: 450,
                              basesPerGraph: 100, offset: 0});
        expect(p).toEqual([0, 100, 200, 300, 400]);
      });

      it('handles a positive offset', function() {
        var p = GC.partition({startBase: 0, endBase: 450,
                              basesPerGraph: 100, offset: 10});
        expect(p).toEqual([ 10, 110, 210, 310, 410]);
      });
      it('handles a negative offset ', function() {
        var p = GC.partition({startBase: 0, endBase: 450,
                              basesPerGraph: 100, offset: -10});
        expect(p).toEqual([ -10, 90, 190, 290, 390]);
      });
      it('handles a stupid negative offset', function() {
        var p = GC.partition({startBase: 0, endBase: 450,
                              basesPerGraph: 100, offset: -110});
        expect(p).toEqual([ -110, -10, 90, 190, 290, 390]);
      });
    });

    describe('trackSizeCalc', function() {
      var testExtractsTrack = {text: 'test', type: 'extracts', active: true, height: 30};
      var testHitsTrack = {text: 'test3', type: 'hits', active: true, height: 20};
      it('with one extract', function() {
        var tracks = [testExtractsTrack];
        expect(GC.trackSizeCalc([testExtractsTrack])).toEqual(30);
        expect(tracks[0].y).toEqual(0);
      });
      it('with many', function() {
        var tracks = [testExtractsTrack, _.clone(testExtractsTrack), testHitsTrack];
        expect(GC.trackSizeCalc(tracks)).toEqual(80);
        expect(tracks[0].y).toEqual(0);
        expect(tracks[1].y).toEqual(30);
        expect(tracks[2].y).toEqual(60);
      });
    });

    describe('metrics', function() {
      var opts = {
        width: 250,
        axisTitle: '% GC',
        profileTicks: 5,
        startBase: 0, endBase: 100,
        tracks: [{height: 5}]
      };
      beforeEach(inject(function(npactConstants) {
        npactConstants.graphStyle.paddingUnit = 5;
        npactConstants.graphStyle.leftPadding = 50;
        npactConstants.graphStyle.tickLength = 5;
      }));

      it('calculates basic metrics', function() {
        var m = GC.chart(opts);
        expect(m).toBeDefined();
        expect(m.height).toEqual(106);
        expect(m.graph.x).toEqual(50);
        expect(m.graph.y).toEqual(10);
      });
    });

    it('can align rectangles', function(){
      var pos = GC.alignRectangles(
        {x: 0, y: 0, width: 100, height: 100},
        {width: 10, height: 10});
      expect(pos).toEqual({x: 45, y: 45});
    });

    it('calculates stops for axes', function(){
      expect(GC.stops(0, 10000)).toEqual({
        interval:1000,
        stops:[0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
      });

      expect(GC.stops(0, 1000)).toEqual({
        interval:100,
        stops: [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
      });

      expect(GC.stops(1000, 2000)).toEqual({
        interval:100,
        stops: [1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000]
      });
    });

    it('calculates stops for offset axes', function(){
      expect(GC.stops(50, 10050)).toEqual({
        interval:1000,
        stops:[0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
      });
    });
    it('calculates stops for given length', function(){
      expect(GC.stops(0, 3000, 10000)).toEqual({
        interval: 1000,
        stops:[0, 1000, 2000, 3000]
      });
    });

    function expectZoom(args, offset, basesPerGraph){
      expect(GC.zoom({startBase: args[0],
                      zoomOnPct: args[1],
                      basesPerGraph: args[2],
                      offset: args[3],
                      zoomingOut: args[4]
                     }))
        .toEqual({offset: offset,
                  basesPerGraph: basesPerGraph});
    }

    it('can zoom on the first row', function(){
      expectZoom([0, 0.5, 10000, 0], 2500, 5000);
      expectZoom([0, 0.1, 10000, 0], 500, 5000);
      expectZoom([900, 0.1, 10000, 0], 1400, 5000);
      expectZoom([0, 0.9, 10000, 0], 4500, 5000);
    });

    it('can zoom on the higher rows', function(){
      expectZoom([10000, 0.5, 10000, 0], 7500, 5000);
      expectZoom([20000, 0.5, 10000, 0], 12500, 5000);
    });
  });

  describe('ExtractParser', function(){
    var EXTRACT = ['H-51*G complement(57104..57904)',
                   'H-53*A complement(58013..59380)',
                   'H-64-C 71945..72100',
                   'G-125-G 88111..88275',
                   'G-124*t 88544..88906',
                   'H-102-C 112734..112931',
                   'H-103*a 112864..113175',
                   'H-115-a complement(129294..>129599)',
                   'XAR_RS17080 complement(<140014..140226)']
          .join('\n'),
        EP, result;

    beforeEach(inject(function(ExtractParser){
      EP = ExtractParser;
      result = EP.parse(EXTRACT);
    }));

    it('can parse many', function(){
      expect(result.length).toBe(9);
    });

    it('can parse a line', function(){
      expect(result[2])
        .toEqual({
          name: 'H-64-C',
          start: 71945, end: 72100,
          complement:0, approximate: false,
          phase: (72100 - 1) % 3
        });
    });
    it('can parse a complement', function(){
      expect(result[0])
        .toEqual({
          name: 'H-51*G',
          start: 57104, end: 57904,
          complement:1, approximate: false,
          phase: (57104 - 1) % 3
        });
    });

    it('can parse a approximates', function(){
      expect(result[7].approximate).toBe(true);
      expect(result[8].approximate).toBe(true);

    });

    it('double parsing is a NOOP', function(){
      expect(EP.parse(result)).toEqual(result);
    });
  });


  describe('GraphConfig', function() {
    var G;
    beforeEach(inject(function(GraphConfig){
      G = GraphConfig;
    }));

    it('calculates profile title', function() {
      G.nucleotides = ['C', 'G'];
      expect(G.profileTitle()).toBe('% CG');
      G.nucleotides = ['A', 'C', 'G'];
      expect(G.profileTitle()).toBe('% ACG');
    });

    describe('.loadTrack', function() {
      it('loads', function() {
        expect(G.activeTracks().length).toBe(0);
        expect(G.findTrack('test')).toBeFalsy();
        G.loadTrack({name: 'test', weight: 10});
        expect(G.findTrack('test')).toBeTruthy();
      });
      it('replaces on duplicate', function() {
        G.loadTrack({name: 'test', type:'extracts', weight: 5, active: false});
        expect(G.findTrack('test')).toBeTruthy();
        G.loadTrack({name: 'test', weight: 10, active: true});
        expect(G.tracks.length).toBe(1);
        expect(G.findTrack('test').weight).toBe(10);
      });
      it('keeps "hits" as the last track', function() {
        G.loadTrack({name: 'test', type: 'extracts', weight: 0, active: true});
        G.loadTrack({name: 'testh', type: 'hits', weight: 100, active: true});
        G.loadTrack({name: 'inactive', active: false, weight: 1000});
        expect(G.activeTracks().length).toBe(2);
        expect(_.last(G.activeTracks()).name).toBe('testh');
        G.loadTrack({name: 'test3', type: 'extracts', weight: 10});
        expect(_.last(G.activeTracks()).name).toBe('testh');
      });
    });

  });
});
