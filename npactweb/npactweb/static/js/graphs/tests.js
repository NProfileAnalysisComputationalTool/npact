describe('Graphs', function(){
  'use strict';

  beforeEach(module('npact'));
  beforeEach(module('templates-main'));

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

  describe('GraphingCalculator', function(){
    var GC, opts = {
      height:150,
      width:255,
      leftPadding: 50,
      rightPadding: 5,
      axisLabelFontsize:10,
      axisTitleFontsize: 20,
      axisTitle: '% GC',
      profileTicks:5,
      startBase: 0, endBase: 100,
      headerY: 5
    };
    beforeEach(inject(function(GraphingCalculator){
      GC = GraphingCalculator;
    }));

    it('calculates graph area', function(){
      var m = GC.chart(opts);
      expect(m.graph).toEqual({
        x: 50,
        y: 10,
        h: 120,
        w: 200
      });
    });

    it('calculates y-axis title', function(){
      var m = GC.chart(opts);
      expect(m.yaxis.titleBox).toEqual(
        {x: 0, y: 10, width: 40, height: 120}
      );
    });

    it('calculates x-axis ticks', function(){

      var ax = GC.xaxis(opts);
      expect(ax.ticks).toEqual([
        {x: 0, y: 0, x2: 0, y2: 5},
        {x: 10, y: 0, x2: 10, y2: 5},
        {x: 20, y: 0, x2: 20, y2: 5},
        {x: 30, y: 0, x2: 30, y2: 5},
        {x: 40, y: 0, x2: 40, y2: 5},
        {x: 50, y: 0, x2: 50, y2: 5},
        {x: 60, y: 0, x2: 60, y2: 5},
        {x: 70, y: 0, x2: 70, y2: 5},
        {x: 80, y: 0, x2: 80, y2: 5},
        {x: 90, y: 0, x2: 90, y2: 5},
        {x: 100, y: 0, x2: 100, y2: 5},
        {x: 110, y: 0, x2: 110, y2: 5}

      ]);
    });

    it('can align rectangles', function(){
      var pos = GC.alignRectangles(
        {x: 0, y: 0, width: 100, height: 100},
        {width: 10, height: 10});
      expect(pos).toEqual({x: 45, y: 45});
    });

    it('calculates stops for axes', function(){
      expect(GC.stops(0,10000)).toEqual({
        interval:1000,
        stops:[0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000]
      });

      expect(GC.stops(0,1000)).toEqual({
        interval:100,
        stops: [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100]
      });

      expect(GC.stops(1000,2000)).toEqual({
        interval:100,
        stops: [900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100]
      });
    });

    it('calculates stops for offset axes', function(){
      expect(GC.stops(50,10050)).toEqual({
        interval:1000,
        stops:[0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000]
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


  describe('TrackReader', function() {
    var T,
        extract = ['H-51*G complement(57104..57904)',
                   'H-53*A complement(58013..59380)',
                   'H-64-C 71945..72100',
                   'G-125-G 88111..88275',
                   'G-124*t 88544..88906',
                   'H-102-C 112734..112931',
                   'H-103*a 112864..113175',
                   'H-115-a complement(129294..129599)',
                   'H-116-G complement(131053..131163)']
          .join('\n')
    ;

    beforeEach(inject(function(TrackReader){
      T = TrackReader;
    }));

    describe('.load',function() {
      it('loads', function() {
        T.load('test', extract).then(function() {
          expect(T.tracks.test).toBeDefined();
        });
        ng.$timeout.flush();
      });
      it('throws on duplicate names', function() {
        T.load('test', extract);
        ng.$timeout.flush();
        expect(function() {
          T.load('test', '');
        }).toThrow(new Err.TrackAlreadyDefined());
      });
    });

    describe('.slice',function() {
      beforeEach(function() {
        T.load('test', extract);
        ng.$timeout.flush();
      });

      it('slices', function() {
        T.slice({name:'test', startBase:72000, endBase:88600})
          .then(function(slice) {
            expect(slice)
              .toContain({start: 71945, end: 72100, complement: 0, name: 'H-64-C', phase: 0, approximate: false});
            expect(slice)
              .toContain({start: 88544, end: 88906, complement: 0, name: 'G-124*t', phase: 0, approximate: false});
            expect(slice)
              .toContain({start: 88111, end: 88275, complement: 0, name: 'G-125-G', phase: 2, approximate: false});
          });
        ng.$timeout.flush();
      });
      it('throws on bad names', function() {
        expect(function() {
          T.slice({name:'test2', startBase:0, endBase:0});
        }).toThrow(new Err.TrackNotFound());
      });
      it('returns empty array on no matches', function() {
        T.slice({name:'test', startBase:0, endBase:100})
          .then(function(slice) {
            expect(slice).toEqual([]);
          });
        ng.$timeout.flush();
      });
    });

  });

  describe('GraphConfig', function() {
    var G;
    beforeEach(inject(function(GraphConfig){
      G = GraphConfig;
    }));

    describe('.loadTrack', function() {
      it('loads', function() {
        expect(G.tracks).toEqual([]);
        expect(G.hasTrack('test')).toBe(false);
        G.loadTrack('test', 'extracts');
        expect(G.hasTrack('test')).toBe(true);
      });
      it('throws on duplicate', function() {
        G.loadTrack('test', 'extracts');
        expect(function() {
          G.loadTrack('test', 'whatever');
        }).toThrow(new Err.TrackAlreadyDefined());
      });
    });

    describe('.headerSpec', function() {
      it('with one extract', function() {
        G.loadTrack('test', 'extracts');
        expect(G.headerSpec()).toEqual({
          headerY: 35,
          headers: [ { text: 'test', lineType: 'extracts', y: 5, height: 30 }]
        });
      });
      it('with many', function() {
        G.loadTrack('test', 'extracts');
        G.loadTrack('test2', 'extracts');
        G.loadTrack('test3', 'hits');
        expect(G.headerSpec()).toEqual({
          headerY: 85,
          headers: [
            { text: 'test', lineType: 'extracts', y: 5, height: 30 },
            { text: 'test2', lineType: 'extracts', y: 35, height: 30 },
            { text: 'test3', lineType: 'hits', y: 65, height: 20 },
          ]});
      });

      it('ignores disabled tracks', function() {
        G.loadTrack('test', 'extracts');
        G.loadTrack('test2', 'extracts');
        G.loadTrack('test3', 'hits');
        _.find(G.tracks, {text:'test3'}).active = false;
        expect(G.headerSpec()).toEqual({
          headerY: 65,
          headers: [
            { text: 'test', lineType: 'extracts', y: 5, height: 30 },
            { text: 'test2', lineType: 'extracts', y: 35, height: 30 },
          ]});
      });
    });
    describe('.partition', function() {
      beforeEach(function() {
        angular.extend(G, {startBase: 0, endBase: 450,
                           basesPerGraph: 100, offset: 0});
      });
      it('partitions', function(){
        var p = G.partition();
        expect(p).toEqual([
          {startBase:0, endBase:99},
          {startBase:100, endBase:199},
          {startBase:200, endBase:299},
          {startBase:300, endBase:399},
          {startBase:400, endBase:499}
        ]);
      });

      it('handles a positive offset', function(){
        G.offset = 10;
        var p = G.partition();

        expect(p).toEqual([
          {startBase:10, endBase:109},
          {startBase:110, endBase:209},
          {startBase:210, endBase:309},
          {startBase:310, endBase:409},
          {startBase:410, endBase:509}
        ]);
      });
    });
  });
});
