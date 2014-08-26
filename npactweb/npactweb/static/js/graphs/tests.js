'use strict';

describe('Graphs', function(){
  beforeEach(module('npact'));
  beforeEach(module('templates-main'));

  var $scope, $compile, $q, $httpBackend;
  beforeEach(inject(function (_$rootScope_, _$compile_, _$q_, _$httpBackend_) {
    $scope = _$rootScope_;
    $compile = _$compile_;
    $q = _$q_;
    $httpBackend = _$httpBackend_;
  }));

  function make(html){
    var el = $compile(html)($scope);
    $scope.$digest();
    return el;
  }

  describe('GraphDealer', function(){
    var GD;

    beforeEach(inject(function(GraphDealer){
      GD = GraphDealer;
    }));

    it('adjusts basesPerGraph based on profile data', function(){
      var p = [{coordinate:0}, {coordinate:5000}];
      expect(GD.opts.basesPerGraph).toBe(10000);
      expect(GD.opts.length).toBe(0);
      GD.setProfile(p);
      $scope.$apply(); // process promises
      expect(GD.opts.length).toBe(5000);
      expect(GD.opts.basesPerGraph).toBe(1000);
    });
  });

  describe('Utils', function(){
    var U;

    beforeEach(inject(function(Utils){
      U = Utils;
    }));

    it('calculates order of magnitude', function(){
      expect(U.orderOfMagnitude(1012)).toBe(1000);
      expect(U.orderOfMagnitude(1012, -1)).toBe(100);
    });
    it('calculates order of magnitude with negative offsets', function(){
      expect(U.orderOfMagnitude(1012, -1)).toBe(100);
    });
    it('calculates order of magnitude with positive offsets', function(){
      expect(U.orderOfMagnitude(1012, 1)).toBe(10000);
    });

  });

  describe('GraphingCalculator', function(){
    var GC, opts = {
      height:150,
      width:255,
      leftPadding: 50,
      rightPadding: 5,
      profileHeight:100,
      axisLabelFontsize:10,
      axisTitleFontsize: 20,
      axisTitle: '% GC',
      profileTicks:5,
      range:[0,100]
    };
    beforeEach(inject(function(GraphingCalculator){
      GC = GraphingCalculator;
    }));
    it('calculates graph area', function(){
      var m = GC.chart(opts);
      expect(m.graph).toEqual({
	x: 50,
	y: 30,
	h: 100,
	w: 200
      });
    });
    it('calculates y-axis ticks', function(){
      var m = GC.chart(opts);
      expect(m.yaxis.ticks).toEqual(
	[
	  {x: 45, y: 30, x2: 50, y2: 30},
	  {x: 45, y: 50, x2: 50, y2: 50},
	  {x: 45, y: 70, x2: 50, y2: 70},
	  {x: 45, y: 90, x2: 50, y2: 90},
	  {x: 45, y: 110, x2: 50, y2: 110},
	  {x: 45, y: 130, x2: 50, y2: 130}
	]);
    });

    it('calculates y-axis labels', function(){
      var m = GC.chart(opts);
      expect(m.yaxis.labels).toEqual([
	{x: 0, y: 25, width: 40, text: 100},
	{x: 0, y: 45, width: 40, text: 80},
	{x: 0, y: 65, width: 40, text: 60},
	{x: 0, y: 85, width: 40, text: 40},
	{x: 0, y: 105, width: 40, text: 20},
	{x: 0, y: 125, width: 40, text: 0}
      ]);
    });

    it('calculates y-axis title', function(){
      var m = GC.chart(opts);
      expect(m.yaxis.titleBox).toEqual(
	{x: 0, y: 30, width: 40, height: 100}
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
      expect(GC.zoom.apply(null, args))
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

});
