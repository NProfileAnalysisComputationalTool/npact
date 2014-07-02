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

  // These are useless right now, delete this if it's still commented out
  // as of 2014-07-01
  
  // describe('npactGraphPage', function(){
  //   var testProfile = {};
  //   it('can render', function(){
  //     $httpBackend.expectGET('../npactweb/static/js/nprofile.json')
  // 	.respond(200, testProfile);
  //     make('<div npact-graph-page>');
  //   });
  // });

  // describe('npactGraph', function(){
  //   it('can render', function(){
  //     make('<div npact-graph>');
  //   });
  // });

  // describe('Grapher', function(){
  //   var g;
  //   beforeEach(inject(function(Grapher){
  //     var el = angular.element('<div>');
  //     g = Grapher({element:el, width:600, height:200});
  //   }));

  //   it('can create', function(){
  //     g.redraw();
  //   });
  // });

  describe('GraphingCalculator', function(){
    var GC, opts = {
      stageHeight:150,
      stageWidth:255,
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


    describe('can zoom', function(){

      var input;

      beforeEach(function(){
	input = {
	  oldZoom: 1,
	  zoom:2,
	  baseScaleX: 1, // px/gn
	  start_gn: 0,
	  length_gn: 100,
	  offsetX:0,
	  centerOn_px: 50
	};
      });
     
      it('after panning left', function(){
	input.offsetX = 20;
	input.centerOn_px = 20; // graph coords
	expect(GC.zoom(input)).toEqual({
	  scaleX:2,
	  offsetX:60,
	  textScaleX:1/2
	});
      });

      it('after panning right', function(){
	input.offsetX = -20;
	input.centerOn_px = 20; // graph coords
	expect(GC.zoom(input)).toEqual({
	  scaleX:2,
	  offsetX:-20,
	  textScaleX:1/2
	});
      });

      it('with no panning', function(){
	input.centerOn_px = 20;
	expect(GC.zoom(input)).toEqual({
	  scaleX:2,
	  offsetX:20,
	  textScaleX:1/2
	});
      });
      
    });
  });

});
