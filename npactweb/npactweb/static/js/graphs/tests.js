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
	{x: 0, y: 30, w: 40, h: 100}
      );
    });

    it('calculates x-axis', function(){
      
      var ax = GC.xaxis(opts);
      expect(ax).toEqual({
	start: 0, end: 100,
	length: 100,
	x: 50, y: 135, scaleX: 2
      });
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
	{x: 100, y: 0, x2: 100, y2: 5} 
      ]);
    });

    it('can align rectangles', function(){
      var pos = GC.alignRectangles(
	{x: 0, y: 0, w: 100, h: 100}, 
	{w: 10, h: 10});
      expect(pos).toEqual({x: 45, y: 45});
    });
    
  });

});
