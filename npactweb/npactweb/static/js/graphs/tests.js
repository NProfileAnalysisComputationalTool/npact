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
      profileHeight:100,      
      axisLabelFontsize:10,
      axisTitleFontsize: 20,
      axisTitle: '% GC',
      profileTicks:5
    };
    beforeEach(inject(function(GraphingCalculator){
      GC = GraphingCalculator;
    }));
    it('calculates graph area', function(){
      var m = GC.chart(opts);
      expect(m.graph).toEqual({
	x: 50,
	y: 35,
	h: 100,
	w: 200
      });
    });
    it('calculates y-axis ticks', function(){
      var m = GC.chart(opts);
      expect(m.yaxis.ticks).toEqual([
	{x:45, y:35, x2:50, y2:35},
	{x:45, y:55, x2:50, y2:55},
	{x:45, y:75, x2:50, y2:75},
	{x:45, y:95, x2:50, y2:95},
	{x:45, y:115, x2:50, y2:115},
	{x:45, y:135, x2:50, y2:135}	
      ]);
    });
    
    it('calculates y-axis labels', function(){
      var m = GC.chart(opts);
      expect(m.yaxis.labels).toEqual([
	{x:0, width: 40, y:30, text:100},
	{x:0, width: 40, y:50, text:80},
	{x:0, width: 40, y:70, text:60},
	{x:0, width: 40, y:90, text:40},
	{x:0, width: 40, y:110, text:20},
	{x:0, width: 40, y:130, text:0}	
      ]);
    });

    it('calculates y-axis title', function(){
      var m = GC.chart(opts);
      expect(m.yaxis.titleBox).toEqual(
	{x:0, y:35, w:40, h:100}
      );
    });

    it('calculates x-axis ticks', function(){
      var m = GC.chart(opts);
      expect(m.yaxis.ticks).toEqual([
	{x:45, y:35, x2:50, y2:35},
	{x:45, y:55, x2:50, y2:55},
	{x:45, y:75, x2:50, y2:75},
	{x:45, y:95, x2:50, y2:95},
	{x:45, y:115, x2:50, y2:115},
	{x:45, y:135, x2:50, y2:135}	
      ]);
    });
    
    it('can align rectangles', function(){
      var pos = GC.alignRectangles({x:0, y:0, w:100, h:100}, {w:10, h:10});
      expect(pos).toEqual({x: 45, y: 45 });
    });
    
  });

});
