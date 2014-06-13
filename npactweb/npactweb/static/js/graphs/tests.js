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
    var GC;
    beforeEach(inject(function(GraphingCalculator){
      GC = GraphingCalculator;
    }));
    it('calculates chart measurements', function(){
      var m = GC.chart({
	stageHeight:150,
	profileHeight:100,
	axisLabelFontsize:10,
	profileTicks:5
      });
      expect(m).toEqual({
	y:35,
	tickY:20,
	labelY:15
      });
    });

    it('can align rectangles', function(){
      var pos = GC.alignRectangles({x:0, y:0, w:100, h:100}, {w:10, h:10});
      expect(pos).toEqual({x: 45, y: 45 });
    });
    
  });

});
