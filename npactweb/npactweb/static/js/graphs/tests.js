'use strict';

describe('Graphs', function(){
  beforeEach(module('npact'));
  beforeEach(module('templates-main'));

  var $scope, $compile, $q;
  beforeEach(inject(function (_$rootScope_, _$compile_, _$q_) {
    $scope = _$rootScope_;
    $compile = _$compile_;
    $q = _$q_;
  }));

  function make(html){
    var el = $compile(html)($scope);
    $scope.$digest();
    return el;
  }

  
  describe('page', function(){
    it('can render', function(){
      make('<div npact-graph-page>');
    });
  });

});
