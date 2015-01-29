describe('Grapher', function() {
  'use strict';
  beforeEach(module('assets'));
  beforeEach(module('npact', function($provide) {
//    $provide.value('$log', console);
    $provide.service('Fetcher', function() { });
  }));

  var Grapher, $digest, $timeout, element, npactConstants,
      makeGrapher,
      sampleOpts = {
        startBase: 0,
        endBase: 10000,
        margin: 1000,
        offset: 0,
        width: 600,
        axisTitle: 'C+G'
      };

  beforeEach(function() { element = angular.element('<div id="testGraph"></div>'); });
  afterEach(function() { angular.element(element).remove(); });

  beforeEach(inject(function(_Grapher_, $rootScope, _$timeout_, _npactConstants_,
                      GraphingCalculator) {
    Grapher = _Grapher_;
    npactConstants = _npactConstants_;
    $digest = $rootScope.$digest;
    $timeout = _$timeout_;

    makeGrapher = function() {
      expect(Grapher).toEqual(jasmine.any(Function));
      var opts = _.clone(sampleOpts);
      opts.colors = npactConstants.lineColors;
      opts = angular.extend(opts,
                            GraphingCalculator.trackSizeCalc(
                              [{text: 'Hits', lineType: 'hits'},
                               {text: 'Extracts', lineType: 'extracts'}]));
      opts.m = GraphingCalculator.chart(opts);
      expect(opts.m.graph).toBeDefined();
      expect(opts.m.graph.w).toBeDefined();
      var g = new Grapher(element[0], opts);
      expect(g).toBeDefined();
      expect(g).toEqual(jasmine.any(Grapher));
      expect(g.xaxis).toBeDefined();
      expect(g.onDragEnd).toEqual(jasmine.any(Function));
      return g;
    };
  }));

  beforeEach(inject(function(NProfiler, TrackReader, $templateCache, $q, $log, $timeout) {
    NProfiler.ddna = $templateCache.get('/js/test-data/sampleDdnaFile.ddna');
    NProfiler.fetching = $q.when(NProfiler.ddna);

    TrackReader.load('Extracts',
                     $templateCache.get('/js/test-data/NC_007760.genes'))
      .then(function() { $log.log("Finished loading extracts"); });
    TrackReader.load('Hits',
                     $templateCache.get('/js/test-data/NC_007760.profiles'))
      .then(function() { $log.log("Finished loading hits"); });
    $timeout.flush();
  }));

  it('should have an element to work with', function() {
    expect(element).toBeDefined();
  });

  it('should be valid', function(done) {
    //This is mostly just a smoke test to run through all the code.

    var g = makeGrapher();
    var redraw = function() {
      expect(g).toBeDefined();
      expect(g).toEqual(jasmine.any(Grapher));
      expect(g.xaxis).toBeDefined();
      return g.redraw(sampleOpts);
    };
    var p = redraw();
    for(var i = 20; i > 0; --i) {
      p = p.then(redraw);
    }
    p.then(done);
    $timeout.flush();
    //KineticJs uses window timeouts, so we need to let that happen
    //and then flush angular's fake $timeout again
    window.setTimeout($timeout.flush, 0);
  });
});
