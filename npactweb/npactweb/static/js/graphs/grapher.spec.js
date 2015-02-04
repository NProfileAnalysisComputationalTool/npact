describe('Grapher', function() {
  'use strict';
  beforeEach(module('assets'));
  beforeEach(module('npact', function($provide) {
//    $provide.value('$log', console);
    $provide.service('Fetcher', function() { });
  }));

  var tracks;
  // Fill out some dependency data
  beforeEach(inject(function(NProfiler, $templateCache, $q, $timeout, Track) {
    NProfiler.ddna = $templateCache.get('/js/test-data/sampleDdnaFile.ddna');
    NProfiler.fetching = $q.when(NProfiler.ddna);
    tracks = [new Track('Extracts',
                        $templateCache.get('/js/test-data/NC_007760.genes'),
                        'extracts'),
              new Track('Hits',
                        $templateCache.get('/js/test-data/NC_007760.profiles'),
                        'hits')];
    try { $timeout.flush(); } catch(e) {}
  }));


  var element;
  //Setup a test dom node we can render to.
  beforeEach(function() {
    element = angular.element('<div id="testGraph"></div>');
    angular.element('body').append(element);
  });
  afterEach(function() { angular.element(element).remove(); });


  var Grapher, $timeout, makeGrapher, $log;

  beforeEach(inject(function(GraphConfig, _Grapher_, $rootScope, _$timeout_, npactConstants,
                      GraphingCalculator, _$log_) {
    Grapher = _Grapher_;
    $timeout = _$timeout_;
    $log = _$log_;
    GraphConfig.endBase = 11140;
    GraphConfig.basesPerGraph = 5000;

    makeGrapher = function(sampleOpts) {
      expect(Grapher).toEqual(jasmine.any(Function));
      var opts = _.clone(sampleOpts);
      opts.tracks = tracks;
      opts.colors = npactConstants.lineColors;
      opts.m = GraphingCalculator.chart(opts);
      expect(opts.m.graph).toBeDefined();
      expect(opts.m.graph.w).toBeDefined();
      var g = new Grapher(element[0], opts);
      expect(g).toBeDefined();
      expect(g).toEqual(jasmine.any(Grapher));
      expect(g.xaxis).toBeDefined();
      expect(g.onPan).toEqual(jasmine.any(Function));
      return g;
    };
  }));

  it('should have an element to work with', function() {
    $log.log("In element test", element);
    expect(element).toBeDefined();
  });

  it('should be valid', function(done) {
    //This is mostly just a smoke test to run through all the code.
    var sampleOpts = {
        startBase: 0,
        margin: 1000,
        offset: 0,
        width: 600,
        axisTitle: 'C+G'
    };
    var g = makeGrapher(sampleOpts);
    var redraw = function() {
      expect(g).toBeDefined();
      expect(g).toEqual(jasmine.any(Grapher));
      expect(g.xaxis).toBeDefined();
      return g.redraw(sampleOpts);
    };
    var p = redraw();
    _.times(2, function() {
      p = p.then(redraw);
    });
    p.then(done);
    $timeout.flush();
    //KineticJs uses window timeouts, so we need to let that happen
    //and then flush angular's fake $timeout again
    window.setTimeout($timeout.flush, 0);
  });
});
