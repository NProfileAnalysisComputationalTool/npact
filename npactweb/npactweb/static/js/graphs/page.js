angular.module('npact')
  // the `GraphDealer` hands out graph data and processes events
  .service('GraphDealer', function($log){

    var opts = {
      basesPerGraph:10000,
      colorBlindFriendly:false,
      page:0
    };

    return {
      opts:opts,
      setColors:function(colorBlindFriendly){
	$log.log('setColors', arguments);
	// TODO: instruct graphs to redraw themselves with the right colors
	opts.colorBlindFriendly = colorBlindFriendly;
      },
      setZoom:function(basesPerGraph){
	$log.log('setZoom', arguments);
	// TODO: repartition data into groups of `basesPerGraph`
	// TODO: redraw graphs with new data
	// TODO: update pages
	opts.basesPerGraph = basesPerGraph;
      },
      zoomTo: function(){
	// TODO: figure out a good API for this
      },
      panTo: function(oldStartBase, newStartBase){
	// TODO: find the graph with the old start base
	// TODO: recalculate graph start/ends
	// TODO: repartition data
	// TODO: redraw graphs with new data
      },
      hasNextPage: function(){
	return true; // TODO: calculate a "max pages" based on zoom and data
      },
      nextPage: function(){
	$log.log('nextPage', arguments);
	// TODO: repartition data
	// TODO: redraw graphs
	opts.page++;
      },
      hasPreviousPage: function(){ return opts.page > 0;},
      previousPage: function(){
	$log.log('previousPage', arguments);
	// TODO: repartition data
	// TODO: redraw graphs
	opts.page--;
      }
    };
  })
  .directive('npactGraphPage', function(STATIC_BASE_URL, $http){
    return {
      restrict:'A',
      templateUrl:STATIC_BASE_URL+'js/graphs/page.html',
      link:function(scope, element, attrs){
	scope.title = 'Moorella thermoacetica Y72, whole genome shotgun sequence';
	scope.ranges = [
	  [0, 10000],
	  [10000, 20000],
	  [20000, 30000],
	  [30000, 40000],
	  [40000, 50000]
	];
	scope.d = {};

	$http.get(STATIC_BASE_URL+'js/nprofile.json').then(function(res){
	  scope.d.profile = res.data;
	});
	$http.get(STATIC_BASE_URL+'js/extract.json').then(function(res){
	  scope.d.cds = res.data;
	});
      }
    };
  });
