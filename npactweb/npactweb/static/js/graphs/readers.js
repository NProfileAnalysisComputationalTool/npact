angular.module('npact')

  .service('TrackReaderRTreeIndex', function(ExtractParser, $log, Utils, Err, $q, $timeout){
    this.tracks = {};
    /**
     * load the given track name and data
     */
    this.load = function(name, data) {
      return (
        this.tracks[name] =
          ExtractParser.parseAsync(data)
          .then(
            function(data) {
              var index = RTree();
              return Utils.forEachAsync(data, function(d) {
                index.insert({ x: d.start, y: 0,
                               w: d.end - d.start, h: 5},
                             d);
              }).then(function() { return index; });
            },
            function(e) { $log.log('TrackReader.load failed', name, e); throw e; }
          ));
    };

    this.slice = function(opts) {
      var track = this.tracks[opts.name];
      if (!track) { throw new Err.TrackNotFound(); }
      return track.then(function(index) {
        return $timeout(function() {
          return index.bbox(opts.startBase, 0, opts.endBase, 5);
        });
      });
    };
  })
  .service('TrackReaderNaive', function(ExtractParser, $log, Utils, Err, $q, $timeout){
    this.tracks = {};
    /**
     * load the given track name and data
     */
    this.load = function(name, data) {
      return (
        this.tracks[name] = ExtractParser.parseAsync(data)
          .catch(
            function(e) { $log.log('TrackReader.load failed', name, e); throw e; }
          ));
    };

    this.slice = function(opts) {
      if (!_.has(this.tracks, opts.name)) { throw new Err.TrackNotFound(); }
      $log.log('TrackReader.slice', opts);

      var minIdx = null, maxIdx, data = this.tracks[opts.name];
      return this.tracks[opts.name].then(function(data) {
        // TODO: use _.sortedIndex to binary search, or some other
        // better index structure/algorithm
        return Utils.forEachAsync(data, function(dataPoint, idx) {
          // extract starts in this range?
          var startsInRange = dataPoint.start >= opts.startBase &&
                dataPoint.start <= opts.endBase,
              // extract ends in this range?
              endsInRange = dataPoint.end >= opts.startBase &&
                dataPoint.end <= opts.endBase;
          if(startsInRange || endsInRange){
            if(minIdx === null) { minIdx = idx; }
            maxIdx = idx;
          }
        }).then(function() {
          return minIdx === null ? [] : data.slice(minIdx, maxIdx+1);
        });
      });
    };
  })
  .service('TrackBinarySearch', function(ExtractParser, $log, Utils, Err, $q, $timeout) {
    this.tracks = {};

    /**
     * load the given track name and data
     */
    this.load = function(name, data) {
      return (this.tracks[name] = ExtractParser.parseAsync(data)
        .then(
          function(data) {
            var byStart = _.sortBy(data, 'start');
            return {
              byStart: byStart,
              byEnd: _(byStart).map(function(d, idx) {
                return {end:d.end, idx:idx};
              }).sortBy('end').value()
            };
          },
          function(e) { $log.log('TrackReader.load failed', name, e); throw e; }
        ));
    };

    this.slice = function(opts) {
      if (!_.has(this.tracks, opts.name)) {
        $log.log('Missing', opts.name);
        throw new Err.TrackNotFound();
      }
      return this.tracks[opts.name].then(function(data) {
        var min = opts.startBase,
            max = opts.endBase,
            rightIdx = _.sortedIndex(data.byStart, {start: max}, 'start'),
            leftPointerIdx = _.sortedIndex(data.byEnd, {end: min}, 'end'),
            leftIdx = data.byStart.length -1;
        while(leftPointerIdx < data.byEnd.length &&
              data.byEnd[leftPointerIdx].end <= max) {
          leftIdx = Math.min(leftIdx, data.byEnd[leftPointerIdx].idx);
          leftPointerIdx++;
        }
        return data.byStart.slice(leftIdx, rightIdx);
      });
    };
  })
  .factory('TrackReader', function(TrackReaderRTreeIndex, TrackReaderNaive, TrackBinarySearch) {

    return TrackBinarySearch; // 906 913 907
    // return TrackBinarySearch; // 958 957 979

    return TrackReaderNaive; // 938 961 936 941

    return TrackReaderRTreeIndex; // 967 917 916
  })
;
