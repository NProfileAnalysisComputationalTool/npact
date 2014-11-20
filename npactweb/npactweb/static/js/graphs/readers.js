angular.module('npact')
  .service('ProfileReader', function($q, Err) {
    'use strict';
    var self = this;

    /**
     * set the current loaded profile, returns the summary
     */
    this.load = function(profile) {
      self.profileData = profile;
      return $q.when(self.summary(profile));
    };

    /**
     * get a shallow copy of the profile for the given range
     */
    this.slice = function(range) {
      if (!self.profileData) { throw new Err.ProfileNotFound(); }

      var sortedIdx = function(coord) {
        return _.sortedIndex(self.profileData, {coordinate: coord}, 'coordinate');
      },
          startIdx = sortedIdx(range.startBase),
          endIdx = sortedIdx(range.endBase) + 1,
          minEndIdx = self.profileData.length - 1
          ;

      // shallow copy of the relevant data, guarding against
      // out-of-bounds indexes
      return self.profileData
        .slice(Math.max(startIdx,0), Math.min(endIdx, minEndIdx));

    };

    /**
     * Calculate summary stats about a profile
     */
    this.summary = function(profile) {
      var p = profile || self.profileData;
      if (!p) { throw new Err.ProfileNotFound(); }

      var start = p[0],
          end = _.last(p);

      return {
        startBase: start.coordinate,
        endBase: end.coordinate,
        length: end.coordinate - start.coordinate
      };
    };

    /**
     * Partition the profile for graphing
     */
    this.partition = function(opts) {
      var p = opts.profileSummary || self.summary(),
          offset = opts.offset || 0,
          startBase = Math.max(p.startBase + offset, 0),
          endBase = p.endBase,
          bpg = opts.basesPerGraph,
          g = new Array(Math.ceil((endBase - startBase) / bpg)),
          idx = 0;
      while(startBase < endBase) {
        g[idx] = {
          startBase: startBase,
          endBase: startBase + bpg - 1
        };
        idx++;
        startBase = startBase + bpg;
      }
      return g;
    };
  })
  .service('TrackReaderRTreeIndex', function(ExtractParser, $log, Utils, Err, $q, $timeout){
    var self = this;

    self.tracks = {};
    /**
     * load the given track name and data
     */
    this.load = function(name, data) {
      if (_.has(self.tracks, name)) { throw new Err.TrackAlreadyDefined(); }
      return ExtractParser.parseAsync(data)
        .then(
          function(data) {
            var index = RTree();
            return Utils.forEachAsync(data, function(d) {
              index.insert({ x: d.start, y: 0,
                             w: d.end - d.start, h: 5},
                           d);
            }).then(function() { self.tracks[name] = index; });
          },
          function(e) { $log.log('TrackReader.load failed', name, e); throw e; }
        );
    };

    this.slice = function(opts) {
      if (!_.has(self.tracks, opts.name)) { throw new Err.TrackNotFound(); }
      return $q(function(resolve, reject) {
        $timeout(function() {
          resolve(self.tracks[opts.name].bbox(opts.startBase, 0, opts.endBase, 5));
        });
      });
    };
  })
  .service('TrackReaderNaive', function(ExtractParser, $log, Utils, Err, $q, $timeout){
    var self = this;

    self.tracks = {};
    /**
     * load the given track name and data
     */
    this.load = function(name, data) {
      if (_.has(self.tracks, name)) { throw new Err.TrackAlreadyDefined(); }
      return ExtractParser.parseAsync(data)
        .then(
          function(data) { return (self.tracks[name] = data); },
          function(e) { $log.log('TrackReader.load failed', name, e); throw e; }
        );
    };

    this.slice = function(opts) {
      if (!_.has(self.tracks, opts.name)) { throw new Err.TrackNotFound(); }
      $log.log('TrackReader.slice', opts);

      var minIdx = null, maxIdx, data = self.tracks[opts.name];
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
    };
  })
  .service('TrackBinarySearch', function(ExtractParser, $log, Utils, Err, $q, $timeout) {
    var self = this;
    self.tracks = {};

        /**
     * load the given track name and data
     */
    this.load = function(name, data) {
      if (_.has(self.tracks, name)) { throw new Err.TrackAlreadyDefined(); }
      return ExtractParser.parseAsync(data)
        .then(
          function(data) {
            var byStart = _.sortBy(data, 'start');
            return (self.tracks[name] = {
              byStart: byStart,
              byEnd: _(byStart).map(function(d, idx) {
                return {end:d.end, idx:idx};
              }).sortBy('end').value()
            });
          },
          function(e) { $log.log('TrackReader.load failed', name, e); throw e; }
        );
    };

    this.slice = function(opts) {
      if (!_.has(self.tracks, opts.name)) { throw new Err.TrackNotFound(); }

      var data = self.tracks[opts.name],
          min = opts.startBase,
          max = opts.endBase,
          rightIdx = _.sortedIndex(data.byStart, {start: max}, 'start'),
          leftPointerIdx = _.sortedIndex(data.byEnd, {end: min}, 'end'),
          leftIdx = data.byStart.length -1;

      while(leftPointerIdx < data.byEnd.length && data.byEnd[leftPointerIdx].end <= max){
        leftIdx = Math.min(leftIdx, data.byEnd[leftPointerIdx].idx);
        leftPointerIdx++;
      }
      return $q.when(data.byStart.slice(leftIdx, rightIdx));
    };
  })
  .factory('TrackReader', function(TrackReaderRTreeIndex, TrackReaderNaive, TrackBinarySearch) {

    return TrackBinarySearch; // 906 913 907
    // return TrackBinarySearch; // 958 957 979

    return TrackReaderNaive; // 938 961 936 941

    return TrackReaderRTreeIndex; // 967 917 916
  })
;
