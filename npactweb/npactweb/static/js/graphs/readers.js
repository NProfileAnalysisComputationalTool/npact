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
      var p = opts.summary || self.summary(),
          offset = opts.offset || 0,
          startBase = p.startBase + offset,
          g = [];

      for(var i = startBase; i < p.endBase; i+= opts.basesPerGraph){
        g.push({
          startBase: Math.max(i, 0),
          endBase: i + opts.basesPerGraph
        });
      }
      return g;
    };
  })
  .service('TrackReader', function(ExtractParser, $log, Utils, Err){
    var self = this;

    self.tracks = {};
    /**
     * load the given track name and data
     */
    this.load = function(name, data) {
      if (_.has(self.tracks, name)) { throw new Err.TrackAlreadyDefined(); }
      return ExtractParser.parseAsync(data)
        .then(function(data) {
          return (self.tracks[name] = data);
        }, function() {
          $log.log('failed to parse data for track', name);
        });
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
;
