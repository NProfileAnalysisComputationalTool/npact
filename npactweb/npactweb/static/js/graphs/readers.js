angular.module('npact')
  .service('ProfileReader', function($q) {
    'use strict';
    var self = this;

    self.ProfileNotFound = new Error('profile not found');

    /**
     * set the current loaded profile, returns the summary
     */
    this.load = function(profile) {
      self.profileData = profile;
      var d = $q.defer();
      d.resolve(self.summary(profile));
      return d.promise;
    };

    /**
     * get a shallow copy of the profile for the given range
     */
    this.slice = function(range) {
      if (!self.profileData) { throw self.ProfileNotFound; }

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
      if (!p) { throw self.ProfileNotFound; }

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
      var p = opts.summary,
          margin = opts.margin || 0,
          offset = opts.offset || 0,
          startBase = p.startBase + offset,
          g = [];

      for(var i = startBase; i < p.endBase; i+= opts.basesPerGraph){
        g.push({
          startBase: Math.max(i - margin, p.startBase),
          endBase: Math.min(i + opts.basesPerGraph + margin, p.endBase)
        });
      }
      return g;
    };
  })
;
