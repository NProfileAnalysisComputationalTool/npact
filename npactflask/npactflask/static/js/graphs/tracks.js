angular.module('npact')
  .factory('Track', function($log, ExtractParser, TrackBinSearchIndex, npactConstants) {
    'use strict';
    var ITrackIndex = TrackBinSearchIndex;
    function Track(name, data, type, weight) {
      this.active = true;
      this.name = name;
      this.type = type;
      if(!_.includes(['extracts', 'hits', 'modified', 'neworfs'], type)) {
        throw new Error("Unknown track type: +", type);
      }
      var trackStyle = npactConstants.trackStyle[type] || {};
      this.style = _.defaults(trackStyle, npactConstants.trackStyle['default']);
      this.weight = weight || 0;
      this.data = this.loadData(data);
      this.index = this.indexData(this.data);
    }
    Track.prototype.loadData = function(data) {
      return ExtractParser.parseAsync(data)
        .catch(function(e) { $log.log('Track.loadData failed', name, e); throw e; });
    };
    Track.prototype.indexData = function(data) {
      return data
        .then(function(data) { return new ITrackIndex(data); })
        .catch(function(e) { $log.log('Track.indexData failed', name, e); throw e; });
    };
    Track.prototype.slice = function(opts) {
      return this.index.then(_.bind(function(index) {
        this.slice = index;
        return index(opts);
      }, this));
    };
    Track.prototype.findByName = function(substr) {
      return this.data.then(function(parsedData) {
        if(!parsedData) return [];
        var searcher = new RegExp(substr.replace(/[.*+?^${}()|[\]\\]/g, "\\$&"), "i");
        return _.filter(parsedData, function(orf) {
          return searcher.test(orf.name);
        });
      });
    };
    return Track;
  })

  .factory('TrackRTreeIndex', function($log, Utils, Err, $q, $timeout) {
    'use strict';
    function TrackRTreeIndex(data) {
      var index = RTree();
      var indexReady = Utils.forEachAsync(data, function(d) {
          index.insert({ x: d.start, y: 0,
                         w: d.end - d.start, h: 5},
                       d);
      });
      return function(opts) {
        return $timeout(function() {
          return indexReady.then(function() {
            return index.bbox(opts.startBase, 0, opts.endBase, 5);
          });
        });
      };
    }
    return TrackRTreeIndex;
  })

  .factory('TrackNaiveIndex', function(Utils, $log, Err, $q, $timeout) {
    'use strict';
    return function(data) {
      return function(opts) {
        var minIdx = null, maxIdx;
        return Utils.forEachAsync(data, function(dataPoint, idx) {
          // extract starts in this range?
          var startsInRange = dataPoint.start >= opts.startBase &&
                dataPoint.start < opts.endBase,
              // extract ends in this range?
              endsInRange = dataPoint.end > opts.startBase &&
                dataPoint.end <= opts.endBase;
          if(startsInRange || endsInRange){
            if(minIdx === null) { minIdx = idx; }
            maxIdx = idx;
          }
        }).then(function() {
          return minIdx === null ? [] : data.slice(minIdx, maxIdx+1);
        });
      };
    };
  })
  .factory('TrackBinSearchIndex', function($log, Utils, Err, $q, $timeout) {
    'use strict';
    return function(data) {
      var byStart = _.sortBy(data, 'start');
      var byEnd = _(byStart).map(function(d, idx) {
                return {end: d.end, idx: idx};
      }).sortBy('end').value();

      return function(opts) {
        return $timeout(function() {
          var min = opts.startBase, max = opts.endBase;
          if(isNaN(min) || isNaN(max)) {
            throw new Error("Invalid start or end coordinate");
          }
          //Is the max coord larger than the start of any?
          var rightIdx = _.sortedIndex(byStart, {start: max}, 'start');
          //Is the min coord smaller than the end of any (with reference back to item)
          var leftPointerIdx = _.sortedLastIndex(byEnd, {end: min}, 'end'),
              leftIdx = rightIdx;
          //Because the byEnd array isn't necessarilty sorted by start
          //walk it to the right and if any of the entries reference
          //one that is less than our current left idx then update the
          //lefidx
          while(leftPointerIdx < byEnd.length) {
            leftIdx = Math.min(leftIdx, byEnd[leftPointerIdx].idx);
            //if the end coordinate is past the max then we're done
            if(byEnd[leftPointerIdx].end >= max) break;
            leftPointerIdx++;
          }
          return _.slice(byStart, leftIdx, rightIdx);
        });
      };
    };
  })
;
