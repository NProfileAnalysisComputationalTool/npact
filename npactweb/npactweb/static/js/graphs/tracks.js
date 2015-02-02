angular.module('npact')
  .factory('Track', function($log, ExtractParser, TrackBinSearchIndex, npactConstants) {
    'use strict';
    var ITrackIndex = TrackBinSearchIndex;
    function Track(name, data, type, weight) {
      this.active = true;
      this.name = name;
      this.type = type;
      this.height = npactConstants.trackHeights[type];
      if(this.height === undefined) {
        throw new Error("Unknown track type: +", type);
      }
      this.weight = weight || 0;
      this.index = this.load(data);
    }
    Track.prototype.load = function(data) {
      return ExtractParser.parseAsync(data)
        .then(function(data) { return new ITrackIndex(data); })
        .catch(function(e) { $log.log('Track.load failed', name, e); throw e; });
    };
    Track.prototype.slice = function(opts) {
      return this.index.then(_.bind(function(index) {
        this.slice = index;
        return index(opts);
      }, this));
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
        $log.log('TrackIndex.slice', opts);
        var minIdx = null, maxIdx;
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
    };
  })
  .factory('TrackBinSearchIndex', function($log, Utils, Err, $q, $timeout) {
    'use strict';
    return function(data) {
      var byStart = _.sortBy(data, 'start');
      var byEnd = _(byStart).map(function(d, idx) {
                return {end:d.end, idx:idx};
      }).sortBy('end').value();

      return function(opts) {
        $log.log('TrackIndex.slice', opts);
        return $timeout(function() {
          var min = opts.startBase,
              max = opts.endBase,
              rightIdx = _.sortedIndex(byStart, {start: max}, 'start'),
              leftPointerIdx = _.sortedIndex(byEnd, {end: min}, 'end'),
              leftIdx = byStart.length -1;
          while(leftPointerIdx < byEnd.length &&
                byEnd[leftPointerIdx].end <= max) {
            leftIdx = Math.min(leftIdx, byEnd[leftPointerIdx].idx);
            leftPointerIdx++;
          }
          return byStart.slice(leftIdx, rightIdx);
        });
      };
    };
  })

;
