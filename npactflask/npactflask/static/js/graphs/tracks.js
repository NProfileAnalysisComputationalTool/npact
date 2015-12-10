angular.module('npact')
  .service('TrackSyncer', function(Fetcher, Track, $http, SAVE_TRACK_URL){
    'use strict';
    var self = this;
    self.fetchTrack= function(filename, name, type, fetcher){
      if(!type) type = "custom";
      if(!name) name = filename;
      if(!fetcher) fetcher = 'fetchFile';
      //console.log('fetchTrack ',name, type , filename);
      config[type] = Fetcher[fetcher](filename).then(function(data){
        var track = new Track(name, data, type, 100, filename);
        GraphConfig.loadTrack(track);
        track.saveTrack = function(){
          self.saveTrack(track);
        };
        return track;
      });
    };
    self.fetchHits = function(HitsFile) {
      return self.fetchTrack(HitsFile, 'Hits', 'hits');
    };

    self.fetchNewOrfs = function(NewOrfsFile) {
      return self.fetchTrack(NewOrfsFile, 'New ORFs', 'neworfs');
    };

    self.fetchModifiedOrfs = function(ModifiedOrfsFile) {
      return self.fetchTrack(ModifiedOrfsFile, 'Modified ORFs', 'modified');
    };

    self.saveTrack = function(track){
      return $http({
        method: 'POST',
        url: SAVE_TRACK_URL,
        data: { track: track },
        headers: {'Content-Type': 'application/json'}
      })
        .then(function (res) {
          return res.data;
        });
    };

  })

  .factory('Track', function($log, ExtractParser, TrackBinSearchIndex, npactConstants, $timeout) {
    'use strict';
    var ITrackIndex = TrackBinSearchIndex;
    function Track(name, data, type, weight, filename) {
      // console.log('new Track ',name, data, type ,weight, filename);
      this.active = true;
      this.name = name;
      this.type = type;
      this.filename = filename;
      if(!_.includes(['extracts', 'hits', 'modified', 'neworfs', 'custom'], type)) {
        throw new Error("Unknown track type: ", type);
        //console.log("Unknown track type: ", type);
      }
      var trackStyle = npactConstants.trackStyle[type] || {};
      this.style = _.defaults(trackStyle, npactConstants.trackStyle['default']);
      this.weight = weight || 0;
      this.data = null;
      this.index = null;
      this.loading = this.loadData(data);
      this.indexing = this.loading.then(_.bind(this.indexData, this));
    }

    Track.prototype.loadData = function(data) {
      return (
        this.loading = ExtractParser.parseAsync(data)
          .then(_.bind(function (data) {
            $log.log("Finished parsing", this.name, ", found ", data.length);
            this.data = data;
            this.loading = false;
            return data;
          }, this))
          .catch(function(e) { $log.log('Track.loadData failed', name, e); throw e; }));
    };

    Track.prototype.indexData = function(data) {
      return (
        this.indexing = $timeout(_.bind(function () {
          this.indexing = false;
          return (this.index = ITrackIndex(data));
        }, this))
          .catch(function(e) { $log.log('Track.indexData failed', name, e); throw e; }));
    };

    Track.prototype.reindex = function () {
      if(this.loading) {
        //shouldn't be able to get here but if this is true then we
        //haven't yet started indexing and nothing to do.
        return;
      }
      else if(this.indexing) {
        this.indexing.then(_.bind(function () {
          this.indexData(this.data);
        }, this));
      }
      else {
        this.indexData(this.data);
      }
   };
    Track.prototype.slice = function(opts) {
      if(this.indexing) {
        return this.indexing.then(function (index) {
          return index(opts);
        });
      }
      else {
        return this.index(opts);
      }
    };

    Track.prototype.findByName = function(substr) {
      var doFind = function () {
        if(!this.data) return [];
        var searcher = new RegExp(substr.replace(/[.*+?^${}()|[\]\\]/g, "\\$&"), "i");
        return _.filter(this.data, function(orf) {
          return searcher.test(orf.name);
        });
      };
      return (this.loading ? this.loading.then : $timeout)(_.bind(doFind, this));
    };

    Track.prototype.remove = function (entry) {
      if(this.loading) {
        throw new Error("Shouldn't be able to remove from a track before the track exists");
      }
      this.data = _.without(this.data, entry);
      this.reindex();
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
