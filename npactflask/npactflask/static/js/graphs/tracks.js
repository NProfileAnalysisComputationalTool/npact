angular.module('npact')
  .factory('Track', function($log,  $q, $timeout, $http,
                      Fetcher, ExtractParser, TrackBinSearchIndex, npactConstants) {
    'use strict';

    var ITrackIndex = TrackBinSearchIndex;
    function Track(filename, data, type, name) {
      // console.log('new Track ',name, data, type ,filename);
      if(!name && !type){
        if(filename.search('profiles')>=0){
          name="Hits";
          type='hits';
        }
        else if(filename.search('modified')>=0){
          name="Modified ORFs";
          type='extracts';
        }
        else if(filename.search('newcds')>=0){
          name="New ORFs";
          type='extracts';
        }
        else if(filename.search('genes')>=0){
          name="Input CDSs";
          type='extracts';
        }
        else{
          name = filename;
          type = 'extracts';
        }
      }
      this.active = true;
      this.filename = filename;
      this.name = name;
      this.type = type;
      if(!_.includes(['extracts', 'hits'], type)) {
        throw new Error("Unknown track type: ", type);
        //console.log("Unknown track type: ", type);
      }
      var trackStyle = npactConstants.trackStyle[type] || {};
      this.style = _.defaults(trackStyle, npactConstants.trackStyle['default']);
      this.data = null;
      this.index = null;
      this.loading = this.loadData(data);
      this.reindex();
    }

    Track.prototype.id = function () {
      var fn = _.last(this.filename.split('/'));
      return _.takeRight(fn.split('.'),2)[0];
    };

    Track.prototype.loadData = function(data) {
      var self = this;
      // console.log(data);
      if(typeof(data) === 'string') { // trackfile to parse
        return (
          self.loading = ExtractParser.parseAsync(data)
            .then(function (data) {
              $log.log("Finished parsing", self.name, ", found ", data.length);
              self.data = [];
              self.metadata = {};
              var cdsidx = 0;
              _.each(data, function(v, k) {
                if(v.type === 'CDS'){
                  self.data.push(v);
                  v.cdsidx = cdsidx++;
                }
                else if(v.type === 'META') self.metadata[k]=v;
              });
              self.loading = false;
              return data;
            })
            .catch(function(e) { $log.log('Track.loadData failed', name, e); throw e; }));
      }
      else { // data must be parsed (probably json version of this object)
        if(Array.isArray(data))
          self.data = data;
        else if(data.data){
          self.data = data.data;
          if(data.name) self.name = data.name;
        }
        return null;
      }
    };

    Track.prototype.indexData = function(data) {
      return (
        this.indexing = $timeout(_.bind(function () {
          _.each(this.data, function(orf, k) {
            orf.cdsidx = k;
          });
          this.indexing = false;
          return (this.index = ITrackIndex(data));
        }, this))
          .catch(function(e) { $log.log('Track.indexData failed', name, e); throw e; }));
    };

    Track.prototype.reindex = function () {
      if(this.loading) {
        this.loading.then(_.bind(this.indexData, this));
      }
      else if(this.indexing) {
        this.indexing.then(_.bind(this.indexData, this, this.data));
      }
      else {
        this.indexData(this.data);
      }
   };
    Track.prototype.slice = function(opts) {
      if(this.loading) {
        return this.loading.then(_.bind(function () {
          return this.slice(opts);
        }, this));
      }
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

    Track.prototype.add = function (entry) {
      if(this.loading) {
        throw new Error("Shouldn't be able to add to a track that doesn't exist");
      }
      this.data.push(entry);
      this.reindex();
    };

    Track.prototype.remove = function (entry) {
      if(!entry) return;
      if(this.loading) {
        throw new Error("Shouldn't be able to remove from a track before the track exists");
      }
      this.data = _.without(this.data, entry);
      this.reindex();
    };

    Track.prototype.save = function () {
      var track = this;
      delete Track.loadedTracks[track.filename];  //this is no longer a valid repr of the url.

      return $http({
        method: 'POST',
        url: Fetcher.buildUrl('save_track'),
        data: { track: _.pick(this, ['filename', 'name', 'type', 'active', 'data', 'metadata']) },
        headers: {'Content-Type': 'application/json'}
      })
        .then(function (res) {
          var tr = _.create(Track.prototype, track);
          tr.filename = res.data.filename;
          Track.loadedTracks[tr.filename] = $q.when(tr);
          return tr;
        });
    };

    /// Static methods on Track.
    Track.loadedTracks = {};

    Track.fetchAllTracks = function(paths){
      return $q.all(_.map(paths, Track.fetchTrack));
    };

    Track.fetchTrack = function(filename){
      console.log('fetchTrack ', filename);
      var tr = Track.loadedTracks[filename];
      if(!tr) {
        tr = Track.loadedTracks[filename] = Fetcher.pollThenFetch(filename)
          .then(function(data) { return new Track(filename, data); });
      }
      return tr;
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
