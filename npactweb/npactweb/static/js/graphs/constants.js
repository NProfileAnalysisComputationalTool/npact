angular.module('npact')
  .constant('npactConstants', {
    graphSpecDefaults: {
      // TODO: determine me dynamically
      leftPadding: 120,

      axisLabelFontsize: 11,
      axisFontcolor: '#444',
      axisTitleFontsize: 20,
      borderColor: '#444',
      rightPadding: 25,
      height: 160,
      profileTicks: 5,


      // header labels and arrows
      headerY: 5,
      headerLabelPadding: 10,
      headerLabelFontcolor: '#444',
      headerLabelFontsize: 11,
      headerArrowHeight: 12,
      headerArrowWidth: 6,
      headerArrowFontsize: 9,
      axisTitle: '% GC'

    },
    lineColors : {
      r: 'red',
      g: 'blue',
      b: 'green'
    },
    colorBlindLineColors : {
      r: 'rgb(213, 94, 0)',
      g: 'rgb(204, 121, 167)',
      b: 'rgb(0, 114, 178)'
    },
    // how much vertical space to leave for different kinds of headers
    headerSizes:{'extracts': 30, 'hits': 20}
  })
  .service('Err',function() {

    var self = this,
        makeError = function(message, name) {
          self[name] = _.partial(Error, message);
        };

    _.map({TrackAlreadyDefined:'Track with this name already defined, each track name must be unique',
           TrackNotFound: 'track not found',
           ProfileNotFound: 'profile not found'
          },
          makeError);

  })
  .constant('Evt', {
    REDRAW: 'redraw',
    NOOP: 'noop',
    REBUILD: 'rebuild',
    PAN: 'pan',
    ZOOM: 'zoom',
    GRAPH_REDRAW_COMPLETE: 'graph redraw complete'
  })
// these are expected to match strings specified in pynpact
  .constant('Pynpact', {
    NEW_CDS: 'File_of_new_CDSs',
    HITS: 'File_of_G+C_coding_potential_regions',
    PDF: 'pdf_filename',
    NPROFILE: 'nprofileData',
    CDS: 'File_of_published_accepted_CDSs',
    ACGT_GAMMA_FILES: 'acgt_gamma_output',
    TITLE: 'first_page_title',
    HAS_CDS: 'isgbk',
    EMAIL: 'email',
    CONFIGURE_URL: 'configure_url'
  })

;
