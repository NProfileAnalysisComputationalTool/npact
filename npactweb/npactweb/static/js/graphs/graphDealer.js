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
      height: 200,
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
;
