'use strict';

module.exports = function(grunt){
  require('load-grunt-tasks')(grunt);

  grunt.initConfig({
    karma: {
      unit: {
	configFile: 'karma.conf.js'
      }
    },
    html2js: {
      options:{
	singleModule:true
      },
      main: {
	src: ['npactweb/static/js/**/*.html'],
	dest: '.tmp/templates.js'
      }
    }
  });

  grunt.registerTask('test', ['newer:html2js:main', 'karma']);
};
