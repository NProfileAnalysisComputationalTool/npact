module.exports = function(grunt){
  'use strict';
  require('load-grunt-tasks')(grunt);
  grunt.loadNpmTasks('grunt-contrib-jshint');
  var jsFiles = ['npactweb/static/js/**/*.js'];

  grunt.initConfig({
    jshint: {
      options: {
        jshintrc: true
      },
      all: jsFiles,
      jenkins: {
        files: {src: jsFiles},
        options: {
          reporter: 'checkstyle',
          reporterOutput: 'jshint.xml',
          force: true
        }
      }
    },
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
