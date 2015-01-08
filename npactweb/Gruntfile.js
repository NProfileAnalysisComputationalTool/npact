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
      options: {
	configFile: 'karma.conf.js'
      },
      unit: {
        reporters: 'dots'
      },
      jenkins: {
        reporters: 'junit',
        outputFile: 'test-results.xml',
        singleRun: true
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

  grunt.registerTask('test', ['karma']);
  grunt.registerTask('jenkins', ['jshint:jenkins',
                                 'html2js:main',
                                 'karma:jenkins']);
};
