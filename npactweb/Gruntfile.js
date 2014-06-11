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
	singleModule:true,
	rename: function(n){ return n.replace('../', '/'); }
      },
      main: {
	src: ['npactweb/static/js/**/*.html'],
	dest: '.tmp/templates.js'
      }
    }
  });

  grunt.registerTask('test', ['newer:html2js:main', 'karma']);
};
