'use strict';

module.exports = function(config){
  config.set({
    files:[
      'npactweb/static/bower_components/jquery/jquery.min.js',
      'npactweb/static/bower_components/angular/angular.js',
      'npactweb/static/bower_components/angular-mocks/angular-mocks.js',
      'npactweb/static/bower_components/kineticjs/kinetic.js',
      '.tmp/templates.js',
      'npactweb/static/js/*.js',
      'npactweb/static/js/**/*.js'
    ],
    exclude: ['npactweb/static/js/*custom.min.js'],
    frameworks:['jasmine'],
    browsers: ['PhantomJS'],
    reporters: ['spec'],
    colors:true,
    singleRun: true
  });
};
