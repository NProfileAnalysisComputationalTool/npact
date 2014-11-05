'use strict';

module.exports = function(config){
  config.set({
    files:[
      'npactweb/static/bower_components/jquery/dist/jquery.min.js',
      'npactweb/static/bower_components/angular/angular.min.js',
      'npactweb/static/bower_components/angular-messages/angular-messages.min.js',
      'npactweb/static/bower_components/angular-mocks/angular-mocks.js',
      'npactweb/static/bower_components/ngInfiniteScroll/build/ng-infinite-scroll.min.js',
      'npactweb/static/bower_components/kineticjs/kinetic.min.js',
      'npactweb/static/bower_components/lodash/dist/lodash.min.js',
      '.tmp/templates.js',
      'npactweb/static/js/*.js',
      'npactweb/static/js/**/*.js'
    ],
    exclude: ['npactweb/static/js/*custom.min.js'],
    frameworks:['jasmine'],
    browsers: ['PhantomJS'],
    reporters: ['spec'],
    colors:true,
    singleRun: false
  });
};
