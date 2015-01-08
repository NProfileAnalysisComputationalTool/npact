module.exports = function(config){
  'use strict';
  config.set({
    files:[
      'npactweb/static/bower_components/jquery/dist/jquery.min.js',
      'npactweb/static/bower_components/angular/angular.min.js',
      'npactweb/static/bower_components/angular-messages/angular-messages.min.js',
      'npactweb/static/bower_components/angular-mocks/angular-mocks.js',
      'npactweb/static/bower_components/kineticjs/kinetic.min.js',
      'npactweb/static/bower_components/lodash/dist/lodash.min.js',
      'npactweb/static/bower_components/rtree/dist/rtree.min.js',
      'npactweb/static/bower_components/ngSticky/sticky.min.js',
      'npactweb/static/js/*.js',
      'npactweb/static/js/**/*.js',
      'npactweb/static/js/**/*.html',
      'npactweb/static/js/**/*.txt',
      'npactweb/static/js/**/*.json',
      'npactweb/static/js/**/*.ddna',

    ],
    exclude: ['npactweb/static/js/*custom.min.js'],
    preprocessors: {
      'npactweb/static/js/**/*.html': ['ng-html2js'],
      'npactweb/static/js/**/*.txt': ['ng-html2js'],
      'npactweb/static/js/**/*.json': ['ng-html2js'],
      'npactweb/static/js/**/*.ddna': ['ng-html2js']
    },
    frameworks:['jasmine'],
    browsers: ['PhantomJS'],
    reporters: ['spec'],
    colors:true,
    singleRun: false,
    logLevel: config.LOG_DEBUG,
    ngHtml2JsPreprocessor: {
      stripPrefix: 'npactweb/static',
      moduleName: 'assets'
    }
  });
};
