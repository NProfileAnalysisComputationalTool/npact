module.exports = function(config){
  'use strict';
  config.set({
    files: [
      'npactweb/static/bower_components/jquery/dist/jquery.min.js',
      'npactweb/static/bower_components/angular/angular.min.js',
      'npactweb/static/bower_components/angular-messages/angular-messages.min.js',
      'npactweb/static/bower_components/angular-sanitize/angular-sanitize.min.js',
      'npactweb/static/bower_components/angular-cookies/angular-cookies.min.js',
      'npactweb/static/bower_components/angular-mocks/angular-mocks.js',
      'npactweb/static/bower_components/angular-jquery-dialog-service/dialog-service.js',
      'npactweb/static/bower_components/konva/konva.min.js',
      'npactweb/static/bower_components/lodash/lodash.min.js',
      'npactweb/static/bower_components/rtree/dist/rtree.min.js',
      'npactweb/static/bower_components/ngSticky/dist/sticky.min.js',
      'npactweb/static/js/*.js',
      'npactweb/static/js/**/*.js',
      'npactweb/static/js/**/*.html',
      'npactweb/static/js/test-data/*'

    ],
    exclude: ['npactweb/static/js/*custom.min.js'],
    preprocessors: {
      'npactweb/static/js/**/*.html': ['ng-html2js'],
      'npactweb/static/js/test-data/*': ['ng-html2js']
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
