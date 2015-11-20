module.exports = function(config){
  'use strict';
  config.set({
    files: [
      'npactflask/static/bower_components/jquery/dist/jquery.min.js',
      'npactflask/static/bower_components/angular/angular.min.js',
      'npactflask/static/bower_components/angular-messages/angular-messages.min.js',
      'npactflask/static/bower_components/angular-sanitize/angular-sanitize.min.js',
      'npactflask/static/bower_components/angular-cookies/angular-cookies.min.js',
      'npactflask/static/bower_components/angular-mocks/angular-mocks.js',
      'npactflask/static/bower_components/konva/konva.min.js',
      'npactflask/static/bower_components/lodash/lodash.min.js',
      'npactflask/static/bower_components/rtree/dist/rtree.min.js',
      'npactflask/static/bower_components/ngSticky/dist/sticky.min.js',
      'http://cdnjs.cloudflare.com/ajax/libs/angular-ui-bootstrap/0.13.1/ui-bootstrap.min.js',
      'http://cdnjs.cloudflare.com/ajax/libs/angular-ui-bootstrap/0.13.1/ui-bootstrap-tpls.min.js',
      'npactflask/static/js/*.js',
      'npactflask/static/js/**/*.js',
      'npactflask/static/js/**/*.html',
      'npactflask/static/js/test-data/*'

    ],
    exclude: ['npactflask/static/js/*custom.min.js'],
    preprocessors: {
      'npactflask/static/js/**/*.html': ['ng-html2js'],
      'npactflask/static/js/test-data/*': ['ng-html2js']
    },
    frameworks:['jasmine'],
    browsers: ['PhantomJS'],
    reporters: ['spec'],
    colors:true,
    singleRun: false,
    logLevel: config.LOG_DEBUG,
    ngHtml2JsPreprocessor: {
      stripPrefix: 'npactflask/static',
      moduleName: 'assets'
    }
  });
};
