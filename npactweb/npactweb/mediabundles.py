# This defines the set of JS/CSS stuff that will be available with
# include_media

# The format for this file is described at
# http://www.allbuttonspressed.com/projects/django-mediagenerator#media-bundles

MEDIA_BUNDLES = (
    ('main.css',
     'css/basic.css',
     'css/custom-theme/jquery-ui-1.10.4.custom.css',
     'css/style.css',
     'bower_components/qtip2/jquery.qtip.min.css'),

    ('print.css',
     'css/print.css'),

    ('jquery.js',
     'bower_components/jquery/dist/jquery.min.js',
     # customized to only load UI elements we want
     'js/jquery-ui-1.10.4.custom.min.js',
     'bower_components/qtip2/jquery.qtip.min.js'),

    ('processing.js',
     'bower_components/kineticjs/kinetic.min.js',
     'bower_components/angular/angular.min.js',
     'bower_components/angular-messages/angular-messages.min.js',
     'bower_components/angular-sanitize/angular-sanitize.min.js',
     'bower_components/lodash/dist/lodash.min.js',
     'bower_components/qtip2/jquery.qtip.min.js',
     'bower_components/rtree/dist/rtree.min.js',
     'bower_components/ngSticky/dist/sticky.min.js',
     'js/app.js',
     'js/graphs/constants.js',
     'js/graphs/utils.js',
     'js/graphs/nprofiler.js',
     'js/graphs/readers.js',
     'js/graphs/page.js',
     'js/graphs/calculator.js',
     'js/graphs/grapher.js',
     'js/graphs/graph.js',
     'js/graphs/config.js')
    )
