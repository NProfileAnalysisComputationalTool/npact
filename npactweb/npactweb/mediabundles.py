# This defines the set of JS/CSS stuff that will be available with
# include_media

# The format for this file is described at
# http://www.allbuttonspressed.com/projects/django-mediagenerator#media-bundles

MEDIA_BUNDLES= (
    ('main.css',
     'css/basic.css',
     'css/custom-theme/jquery-ui-1.10.4.custom.css',
     'css/style.css',
     'qtip2/jquery.qtip.css'),

    ('print.css',
     'css/print.css'),

    ('jquery.js',
     'js/jquery-1.10.2.js',
     'js/jquery-ui-1.10.4.custom.min.js',
     'qtip2/jquery.qtip.js'),

    ('processing.js',
     'js/processing.js')
    )
