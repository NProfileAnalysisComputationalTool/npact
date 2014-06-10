
# About
This is the code for the website.

# Important files:

(in this folder)

* `settings.py` is the configuration for this django
  website. Tells what modules and applications to load; configures
  where paths for different content are.

* `npactweb/static/` holds the static images, CSS, and JS used on the
  website.

* `npactweb/templates/` holds the templates that generate html on the
  website. Go here to change text on the site, rearrange blocks, or
  add new content.

* `npactweb/urls.py` controls the web sites url space; django uses
  this to dispatch incoming requests to views.

* `npactweb/views/` holds the code that gets run on an
  incoming request.

* `cleanup.py` used to clean out old uploaded files so they
  don't take too much space on the server.

* `django.fcgi` is the entry point for fcgi webservers.
