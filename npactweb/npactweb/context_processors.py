"""Add context processors for template processing"""

from django.core.urlresolvers import resolve
from npactweb import getrelpath


def resolvermatch(request):
    """Add `resolved`-- the current url's ResolverMatch object to context

    Example use: Set class=active to a menu item, based on the
    namespace of the currently matched url.

        <li {% if resolved.namespace == "home" %}class="active"{% endif %}>
            home
        </li>

    or more specifically:

        <li {% if resolved.view_name == "contact" %}class="active"{% endif %}>
            contact
        </li>

    See https://docs.djangoproject.com/en/1.3/topics/http/urls/#django.core.urlresolvers.ResolverMatch
    """
    return {'resolved': resolve(request.path)}


def addgetrelpath(request):
    return {'getrelpath': lambda p: p and getrelpath(p)}
