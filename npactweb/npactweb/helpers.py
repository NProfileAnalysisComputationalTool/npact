from django.utils.http import urlencode


def dict_to_querystring(dict):
    if len(dict):
        return "?" + urlencode(dict, True)
    else:
        return ""
