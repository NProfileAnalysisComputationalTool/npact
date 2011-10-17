import logging
import os.path

from django.conf import settings

logger = logging.getLogger(__name__)

logger.info("Starting the main spat site.")



def library_root() :
    abspath = getabspath("library",False)
    if not os.path.exists(abspath):
        logger.info("Creating library root at %s", abspath)
        os.makedirs(abspath)
    return abspath


def getrelpath(abspath):
    return os.path.relpath(abspath, settings.MEDIA_ROOT)

def is_clean_path(path) :
    return path.find('..') < 0

def getabspath(relpath, raise_on_missing=True):
    if not is_clean_path(relpath):
        logger.error("Illegal path submitted", relpath)
        raise Exception("Illegal path")
    abspath = os.path.join(settings.MEDIA_ROOT, relpath)
    if raise_on_missing and not os.path.exists(abspath):
        raise IOError("Path not found: " + relpath)
    return abspath

