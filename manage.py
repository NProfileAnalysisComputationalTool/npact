#!ve/bin/python
import os

os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'  # or whatever

from django.core import management

if __name__ == "__main__":
    management.execute_from_command_line()
