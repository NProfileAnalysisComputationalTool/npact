from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

import taskqueue
import taskqueue.tqdaemon

class Command(BaseCommand):
    args = '<run|start|stop|restart|status|kill>'
    help = 'Runs the background processing taskqueue daemon'

    def handle(self, *args, **options):
        taskqueue.BASE_DIR = settings.TQ_DIR
        if len(args) != 1:
            raise ArgumentError("Invalid command to tqdaemon.")
        else:
            cmd = getattr(taskqueue.tqdaemon, args[0], False)
            if cmd:
                return str(cmd()) + '\n'
