from django.core.management.base import BaseCommand, CommandError
import taskqueue

class Command(BaseCommand):
    args = ''
    help = 'Runs the background processing taskqueue daemon'

    def handle(self, *args, **options):
        taskqueue.daemonize()
