import taskqueue
import taskqueue.tqdaemon
from npactflask import app

class Command(BaseCommand):
    args = '<run|start|stop|restart|status|kill>'
    help = 'Runs the background processing taskqueue daemon'

    def handle(self, *args, **options):
        taskqueue.BASE_DIR = app.config['TQ_DIR']
        if len(args) != 1:
            raise ArgumentError("Invalid command to tqdaemon.")
        else:
            cmd = getattr(taskqueue.tqdaemon, args[0], False)
            if cmd:
                return str(cmd()) + '\n'
