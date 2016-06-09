import cherrypy
from cherrypy.process import plugins
import multiprocessing as mp
import os.path
import signal
import sys
import threading
import webbrowser
from myPhyloDB.wsgi import application

from database.queue import process
from database.utils import threads
from database.dataqueue import dataprocess


cherrypy.tree.graft(application)


class Server(object):
    def __init__(self):
        self.base_dir = os.path.join(os.path.abspath(os.getcwd()), "myPhyloDB")

        cherrypy.config.update("config/server.cfg")
        DjangoAppPlugin(cherrypy.engine, self.base_dir).subscribe()

    def browse(self):
        url = ''
        f = open("config/server.cfg")
        lines = f.readlines()
        for line in lines:
            if "server.socket_port: " in line:
                port = line.split(' ')[1]
                url = "http://127.0.0.1:" + str(port.rstrip('\n')) + "/myPhyloDB/home/"
        webbrowser.open_new(url)

    def run(self):
        engine = cherrypy.engine

        if hasattr(engine, "signal_handler"):
            engine.signal_handler.subscribe()
        if hasattr(engine, "console_control_handler"):
            engine.console_control_handler.subscribe()

        cherrypy.engine.subscribe('engine.start', Server.browse(self), priority=90)

        engine.start()
        engine.block()


class DjangoAppPlugin(plugins.SimplePlugin):
    def __init__(self, bus, base_dir):
        plugins.SimplePlugin.__init__(self, bus)
        self.base_dir = base_dir

    def start(self):
        import django
        django.setup()

        from config.local_cfg import update
        update()

        staticpath = os.path.abspath(self.base_dir)
        staticpath = os.path.split(staticpath)[0]
        staticpath = os.path.join(staticpath, 'templates')
        static_handler = cherrypy.tools.staticdir.handler(section="/", dir=staticpath, root='')
        cherrypy.tree.mount(static_handler, '/templates')

        staticpath = os.path.abspath(self.base_dir)
        staticpath = os.path.split(staticpath)[0]
        staticpath = os.path.join(staticpath, 'media')
        static_handler = cherrypy.tools.staticdir.handler(section="/", dir=staticpath, root='')
        cherrypy.tree.mount(static_handler, '/media')

        staticpath = os.path.abspath(self.base_dir)
        staticpath = os.path.split(staticpath)[0]
        staticpath = os.path.join(staticpath, 'sample_files')
        static_handler = cherrypy.tools.staticdir.handler(section="/", dir=staticpath, root='')
        cherrypy.tree.mount(static_handler, '/sample_files')

        staticpath = os.path.abspath(self.base_dir)
        staticpath = os.path.split(staticpath)[0]
        staticpath = os.path.join(staticpath, 'instructions')
        static_handler = cherrypy.tools.staticdir.handler(section="/", dir=staticpath, root='')
        cherrypy.tree.mount(static_handler, '/instructions')


def signal_handler(signal, frame):
    print 'Exiting...'
    cherrypy.engine.exit()


signal.signal(signal.SIGINT, signal_handler)
num_threads = threads()


if __name__ == '__main__':
    import django.core.management
    #from django.core.management import execute_from_command_line
    #execute_from_command_line(sys.argv)

    mp.freeze_support()

    for pid in xrange(num_threads):
        thread = threading.Thread(target=process, args=(pid, ))
        thread.setDaemon(True)
        thread.start()

    dataThread = threading.Thread(target=dataprocess, args=(0, ))
    dataThread.setDaemon(True)
    dataThread.start()

    Server().run()

