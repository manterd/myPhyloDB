import os
import sys

os.environ['DJANGO_SETTINGS_MODULE'] = 'myPhyloDB.settings'

import django
django.setup()

import cherrypy
import multiprocessing as mp
import signal
import threading
import webbrowser

from myPhyloDB.wsgi import application


class Server(object):
    def __init__(self):
        self.base_dir = os.path.join(os.path.abspath(os.getcwd()), "myPhyloDB")

        cherrypy.config.update("config/server.cfg")
        DjangoAppPlugin(cherrypy.engine, self.base_dir).subscribe()

    def browse(self):
        f = open("config/server.cfg")
        lines = f.readlines()
        ip = ''
        port = ''
        for line in lines:
            if line.startswith('server.socket_host:'):
                ip = line.split('"')[1]
            if line.startswith('server.socket_port:'):
                port = line.split(' ')[1]
        url = "http://" + str(ip.rstrip('\n')) + ":" + str(port.rstrip('\n')) + '/myPhyloDB/home/'
        webbrowser.open_new(url)

    def run(self):
        engine = cherrypy.engine

        if hasattr(engine, "signal_handler"):
            engine.signal_handler.subscribe()
        elif hasattr(engine, "console_control_handler"):
            engine.console_control_handler.subscribe()

        cherrypy.engine.subscribe('engine.start', Server.browse(self), priority=90)

        engine.start()
        engine.block()


class DjangoAppPlugin(cherrypy.process.plugins.SimplePlugin):
    def __init__(self, bus, base_dir):
        cherrypy.process.plugins.SimplePlugin.__init__(self, bus)
        self.base_dir = base_dir

    def start(self):
        staticpath = os.path.abspath(self.base_dir)
        staticpath = os.path.split(staticpath)[0]
        staticpath = os.path.join(staticpath, 'myPhyloDB/media')
        static_handler = cherrypy.tools.staticdir.handler(section="/", dir=staticpath, root='')
        cherrypy.tree.mount(static_handler, '/myPhyloDB/media')

        staticpath = os.path.abspath(self.base_dir)
        staticpath = os.path.split(staticpath)[0]
        staticpath = os.path.join(staticpath, 'uploads')
        static_handler = cherrypy.tools.staticdir.handler(section="/", dir=staticpath, root='')
        cherrypy.tree.mount(static_handler, '/uploads')

        staticpath = os.path.abspath(self.base_dir)
        staticpath = os.path.split(staticpath)[0]
        staticpath = os.path.join(staticpath, 'instructions')
        static_handler = cherrypy.tools.staticdir.handler(section="/", dir=staticpath, root='')
        cherrypy.tree.mount(static_handler, '/instructions')

        cherrypy.tree.graft(application)


def signal_handler(signal, frame):
    print 'Exiting...'
    cherrypy.engine.exit()


if __name__ == '__main__':
    from django.core.management import execute_from_command_line
    execute_from_command_line(sys.argv)

    from functions.utils.utils_df import analysisThreads
    from functions.queues.queue import process
    from functions.queues.dataqueue import dataprocess

    signal.signal(signal.SIGINT, signal_handler)
    num_threads = analysisThreads()

    for pid in xrange(num_threads):
        thread = threading.Thread(target=process, args=(pid, ))
        thread.setDaemon(True)
        thread.start()

    dataThread = threading.Thread(target=dataprocess, args=(0, ))
    dataThread.setDaemon(True)
    dataThread.start()

    mp.freeze_support()
    Server().run()

