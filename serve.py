#!/usr/bin/env python
import os
import os.path
import sys
import cherrypy
from cherrypy.process import plugins
from django.core.handlers.wsgi import WSGIHandler
import webbrowser
import multiprocessing as mp

#CHANGE!!
class Server(object):
    def __init__(self):
        self.base_dir = os.path.join(os.path.abspath(os.getcwd()), "phyloDB")

        conf = {
            'server.socket_host': "0.0.0.0",
            'server.socket_port': 8000,
            'server.thread_pool': 10,
            'checker.on': False,
            'engine.autoreload.on': False,
        }

        cherrypy.config.update(conf)
        DjangoAppPlugin(cherrypy.engine, self.base_dir).subscribe()

    def browse(self):
        webbrowser.open_new("http://127.0.0.1:8000/myPhyloDB/home")

    def run(self):
        engine = cherrypy.engine
        if hasattr(engine, "signal_handler"):
            engine.signal_handler.subscribe()
        if hasattr(engine, "console_control_handler"):
            engine.console_control_handler.subscribe()
        cherrypy.engine.subscribe('engine.start', Server.browse(self), priority=90)
        cherrypy.log.screen = None
        engine.start()
        engine.block()


class DjangoAppPlugin(plugins.SimplePlugin):
    def __init__(self, bus, base_dir):
        plugins.SimplePlugin.__init__(self, bus)
        self.base_dir = base_dir

    def start(self):
        os.environ.setdefault("DJANGO_SETTINGS_MODULE", "phyloDB.settings")
        import django.test
        import HTMLParser
        import Cookie
        import django.contrib.sessions.serializers

        cherrypy.tree.graft(WSGIHandler())
        staticpath = os.path.abspath(self.base_dir)
        staticpath = os.path.split(staticpath)[0]
        staticpath = os.path.join(staticpath, 'media')

        static_handler = cherrypy.tools.staticdir.handler(section="/", dir=staticpath, root='')
        cherrypy.tree.mount(static_handler, '/media')

if __name__ == '__main__':
    mp.freeze_support()
    Server().run()