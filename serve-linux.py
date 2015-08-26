#!/usr/bin/env python
import os
import os.path
import sys
import signal
import logging

import cherrypy
from cherrypy.process import plugins
from cherrypy import _cplogging, _cperror
from django.core.handlers.wsgi import WSGIHandler
from django.http import HttpResponseServerError
import webbrowser
import multiprocessing as mp


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

        cherrypy.tree.graft(HTTPLogger(WSGIHandler()))
        
        static_handler = cherrypy.tools.staticdir.handler(section="/", dir=staticpath, root='')
        cherrypy.tree.mount(static_handler, '/media')


class HTTPLogger(_cplogging.LogManager):
    def __init__(self, app):
        _cplogging.LogManager.__init__(self, id(self), cherrypy.log.logger_root)
        self.app = app

    def __call__(self, environ, start_response):
        try:
            response = self.app(environ, start_response)
            self.access(environ, response)
            return response
        except:
            self.error(traceback=True)
            return HttpResponseServerError(_cperror.format_exc())

    def access(self, environ, response):
        atoms = {'h': environ.get('REMOTE_ADDR', ''),
                 'l': '-',
                 'u': "-",
                 't': self.time(),
                 'r': "%s %s %s" % (environ['REQUEST_METHOD'], environ['REQUEST_URI'], environ['SERVER_PROTOCOL']),
                 's': response.status_code,
                 'b': str(len(response.content)),
                 'f': environ.get('HTTP_REFERER', ''),
                 'a': environ.get('HTTP_USER_AGENT', ''),
                 }
        for k, v in atoms.items():
            if isinstance(v, unicode):
                v = v.encode('utf8')
            elif not isinstance(v, str):
                v = str(v)

            v = repr(v)[1:-1]
            atoms[k] = v.replace('"', '\\"')

        try:
            self.access_log.log(logging.INFO, self.access_log_format % atoms)
        except:
            self.error(traceback=True)


def signal_handler(signal, frame):
    print 'Exiting...'
    cherrypy.engine.exit()
    sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)


if __name__ == '__main__':
    mp.freeze_support()
    Server().run()
