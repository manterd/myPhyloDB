# -*- mode: python -*-

a = Analysis([
        'serve-linux.py'],
        pathex=[
            '/home/manterd/PycharmProjects/myPhyloDB'
            ],
            hiddenimports=[
                'django.conf.urls',
                'django.conf.urls.shortcut',
                'django.contrib.admin',
                'django.contrib.admin.templatetags.admin_list',
                'django.contrib.admin.templatetags.admin_modify',
                'django.contrib.admin.templatetags.log',
                'django.contrib.admin.views.auth',
                'django.contrib.admin.views.doc',
                'django.contrib.admin.views.main',
                'django.contrib.admin.views.template',
                'django.contrib.auth',
                'django.contrib.auth.backends',
                'django.contrib.auth.middleware',
                'django.contrib.auth.views',
                'django.contrib.contenttypes',
                'django.contrib.messages.middleware',
                'django.contrib.sessions',
                'django.contrib.sessions.backends.db',
                'django.contrib.sessions.middleware',
                'django.contrib.sites',
                'django.core.cache.backends',
                'django.core.cache.backends.locmem',
                'django.core.context_processors',
                'django.db.backends.sqlite3.base',
                'django.db.backends.sqlite3.client',
                'django.db.backends.sqlite3.creation',
                'django.db.backends.sqlite3.introspection',
                'django.middleware.common',
                'django.middleware.doc',
                'django.template.defaultfilters',
                'django.template.defaulttags',
                'django.template.loader_tags',
                'django.template.loaders.app_directories',
                'django.template.loaders.filesystem',
                'django.templatetags.i18n',
                'django.views.defaults',
                'django.views.generic.list_detail',
                'django.views.i18n'
                'django.views.static',
                'registration.admin',
                'registration.forms',
                'registration.urls',
            ],
            hookspath=None,
            runtime_hooks=None,
            excludes=[
                    '_gtkagg',
                    '_tkagg',
                    '_agg2',
                    '_cairo',
                    '_cocoaagg',
                    '_fltkagg',
                    '_gtk',
                    '_gtkcairo',
                    'backend_qt',
                    'backend_qt4',
                    'backend_qtagg'
                    'backend_cairo',
                    'backend_cocoagg',
                    'Tkconstants',
                    'Tkinter',
                    'tcl',
                    '_imagingtk',
                    'PIL._imagingtk',
                    'ImageTk',
                    'PIL.ImageTk',
                    'TixTk'
                    ],
                )

def Datafiles(*filenames, **kw):
    import os

    def datafile(path, strip_path=True):
        parts = path.split('/')
        path = name = os.path.join(*parts)
        if strip_path:
            name = os.path.basename(path)
        return name, path, 'DATA'

    strip_path = kw.get('strip_path', True)
    return TOC(
        datafile(filename, strip_path=strip_path)
        for filename in filenames
        if os.path.isfile(filename))


db = Datafiles('db.Microbe')

def extra_datas(mydir):
    def rec_glob(p, files):
        import os
        import glob
        for d in glob.glob(p):
            if os.path.isfile(d):
                files.append(d)
            rec_glob("%s/*" % d, files)
    files = []
    rec_glob("%s/*" % mydir, files)
    extra_datas = []
    for f in files:
        extra_datas.append((f, f, 'DATA'))

    return extra_datas

a.datas += extra_datas('config')
a.datas += extra_datas('instructions/current')
a.datas += extra_datas('media')
a.datas += extra_datas('mothur/mothur-linux')
a.datas += extra_datas('mothur/reference')
a.datas += extra_datas('R/R-Linux')
a.datas += extra_datas('templates')
a.datas += extra_datas('uploads')

pyz = PYZ(a.pure)

exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='myPhyloDB.exe',
          debug=False,
          strip=None,
          upx=True,
          console=True
          )

coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               db,
               strip=None,
               upx=True,
               name='myPhyloDB'
               )

