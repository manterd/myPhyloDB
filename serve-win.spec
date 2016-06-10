# -*- mode: python -*-

a = Analysis(
            ['serve-win.py'],
            pathex=[
                'C:\\Users\\daniel.manter\\Documents\\GitHub\\myPhyloDB',
                'C:\\Users\\daniel.manter\\AppData\\Local\\Continuum\\Anaconda2\\Lib\\site-packages'
                ],
            hiddenimports=[
                'django.apps',
                'django.apps.config',
                'django.apps.registry',
                'django.contrib.admin.templatetags.log',
                'django.contrib.admin.templatetags.__init__',
                'django.contrib.admin.templatetags.admin_list',
                'django.contrib.admin.templatetags.admin_modify',
                'django.contrib.admin.templatetags.admin_static',
                'django.contrib.admin.templatetags.admin_urls',
                'django.contrib.auth',
                'django.contrib.auth.backends',
                'django.contrib.auth.middleware',
                'django.contrib.auth.views',
                'django.contrib.messages.middleware',
                'django.contrib.sessions',
                'django.contrib.sessions.backends.db',
                'django.contrib.sessions.middleware',
                'django.core.cache.backends',
                'django.core.cache.backends.locmem',
                'django.template.defaulttags',
                'django.templatetags.__init__',
                'django.templatetags.cache',
                'django.templatetags.l10n',
                'django.templatetags.static',
                'django.templatetags.tz',
                'django.templatetags.future',
                'django.views.defaults',
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


db1 = Datafiles('db.Microbe')
db2 = Datafiles('db.PICRUSt')

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
a.datas += extra_datas('database/migrations')
a.datas += extra_datas('registration/migrations')
a.datas += extra_datas('instructions/current')
a.datas += extra_datas('media')
a.datas += extra_datas('mothur/mothur-win')
a.datas += extra_datas('mothur/reference/align')
a.datas += extra_datas('mothur/reference/taxonomy')
a.datas += extra_datas('mothur/reference/template')
a.datas += extra_datas('R/R-Portable')
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
               db1,
               db2,
               strip=None,
               upx=True,
               name='myPhyloDB'
               )

