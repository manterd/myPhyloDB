# -*- mode: python -*-
a = Analysis(['serve-win.py'],
             pathex=['C:\\Users\\daniel.manter\\Documents\\GitHub\\myPhyloDB'],
             hiddenimports=[
                'scipy.special._ufuncs_cxx',
                'django.templatetags.future',
                'django.templatetags.i18n',
                'django.templatetags.__init__',
                'django.templatetags.cache',
                'django.templatetags.l10n',
                'django.templatetags.static',
                'django.templatetags.tz',
                'django.contrib.admin.templatetags.log',
                'django.contrib.admin.templatetags.__init__',
                'django.contrib.admin.templatetags.admin_list',
                'django.contrib.admin.templatetags.admin_modify',
                'django.contrib.admin.templatetags.admin_static',
                'django.contrib.admin.templatetags.admin_urls',
            ],
             hookspath=None,
             runtime_hooks=None)


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


db = Datafiles('dbMicrobe')
cfg = Datafiles('server.cfg')

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

a.datas += extra_datas('instructions')
a.datas += extra_datas('media')
a.datas += extra_datas('mothur/mothur-win')
a.datas += extra_datas('mothur/reference')
a.datas += extra_datas('R/R-Portable')
a.datas += extra_datas('sample_files')
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
               cfg,
               strip=None,
               upx=True,
               name='myPhyloDB'
               )

