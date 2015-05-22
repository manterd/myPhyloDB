# -*- mode: python -*-
a = Analysis(['serve.py'],
             pathex=['C:\\Users\\daniel.manter\\Documents\\GitHub\\myPhyloDB'],
             #pathex=['/home/manterd/PycharmProjects/myPhyloDB'],
             hiddenimports=[
                'scipy.special._ufuncs_cxx',
                'django.templatetags.future',
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
a.datas += extra_datas('media/admin/img')
a.datas += extra_datas('media/css')
a.datas += extra_datas('media/datatables/css')
a.datas += extra_datas('media/datatables/images')
a.datas += extra_datas('media/datatables/js')
a.datas += extra_datas('media/dynatree')
a.datas += extra_datas('media/dynatree/src')
a.datas += extra_datas('media/dynatree/src/skin')
a.datas += extra_datas('media/highcharts')
a.datas += extra_datas('media/highcharts/themes')
a.datas += extra_datas('media/images')
a.datas += extra_datas('media/jquery')
a.datas += extra_datas('media/jquery-ui')
a.datas += extra_datas('media/tabletools/css')
a.datas += extra_datas('media/tabletools/images')
a.datas += extra_datas('media/tabletools/images/psd')
a.datas += extra_datas('media/tabletools/js')
a.datas += extra_datas('media/tabletools/swf')
a.datas += extra_datas('mothur/mothur-linux')
a.datas += extra_datas('mothur/mothur-linux/lookupFiles')
a.datas += extra_datas('mothur/mothur-linux/uchime_src')
a.datas += extra_datas('mothur/mothur-win')
a.datas += extra_datas('mothur/mothur-win/lookupFiles')
a.datas += extra_datas('mothur/reference')
a.datas += extra_datas('sample_files')
a.datas += extra_datas('templates')
a.datas += extra_datas('templates/admin')

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

