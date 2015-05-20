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

db = Datafiles('dbMicrobe', strip_path=False)

templates = Datafiles(
    'templates/anova.html',
    'templates/base.html',
    'templates/home.html',
    'templates/login.html',
    'templates/pcoa.html',
    'templates/select.html',
    'templates/taxa.html',
    'templates/upload.html',
    'templates/users.html',
    'templates/admin/base_site.html',
    strip_path=False
    )

media = Datafiles(
    'media/css/style.css',
    'media/datatables/css/jquery.dataTables.css',
    'media/datatables/css/jquery.dataTables_themeroller.css',
    'media/datatables/images/back_disabled.png',
    'media/datatables/images/back_enabled.png',
    'media/datatables/images/back_enabled_hover.png',
    'media/datatables/images/favicon.ico',
    'media/datatables/images/forward_disabled.png',
    'media/datatables/images/forward_enabled.png',
    'media/datatables/images/forward_enabled_hover.png',
    'media/datatables/images/sort_asc.png',
    'media/datatables/images/sort_asc_disabled.png',
    'media/datatables/images/sort_both.png',
    'media/datatables/images/sort_desc.png',
    'media/datatables/images/sort_desc_disabled.png',
    'media/datatables/images/Sorting icons.psd',
    'media/datatables/js/jquery.dataTables.js',
    'media/dynatree/src/skin/icons.gif',
    'media/dynatree/src/skin/icons-rtl.gif',
    'media/dynatree/src/skin/loading.gif',
    'media/dynatree/src/skin/ui.dynatree.css',
    'media/dynatree/src/skin/vline.gif',
    'media/dynatree/src/skin/vline-rtl.gif',
    'media/dynatree/src/jquery.dynatree.js',
    'media/dynatree/jquery.cookie.js',
    'media/highcharts/themes/dark-blue.js',
    'media/highcharts/themes/dark-green.js',
    'media/highcharts/themes/dark-unica.js',
    'media/highcharts/themes/gray.js',
    'media/highcharts/themes/grid.js',
    'media/highcharts/themes/grid-light.js',
    'media/highcharts/themes/sand-signika.js',
    'media/highcharts/themes/skies.js',
    'media/highcharts/exporting.js',
    'media/highcharts/grouped-categories.js',
    'media/highcharts/highcharts.js',
    'media/images/ARS.jpeg',
    'media/images/database_2_48.ico',
    'media/images/database_2_48.png',
    'media/images/dbimage.jpeg',
    'media/images/DNA.jpeg',
    'media/images/helix.jpeg',
    'media/images/sf_menu.png',
    'media/jquery/jquery-1.11.1.min.js',
    'media/jquery-ui/jquery-ui.css',
    'media/jquery-ui/jquery-ui.js',
    'media/sample_files/final.shared',
    'media/sample_files/final.taxonomy',
    'media/sample_files/Project.csv',
    'media/sample_files/Sample.csv',
    'media/tabletools/css/dataTables.tableTools.css',
    'media/tabletools/images/psd/collection.psd',
    'media/tabletools/images/psd/copy document.psd',
    'media/tabletools/images/psd/file_types.psd',
    'media/tabletools/images/psd/printer.psd',
    'media/tabletools/images/background.png',
    'media/tabletools/images/collection.png',
    'media/tabletools/images/collection_hover.png',
    'media/tabletools/images/copy.png',
    'media/tabletools/images/copy_hover.png',
    'media/tabletools/images/csv.png',
    'media/tabletools/images/csv_hover.png',
    'media/tabletools/images/pdf.png',
    'media/tabletools/images/pdf_hover.png',
    'media/tabletools/images/print.png',
    'media/tabletools/images/print_hover.png',
    'media/tabletools/images/xls.png',
    'media/tabletools/images/xls_hover.png',
    'media/tabletools/js/dataTables.tableTools.js',
    'media/tabletools/swf/copy_csv_xls_pdf.swf',
    strip_path=False
    )

sample_files = Datafiles(
    'sample_files/Example1.project.csv',
    'sample_files/Example1.sample.csv',
    'sample_files/Example1.taxonomy',
    'sample_files/Example1.shared',
    'sample_files/Example2.project.csv',
    'sample_files/Example2.sample.csv',
    'sample_files/Example2.batch',
    'sample_files/Example2.oligos',
    'sample_files/Example2.sff',
    strip_path=False
    )

mothur_files = Datafiles(
    'mothur/mothur.exe',
    'mothur/uchime.exe',
    'mothur/lookupFiles/LookUp_GS20.pat',
    'mothur/lookupFiles/LookUp_GSFLX.pat',
    'mothur/lookupFiles/LookUp_Titanium.pat',
    'mothur/reference/gg_13_5_99.fasta',
    'mothur/reference/gg_13_5_99.pds.tax',
    'mothur/reference/silva.archaea.fasta',
    'mothur/reference/silva.bacteria.fasta',
    'mothur/reference/silva.eukarya.fasta',
    'mothur/reference/silva.gold.ng.fasta',
    strip_path=False
    )

instructions = Datafiles(
    'instructions/Manual.pdf',
    strip_path=False
    )

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
			   templates,
			   media,
			   sample_files,
			   instructions,
			   mothur_files,
               strip=None,
               upx=True,
               name='myPhyloDB')

