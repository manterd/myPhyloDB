from django.conf.urls import url

import functions

urlpatterns = [
    ### urls from parser page
    url(r'^status/$', functions.status, name='status'),
    url(r'^reanalyze/$', functions.reanalyze, name='reanalyze'),

    ### urls from trees page
    url(r'^getProjectTree/$', functions.getProjectTree, name='getProjectTree'),
    url(r'^getProjectTreeChildren/$', functions.getProjectTreeChildren, name='getProjectTreeChildren'),
    url(r'^getSampleCatTree/$', functions.getSampleCatTree, name='getSampleCatTree'),
    url(r'^getSampleCatTreeChildren/$', functions.getSampleCatTreeChildren, name='getSampleCatTreeChildren'),
    url(r'^getSampleQuantTree/$', functions.getSampleQuantTree, name='getSampleQuantTree'),
    url(r'^getSampleQuantTreeChildren/$', functions.getSampleQuantTreeChildren, name='getSampleQuantTreeChildren'),
    url(r'^getTaxaTree/$', functions.getTaxaTree, name='getTaxaTree'),
    url(r'^getTaxaTreeChildren/$', functions.getTaxaTreeChildren, name='getTaxaTreeChildren'),
    url(r'^getKEGGTree/$', functions.getKEGGTree, name='getKEGGTree'),
    url(r'^getKEGGTreeChildren/$', functions.getKEGGTreeChildren, name='getKEGGTreeChildren'),
    url(r'^getKEGGTree2/$', functions.getKEGGTree2, name='getKEGGTree2'),
    url(r'^getNZTree/$', functions.getNZTree, name='getNZTree'),
    url(r'^getNZTreeChildren/$', functions.getNZTreeChildren, name='getNZTreeChildren'),
    url(r'^makeUpdateTree/$', functions.makeUpdateTree, name='makeUpdateTree'),
    url(r'^makeFilesTree/$', functions.makeFilesTree, name='makeFilesTree'),
    url(r'^makeReproTree/$', functions.makeReproTree, name='makeReproTree'),
    url(r'^getDownloadTree/$', functions.getDownloadTree, name='getDownloadTree'),
    url(r'^getDownloadTreeChildren/$', functions.getDownloadTreeChildren, name='getDownloadTreeChildren'),
    url(r'^getPermissionTree/$', functions.getPermissionTree, name='getPermissionTree'),
    url(r'^getFilePermTree/$', functions.getFilePermTree, name='getFilePermTree'),
    url(r'^getLocationSamplesTree/$', functions.getLocationSamplesTree, name='getLocationSamplesTree'),

    ### urls from utils page
    url(r'^getRawDataTab/$', functions.getRawDataTab, name='getRawDataTab'),
    url(r'^getRawDataBiom/$', functions.getRawDataBiom, name='getRawDataBiom'),
    url(r'^removeFiles/$', functions.removeFiles, name='removeFiles'),
]
