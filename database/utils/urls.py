from django.conf.urls import url

import database.utils

urlpatterns = [
    ### urls from parser page
    url(r'^status/$', database.utils.status, name='status'),
    url(r'^reanalyze/$', database.utils.reanalyze, name='reanalyze'),

    ### urls from trees page
    url(r'^getProjectTree/$', database.utils.getProjectTree, name='getProjectTree'),
    url(r'^getProjectTreeChildren/$', database.utils.getProjectTreeChildren, name='getProjectTreeChildren'),
    url(r'^getSampleCatTree/$', database.utils.getSampleCatTree, name='getSampleCatTree'),
    url(r'^getSampleCatTreeChildren/$', database.utils.getSampleCatTreeChildren, name='getSampleCatTreeChildren'),
    url(r'^getSampleQuantTree/$', database.utils.getSampleQuantTree, name='getSampleQuantTree'),
    url(r'^getSampleQuantTreeChildren/$', database.utils.getSampleQuantTreeChildren, name='getSampleQuantTreeChildren'),
    url(r'^getTaxaTree/$', database.utils.getTaxaTree, name='getTaxaTree'),
    url(r'^getTaxaTreeChildren/$', database.utils.getTaxaTreeChildren, name='getTaxaTreeChildren'),
    url(r'^getKEGGTree/$', database.utils.getKEGGTree, name='getKEGGTree'),
    url(r'^getKEGGTreeChildren/$', database.utils.getKEGGTreeChildren, name='getKEGGTreeChildren'),
    url(r'^getKEGGTree2/$', database.utils.getKEGGTree2, name='getKEGGTree2'),
    url(r'^getNZTree/$', database.utils.getNZTree, name='getNZTree'),
    url(r'^getNZTreeChildren/$', database.utils.getNZTreeChildren, name='getNZTreeChildren'),
    url(r'^makeUpdateTree/$', database.utils.makeUpdateTree, name='makeUpdateTree'),
    url(r'^makeReproTree/$', database.utils.makeReproTree, name='makeReproTree'),
    url(r'^getDownloadTree/$', database.utils.getDownloadTree, name='getDownloadTree'),
    url(r'^getDownloadTreeChildren/$', database.utils.getDownloadTreeChildren, name='getDownloadTreeChildren'),
    url(r'^getPermissionTree/$', database.utils.getPermissionTree, name='getPermissionTree'),

    ### urls from utils page
    url(r'^getRawData/$', database.utils.getRawData, name='getRawData'),
    url(r'^removeFiles/$', database.utils.removeFiles, name='removeFiles'),
]
