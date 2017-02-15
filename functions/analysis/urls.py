from django.conf.urls import url

import functions

urlpatterns = [
    ### urls from normalization page
    url(r'^getNorm/$', functions.getNorm, name='getNorm'),
    url(r'^getTab/$', functions.getTab, name='getTab'),
    url(r'^getBiom/$', functions.getBiom, name='getBiome'),

    ### urls from analysis pages
    url(r'^getCatUnivData/$', functions.getCatUnivData, name='getCatUnivData'),
    url(r'^getQuantUnivData/$', functions.getQuantUnivData, name='getQuantUnivData'),
    url(r'^getDiffAbund/$', functions.getDiffAbund, name='getDiffAbund'),
    url(r'^getGAGE/$', functions.getGAGE, name='getGAGE'),
    url(r'^getPCA/$', functions.getPCA, name='getPCA'),
    url(r'^getPCoA/$', functions.getPCoA, name='getPCoA'),
    url(r'^statusPyBake/$', functions.statusPyBake, name='statusPyBake'),
    url(r'^getRF/$', functions.getRF, name='getRF'),
    url(r'^getsoil_index/$', functions.getsoil_index, name='getsoil_index'),
    url(r'^getSpAC/$', functions.getSpAC, name='getSpAC'),
    url(r'^getSPLS/$', functions.getSPLS, name='getSPLS'),
    url(r'^getWGCNA/$', functions.getWGCNA, name='getWGCNA')
]