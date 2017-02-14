from django.conf.urls import url

import database.analysis

urlpatterns = [
    ### urls from normalization page
    url(r'^getNorm/$', database.analysis.getNorm, name='getNorm'),
    url(r'^getTab/$', database.analysis.getTab, name='getTab'),
    url(r'^getBiom/$', database.analysis.getBiom, name='getBiome'),

    ### urls from analysis pages
    url(r'^getCatUnivData/$', database.analysis.getCatUnivData, name='getCatUnivData'),
    url(r'^getQuantUnivData/$', database.analysis.getQuantUnivData, name='getQuantUnivData'),
    url(r'^getDiffAbund/$', database.analysis.getDiffAbund, name='getDiffAbund'),
    url(r'^getGAGE/$', database.analysis.getGAGE, name='getGAGE'),
    url(r'^getPCA/$', database.analysis.getPCA, name='getPCA'),
    url(r'^getPCoA/$', database.analysis.getPCoA, name='getPCoA'),
    url(r'^statusPyBake/$', database.analysis.statusPyBake, name='statusPyBake'),
    url(r'^getRF/$', database.analysis.getRF, name='getRF'),
    url(r'^getsoil_index/$', database.analysis.getsoil_index, name='getsoil_index'),
    url(r'^getSpAC/$', database.analysis.getSpAC, name='getSpAC'),
    url(r'^getSPLS/$', database.analysis.getSPLS, name='getSPLS'),
    url(r'^getWGCNA/$', database.analysis.getWGCNA, name='getWGCNA')
]