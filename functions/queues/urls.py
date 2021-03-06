from django.conf.urls import url

import functions

urlpatterns = [
    ### urls from analysis queue
    url(r'^funcCall/$', functions.funcCall, name='funcCall'),
    url(r'^stop/$', functions.stop, name='stop'),
    url(r'^getAnalysisQueue/$', functions.getAnalysisQueue, name='getAnalysisQueue'),
    url(r'^getAnalysisHistory/$', functions.getAnalysisHistory, name='getAnalysisHistory'),

    ### data from upload/reprocessing queue
    url(r'^datfuncCall/$', functions.datfuncCall, name='datfuncCall'),
    url(r'^datstop/$', functions.datstop, name='datstop'),
    url(r'^getDataQueue/$', functions.getDataQueue, name='getDataQueue'),
]
