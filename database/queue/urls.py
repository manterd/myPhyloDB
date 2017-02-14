from django.conf.urls import url

import database.queue

urlpatterns = [
    ### urls from analysis queue
    url(r'^funcCall/$', database.queue.funcCall, name='funcCall'),
    url(r'^stop/$', database.queue.stop, name='stop'),
    url(r'^removeRID/$', database.queue.removeRID, name='removeRID'),

    ### data from upload/reprocessing queue
    url(r'^datfuncCall/$', database.queue.datfuncCall, name='datfuncCall'),
    url(r'^datstop/$', database.queue.datstop, name='datstop'),
]