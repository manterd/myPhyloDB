from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^home/$', views.home, name='home'),
    url(r'^upload/$', views.upload, name='upload'),
    url(r'^select/$', views.select, name='select'),
    url(r'^taxa/$', views.taxa, name='taxa'),
    url(r'^anova/$', views.ANOVA, name='anova'),
    url(r'^norm/$', views.norm, name='norm'),
    url(r'^diffabund/', views.DiffAbund, name='diffabund'),
    url(r'^pcoa/$', views.PCoA, name='pcoa'),
    url(r'^spls/$', views.SPLS, name='spls'),
    url(r'^reprocess/$', views.reprocess, name='reprocess'),
    url(r'^update/$', views.update, name='update'),
]
