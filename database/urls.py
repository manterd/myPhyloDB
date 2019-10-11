from django.conf.urls import url

import views

urlpatterns = [
    url(r'^home/$', views.home, name='home'),
    url(r'^process/$', views.process, name='process'),
    url(r'^download/$', views.download, name='download'),
    url(r'^select/$', views.select, name='select'),
    url(r'^taxa/$', views.taxa, name='taxa'),
    url(r'^kegg_path/$', views.kegg_path, name='kegg_path'),
    url(r'^kegg_enzyme/$', views.kegg_enzyme, name='kegg_enzyme'),
    url(r'^anova/$', views.ANOVA, name='anova'),
    url(r'^corr/$', views.CORR, name='corr'),
    url(r'^spac/$', views.spac, name='spac'),
    url(r'^soil_index/$', views.soil_index, name='soil_index'),
    url(r'^norm/$', views.norm, name='norm'),
    url(r'^diffabund/', views.DiffAbund, name='diffabund'),
    url(r'^gage/', views.GAGE, name='gage'),
    url(r'^pca/$', views.PCA, name='pca'),
    url(r'^pcoa/$', views.PCoA, name='pcoa'),
    url(r'^caret/$', views.RF, name='rf'),
    url(r'^spls/$', views.SPLS, name='spls'),
    url(r'^wgcna/$', views.WGCNA, name='wgcna'),
    url(r'^reprocess/$', views.reprocess, name='reprocess'),
    url(r'^update/$', views.update, name='update'),
    url(r'^pybake/$', views.pybake, name='pybake'),
    url(r'^upload/$', views.files, name='files'),
    url(r'^console/$', views.admin_console, name='console'),
    url(r'^history/$', views.history, name='history'),
    url(r'^core/$', views.core, name='core'),
]
