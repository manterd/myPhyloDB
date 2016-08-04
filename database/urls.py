from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^home/$', views.home, name='home'),
    url(r'^upload/$', views.upload, name='upload'),
    url(r'^select/$', views.select, name='select'),
    url(r'^taxa/$', views.taxa, name='taxa'),
    url(r'^kegg_path/$', views.kegg_path, name='kegg_path'),
    url(r'^kegg_enzyme/$', views.kegg_enzyme, name='kegg_enzyme'),
    url(r'^anova/$', views.ANOVA, name='anova'),
    url(r'^spac/$', views.rich, name='spac'),
    url(r'^soil_index/$', views.soil_index, name='soil_index'),
    url(r'^norm/$', views.norm, name='norm'),
    url(r'^diffabund/', views.DiffAbund, name='diffabund'),
    url(r'^gage/', views.GAGE, name='gage'),
    url(r'^pca/$', views.PCA, name='pca'),
    url(r'^pcoa/$', views.PCoA, name='pcoa'),
    url(r'^spls/$', views.SPLS, name='spls'),
    url(r'^wgcna/$', views.WGCNA, name='wgcna'),
    url(r'^reprocess/$', views.reprocess, name='reprocess'),
    url(r'^update/$', views.update, name='update'),
    url(r'^pybake/$', views.pybake, name='pybake'),
    url(r'^accounts/profile/$', views.profile, name='profile'),
    url(r'^accounts/changeuser/$', views.changeuser, name='changeuser'),
]
