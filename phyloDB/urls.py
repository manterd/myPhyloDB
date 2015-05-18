from django.conf.urls import *
from django.contrib import admin
admin.autodiscover()


urlpatterns = patterns('',
    (r'^myPhyloDB/admin/', include(admin.site.urls)),
    (r'^myPhyloDB/login/$', 'django.contrib.auth.views.login', {'template_name': 'login.html'}),
    url(r'^myPhyloDB/home/$', 'database.views.home', name='home'),
    url(r'^myPhyloDB/upload/$', 'database.views.upload', name='upload'),
    url(r'^myPhyloDB/select/$', 'database.views.select', name='select'),
    url(r'^myPhyloDB/taxa/$', 'database.views.taxa', name='taxa'),
    url(r'^myPhyloDB/ANOVA/$', 'database.views.ANOVA', name='anova'),
    url(r'^myPhyloDB/PCoA/$', 'database.views.PCoA', name='pcoa'),
    url(r'^myPhyloDB/users/$', 'database.views.users', name='users'),

    url(r'^saveCookie/$', 'database.views.saveCookie', name='saveCookie'),
    url(r'^getCookie/$', 'database.views.getCookie', name='getCookie'),

    url(r'^getProjectTree/$', 'database.trees.getProjectTree', name='getProjectTree'),
    url(r'^getProjectTreeChildren/$', 'database.trees.getProjectTreeChildren', name='getProjectTreeChildren'),

    url(r'^getSampleCatTree/$', 'database.trees.getSampleCatTree', name='getSampleCatTree'),
    url(r'^getSampleCatTreeChildren/$', 'database.trees.getSampleCatTreeChildren', name='getSampleCatTreeChildren'),

    url(r'^getSampleQuantTree/$', 'database.trees.getSampleQuantTree', name='getSampleQuantTree'),
    url(r'^getSampleQuantTreeChildren/$', 'database.trees.getSampleQuantTreeChildren', name='getSampleQuantTreeChildren'),

    url(r'^getTaxaTree/$', 'database.trees.getTaxaTree', name='getTaxaTree'),
    url(r'^getTaxaTreeChildren/$', 'database.trees.getTaxaTreeChildren', name='getTaxaTreeChildren'),

    url(r'^getCatUnivData', 'database.anova_graphs.getCatUnivData', name='getCatUnivData'),
    url(r'^getQuantUnivData', 'database.anova_graphs.getQuantUnivData', name='getQuantUnivData'),
    url(r'^getCatPCoAData', 'database.pcoa_graphs.getCatPCoAData', name='getCatPCoAData'),
    url(r'^getQuantPCoAData', 'database.pcoa_graphs.getQuantPCoAData', name='getQuantPCoAData'),

    url(r'^instructions', 'database.views.instructions', name='instructions'),
    url(r'^project_file', 'database.views.project_file', name='project_file'),
    url(r'^sample_file', 'database.views.sample_file', name='sample_file'),
    url(r'^shared_file', 'database.views.shared_file', name='shared_file'),
    url(r'^taxonomy_file', 'database.views.taxonomy_file', name='taxonomy_file'),

    url(r'^statusANOVA', 'database.anova_graphs.statusANOVA', name='statusANOVA'),
    url(r'^statusPCoA', 'database.pcoa_graphs.statusPCoA', name='statusPCoA'),
    url(r'^status', 'database.parsers.status', name='status'),
)


