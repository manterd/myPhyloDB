from django.conf.urls import *
from django.contrib import admin
admin.autodiscover()


urlpatterns = patterns('',
    (r'^myPhyloDB/admin/', include(admin.site.urls)),
    (r'^myPhyloDB/login/$', 'django.contrib.auth.views.login', {'template_name': 'login.html'}),
    url(r'^myPhyloDB/logout/$', 'django.contrib.auth.views.logout', name='logout'),
    url(r'^myPhyloDB/home/$', 'database.views.home', name='home'),
    url(r'^myPhyloDB/upload/$', 'database.views.upload', name='upload'),
    url(r'^myPhyloDB/select/$', 'database.views.select', name='select'),
    url(r'^myPhyloDB/taxa/$', 'database.views.taxa', name='taxa'),
    url(r'^myPhyloDB/ANOVA/$', 'database.views.ANOVA', name='anova'),
    url(r'^myPhyloDB/DiffAbund/$', 'database.views.DiffAbund', name='DiffAbund'),
    url(r'^myPhyloDB/PCoA/$', 'database.views.PCoA', name='pcoa'),
    url(r'^myPhyloDB/users/$', 'database.views.users', name='users'),
    url(r'^myPhyloDB/reprocess/$', 'database.views.reprocess', name='reprocess'),
    url(r'^myPhyloDB/update/$', 'database.views.update', name='update'),

    url(r'^saveCookie/$', 'database.views.saveCookie', name='saveCookie'),
    url(r'^getCookie/$', 'database.views.getCookie', name='getCookie'),

    url(r'^reanalyze/$', 'database.parsers.reanalyze', name='reanalyze'),

    url(r'^getProjectTree/$', 'database.trees.getProjectTree', name='getProjectTree'),
    url(r'^getProjectTreeChildren/$', 'database.trees.getProjectTreeChildren', name='getProjectTreeChildren'),
    url(r'^makeUpdateTree/$', 'database.trees.makeUpdateTree', name='makeUpdateTree'),
    url(r'^makeReproTree/$', 'database.trees.makeReproTree', name='makeReproTree'),

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

    url(r'^statusANOVA', 'database.anova_graphs.statusANOVA', name='statusANOVA'),
    url(r'^statusPCoA', 'database.pcoa_graphs.statusPCoA', name='statusPCoA'),
    url(r'^status', 'database.parsers.status', name='status'),
    url(r'^updateDiffAbund', 'database.diffabund_graphs.updateDiffAbund', name='updateDiffAbund'),

    url('getDiffAbund', 'database.diffabund_graphs.getDiffAbund', name='getDiffAbund'),

    url('removeRIDANOVA', 'database.anova_graphs.removeRIDANOVA', name='removeRIDANOVA'),
    url('removeRIDDIFF', 'database.diffabund_graphs.removeRIDDIFF', name='removeRIDDIFF'),
    url('removeRIDPCOA', 'database.pcoa_graphs.removeRIDPCOA', name='removeRIDPCOA'),

)


