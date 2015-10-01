from django.conf.urls import *
from django.contrib import admin
from registration.backends.simple.views import RegistrationView

admin.autodiscover()


class MyRegistrationView(RegistrationView):
    def get_success_url(self, request=None, user=None):
        return '/myPhyloDB/select/'


urlpatterns = patterns('',
    url(r'^admin/', include(admin.site.urls)),
    url(r'^accounts/', include('registration.backends.simple.urls')),
    url(r'^myPhyloDB/register/$', MyRegistrationView.as_view(), name='registration_register'),

    url(r'^myPhyloDB/home/$', 'database.views.home', name='home'),
    url(r'^myPhyloDB/upload/$', 'database.views.upload', name='upload'),
    url(r'^myPhyloDB/select/$', 'database.views.select', name='select'),
    url(r'^myPhyloDB/taxa/$', 'database.views.taxa', name='taxa'),
    url(r'^myPhyloDB/ANOVA/$', 'database.views.ANOVA', name='anova'),
    url(r'^myPhyloDB/DiffAbund/$', 'database.views.DiffAbund', name='DiffAbund'),
    url(r'^myPhyloDB/PCoA/$', 'database.views.PCoA', name='pcoa'),
    url(r'^myPhyloDB/reprocess/$', 'database.views.reprocess', name='reprocess'),
    url(r'^myPhyloDB/update/$', 'database.views.update', name='update'),
    url(r'^saveCookie/$', 'database.views.saveCookie', name='saveCookie'),
    url(r'^getCookie/$', 'database.views.getCookie', name='getCookie'),

    url(r'^status/$', 'database.parsers.status', name='status'),
    url(r'^reanalyze/$', 'database.parsers.reanalyze', name='reanalyze'),

    url(r'^getProjectTree/$', 'database.trees.getProjectTree', name='getProjectTree'),
    url(r'^getProjectTreeChildren/$', 'database.trees.getProjectTreeChildren', name='getProjectTreeChildren'),
    url(r'^getSampleCatTree/$', 'database.trees.getSampleCatTree', name='getSampleCatTree'),
    url(r'^getSampleCatTreeChildren/$', 'database.trees.getSampleCatTreeChildren', name='getSampleCatTreeChildren'),
    url(r'^getSampleQuantTree/$', 'database.trees.getSampleQuantTree', name='getSampleQuantTree'),
    url(r'^getSampleQuantTreeChildren/$', 'database.trees.getSampleQuantTreeChildren', name='getSampleQuantTreeChildren'),
    url(r'^getTaxaTree/$', 'database.trees.getTaxaTree', name='getTaxaTree'),
    url(r'^getTaxaTreeChildren/$', 'database.trees.getTaxaTreeChildren', name='getTaxaTreeChildren'),
    url(r'^makeUpdateTree/$', 'database.trees.makeUpdateTree', name='makeUpdateTree'),
    url(r'^makeReproTree/$', 'database.trees.makeReproTree', name='makeReproTree'),

    url(r'^getCatUnivData/$', 'database.anova_graphs.getCatUnivData', name='getCatUnivData'),
    url(r'^getQuantUnivData/$', 'database.anova_graphs.getQuantUnivData', name='getQuantUnivData'),
    url(r'^statusANOVA/$', 'database.anova_graphs.statusANOVA', name='statusANOVA'),
    url(r'^removeRIDANOVA/$', 'database.anova_graphs.removeRIDANOVA', name='removeRIDANOVA'),

    url(r'^getCatPCoAData/$', 'database.pcoa_graphs.getCatPCoAData', name='getCatPCoAData'),
    url(r'^statusPCoA/$', 'database.pcoa_graphs.statusPCoA', name='statusPCoA'),
    url(r'^removeRIDPCOA/$', 'database.pcoa_graphs.removeRIDPCOA', name='removeRIDPCOA'),

    url(r'^updateDiffAbund/$', 'database.diffabund_graphs.updateDiffAbund', name='updateDiffAbund'),
    url(r'^getDiffAbund/$', 'database.diffabund_graphs.getDiffAbund', name='getDiffAbund'),
    url(r'^removeRIDDIFF/$', 'database.diffabund_graphs.removeRIDDIFF', name='removeRIDDIFF'),
)


