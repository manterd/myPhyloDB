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
    url(r'^myPhyloDB/export/$', 'database.views.export', name='export'),
    url(r'^myPhyloDB/DiffAbund/$', 'database.views.DiffAbund', name='DiffAbund'),
    url(r'^myPhyloDB/PCoA/$', 'database.views.PCoA', name='pcoa'),
    url(r'^myPhyloDB/sPLS/$', 'database.views.SPLS', name='spls'),
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

    url(r'^getExCatData/$', 'database.export.export_graphs.getExCatData', name='getExCatData'),
    url(r'^statusExport/$', 'database.export.export_graphs.statusExport', name='statusExport'),
    url(r'^statusRIDExport/$', 'database.export.export_graphs.removeRIDExport', name='removeRIDExport'),

    url(r'^getCatUnivData/$', 'database.anova.anova_graphs.getCatUnivData', name='getCatUnivData'),
    url(r'^getQuantUnivData/$', 'database.anova.anova_graphs.getQuantUnivData', name='getQuantUnivData'),
    url(r'^statusANOVA/$', 'database.anova.anova_graphs.statusANOVA', name='statusANOVA'),
    url(r'^removeRIDANOVA/$', 'database.anova.anova_graphs.removeRIDANOVA', name='removeRIDANOVA'),
    url(r'^stopANOVA/$', 'database.anova.anova_graphs.stopANOVA', name='stopANOVA'),

    url(r'^getPCoA/$', 'database.pcoa.pcoa_graphs.getPCoA', name='getPCoA'),
    url(r'^statusPCoA/$', 'database.pcoa.pcoa_graphs.statusPCoA', name='statusPCoA'),
    url(r'^removeRIDPCoA/$', 'database.pcoa.pcoa_graphs.removeRIDPCoA', name='removeRIDPCoA'),
    url(r'^removegraphPCoA/$', 'database.pcoa.pcoa_graphs.removegraphPCoA', name='removegraphPCoA'),

    url(r'^getSPLS/$', 'database.spls.spls_graphs.getSPLS', name='getSPLS'),
    url(r'^statusSPLS/$', 'database.spls.spls_graphs.statusSPLS', name='statusSPLS'),
    url(r'^removeRIDSPLS/$', 'database.spls.spls_graphs.removeRIDSPLS', name='removeRIDSPLS'),
    url(r'^removegraphSPLS/$', 'database.spls.spls_graphs.removegraphSPLS', name='removegraphSPLS'),

    url(r'^updateDiffAbund/$', 'database.diffabund.diffabund_graphs.updateDiffAbund', name='updateDiffAbund'),
    url(r'^getDiffAbund/$', 'database.diffabund.diffabund_graphs.getDiffAbund', name='getDiffAbund'),
    url(r'^removeRIDDIFF/$', 'database.diffabund.diffabund_graphs.removeRIDDIFF', name='removeRIDDIFF'),
)


