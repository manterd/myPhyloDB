from django.conf.urls import include, url
from django.contrib import admin
from registration.backends.simple.views import RegistrationView

from database import views, parsers, trees
from database.anova import anova_graphs
from database.diffabund import diffabund_graphs
from database.norm import norm_graphs
from database.pcoa import pcoa_graphs
from database.spls import spls_graphs


admin.autodiscover()


class MyRegistrationView(RegistrationView):
    def get_success_url(self, request=None, user=None):
        return '/myPhyloDB/select/'


urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^accounts/', include('registration.urls')),
    url(r'^myPhyloDB/register/$', MyRegistrationView.as_view(), name='registration_register'),
    url(r'^myPhyloDB/', include('database.urls')),

    url(r'^saveSampleCookie/$', views.saveSampleCookie, name='saveSampleCookie'),
    url(r'^getSampleCookie/$', views.getSampleCookie, name='getSampleCookie'),
    url(r'^clearNormCookie/$', views.clearNormCookie, name='clearNormCookie'),

    url(r'^status/$', parsers.status, name='status'),
    url(r'^reanalyze/$', parsers.reanalyze, name='reanalyze'),

    url(r'^getProjectTree/$', trees.getProjectTree, name='getProjectTree'),
    url(r'^getProjectTreeChildren/$', trees.getProjectTreeChildren, name='getProjectTreeChildren'),
    url(r'^getSampleCatTree/$', trees.getSampleCatTree, name='getSampleCatTree'),
    url(r'^getSampleCatTreeChildren/$', trees.getSampleCatTreeChildren, name='getSampleCatTreeChildren'),
    url(r'^getSampleQuantTree/$', trees.getSampleQuantTree, name='getSampleQuantTree'),
    url(r'^getSampleQuantTreeChildren/$', trees.getSampleQuantTreeChildren, name='getSampleQuantTreeChildren'),
    url(r'^getTaxaTree/$', trees.getTaxaTree, name='getTaxaTree'),
    url(r'^getTaxaTreeChildren/$', trees.getTaxaTreeChildren, name='getTaxaTreeChildren'),
    url(r'^makeUpdateTree/$', trees.makeUpdateTree, name='makeUpdateTree'),
    url(r'^makeReproTree/$', trees.makeReproTree, name='makeReproTree'),

    url(r'^getNormCatData/$', norm_graphs.getNormCatData, name='getNormCatData'),
    url(r'^statusNorm/$', norm_graphs.statusNorm, name='statusNorm'),
    url(r'^removeRIDNorm/$', norm_graphs.removeRIDNorm, name='removeRIDNorm'),
    url(r'^stopNorm/$', norm_graphs.stopNorm, name='stopNorm'),

    url(r'^getCatUnivData/$', anova_graphs.getCatUnivData, name='getCatUnivData'),
    url(r'^getQuantUnivData/$', anova_graphs.getQuantUnivData, name='getQuantUnivData'),
    url(r'^statusANOVA/$', anova_graphs.statusANOVA, name='statusANOVA'),
    url(r'^removeRIDANOVA/$', anova_graphs.removeRIDANOVA, name='removeRIDANOVA'),
    url(r'^stopANOVA/$', anova_graphs.stopANOVA, name='stopANOVA'),

    url(r'^getPCoA/$', pcoa_graphs.getPCoA, name='getPCoA'),
    url(r'^statusPCoA/$', pcoa_graphs.statusPCoA, name='statusPCoA'),
    url(r'^removeRIDPCoA/$', pcoa_graphs.removeRIDPCoA, name='removeRIDPCoA'),
    url(r'^removegraphPCoA/$', pcoa_graphs.removegraphPCoA, name='removegraphPCoA'),
    url(r'^stopPCoA/$', pcoa_graphs.stopPCoA, name='stopPCoA'),

    url(r'^getSPLS/$', spls_graphs.getSPLS, name='getSPLS'),
    url(r'^statusSPLS/$', spls_graphs.statusSPLS, name='statusSPLS'),
    url(r'^removeRIDSPLS/$', spls_graphs.removeRIDSPLS, name='removeRIDSPLS'),
    url(r'^removegraphSPLS/$', spls_graphs.removegraphSPLS, name='removegraphSPLS'),
    url(r'^stopSPLS/$', spls_graphs.stopSPLS, name='stopSPLS'),

    url(r'^updateDiffAbund/$', diffabund_graphs.updateDiffAbund, name='updateDiffAbund'),
    url(r'^getDiffAbund/$', diffabund_graphs.getDiffAbund, name='getDiffAbund'),
    url(r'^removeRIDDIFF/$', diffabund_graphs.removeRIDDIFF, name='removeRIDDIFF'),
    url(r'^stopDiffAbund/$', diffabund_graphs.stopDiffAbund, name='stopDiffAbund'),

]
