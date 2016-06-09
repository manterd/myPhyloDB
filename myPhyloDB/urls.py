from django.conf.urls import include, url
from django.contrib import admin
from registration.backends.simple.views import RegistrationView

from database import views, parsers, trees, queue
from database.anova import anova_graphs
from database.spac import spac_graphs
from database.diffabund import diffabund_graphs
from database.gage import gage_graphs
from database.norm import norm_graphs
from database.pca import pca_graphs
from database.pcoa import pcoa_graphs
from database.spls import spls_graphs
from database.soil_index import soil_index_graphs
from database.wgcna import wgcna_graphs
from database.pybake import pybake

admin.autodiscover()


class MyRegistrationView(RegistrationView):
    def get_success_url(self, request=None, user=None):
        return '/myPhyloDB/select/'


urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^accounts/', include('registration.urls')),
    url(r'^myPhyloDB/register/$', MyRegistrationView.as_view(), name='registration_register'),
    url(r'^myPhyloDB/', include('database.urls')),

    url(r'^saveSampleList/$', views.saveSampleList, name='saveSampleList'),

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

    url(r'^getKEGGTree/$', trees.getKEGGTree, name='getKEGGTree'),
    url(r'^getKEGGTreeChildren/$', trees.getKEGGTreeChildren, name='getKEGGTreeChildren'),

    url(r'^getKEGGTree2/$', trees.getKEGGTree2, name='getKEGGTree2'),

    url(r'^getNZTree/$', trees.getNZTree, name='getNZTree'),
    url(r'^getNZTreeChildren/$', trees.getNZTreeChildren, name='getNZTreeChildren'),

    url(r'^makeUpdateTree/$', trees.makeUpdateTree, name='makeUpdateTree'),
    url(r'^makeReproTree/$', trees.makeReproTree, name='makeReproTree'),

    url(r'^getNorm/$', norm_graphs.getNorm, name='getNorm'),
    url(r'^statusNorm/$', norm_graphs.statusNorm, name='statusNorm'),
    url(r'^getTab/$', norm_graphs.getTab, name='getTab'),
    url(r'^getBiom/$', norm_graphs.getBiom, name='getBiome'),

    url(r'^getCatUnivData/$', anova_graphs.getCatUnivData, name='getCatUnivData'),
    url(r'^getQuantUnivData/$', anova_graphs.getQuantUnivData, name='getQuantUnivData'),
    url(r'^statusANOVA/$', anova_graphs.statusANOVA, name='statusANOVA'),
    url(r'^removeRIDANOVA/$', anova_graphs.removeRIDANOVA, name='removeRIDANOVA'),
    url(r'^getTabANOVA/$', anova_graphs.getTabANOVA, name='getTabANOVA'),
    url(r'^removeANOVAFiles/$', anova_graphs.removeANOVAFiles, name='removeANOVAFiles'),

    url(r'^getSpAC/$', spac_graphs.getSpAC, name='getSpAC'),
    url(r'^statusSpAC/$', spac_graphs.statusSpAC, name='statusSpAC'),
    url(r'^removeRIDSpAC/$', spac_graphs.removeRIDSpAC, name='removeRIDSpAC'),
    url(r'^removeSpACFiles/$', spac_graphs.removeSpACFiles, name='removeSpACFiles'),
    url(r'^getTabSpAC/$', spac_graphs.getTabSpAC, name='getTabSpAC'),

    url(r'^getsoil_index/$', soil_index_graphs.getsoil_index, name='getsoil_index'),
    url(r'^statussoil_index/$', soil_index_graphs.statussoil_index, name='statussoil_index'),
    url(r'^removeRIDsoil_index/$', soil_index_graphs.removeRIDsoil_index, name='removeRIDsoil_index'),
    url(r'^removesoil_indexFiles/$', soil_index_graphs.removesoil_indexFiles, name='removesoil_indexFiles'),
    url(r'^getTabsoil_index/$', soil_index_graphs.getTabsoil_index, name='getTabsoil_index'),

    url(r'^getPCoA/$', pcoa_graphs.getPCoA, name='getPCoA'),
    url(r'^statusPCoA/$', pcoa_graphs.statusPCoA, name='statusPCoA'),
    url(r'^removeRIDPCoA/$', pcoa_graphs.removeRIDPCoA, name='removeRIDPCoA'),
    url(r'^removePCoAFiles/$', pcoa_graphs.removePCoAFiles, name='removePCoAFiles'),
    url(r'^getTabPCoA/$', pcoa_graphs.getTabPCoA, name='getTabPCoA'),

    url(r'^getWGCNA/$', wgcna_graphs.getWGCNA, name='getWGCNA'),
    url(r'^statusWGCNA/$', wgcna_graphs.statusWGCNA, name='statusWGCNA'),
    url(r'^removeRIDWGCNA/$', wgcna_graphs.removeRIDWGCNA, name='removeRIDWGCNA'),
    url(r'^removeWGCNAFiles/$', wgcna_graphs.removeWGCNAFiles, name='removeWGCNAFiles'),
    url(r'^getTabWGCNA/$', wgcna_graphs.getTabWGCNA, name='getTabWGCNA'),

    url(r'^getPCA/$', pca_graphs.getPCA, name='getPCA'),
    url(r'^statusPCA/$', pca_graphs.statusPCA, name='statusPCA'),
    url(r'^removeRIDPCA/$', pca_graphs.removeRIDPCA, name='removeRIDPCA'),
    url(r'^removePCAFiles/$', pca_graphs.removePCAFiles, name='removePCAFiles'),
    url(r'^getTabPCA/$', pca_graphs.getTabPCA, name='getTabPCA'),

    url(r'^getSPLS/$', spls_graphs.getSPLS, name='getSPLS'),
    url(r'^statusSPLS/$', spls_graphs.statusSPLS, name='statusSPLS'),
    url(r'^removeRIDSPLS/$', spls_graphs.removeRIDSPLS, name='removeRIDSPLS'),
    url(r'^removeSPLSFiles/$', spls_graphs.removeSPLSFiles, name='removeSPLSFiles'),
    url(r'^getTabSPLS/$', spls_graphs.getTabSPLS, name='getTabSPLS'),

    url(r'^updateDiffAbund/$', diffabund_graphs.updateDiffAbund, name='updateDiffAbund'),
    url(r'^getDiffAbund/$', diffabund_graphs.getDiffAbund, name='getDiffAbund'),
    url(r'^removeRIDDIFF/$', diffabund_graphs.removeRIDDIFF, name='removeRIDDIFF'),
    url(r'^getTabDiffAbund/$', diffabund_graphs.getTabDiffAbund, name='getTabDiffAbund'),
    url(r'^removeDiffAbundFiles/$', diffabund_graphs.removeDiffAbundFiles, name='removeDiffAbundFiles'),

    url(r'^updateGAGE/$', gage_graphs.updateGAGE, name='updateGAGE'),
    url(r'^getGAGE/$', gage_graphs.getGAGE, name='getGAGE'),
    url(r'^removeRIDGAGE/$', gage_graphs.removeRIDGAGE, name='removeRIDGAGE'),
    url(r'^getTabGAGE/$', gage_graphs.getTabGAGE, name='getTabGAGE'),
    url(r'^removeGAGEFiles/$', gage_graphs.removeGAGEFiles, name='removeGAGEFiles'),

    url(r'^statusPyBake/$', pybake.statusPyBake, name='statusPyBake'),

    url(r'^taxaJSON/$', views.taxaJSON, name='taxaJSON'),
    url(r'^nzJSON/$', views.nzJSON, name='nzJSON'),
    url(r'^pathJSON/$', views.pathJSON, name='pathJSON'),

    url(r'^myPhyloDB/select$', views.uploadNorm, name='uploadNorm'),

    # queue caller
    url(r'^funcCall/$', queue.funcCall, name='funcCall'),
    url(r'^stop/$', queue.stop, name='stop')
]
