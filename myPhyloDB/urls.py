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
from database import dataqueue

admin.autodiscover()


class MyRegistrationView(RegistrationView):
    def get_success_url(self, request=None, user=None):
        return '/myPhyloDB/select/'


urlpatterns = [
    url(r'^myPhyloDB/admin/', admin.site.urls),
    url(r'^myPhyloDB/accounts/', include('registration.urls')),
    url(r'^myPhyloDB/register/$', MyRegistrationView.as_view(), name='registration_register'),
    url(r'^myPhyloDB/', include('database.urls')),

    url(r'^myPhyloDB/saveSampleList/$', views.saveSampleList, name='saveSampleList'),

    url(r'^myPhyloDB/status/$', parsers.status, name='status'),
    url(r'^myPhyloDB/reanalyze/$', parsers.reanalyze, name='reanalyze'),

    url(r'^myPhyloDB/getProjectTree/$', trees.getProjectTree, name='getProjectTree'),
    url(r'^myPhyloDB/getProjectTreeChildren/$', trees.getProjectTreeChildren, name='getProjectTreeChildren'),

    url(r'^myPhyloDB/getSampleCatTree/$', trees.getSampleCatTree, name='getSampleCatTree'),
    url(r'^myPhyloDB/getSampleCatTreeChildren/$', trees.getSampleCatTreeChildren, name='getSampleCatTreeChildren'),
    url(r'^myPhyloDB/getSampleQuantTree/$', trees.getSampleQuantTree, name='getSampleQuantTree'),
    url(r'^myPhyloDB/getSampleQuantTreeChildren/$', trees.getSampleQuantTreeChildren, name='getSampleQuantTreeChildren'),

    url(r'^myPhyloDB/getTaxaTree/$', trees.getTaxaTree, name='getTaxaTree'),
    url(r'^myPhyloDB/getTaxaTreeChildren/$', trees.getTaxaTreeChildren, name='getTaxaTreeChildren'),

    url(r'^myPhyloDB/getKEGGTree/$', trees.getKEGGTree, name='getKEGGTree'),
    url(r'^myPhyloDB/getKEGGTreeChildren/$', trees.getKEGGTreeChildren, name='getKEGGTreeChildren'),

    url(r'^myPhyloDB/getKEGGTree2/$', trees.getKEGGTree2, name='getKEGGTree2'),

    url(r'^myPhyloDB/getNZTree/$', trees.getNZTree, name='getNZTree'),
    url(r'^myPhyloDB/getNZTreeChildren/$', trees.getNZTreeChildren, name='getNZTreeChildren'),

    url(r'^myPhyloDB/makeUpdateTree/$', trees.makeUpdateTree, name='makeUpdateTree'),
    url(r'^myPhyloDB/makeReproTree/$', trees.makeReproTree, name='makeReproTree'),

    url(r'^myPhyloDB/getNorm/$', norm_graphs.getNorm, name='getNorm'),
    url(r'^myPhyloDB/statusNorm/$', norm_graphs.statusNorm, name='statusNorm'),
    url(r'^myPhyloDB/getTab/$', norm_graphs.getTab, name='getTab'),
    url(r'^myPhyloDB/getBiom/$', norm_graphs.getBiom, name='getBiome'),

    url(r'^myPhyloDB/getCatUnivData/$', anova_graphs.getCatUnivData, name='getCatUnivData'),
    url(r'^myPhyloDB/getQuantUnivData/$', anova_graphs.getQuantUnivData, name='getQuantUnivData'),
    url(r'^myPhyloDB/statusANOVA/$', anova_graphs.statusANOVA, name='statusANOVA'),
    url(r'^myPhyloDB/removeRIDANOVA/$', anova_graphs.removeRIDANOVA, name='removeRIDANOVA'),
    url(r'^myPhyloDB/getTabANOVA/$', anova_graphs.getTabANOVA, name='getTabANOVA'),
    url(r'^myPhyloDB/removeANOVAFiles/$', anova_graphs.removeANOVAFiles, name='removeANOVAFiles'),

    url(r'^myPhyloDB/getSpAC/$', spac_graphs.getSpAC, name='getSpAC'),
    url(r'^myPhyloDB/statusSpAC/$', spac_graphs.statusSpAC, name='statusSpAC'),
    url(r'^myPhyloDB/removeRIDSpAC/$', spac_graphs.removeRIDSpAC, name='removeRIDSpAC'),
    url(r'^myPhyloDB/removeSpACFiles/$', spac_graphs.removeSpACFiles, name='removeSpACFiles'),
    url(r'^myPhyloDB/getTabSpAC/$', spac_graphs.getTabSpAC, name='getTabSpAC'),

    url(r'^myPhyloDB/getsoil_index/$', soil_index_graphs.getsoil_index, name='getsoil_index'),
    url(r'^myPhyloDB/statussoil_index/$', soil_index_graphs.statussoil_index, name='statussoil_index'),
    url(r'^myPhyloDB/removeRIDsoil_index/$', soil_index_graphs.removeRIDsoil_index, name='removeRIDsoil_index'),
    url(r'^myPhyloDB/removesoil_indexFiles/$', soil_index_graphs.removesoil_indexFiles, name='removesoil_indexFiles'),
    url(r'^myPhyloDB/getTabsoil_index/$', soil_index_graphs.getTabsoil_index, name='getTabsoil_index'),

    url(r'^myPhyloDB/getPCoA/$', pcoa_graphs.getPCoA, name='getPCoA'),
    url(r'^myPhyloDB/statusPCoA/$', pcoa_graphs.statusPCoA, name='statusPCoA'),
    url(r'^myPhyloDB/removeRIDPCoA/$', pcoa_graphs.removeRIDPCoA, name='removeRIDPCoA'),
    url(r'^myPhyloDB/removePCoAFiles/$', pcoa_graphs.removePCoAFiles, name='removePCoAFiles'),
    url(r'^myPhyloDB/getTabPCoA/$', pcoa_graphs.getTabPCoA, name='getTabPCoA'),

    url(r'^myPhyloDB/getWGCNA/$', wgcna_graphs.getWGCNA, name='getWGCNA'),
    url(r'^myPhyloDB/statusWGCNA/$', wgcna_graphs.statusWGCNA, name='statusWGCNA'),
    url(r'^myPhyloDB/removeRIDWGCNA/$', wgcna_graphs.removeRIDWGCNA, name='removeRIDWGCNA'),
    url(r'^myPhyloDB/removeWGCNAFiles/$', wgcna_graphs.removeWGCNAFiles, name='removeWGCNAFiles'),
    url(r'^myPhyloDB/getTabWGCNA/$', wgcna_graphs.getTabWGCNA, name='getTabWGCNA'),

    url(r'^myPhyloDB/getPCA/$', pca_graphs.getPCA, name='getPCA'),
    url(r'^myPhyloDB/statusPCA/$', pca_graphs.statusPCA, name='statusPCA'),
    url(r'^myPhyloDB/removeRIDPCA/$', pca_graphs.removeRIDPCA, name='removeRIDPCA'),
    url(r'^myPhyloDB/removePCAFiles/$', pca_graphs.removePCAFiles, name='removePCAFiles'),
    url(r'^myPhyloDB/getTabPCA/$', pca_graphs.getTabPCA, name='getTabPCA'),

    url(r'^myPhyloDB/getSPLS/$', spls_graphs.getSPLS, name='getSPLS'),
    url(r'^myPhyloDB/statusSPLS/$', spls_graphs.statusSPLS, name='statusSPLS'),
    url(r'^myPhyloDB/removeRIDSPLS/$', spls_graphs.removeRIDSPLS, name='removeRIDSPLS'),
    url(r'^myPhyloDB/removeSPLSFiles/$', spls_graphs.removeSPLSFiles, name='removeSPLSFiles'),
    url(r'^myPhyloDB/getTabSPLS/$', spls_graphs.getTabSPLS, name='getTabSPLS'),

    url(r'^myPhyloDB/updateDiffAbund/$', diffabund_graphs.updateDiffAbund, name='updateDiffAbund'),
    url(r'^myPhyloDB/getDiffAbund/$', diffabund_graphs.getDiffAbund, name='getDiffAbund'),
    url(r'^myPhyloDB/removeRIDDIFF/$', diffabund_graphs.removeRIDDIFF, name='removeRIDDIFF'),
    url(r'^myPhyloDB/getTabDiffAbund/$', diffabund_graphs.getTabDiffAbund, name='getTabDiffAbund'),
    url(r'^myPhyloDB/removeDiffAbundFiles/$', diffabund_graphs.removeDiffAbundFiles, name='removeDiffAbundFiles'),

    url(r'^myPhyloDB/updateGAGE/$', gage_graphs.updateGAGE, name='updateGAGE'),
    url(r'^myPhyloDB/getGAGE/$', gage_graphs.getGAGE, name='getGAGE'),
    url(r'^myPhyloDB/removeRIDGAGE/$', gage_graphs.removeRIDGAGE, name='removeRIDGAGE'),
    url(r'^myPhyloDB/getTabGAGE/$', gage_graphs.getTabGAGE, name='getTabGAGE'),
    url(r'^myPhyloDB/removeGAGEFiles/$', gage_graphs.removeGAGEFiles, name='removeGAGEFiles'),

    url(r'^myPhyloDB/statusPyBake/$', pybake.statusPyBake, name='statusPyBake'),

    url(r'^myPhyloDB/taxaJSON/$', views.taxaJSON, name='taxaJSON'),
    url(r'^myPhyloDB/nzJSON/$', views.nzJSON, name='nzJSON'),
    url(r'^myPhyloDB/pathJSON/$', views.pathJSON, name='pathJSON'),

    # For this url, preference would be for /myPhyloDB/select/
    # However, page names must be unique to render correctly
    url(r'^myPhyloDB/select_upload/$', views.uploadNorm, name='uploadNorm'),

    # queue caller
    url(r'^myPhyloDB/funcCall/$', queue.funcCall, name='funcCall'),
    url(r'^myPhyloDB/stop/$', queue.stop, name='stop'),

    # data queue
    url(r'^myPhyloDB/datfuncCall/$', dataqueue.datfuncCall, name='datfuncCall'),
    url(r'^myPhyloDB/datstop/$', dataqueue.datstop, name='datstop'),

    # usrFiles
    url(r'^myPhyloDB/usrFiles/$', views.usrFiles, name='usrFiles'),

    # Json data for select page
    url(r'^myPhyloDB/projectTableJSON/$', views.projectTableJSON, name='projectTableJSON'),
    url(r'^myPhyloDB/sampleTableJSON/$', views.sampleTableJSON, name='sampleTableJSON'),
    url(r'^myPhyloDB/referenceTableJSON/$', views.referenceTableJSON, name='referenceTableJSON'),
    url(r'^myPhyloDB/airTableJSON/$', views.airTableJSON, name='airTableJSON'),
    url(r'^myPhyloDB/associatedTableJSON/$', views.associatedTableJSON, name='associatedTableJSON'),
    url(r'^myPhyloDB/microbialTableJSON/$', views.microbialTableJSON, name='microbialTableJSON'),
    url(r'^myPhyloDB/soilTableJSON/$', views.soilTableJSON, name='soilTableJSON'),
    url(r'^myPhyloDB/waterTableJSON/$', views.waterTableJSON, name='waterTableJSON'),
    url(r'^myPhyloDB/userTableJSON/$', views.userTableJSON, name='userTableJSON'),

]

