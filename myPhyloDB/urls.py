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
from database import utils
from database import dataqueue

admin.autodiscover()


class MyRegistrationView(RegistrationView):
    def get_success_url(self, request=None, user=None):
        return '/myPhyloDB/select/'


urlpatterns = [
    ### administration, registration, and main myPhyloDB pages
    url(r'^myPhyloDB/admin/', admin.site.urls),
    url(r'^myPhyloDB/accounts/register/$', MyRegistrationView.as_view(), name='registration_register'),
    url(r'^myPhyloDB/accounts/', include('registration.backends.simple.urls')),
    url(r'^myPhyloDB/', include('database.urls')),

    ### urls from views page
    url(r'^myPhyloDB/saveSampleList/$', views.saveSampleList, name='saveSampleList'),
    url(r'^myPhyloDB/taxaJSON/$', views.taxaJSON, name='taxaJSON'),
    url(r'^myPhyloDB/nzJSON/$', views.nzJSON, name='nzJSON'),
    url(r'^myPhyloDB/pathJSON/$', views.pathJSON, name='pathJSON'),
    url(r'^myPhyloDB/projectTableJSON/$', views.projectTableJSON, name='projectTableJSON'),
    url(r'^myPhyloDB/sampleTableJSON/$', views.sampleTableJSON, name='sampleTableJSON'),
    url(r'^myPhyloDB/referenceTableJSON/$', views.referenceTableJSON, name='referenceTableJSON'),
    url(r'^myPhyloDB/airTableJSON/$', views.airTableJSON, name='airTableJSON'),
    url(r'^myPhyloDB/associatedTableJSON/$', views.associatedTableJSON, name='associatedTableJSON'),
    url(r'^myPhyloDB/microbialTableJSON/$', views.microbialTableJSON, name='microbialTableJSON'),
    url(r'^myPhyloDB/soilTableJSON/$', views.soilTableJSON, name='soilTableJSON'),
    url(r'^myPhyloDB/waterTableJSON/$', views.waterTableJSON, name='waterTableJSON'),
    url(r'^myPhyloDB/userTableJSON/$', views.userTableJSON, name='userTableJSON'),
    url(r'^myPhyloDB/usrFiles/$', views.usrFiles, name='usrFiles'),

    ### urls from parser page
    url(r'^myPhyloDB/status/$', parsers.status, name='status'),
    url(r'^myPhyloDB/reanalyze/$', parsers.reanalyze, name='reanalyze'),

    ### urls from trees page
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

    ### urls from normalization page
    url(r'^myPhyloDB/getNorm/$', norm_graphs.getNorm, name='getNorm'),
    url(r'^myPhyloDB/getTab/$', norm_graphs.getTab, name='getTab'),
    url(r'^myPhyloDB/getBiom/$', norm_graphs.getBiom, name='getBiome'),

    ### urls from analysis pages
    url(r'^myPhyloDB/getCatUnivData/$', anova_graphs.getCatUnivData, name='getCatUnivData'),
    url(r'^myPhyloDB/getQuantUnivData/$', anova_graphs.getQuantUnivData, name='getQuantUnivData'),
    url(r'^myPhyloDB/getDiffAbund/$', diffabund_graphs.getDiffAbund, name='getDiffAbund'),
    url(r'^myPhyloDB/getGAGE/$', gage_graphs.getGAGE, name='getGAGE'),
    url(r'^myPhyloDB/getPCA/$', pca_graphs.getPCA, name='getPCA'),
    url(r'^myPhyloDB/getPCoA/$', pcoa_graphs.getPCoA, name='getPCoA'),
    url(r'^myPhyloDB/getsoil_index/$', soil_index_graphs.getsoil_index, name='getsoil_index'),
    url(r'^myPhyloDB/getSpAC/$', spac_graphs.getSpAC, name='getSpAC'),
    url(r'^myPhyloDB/getSPLS/$', spls_graphs.getSPLS, name='getSPLS'),
    url(r'^myPhyloDB/getWGCNA/$', wgcna_graphs.getWGCNA, name='getWGCNA'),

    ### urls from analysis queue
    url(r'^myPhyloDB/funcCall/$', queue.funcCall, name='funcCall'),
    url(r'^myPhyloDB/stop/$', queue.stop, name='stop'),
    url(r'^myPhyloDB/removeRID/$', queue.removeRID, name='removeRID'),

    ### data from upload/reprocessing queue
    url(r'^myPhyloDB/datfuncCall/$', dataqueue.datfuncCall, name='datfuncCall'),
    url(r'^myPhyloDB/datstop/$', dataqueue.datstop, name='datstop'),

    ### urls from pybake
    url(r'^myPhyloDB/statusPyBake/$', pybake.statusPyBake, name='statusPyBake'),

    ### urls from utils file
    url(r'^myPhyloDB/getRawData/$', utils.getRawData, name='getRawData'),
    url(r'^myPhyloDB/removeFiles/$', utils.removeFiles, name='removeFiles'),
]

