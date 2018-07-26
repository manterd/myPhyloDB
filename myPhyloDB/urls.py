from django.conf.urls import include, url
from django.contrib import admin

from database import views

admin.site.site_url = '/myPhyloDB/'
admin.autodiscover()

urlpatterns = [
    # allauth Urls
    url(r'^myPhyloDB/accounts/', include('allauth.urls')),

    # Admin Urls
    url(r'^myPhyloDB/admin/', admin.site.urls),

    url(r'^myPhyloDB/', include('database.urls')),
    url(r'^myPhyloDB/', include('functions.analysis.urls')),
    url(r'^myPhyloDB/', include('functions.queues.urls')),
    url(r'^myPhyloDB/', include('functions.utils.urls')),

    ### urls from views page
    url(r'^myPhyloDB/saveSampleList/$', views.saveSampleList, name='saveSampleList'),
    url(r'^myPhyloDB/taxaJSON/$', views.taxaJSON, name='taxaJSON'),
    url(r'^myPhyloDB/nzJSON/$', views.nzJSON, name='nzJSON'),
    url(r'^myPhyloDB/nzTaxaJSON/$', views.nzTaxaJSON, name='nzTaxaJSON'),
    url(r'^myPhyloDB/pathJSON/$', views.pathJSON, name='pathJSON'),
    url(r'^myPhyloDB/pathTaxaJSON/$', views.pathTaxaJSON, name='pathTaxaJSON'),
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
    url(r'^myPhyloDB/getProjectFiles/$', views.getProjectFiles, name='getProjectFiles'),
    url(r'^myPhyloDB/remProjectFiles/$', views.remProjectFiles, name='remProjectFiles'),
    url(r'^myPhyloDB/removeUploads/$', views.removeUploads, name='removeUploads'),
    url(r'^myPhyloDB/addPerms/$', views.addPerms, name='addPerms'),
    url(r'^myPhyloDB/remPerms/$', views.remPerms, name='remPerms'),
    url(r'^myPhyloDB/updateFilePerms/$', views.updateFilePerms, name='updateFilePerms'),


    url(r'^myPhyloDB/accounts/profile/$', views.profile, name='profile'),
    url(r'^myPhyloDB/accounts/changeuser/$', views.changeuser, name='changeuser'),
    url(r'^myPhyloDB/accounts/updateInfo/$', views.updateInfo, name='updateInfo'),
]
