from django.conf.urls import patterns, include
from django.conf import settings
from django.conf.urls.static import static


urlpatterns = patterns('',
    (r'^', include('phyloDB.urls')),
) + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)