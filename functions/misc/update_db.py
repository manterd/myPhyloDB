# Do not change order of imports...
import os
from django.core.wsgi import get_wsgi_application
os.environ["DJANGO_SETTINGS_MODULE"] = "myPhyloDB.settings"
application = get_wsgi_application()

# Do not change order of imports...
from django.core.management import call_command
from django.db.models import Sum
from database.models import Sample, Profile


def update_reads():
    sampleList = Sample.objects.values_list('sampleid', flat=True)
    for ID in sampleList:
        sample = Sample.objects.get(sampleid=ID)
        reads = Profile.objects.filter(sampleid=ID).aggregate(count=Sum('count'))
        print ID, 'has', reads['count'], 'reads'
        sample.reads = reads['count']
        sample.save()

if __name__ == "__main__":
    call_command('migrate')
    call_command('migrate', '--database=picrust')
    call_command('migrate', '--run-syncdb')
    update_reads()

