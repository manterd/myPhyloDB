import os

os.environ['DJANGO_SETTINGS_MODULE'] = 'myPhyloDB.settings'

import django
django.setup()

from database.views import populateKoOtuList

print "Running popKegg"

populateKoOtuList()
