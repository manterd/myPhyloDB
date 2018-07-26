import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'myPhyloDB.settings'
import django
django.setup()

from database.views import getAdminLog

print "Getting log"
getAdminLog()
print "Complete, check admin_log.txt"
