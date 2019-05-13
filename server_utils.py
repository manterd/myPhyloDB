import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'myPhyloDB.settings'
import django
django.setup()

# this py script is for running functions outside the normal code which require django to be active


def createPublicList():
    from database.models import PublicProjects
    PublicProjects.objects.all().delete()
    theList = PublicProjects.objects.create()
    print "Object made"
    theList.update_list()   # no args to make the list from scratch
    theList.save()
    print "List made"


def createPrivateLists():
    from database.models import UserProfile
    allUserProfs = UserProfile.objects.all()
    print "Got profs"
    for prof in allUserProfs:
        prof.update_list()  # creating initial lists for all users presently in database
        prof.save()
        # need to run this when a new user is created???
    print "Made all lists"


def adminLog():
    from database.views import getAdminLog

    print "Getting log"
    getAdminLog()
    print "Complete, check admin_log.txt"



def fixISV():
    from database.models import OTU_99
    # need to go through existing database for OTU's with isv__ names, ensure sequential numbers and no duplicates
    # get all the current ISV's
    oldISV = OTU_99.objects.filter(otuName__startswith='isv_')
    countISV = oldISV.count()
    print("Cleaning up", countISV, "ISVs")
    curISV = 0
    for isv in oldISV:
        myName = isv.otuName
        myISV = int(myName.split("_")[1])
        if myISV != curISV:
            # we skipped! OR we ran through a duplicate
            isv.otuName = "isv_"+str(curISV)    # does it matter if we just write over it? is this referenced somewhere?
        curISV += 1
    pass