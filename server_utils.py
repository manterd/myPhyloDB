import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'myPhyloDB.settings'
import django
django.setup()
import sys

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
        oldName = isv.otuName
        isv.otuName = "isv_"+str(curISV)    # strictly used as an ID with no references in the database, so overwrite
        #print(oldName, "-->", isv.otuName)
        curISV += 1
        isv.save()
    print("Overwrote old ISV names, saving")


def fixBlankTaxaNames():
    # new problem: this function can wind up creating duplicate name sets
    # (multiple unclassified entries with same parents)
    # need to merge duplicates into the same entry (only for unclassified in this context, ignore sequence)
    from database.models import Kingdom, Phyla, Class, Order, Family, Genus, Species
    # for all levels of taxa that aren't OTU, find blanks and replace with 'unclassified'
    print "Starting blank taxa fix..."
    # kingdom
    blanks = Kingdom.objects.filter(kingdomName="")
    for kingdom in blanks:
        kingdom.kingdomName = "unclassified"
        kingdom.save()
    print "Kingdoms complete"

    # phyla
    blanks = Phyla.objects.filter(phylaName="")
    for phyla in blanks:
        phyla.phylaName = "unclassified"
        phyla.save()
    print "Phyla complete"

    # class
    blanks = Class.objects.filter(className="")
    for clas in blanks:
        clas.className = "unclassified"
        clas.save()
    print "Class complete"

    # order
    blanks = Order.objects.filter(orderName="")
    for order in blanks:
        order.orderName = "unclassified"
        order.save()
    print "Order complete"

    # family
    blanks = Family.objects.filter(familyName="")
    for family in blanks:
        family.familyName = "unclassified"
        family.save()
    print "Family complete"

    # genus
    blanks = Genus.objects.filter(genusName="")
    for genus in blanks:
        genus.genusName = "unclassified"
        genus.save()
    print "Genus complete"

    # species
    blanks = Species.objects.filter(speciesName="")
    for species in blanks:
        # this needs to check for existing versions of this, profile use, etc
        species.speciesName = "unclassified"
        species.save()
    print "Species complete"

    # could do a step for otu, bear in mind that sequences are potentially available here

    print "Blank taxa names fixed"


# this is main, check arg[1] for command and run it
if len(sys.argv) != 2:
    print("Incorrect number of arguments")
else:
    myCommand = sys.argv[1]
    if myCommand == "adminLog":
        adminLog()
    elif myCommand == "fixISV":
        fixISV()
    elif myCommand == "createPublicList":
        createPublicList()
    elif myCommand == "createPrivateLists":
        createPrivateLists()
    elif myCommand == "fixBlankTaxaNames":
        print "Function currently disabled" #fixBlankTaxaNames()
    else:
        print("Incorrect command: Options are adminLog fixISV fixBlankTaxaNames createPublicList createPrivateLists")
