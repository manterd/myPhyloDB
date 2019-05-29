from django.contrib.auth.models import User
from database.models import UserProfile, Project, PublicProjects
from django.http import HttpResponse
import json
import functions
import numpy as np

"""
File intended to contain all functions related to permissions alteration, for organization and ease of design

Current functions:
    pubToPriv: Given a project id (p_uuid), adjust permissions related to the project with matching id as a public
                project being switched to a private project. This involves removal of p_uuid from public list, and
                searching through the owner's permissions list on top of the project's to populate the private lists

    privToPub: Given a project id (p_uuid), adjust permissions related to the project with matching id as a private
                project being switched to a public project. This involves removal of the project from zero to many users
                private permissions lists, and adding the project to the global public project list
"""


# Add p_uuid to public project list, assume it does not exist already
def newPub(p_uuid):
    # get the single instance (at least there should only be one) of the public projects list
    publicList = PublicProjects.objects.all().first()
    # add this p_uuid to the list, using method implemented in the model
    publicList.update_list(add=p_uuid)
    publicList.save()
    # done
    return


# Add p_uuid to relevant private lists, assume it did not exist before
def newPriv(p_uuid):
    # get actual project object
    myProj = Project.objects.get(projectid=p_uuid)
    # get owner of project
    myOwner = myProj.owner
    # get profile of owner
    myProfile = UserProfile.objects.get(user=myOwner)
    # get owner's list of shared user ids
    myGivenPerms = myProfile.gavePermsTo.split(";")
    # get list of users that corresponds to ids
    myGivenPermsUsers = User.objects.filter(id__in=myGivenPerms)
    # get profiles of shared users
    mySharedUsers = UserProfile.objects.filter(user__in=myGivenPermsUsers)
    # add this project to shared users lists
    for share in mySharedUsers:
        share.update_list(add=p_uuid)
        share.save()
    # done
    return


# Remove p_uuid from public
def remPub(p_uuid):
    # get the single instance (at least there should only be one) of the public projects list
    publicList = PublicProjects.objects.all().first()
    # remove p_uuid from list
    publicList.update_list(remove=p_uuid)
    publicList.save()
    # done
    return


# Remove p_uuid from private lists, remember to check the project's whitelist as well as the owner's
def remPriv(p_uuid):
    # get actual project object
    myProj = Project.objects.get(projectid=p_uuid)
    # get owner of project
    myOwner = myProj.owner
    # get profile of owner
    myProfile = UserProfile.objects.get(user=myOwner)
    # get owner's list of shared user ids
    myGivenPerms = myProfile.gavePermsTo.split(";")
    # get actual user objects from ids
    myGivenPermsUsers = User.objects.filter(id__in=myGivenPerms)
    # get profiles of shared users
    mySharedUsers = UserProfile.objects.filter(user__in=myGivenPermsUsers)
    # add this project to shared users lists
    for share in mySharedUsers:
        share.update_list(remove=p_uuid)
        share.save()
    # check the whitelists for this project, cleanup whitelisted users
    # we only need the viewlist since anyone added to the edit list is automatically added to viewlist
    myViewList = myProj.whitelist_view.split(";")
    myViewUsers = User.objects.filter(id__in=myViewList)
    whitelistedUsers = UserProfile.objects.filter(user__in=myViewUsers)
    for whitelisted in whitelistedUsers:
        whitelisted.update_list(remove=p_uuid)
        whitelisted.save()
    # wipe out the whitelists
    myProj.whitelist_view = ""
    myProj.whitelist_edit = ""
    myProj.save()
    # done
    return


# We cannot make any assumptions as to the permanence of any given private permission, as the owner's permissions list
# can change between permissions changes for this project. As such, we must use only the data we know is presently valid
# which means that the owner of this project now wants a specific permissions mode, and AT PRESENT has certain users
# included in their whitelist

def pubToPriv(p_uuid):
    remPub(p_uuid)
    newPriv(p_uuid)
    return


def privToPub(p_uuid):
    remPriv(p_uuid)
    newPub(p_uuid)
    return


# Entering the realm of older code: beware of nonsense
def updateAccPerms(request):
    # gets called from ajax, given list of names to remove or add
    functions.log(request, "FUNCTION", "ACCPERMS")
    if request.is_ajax():
        allJson = json.loads(request.GET["all"])
        # get selected names from list, to be removed
        remList = allJson['keys']
        # get list of names to add
        nameList = allJson['names']
        nameList = nameList.split(';')
        # get the user performing changes
        thisUser = User.objects.get(username=request.user.username)
        # give permission from thisUser to users in nameList
        giveFilePerms(thisUser, nameList)
        # revoke permission from thisUser to users in remList
        removeFilePerms(thisUser, remList)

    else:
        print "Request was not ajax!"
        functions.log(request, "NON_AJAX", "ACCPERMS")

    retDict = {"error": "none"}
    res = json.dumps(retDict)
    return HttpResponse(res, content_type='application/json')


# called by updateAccPerms, but can stand on its own
# given a user (owner) and a list of usernames to give permissions to (addList), update info for all users involved
# such that all of owner's projects are visible to all users in addList (and previous list members)
def giveFilePerms(owner, addList):
    # given a user and a list of usernames, add username list to user's whitelist
    ownProfile = UserProfile.objects.get(user=owner)  # owner is a user object, ie request.user
    currentList = ownProfile.gavePermsTo.split(';')
    addUsers = []
    for addName in addList:
        try:
            curUser = User.objects.get(username=addName)
            if addName != "":
                addUsers.append(curUser)
                ownProfile.update_gavePerms(add=addName)
        except:
            pass
    # save changes
    ownProfile.save()
    # with gavePermsTo updated, we now update hasPermsFrom and privateProjectList fields
    # get all projects from the owner for later use
    ownerProjects = Project.objects.filter(owner=owner)
    ownerProjIDs = []
    for proj in ownerProjects:
        ownerProjIDs.append(proj.projectid)
    # for each user in the list of users, get their profile and update the lists
    for prof in addUsers:
        if prof.username in currentList:
            # update hasPermsFrom list
            # check if this user has owner on their list, if not add
            thisProf = UserProfile.objects.get(user=prof)
            thisPermsList = thisProf.hasPermsFrom.split(';')
            if owner.username not in thisPermsList:
                thisProf.update_hasPerms(add=owner.username)
            # now update privateProjectList
            for projID in ownerProjIDs:
                thisProf.update_list(add=projID)
            # save changes to this user's profile
            thisProf.save()
    return


def removeFilePerms(owner, remList):
    # use .lower on everything, remove users in remList from owner's whitelist (if they are there already)
    # also, go to each user on remList and remove owner from their "added by" list (again, check if its there already)
    ownProfile = UserProfile.objects.get(user=owner)  # owner is a user object, ie request.user
    currentList = ownProfile.gavePermsTo.split(';')
    remUsers = []
    for curName in currentList:
        if curName != "":
            try:
                checkUser = User.objects.get(username=curName)
                if curName in remList:
                    # user is to be removed, save for profile check later
                    remUsers.append(checkUser)
                    # actually remove from list
                    ownProfile.update_gavePerms(remove=curName)
            except:
                pass
    # save changes
    ownProfile.save()
    ownerProjects = Project.objects.filter(owner=owner)
    # WE HAVE THE LIST, just update those users
    # Also, bear in mind that a project might have whitelisted a user, so check before removing from private list
    for curUser in remUsers:
        if curUser.username not in currentList:
            curProf = UserProfile.objects.get(user=curUser)
            curHasPerms = curProf.hasPermsFrom.split(';')
            if owner.username in curHasPerms:
                # remove it!, save it again after
                curProf.update_hasPerms(remove=owner.username)
                # update privateProjectList
                # check if project has this user whitelisted before removing
                for proj in ownerProjects:
                    if curUser.username not in proj.whitelist_view.split(";"):
                        curProf.update_list(rem=proj.id)
                # save changes
                curProf.save()
    return


def updateProjPerms(request):  # this is the project whitelisting section, a tree of projects and a list of names were sent
    functions.log(request, "FUNCTION", "PROJPERMS")
    if request.is_ajax():
        allJson = json.loads(request.GET["all"])
        # get selected projects list (files? not samples as subsets though)
        selList = allJson['keys']
        # get list of names to remove, in the form of a list of username-projectid pairs
        remKeys = allJson['remKeys']
        # for now, assuming view perms are being revoked, as such, edit perms will be wiped if present
        # get list of names to add
        nameList = allJson['names']
        nameList = nameList.split(';')
        # get permission mode (view vs edit)  * view is redundant if project is public
        permLevel = allJson['mode']
        # permLevel 0 means view only
        # permLevel 1 means editing as well

        thisUser = User.objects.get(username=request.user.username)  # username uniqueness is case incensitive
        # (as in no accounts named ADMIN since admin exists, etc) So we can run everything in .lower for equality checks

        # loop through nameList, verify each exists
        errorList = []
        finalNameList = nameList[:]
        for name in nameList:
            if name != "":
                if not User.objects.filter(username=name).exists():
                    errorList.append(name)
                    finalNameList.remove(name)
                    # remove name from nameList if user does not exist

        # get project by id from selList
        for pid in selList:
            curProj = Project.objects.get(projectid=pid)
            # check if thisUser is owner of project (or superuser)
            if thisUser == curProj.owner or thisUser.is_superuser:
                # loop through nameList, check if name is present on permList already
                for name in finalNameList:
                    curProj.update_viewPerms(add=name)
                    # if being given edit perms, update that list as well
                    if permLevel:
                        curProj.update_editPerms(add=name)
                    # we're modifying project visibility either way, so add this project to each user's private list
                    # already ran a username exists check on the whole list, so can assume the objects exist
                    curUser = UserProfile.objects.get(user=User.objects.get(username=name))
                    curUser.update_list(add=pid)
                    curUser.save()
            else:
                # user does not have permission to make these permissions changes
                print(thisUser.username, "has invalid permissions for permissions change on project owned by", curProj.owner.username)
                functions.log(request, "PERM_ERROR", "PROJPERMS")
            # save changes
            curProj.save()
        # return new profile page (boxes cleared, perhaps success alert message)
        remDict = {}
        for pair in remKeys:
            mySplit = pair.split(";")
            if mySplit[1] not in remDict.keys():
                remDict[mySplit[1]] = []
            remDict[mySplit[1]].append(str(mySplit[0]))
        # now we have a dictionary with key projectID and values are a list of names to revoke perms of
        for projID in remDict.keys():
            myProj = Project.objects.get(projectid=projID)
            for name in remDict[projID]:
                myProj.update_viewPerms(remove=name)
                myProj.update_editPerms(remove=name)
                # now remove this project from the user's whitelist IFF they have not been granted account-based perms
                curUser = UserProfile.objects.get(user=User.objects.get(username=name))
                if request.user.username not in curUser.hasPermsFrom.split(";"):
                    curUser.update_list(remove=projID)
                    curUser.save()
            myProj.save()
        # check if errorList has entries
        if len(errorList) > 0:
            # set return text to be an error reporting failed names
            text = "Name(s) not found: "
            iter = 1
            for name in errorList:
                if iter == len(errorList):
                    text += str(name)
                else:
                    text += str(name) + ", "
                iter += 1
        else:
            text = 'Project permissions updated successfully'

        retDict = {"error": text}
        res = json.dumps(retDict)
        return HttpResponse(res, content_type='application/json')
    else:
        print "Request was not ajax!"
        functions.log(request, "NON_AJAX", "PROJPERMS")


def getViewProjects(request):   # use this function as often as possible for project queries, put all perms stuff here
    # permissions are both project and account based, for improved usability

    # give account perms to co-workers and common collaborators
    # give project based permissions to uncommon joint project groups

    projects = Project.objects.none()
    if request.user.is_superuser:
        projects = Project.objects.order_by('project_name')
        # check if project is WIP, flag for dynatree highlighting?

    elif request.user.is_authenticated():
        # run through list of projects, when valid project is found, append filterIDS with ID
        # projects will be a queryset set to all projects, then filtered by ids in filterIDS
        # queryStartTime = time.time()

        publicIDs = PublicProjects.objects.all().first().List.split(",")    # got the public
        # print "publicIDs:", publicIDs
        # get private IDs here
        privateIDs = UserProfile.objects.get(user=request.user).privateProjectList.split(",")
        # print "privateIDs:", privateIDs
        filterIDS = np.unique(privateIDs+publicIDs)

        # queryTime = time.time() - queryStartTime
        # print "Select query time:", queryTime
        projects = Project.objects.filter(projectid__in=filterIDS).order_by('project_name')

    if not request.user.is_superuser and not request.user.is_authenticated():
        # impossible to have guest user be on whitelist (hopefully), so public only check
        projects = Project.objects.filter(status='public').order_by('project_name')

    return projects


def getEditProjects(request):   # TODO check all permissions required trees are verified on backend afterwards
    # TODO implement privateProjectList equivalent for editing permissions, to improve scalability
    projects = Project.objects.none()
    if request.user.is_superuser:
        projects = Project.objects.order_by('project_name')

    elif request.user.is_authenticated():
        # run through list of projects, when valid project is found, append filterIDS with ID
        # projects will be a queryset set to all projects, then filtered by ids in filterIDS
        filterIDS = []
        for proj in Project.objects.all():
            good = False  # good to add to list
            if proj.owner == request.user:
                good = True
            checkList = proj.whitelist_edit.split(';')
            for cid in checkList:
                if cid == request.user.id:
                    good = True
            if good:
                filterIDS.append(proj.projectid)
        projects = Project.objects.filter(projectid__in=filterIDS).order_by('project_name')

    return projects