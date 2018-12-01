from django.http import HttpResponse, HttpRequest
from django.shortcuts import render
import json
from Queue import Queue
from time import sleep
import datetime
import functions
import database.views
from database.forms import UploadForm1, UploadForm2, UploadForm5
from database.models import Reference


datQ = Queue(maxsize=0)
datActiveList = [0]
datQueueList = {}
datQueueFuncs = {}
datQueueTimes = {}
datQueueUsers = {}
datStopList = [0]
datRecent = {}
datQList = []
datStatDict = {}
datStopDict = {}
datStopped = 0


def datstop(request):
    global datActiveList, datQueueList, datStopList, datStopDict, datStopped, datRecent
    # datStopped is number of queued requests which have been stopped in advance. Was intended for progress tracking
    # BUT datStopped fails to account for queue removals which occur AFTER a request (ie can make it look like your
    # request is earlier in the queue because someone after you cancelled. Need to make datStopped into a dictionary
    # basically everything in this request manager needs a dictionary in order to be practical
    RID = str(request.GET['all'])
    datStopDict[RID] = RID  # this is for queue skipping
    try:
        pid = datActiveList.index(RID)  # see if RID is an active process, we'll need to update stopList if so
        datStopList[pid] = RID
        datActiveList[pid] = 0
        functions.termP()   # this ends mothur subprocess(es) directly, this assumes only one mothur can run at a time
    except Exception:   # not in activelist, either done/stopped already or its in the queue
        try:
            # check queuelist, if its there flag its RID so the main process skips it, return stop msg to user page
            # wait for status, then cleanup
            # get function from funccall, return correct page rendering
            print "Stopping queued request"
            pid = datQueueList[RID]
            thisFunc = datQueueFuncs[pid]

            # available functions: uploadFunc reanalyze updateFunc pybake
            # log the stop request
            functions.log(request, "QSTOP", thisFunc)

            print "Getting page returns"
            retStuff = None

            if thisFunc == "uploadFunc":    # just directly rendering page because request chaining is clunky
                projects = Reference.objects.none()
                if request.user.is_superuser:
                    projects = Reference.objects.all().order_by('projectid__project_name', 'path')
                elif request.user.is_authenticated():
                    projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(
                        author=request.user)
                retStuff = render(
                    request,
                    'upload.html',
                    {'projects': projects,
                     'form1': UploadForm1,
                     'form2': UploadForm2,
                     'error': "Upload stopped"
                     }
                )

            elif thisFunc == "reanalyze":
                # reprocess and pybake are fine without stoplists since they are page requests
                retStuff = database.views.reprocess(request)

            elif thisFunc == "updateFunc":
                functions.log(request, "STOP", "UPDATE")
                state = "Update stopped"
                retStuff = render(
                    request,
                    'update.html',
                    {'form5': UploadForm5,
                     'state': state}
                )

            elif thisFunc == "pybake":
                retStuff = database.views.pybake(request)

            print "Ret!"
            if retStuff is None:
                print "This is a problem. Hacked?"  # or someone added a new dataqueue function without updating here
                # TODO more modular code, make it easier to add new dataqueue functions (fewer sections to update)
            datRecent[pid] = retStuff

        except Exception as er:   # not in queuelist either. actual error
            print "Error with datStop:", er
            myDict = {'error': 'Analysis not running'}
            stop = json.dumps(myDict)
            return HttpResponse(stop, content_type='application/json')

    cleanup(RID)
    myDict = {'error': 'none', 'message': 'Your database changes have been cancelled!'}
    stop = json.dumps(myDict)
    return HttpResponse(stop, content_type='application/json')


# TODO put general function description (going over parameters and usage) at the start of each function definition
def getDataQueue(request):
    if not request.user.is_superuser or not request.user.is_authenticated:
        output = json.dumps({'display': "Invalid Permissions"})
        return HttpResponse(output, content_type='application/json')
    stringDict = {}
    queueString = ""

    # Dictionaries, all keyed by RID. We get these from the global list, not declared here since we aren't mutating
    '''
    datQueueList = {}       rid
    datQueueFuncs = {}      function called
    datQueueTimes = {}      when did they call it
    datQueueUsers = {}      who made the call
    '''

    # Get queued processes, format string for display
    #print "Queued"
    for dataReq in datQueueList:
        stringDict[str(datQueueTimes[dataReq])] = str(dataReq) + ";" + str(datQueueTimes[dataReq]) + ";" + \
                       str(datQueueFuncs[dataReq]) + ";" + str(datQueueUsers[dataReq]) + ";QUEUED\n"

    # Get active processes
    #print "Active"
    try:
        for dataProc in datActiveList:
            if len(str(dataProc)) > 1:
                stringDict[str(datQueueTimes[dataProc])] = str(dataProc) + ";" + str(datQueueTimes[dataProc]) + ";" + \
                               str(datQueueFuncs[dataProc]) + ";" + str(datQueueUsers[dataProc]) + ";ACTIVE\n"
    except Exception as ex:
        print "Problem with dataqueue console:", ex
    # sort by timestamp, reversed order so the "biggest" (newest) times are first
    # the idea being to display the queue in a sequence that "falls" (both more intuitive and follows the pattern set
    # by the rest of our output texts thus far)
    #print "Sort"
    for key in sorted(stringDict.keys(), reverse=True):
        queueString += stringDict[key]

    # debug string, could put these under a flag. Actually, TODO make global debug flag for print spam / verbose output
    #print "Strings:", stringDict, "Keys:", stringDict.keys(), "Sorted:", sorted(stringDict.keys()), "Queue:", queueString

    #print "Display"
    output = json.dumps({'display': queueString})
    return HttpResponse(output, content_type='application/json')


def decremQ():
    for thing in datStatDict:
        datStatDict[thing] -= 1


def dataprocess(pid):
    global datActiveList, datQueueList, datQueueFuncs, datRecent
    while True:
        data = datQ.get(block=True, timeout=None)   # .get will block until there's something in the queue (thread-safe)
        decremQ()
        RID = data['RID']
        if RID in datStopDict:
            datStopDict.pop(RID, 0)
        else:
            funcName = data['funcName']
            request = data['request']
            stopList = data['stop']  # not sure what this is for, this whole file has a lot of cleaning up to be done
            datQueueList.pop(RID, 0)
            datActiveList[pid] = RID
            functions.log(request, "QSTART", funcName)
            if datActiveList[pid] == RID:
                if funcName == "uploadFunc":
                    datRecent[RID] = database.views.uploadFunc(request, datStopList)
                if funcName == "fileUpFunc":
                    datRecent[RID] = database.views.fileUpFunc(request, datStopList)
                if funcName == "reanalyze":
                    resp = functions.reanalyze(request, datStopList)
                    if resp is None:
                        resp = database.views.reprocess(request)
                    if resp == "Stopped":
                        # why does this get a normal call first?
                        resp = database.views.reprocess(request)
                        resp['error'] = "Reprocessing stopped"
                    datRecent[RID] = resp
                if funcName == "updateFunc":
                    datRecent[RID] = database.views.updateFunc(request, datStopList)
                if funcName == "geneParse":
                    #print 'starting'
                    datRecent[RID] = functions.geneParse(request)
                if funcName == "koParse":
                    datRecent[RID] = functions.koParse(request)
                if funcName == "nzParse":
                    datRecent[RID] = functions.nzParse(request)
            datActiveList[pid] = ''
            functions.log(request, "QFINISH", funcName)
        cleanup(RID)


def cleanup(RID):
    global datQueueFuncs, datQueueTimes, datQueueUsers, datQueueList
    # could potentially have some gone before others, for now just putting all in try's
    try:
        datQueueList.pop(RID, 0)
    except:
        pass
    try:
        datQueueFuncs.pop(RID, 0)
    except:
        pass
    try:
        datQueueTimes.pop(RID, 0)
    except:
        pass
    try:
        datQueueUsers.pop(RID, 0)
    except:
        pass


def datfuncCall(request):
    global datActiveList, datQueueList, datQueueFuncs, datStopList, datStopDict, datStatDict, datQList, datQueueTimes, datQueueUsers
    RID = request.POST['RID']
    funcName = request.POST['funcName']
    datQueueList[RID] = RID  # add to queuelist, remove when processed or stopped
    datQueueFuncs[RID] = funcName
    datQueueTimes[RID] = str(datetime.datetime.now())
    datQueueUsers[RID] = request.user.username
    # add a timestamp dict as well, to sort the others for logging purposes
    qDict = {'RID': RID, 'funcName': funcName, 'request': request, 'stop': datStopList}
    datQList.append(qDict)
    datQ.put(qDict, True)
    datStatDict[RID] = int(datQ.qsize())

    # print log info, need to write this to a file somewhere
    functions.log(request, "QADD", funcName)
    while True:
        try:
            results = datRecent[RID]
            datRecent.pop(RID, 0)
            datStatDict.pop(RID, 0)
            # this return resets data like mothurStatSave
            return results
        except KeyError:
            if RID in datStopList:
                datStatDict.pop(RID, 0)
        except Exception as e:
            print "Unexpected exception: "+str(e)
        sleep(1)


def datstat(RID):
    try:
        datStopDict[RID] = datStopDict[RID]
        if RID in datQueueList:
            return -512
        else:
            return -1024
    except Exception:
        # if not in stopDict (ie not called to stop) queuepos is current (minus number of pre stops)
        return datStatDict[RID] - datStopped


def getStops():
    return datStopList

