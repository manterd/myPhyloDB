from django.http import HttpResponse, HttpRequest
from django.shortcuts import render
import json
from Queue import Queue
from time import sleep
import datetime
import functions
from functions.utils.debug import debug
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


def datstop(request):
    global datActiveList, datQueueList, datStopList, datStopDict, datRecent
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
            pid = datQueueList[RID]
            thisFunc = datQueueFuncs[pid]

            # available functions: processFunc reanalyze updateFunc pybake
            # log the stop request
            functions.log(request, "QSTOP", thisFunc)

            retStuff = None

            if thisFunc == "processFunc":
                retStuff = database.views.process(request, errorText="Processing stopped")

            elif thisFunc == "reanalyze":
                # do these pages support errorText flag? TODO 1.3 they need to
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


            if retStuff is None:
                print "This is a problem. Hacked?"  # or someone added a new dataqueue function without updating here
                # TODO 1.4 more modular code, make it easier to add new dataqueue functions (fewer sections to update)
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


# TODO 1.4 put general function description (going over parameters and usage) at the start of each function definition
def getDataQueue(request):  # TODO 1.4 console display items getting stuck when errored (not urgent, self fixes with new queued item)
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

    #debug("Strings:", stringDict, "Keys:", stringDict.keys(), "Sorted:", sorted(stringDict.keys()), "Queue:", queueString)

    #print "Display"
    output = json.dumps({'display': queueString})
    return HttpResponse(output, content_type='application/json')


def decremQ(start=0):  # Reduce the waitlist stat for each index after and including given index
    # this is used when a request is cancelled before entering processing, since we will skip it once its turn arrives
    curIndex = 0
    for thing in datStatDict:
        if curIndex >= start:
            datStatDict[thing] -= 1
        curIndex += 1


def dataprocess(pid):   # TODO 1.3 uploads end with "{"error": "Exception: invalid literal for int() with base 10: ''"}"
    # error page (instead of process page), even when successful (not biom upload, fastq upload does this only?)
    global datActiveList, datQueueList, datQueueFuncs, datRecent
    RID = "NULL_RID"  # only if exception occurs before first process
    request = None
    while True:
        try:
            data = datQ.get(block=True, timeout=None)   # .get will block until there's something in the queue (thread-safe)
            decremQ()
            RID = data['RID']
            if RID in datStopDict:
                datStopDict.pop(RID, 0)
            else:
                funcName = data['funcName']
                request = data['request']
                # TODO 1.4 cleanup
                stopList = data['stop']  # not sure what this is for, this whole file has a lot of cleaning up to be done
                datQueueList.pop(RID, 0)
                datActiveList[pid] = RID    # TODO 1.4 cleanup datActiveList \/
                functions.log(request, "QSTART", funcName)
                if datActiveList[pid] == RID:   # TODO 1.4 why do we assign and immediately check?
                    debug("About to run", funcName, "for", request.user.username)
                    if funcName == "processFunc":
                        datRecent[RID] = database.views.processFunc(request, datStopList)
                    elif funcName == "reanalyze":
                        resp = functions.reanalyze(request, datStopList)    # run reanalyze, check results for errors
                        if resp is None:    # if no errors, return a standard reprocess page
                            resp = database.views.reprocess(request)
                        if resp == "Stopped":   # if an error occurred, return a reprocess page with an error message
                            resp = database.views.reprocess(request)
                            resp['error'] = "Reprocessing stopped"
                        datRecent[RID] = resp
                    elif funcName == "updateFunc":
                        datRecent[RID] = database.views.updateFunc(request, datStopList)
                    elif funcName == "geneParse":
                        #print 'starting'
                        datRecent[RID] = functions.geneParse(request)
                    elif funcName == "koParse":
                        datRecent[RID] = functions.koParse(request)
                    elif funcName == "nzParse":
                        datRecent[RID] = functions.nzParse(request)

                    else:
                        # received a function name that we don't support
                        # either a developer typo exists (here or on a web page) OR someone is trying to break things
                        # either way, this is a notable event, so we print to console and log
                        print "Security check:", request.user.username, "attempting to call function", str(funcName), "in DQ"
                        functions.log(request, "INVALID_FUNCTION_NAME_DQ", str(funcName))
                        # should probably return something to the user page so they don't get stuck
                        # return is awkward here since the page is implied by funcName, which is currently invalid
                        # so we'll just send 'em to the home page
                        datRecent[RID] = database.views.home(request, errorText="Invalid function name")
                datActiveList[pid] = ''
                functions.log(request, "QFINISH", funcName)
            cleanup(RID)
        except Exception as e:
            print "Error during data queue:", e
            if request is not None:
                functions.log(request, "ERROR_DQ", str(e)+"\n")
            pid = datQueueList[RID]
            thisFunc = datQueueFuncs[pid]
            if thisFunc == "processFunc":
                retStuff = database.views.process(request, errorText=str(e))

            elif thisFunc == "reanalyze":
                # do these pages support errorText flag? TODO 1.3 they need to
                retStuff = database.views.reprocess(request)

            elif thisFunc == "updateFunc":
                functions.log(request, "ERROR", "UPDATE")
                state = str(e)
                retStuff = render(
                    request,
                    'update.html',
                    {'form5': UploadForm5,
                     'state': state}
                )

            elif thisFunc == "pybake":
                retStuff = database.views.pybake(request)

            if retStuff is None:
                print "DataQueue error has no valid func?"  # "unreachable state" haha
            datRecent[pid] = retStuff


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
    debug("DatFuncCall called by", request.user.username)
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
    debug("Added statDict for RID", RID)

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
        try:
            return datStatDict[RID]  # TODO 1.4 some way to track preemptive stops, decremQ with arg?
        except Exception as ex:
            print "Exception with datstat:", ex


def getStops():
    return datStopList
