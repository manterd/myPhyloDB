from django.http import HttpResponse
from Queue import Queue
import json
import multiprocessing as mp
from time import sleep, time
import gc
import datetime
import config.local_cfg
from django.contrib.auth.models import User
from database.models import UserProfile
import os
import shutil

import functions
from functions.analysis import analysis
from functions.utils.debug import debug


def analysisThreads():
    try:
        usr_threads = int(config.local_cfg.usr_num_threads)
    except Exception:
        usr_threads = 1

    if usr_threads <= 0:
        num_threads = 1
    elif usr_threads > mp.cpu_count():
        num_threads = int(mp.cpu_count())
    else:
        num_threads = usr_threads
    return int(num_threads)


q = Queue(maxsize=0)
activeList = [0] * analysisThreads()
stopList = [0] * analysisThreads()

queueList = {}
queueFuncs = {}
queueTimes = {}
queueUsers = {}

recent = {}
qList = []
statDict = {}
stopDict = {}
stopped = 0

cleanupQueue = Queue(maxsize=0)

# merged status:
base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}
complete = {}


def setBase(RID, val):
    global base
    base[RID] = val


def getBase(RID):
    return base[RID]


def getAnalysisQueue(request):  # main queue equivalent of working dataqueue tracker
    if not request.user.is_superuser or not request.user.is_authenticated:
        output = json.dumps({'display': "Invalid Permissions"})
        return HttpResponse(output, content_type='application/json')
    stringDict = {}
    queueString = ""

    # Dictionaries, all keyed by RID. We get these from the global list, not declared here since we aren't mutating
    '''
    queueList = {}       rid
    queueFuncs = {}      function called
    queueTimes = {}      when did they call it
    queueUsers = {}      who made the call
    '''

    # Get queued processes, format string for display
    #print "Queued"
    for dataReq in queueList:
        print dataReq
        stringDict[str(queueTimes[dataReq])] = str(dataReq) + ";" + str(queueTimes[dataReq]) + ";" + \
                       str(queueFuncs[dataReq]) + ";" + str(queueUsers[dataReq]) + ";QUEUED\n"

    # Get active processes
    #print "Active"
    try:
        for dataProc in activeList:
            if len(str(dataProc)) > 1:
                stringDict[str(queueTimes[dataProc])] = str(dataProc) + ";" + str(queueTimes[dataProc]) + ";" + \
                               str(queueFuncs[dataProc]) + ";" + str(queueUsers[dataProc]) + ";ACTIVE\n"
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
    queueDict = {'display': queueString}
    output = json.dumps(queueDict)
    return HttpResponse(output, content_type='application/json')


def getAnalysisHistory(request):
    if not request.user.is_authenticated:
        output = json.dumps({'display': "Invalid Permissions"})
        return HttpResponse(output, content_type='application/json')
    stringDict = {}
    queueString = ""
    recentList = UserProfile.objects.get(user=request.user).recentRIDs.split(";")
    for rid in recentList:
        myRid = str(rid)
        try:
            if myRid != "":
                stringDict[str(queueTimes[myRid])] = str(myRid) + ";" + str(queueTimes[myRid]) + ";" + \
                           str(queueFuncs[myRid]) + ";" + str(queueUsers[myRid]) + "\n"
        except:
            pass  # probably restarted the server since this particular analysis, the RID is not available
    # Get active processes
    '''
    try:
        for dataProc in activeList:
            if len(str(dataProc)) > 1:
                if queueUsers[dataProc] == request.user.username:
                    stringDict[dataProc] = str(dataProc) + ";" + str(queueTimes[dataProc]) + ";" + \
                                            str(queueFuncs[dataProc]) + ";" + str(queueUsers[dataProc]) + ";ACTIVE\n"
    except Exception as ex:
        print "Problem with dataqueue console:", ex'''

    for key in sorted(stringDict.keys()):
        queueString += stringDict[key]
    queueDict = {'display': queueString}
    output = json.dumps(queueDict)
    return HttpResponse(output, content_type='application/json')


def stop(request):
    global activeList, stopList, stopped, stopDict, queueList, queueFuncs, queueTimes, queueUsers
    RID = str(request.GET['all'])
    # put RID in stopDict and the request will not run (RID:True key value pair, we just check if the key exists so a
    # small value makes this a smaller structure if it gets upscaled for some reason (unlikely outside of DoS attempts)
    stopDict[RID] = True
    try:
        try:
            # this section will error if request is being processed already, ie its for removing queued requests only
            # the section below should handle active processes
            qid = queueList[RID]
            thisFunc = queueFuncs[qid]
            functions.log(request, "QSTOP", thisFunc)
            queueList.pop(RID, 0)   # try to remove from queuelist before processing
        except:
            # if we are here then most likely the process was ongoing when stop came through
            pass
        # find request in activeList, if its there we can tell the function to stop it ASAP and return user success msg
        pid = activeList.index(RID)
        stopList[pid] = RID
        activeList[pid] = 0
        functions.log(request, "QSTOP", request.user)
        myDict = {'error': 'none', 'message': 'Your analysis has been stopped!'}
        stop = json.dumps(myDict)
        return HttpResponse(stop, content_type='application/json')
    except Exception as e:
        myDict = {'error': 'Analysis not running'}
        stop = json.dumps(myDict)
        return HttpResponse(stop, content_type='application/json')


def decremQ():
    for thing in statDict:
        statDict[thing] -= 1


def process(pid):
    RID = "NULL_RID"    # this should only come up if the queue errors BEFORE reading request RID, an odd occurrence
    request = None
    while True:
        try:
            global activeList, stopList, stopDict

            # get next entry from queue
            data = q.get(block=True, timeout=None)
            decremQ()
            RID = data['RID']
            queueList.pop(RID, 0)   # remove from queue position tracker
            if RID in stopDict.keys():
                stopDict.pop(RID, 0)
            else:
                funcName = data['funcName']
                request = data['request']
                activeList[pid] = RID
                functions.log(request, "QSTART", funcName)
                if activeList[pid] == RID:
                    # TODO 1.3 finish moving analyses into analysis.py classes
                    # TODO 1.3 cleanup and history are treating getNorm like a proper analysis (no point in norm history view)
                    # TODO 1.3 put a limit on how many requests can be queued from the same user at one time (queueUsers has the needed info)
                    if funcName == "getNorm":   # at present likely not worthwhile to port norm to analysis class
                        recent[RID] = functions.getNorm(request, RID, stopList, pid)
                    elif funcName == "getCatUnivData":
                        myAnalysis = analysis.Anova(request, RID, stopList, pid)
                        recent[RID] = myAnalysis.run()  # can make this line shared once all analyses are ported over
                    elif funcName == "getQuantUnivData":
                        myAnalysis = analysis.Anova(request, RID, stopList, pid)
                        recent[RID] = myAnalysis.run(quant=True)
                    elif funcName == "getCorr":
                        myAnalysis = analysis.Corr(request, RID, stopList, pid)
                        recent[RID] = myAnalysis.run()
                    elif funcName == "getPCA":
                        myAnalysis = analysis.PCA(request, RID, stopList, pid)
                        recent[RID] = myAnalysis.run()
                    elif funcName == "getPCoA":
                        myAnalysis = analysis.PCoA(request, RID, stopList, pid)
                        recent[RID] = myAnalysis.run()
                    elif funcName == "getRF":
                        myAnalysis = analysis.Caret(request, RID, stopList, pid)
                        recent[RID] = myAnalysis.run()
                    elif funcName == "getDiffAbund":
                        myAnalysis = analysis.diffAbund(request, RID, stopList, pid)
                        recent[RID] = myAnalysis.run()
                    elif funcName == "getGAGE":
                        myAnalysis = analysis.Gage(request, RID, stopList, pid)
                        recent[RID] = myAnalysis.run()
                    elif funcName == "getCore":
                        myAnalysis = analysis.Core(request, RID, stopList, pid)
                        recent[RID] = myAnalysis.run()
                    elif funcName == "getSPLS":
                        recent[RID] = functions.getSPLS(request, stopList, RID, pid)
                    elif funcName == "getWGCNA":
                        recent[RID] = functions.getWGCNA(request, stopList, RID, pid)
                    elif funcName == "getSpAC":
                        recent[RID] = functions.getSpAC(request, stopList, RID, pid)
                    elif funcName == "getsoil_index":
                        recent[RID] = functions.getsoil_index(request, stopList, RID, pid)
                    else:
                        # received a function name that we don't support
                        # either a developer typo exists (here or on a web page) OR someone is trying to break things
                        # either way, this is a notable event, so we print to console and log
                        print "Security check:", request.user.username, "attempting to call function", funcName, "in AQ"
                        functions.log(request, "INVALID_FUNCTION_NAME_AQ", funcName)
                        # should probably return something to the user page so they don't get stuck
                        myDict = {'error': 'Invalid function name'}
                        stop = json.dumps(myDict)
                        recent[RID] = HttpResponse(stop, content_type='application/json')
                debug("Finished an analysis iteration")
                activeList[pid] = 0
                stopDict.pop(RID, 0)
                if "{\"error\"" in str(recent[RID]).split(":")[1]:  # This assumes error messages will always put error in front (its the ony section)
                    # while a proper results set will have 'error: none' much later
                    functions.log(request, "ERROR_AQ", funcName)
                else:
                    functions.log(request, "QFINISH", funcName)
            cleanup(RID, queueUsers[RID])
        except Exception as e:
            print "Error during analysis queue:", e
            if request is not None:
                functions.log(request, "ERROR_AQ", str(e)+"\n")  # atm this triggers when an UNHANDLED exception occurs. Handled check is the other ERROR_AQ line
            myDict = {'error': "Exception: "+str(e.message)}
            # this depends on the page in question directly, as each is responsible for handling 'error' in response
            stop = json.dumps(myDict)
            recent[RID] = HttpResponse(stop, content_type='application/json')


def cleanup(RID, username):   # cleanup and removeRID two parts of the same concept, group and add recent and such cleanup
    global queueFuncs, queueTimes, queueUsers, queueList, cleanupQueue, base, stage, time1, time2, TimeDiff
    # TODO 1.3 cleanup rplots and such (find directories named RID and remove them and their content)
    # get the user's recent RID list, add the new RID, and potentially receive an older RID to add to the queue
    myProf = UserProfile.objects.get(user=User.objects.get(username=username))
    toClean = myProf.addRecentRID(RID)
    if toClean is not None:
        cleanupQueue.put(toClean, True)
    myProf.save()
    #print "CLEANUP CALL FROM", RID, ", QUEUEING", toClean, "FOR CLEANUP
    cleanRID = ''
    try:
        # from here, see if something's in the queue, try to clean it up
        cleanRID = cleanupQueue.get(block=False, timeout=1)
    except:
        pass
    try:
        queueList.pop(cleanRID, 0)
    except:
        pass

    try:
        queueFuncs.pop(cleanRID, 0)
    except:
        pass

    try:
        # cleaning up old files
        delPath = "myPhyloDB/media/temp"
        for fname in os.walk(delPath):
            for sub in fname[2]:
                try:
                    if sub.split('.')[0] == str(cleanRID):
                        fullFname = str(fname[0]) + "/" + str(sub)
                        #print "Removing file:", fullFname
                        os.remove(fullFname)
                except Exception as fileEx:
                    print "Error during tempfile deletion:", fileEx
                    pass
            for sub in fname[1]:        # atm these check filenames only, not directories
                try:
                    if sub.split('.')[0] == str(cleanRID):
                        fullFname = str(fname[0]) + "/" + str(sub)
                        #print "Removing directory:", fullFname
                        shutil.rmtree(fullFname)
                except Exception as fileEx:
                    print "Error during tempfile deletion:", fileEx
                    pass

    except Exception as ex:
        print "Error during file cleanup:", ex
        pass

    try:
        queueTimes.pop(cleanRID, 0)
    except:
        pass

    try:
        queueUsers.pop(cleanRID, 0)
    except:
        pass

    try:
        recent.pop(cleanRID, 0)
    except:
        pass

    try:
        statDict.pop(cleanRID, 0)
    except:
        pass

    try:
        base.pop(cleanRID, None)
    except:
        pass

    try:
        stage.pop(cleanRID, None)
    except:
        pass

    try:
        time1.pop(cleanRID, None)
    except:
        pass

    try:
        time2.pop(cleanRID, None)
    except:
        pass

    try:
        TimeDiff.pop(cleanRID, None)
    except:
        pass

    try:
        complete.pop(cleanRID, None)
    except:
        pass

    gc.collect()  # attempting to cleanup memory leaks since processes technically are still there


def funcCall(request):
    # new hybrid call
    global activeList, stopList, statDict, qList, base, stage, time1, time2, TimeDiff, complete
    allJson = request.body.split('&')[0]
    data = json.loads(allJson)
    reqType = data['reqType']
    RID = data['RID']
    if reqType == "call":
        funcName = data['funcName']
        try:
            dataID = data['dataID']
        except Exception:
            myDict = {}
            myDict['error'] = "Error: Missing dataID"
            json_data = json.dumps(myDict)
            return HttpResponse(json_data,  content_type='application/json')

        if dataID == UserProfile.objects.get(user=request.user).dataID:
            time1[RID] = time()
            qDict = {'RID': RID, 'funcName': funcName, 'request': request}
            qList.append(qDict)  # what on earth does this do?
            # add to queue
            queueList[RID] = RID
            queueFuncs[RID] = funcName  # could change funcName here to match page name, for readability
            queueTimes[RID] = str(datetime.datetime.now())
            queueUsers[RID] = request.user.username
            q.put(qDict, True)
            statDict[RID] = int(q.qsize())
            complete[RID] = False

            functions.log(request, "QADD", funcName)

            myDict = {}
            myDict['resType'] = "status"
            myDict['error'] = "none"
            json_data = json.dumps(myDict)
            return HttpResponse(json_data, content_type='application/json')
        else:
            myDict = {}
            myDict['error'] = "Error: Selected data has changed, please refresh the page"
            json_data = json.dumps(myDict)
            return HttpResponse(json_data, content_type='application/json')

    if reqType == "status":
        try:
            results = recent[RID]   # recent[RID] is only initialized if function has returned with data
            # return results to client
            #print "Returning to client for RID", RID
            return results  # WE DID IT
        except KeyError:
            time2[RID] = time()
            try:
                TimeDiff[RID] = time2[RID] - time1[RID]
            except Exception:
                TimeDiff[RID] = 0

            try:
                if RID in queueList:  # has not started running yet but data is valid
                    stage[RID] = 'Analysis has been placed in queue, there are ' + str(stat(RID)) + ' others in front of you.'
                else:   # analysis has started, data is valid
                    stage[RID] = str(base[RID]) + '\n<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
            except Exception:   # data is not valid
                if TimeDiff[RID] == 0:
                    try:
                        if not complete[RID]:   # complete is valid, but false. Bug if encountered?
                            stage[RID] = 'In queue'
                        else:   # complete is valid and true, results are probably on the way
                            stage[RID] = 'Analysis complete, preparing results'  # stuck on one of these occasionally
                            #print "Finished"
                    except Exception:   # complete is not valid, implying variables have been cleaned up (or not created yet?)
                        stage[RID] = 'Request already completed'  # really shouldn't see this, if so BUG!!
                else:   # timediff not zero, so analysis has started to run despite data not valid (?)
                    stage[RID] = 'Analysis starting'    # getting this when queued
                    #print "Start"

            myDict = {'stage': stage[RID], 'resType': 'status'}
            json_data = json.dumps(myDict)
            return HttpResponse(json_data, content_type='application/json')
        except Exception as e:
            print "Error: ", e


def stat(RID):
    return statDict[RID] - stopped

