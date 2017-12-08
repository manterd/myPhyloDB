from django.http import HttpResponse
from Queue import Queue
import json
import multiprocessing as mp
from time import sleep, time
import threading
import gc

import config.local_cfg
from database.models import UserProfile

import functions

from functions.analysis import analysis

#from pycallgraph import PyCallGraph
#from pycallgraph.output import GephiOutput
#from pycallgraph import Config

import datetime


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
recent = {}
qList = []
threadDict = {}
statDict = {}
stopDict = {}
stopped = 0

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


def removeRID(RID):
    global base, stage, time1, time2, TimeDiff
    try:
        base.pop(RID, None)
        stage.pop(RID, None)
        time1.pop(RID, None)
        time2.pop(RID, None)
        TimeDiff.pop(RID, None)
        complete.pop(RID, None)
        return True
    except Exception:
        return False


def stop(request):
    global activeList, stopList, threadDict, stopDict, queueList, queueFuncs   # needs queuelist (maybe queuefuncs) from dataqueue, currently not at the same level
    RID = str(request.GET['all'])
    stopDict[RID] = RID
    try:
        pid = activeList.index(RID)
        stopList[pid] = RID
        activeList[pid] = 0
        try:
            pid = queueList[RID]    # hi
            thisFunc = queueFuncs[pid]
            functions.log(request, "QSTOP", thisFunc)
            queueList.pop(RID, 0)   # try to remove from queuelist
        except:
            pass    # already removed probably, moving on

        threadName = threadDict[RID]
        threads = threading.enumerate()
        for thread in threads:
            if thread.name == threadName:
                thread.terminate()      # is this necessary? can simplify this code a lot if not
                stopList[pid] = RID
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
    global activeList, threadDict, stopped
    #count = 0   # counts number of functions cycled through
    #graphConfig = Config(max_depth=3)
    while True:
        try:
            #graphOut = GephiOutput(output_file='funcCall'+str(count)+'.gdf')
            #count += 1
            #with PyCallGraph(output=graphOut, config=graphConfig):  # when not debugging, tab back under here
            data = q.get(block=True, timeout=None)
            decremQ()
            RID = data['RID']
            queueList.pop(RID, 0)   # remove from queue position tracker
            if RID in stopDict:
                stopDict.pop(RID, 0)
            else:
                funcName = data['funcName']
                request = data['request']
                activeList[pid] = RID
                thread = threading.current_thread()
                threadDict[RID] = thread.name
                functions.log(request, "QSTART", funcName)
                if activeList[pid] == RID:
                    if funcName == "getNorm":
                        recent[RID] = functions.getNorm(request, RID, stopList, pid)
                    if funcName == "getCatUnivData":
                        myAnalysis = analysis.Anova(request, RID, stopList, pid)
                        recent[RID] = myAnalysis.run()
                    if funcName == "getQuantUnivData":
                        recent[RID] = functions.getQuantUnivData(request, RID, stopList, pid)
                    if funcName == "getCorr":
                        myAnalysis = analysis.Corr(request, RID, stopList, pid)
                        recent[RID] = myAnalysis.run()
                    if funcName == "getPCA":
                        myAnalysis = analysis.PCA(request, RID, stopList, pid)
                        recent[RID] = myAnalysis.run()
                    if funcName == "getPCoA":
                        myAnalysis = analysis.PCoA(request, RID, stopList, pid)
                        recent[RID] = myAnalysis.run()
                    if funcName == "getRF":
                        recent[RID] = functions.getRF(request, stopList, RID, pid)
                    if funcName == "getDiffAbund":
                        recent[RID] = functions.getDiffAbund(request, stopList, RID, pid)
                    if funcName == "getGAGE":
                        recent[RID] = functions.getGAGE(request, stopList, RID, pid)
                    if funcName == "getSPLS":
                        recent[RID] = functions.getSPLS(request, stopList, RID, pid)
                    if funcName == "getWGCNA":
                        recent[RID] = functions.getWGCNA(request, stopList, RID, pid)
                    if funcName == "getSpAC":
                        recent[RID] = functions.getSpAC(request, stopList, RID, pid)
                    if funcName == "getsoil_index":
                        recent[RID] = functions.getsoil_index(request, stopList, RID, pid)
                activeList[pid] = ''
                threadDict.pop(RID, 0)
                stopDict.pop(RID, 0)
                functions.log(request, "QFINISH", funcName)
            sleep(1)
        except Exception as e:
            print "Error during primary queue:", e


def funcCall(request):
    # new hybrid call
    global activeList, stopList, stopDict, statDict, qList, base, stage, time1, time2, TimeDiff
    allJson = request.body.split('&')[0]
    data = json.loads(allJson)
    reqType = data['reqType']
    RID = data['RID']
    funcName = data['funcName']
    if reqType == "call":
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
            queueFuncs[RID] = funcName
            q.put(qDict, True)
            statDict[RID] = int(q.qsize())
            complete[RID] = False

            # print log info, need to write this to a file somewhere
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
            results = recent[RID]
            recent.pop(RID, 0)
            statDict.pop(RID, 0)
            removeRID(RID)
            gc.collect()  # attempting to cleanup memory leaks since processes technically are still there
            # results on anova quant not working occasionally, depends on categorical selection
            # return results to client, verify timing?
            #print "Returning a thing!"
            return results
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
                        stage[RID] = 'Downloading results'  # bug with timing if this sticks ?
                        # bug is more likely on the page side
                        #print "Done?"
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

