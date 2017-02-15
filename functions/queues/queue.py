from django.http import HttpResponse
from Queue import Queue
import json
import multiprocessing as mp
from time import sleep, time
import threading

import config.local_cfg
from database.models import UserProfile

import functions


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
    global activeList, stopList, threadDict, stopDict
    RID = str(request.GET['all'])
    stopDict[RID] = True
    try:
        pid = activeList.index(RID)
        stopList[pid] = RID
        activeList[pid] = 0
        threadName = threadDict[RID]
        threads = threading.enumerate()
        for thread in threads:
            if thread.name == threadName:
                thread.terminate()
                stopList[pid] = RID
                myDict = {'error': 'none', 'message': 'Your analysis has been stopped!'}
                stop = json.dumps(myDict)
                return HttpResponse(stop, content_type='application/json')
    except Exception:
        myDict = {'error': 'Analysis not running'}
        stop = json.dumps(myDict)
        return HttpResponse(stop, content_type='application/json')


def decremQ():
    for thing in statDict:
        statDict[thing] -= 1


def process(pid):
    global activeList, threadDict, stopped
    while True:
        data = q.get(block=True, timeout=None)
        decremQ()
        RID = data['RID']
        if RID in stopDict:
            stopDict.pop(RID, 0)
        else:
            funcName = data['funcName']
            request = data['request']
            activeList[pid] = RID
            thread = threading.current_thread()
            threadDict[RID] = thread.name
            if activeList[pid] == RID:
                if funcName == "getNorm":
                    recent[RID] = functions.getNorm(request, RID, stopList, pid)
                if funcName == "getCatUnivData":
                    recent[RID] = functions.getCatUnivData(request, RID, stopList, pid)
                if funcName == "getQuantUnivData":
                    recent[RID] = functions.getQuantUnivData(request, RID, stopList, pid)
                if funcName == "getPCA":
                    recent[RID] = functions.getPCA(request, stopList, RID, pid)
                if funcName == "getPCoA":
                    recent[RID] = functions.pcoa_graphs.getPCoA(request, stopList, RID, pid)
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
        sleep(1)


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
            print "Missing dataID"
            myDict = {}
            myDict['error'] = "Error: Dev done goofed!"
            json_data = json.dumps(myDict)
            return HttpResponse(json_data,  content_type='application/json')

        if dataID == UserProfile.objects.get(user=request.user).dataID:
            time1[RID] = time()
            qDict = {'RID': RID, 'funcName': funcName, 'request': request}
            qList.append(qDict)
            q.put(qDict, True)
            statDict[RID] = int(q.qsize())
            complete[RID] = False

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
            stopDict.pop(RID, 0)
            removeRID(RID)
            # results on anova quant not working occasionally, depends on categorical selection
            return results
        except KeyError:
            time2[RID] = time()
            try:
                TimeDiff[RID] = time2[RID] - time1[RID]
            except Exception:
                TimeDiff[RID] = 0

            try:
                if TimeDiff[RID] == 0:
                    stage[RID] = 'Analysis has been placed in queue, there are ' + str(stat(RID)) + ' others in front of you.'
                else:
                    stage[RID] = str(base[RID]) + '\n<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
            except Exception:
                if TimeDiff[RID] == 0:
                    try:
                        if not complete[RID]:
                            stage[RID] = 'In queue'
                        else:
                            stage[RID] = 'Analysis complete, preparing results'  # stuck on one of these occasionally
                    except Exception:
                        stage[RID] = 'Downloading results'  # maybe?
                else:
                    stage[RID] = 'Analysis starting'

            myDict = {'stage': stage[RID], 'resType': 'status'}
            json_data = json.dumps(myDict)
            return HttpResponse(json_data, content_type='application/json')
        except Exception as e:
            print "Error: ", e


def stat(RID):
    return statDict[RID] - stopped

