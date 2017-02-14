from django.http import HttpResponse
import json
from Queue import Queue
from time import sleep
import threading

import database.utils
import database.views


q = Queue(maxsize=0)
activeList = [0]
stopList = [0]
recent = {}
qList = []
threadDict = {}
statDict = {}
stopDict = {}
stopped = 0


def datstop(request):
    global activeList, stopList, threadDict, stopDict, stopped
    RID = str(request.GET['all'])
    stopDict[RID] = True
    try:
        pid = activeList.index(RID)
        stopList[pid] = RID
        activeList[pid] = 0
        threadName = threadDict[RID]
        database.utils.termP()
        threads = threading.enumerate()
        for thread in threads:
            if thread.name == threadName:
                thread.terminate()
                stopList[pid] = RID
                myDict = {'error': 'none', 'message': 'Your analysis has been stopped!'}
                stop = json.dumps(myDict)
                stopped += 1
                return HttpResponse(stop, content_type='application/json')
    except Exception:
        myDict = {'error': 'Analysis not running'}
        stop = json.dumps(myDict)
        return HttpResponse(stop, content_type='application/json')


def decremQ():
    for thing in statDict:
        statDict[thing] -= 1


def dataprocess(pid):
    global activeList, threadDict, stopped, recent
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
                if funcName == "uploadFunc":
                    recent[RID] = database.views.uploadFunc(request, stopList)
                if funcName == "reanalyze":
                    resp = database.utils.reanalyze(request, stopList)
                    if resp is None:
                        resp = database.views.reprocess(request)
                    if resp == "Stopped":
                        resp = database.views.reprocess(request)
                        resp['error'] = "Reprocessing stopped"
                    recent[RID] = resp
                if funcName == "updateFunc":
                    recent[RID] = database.views.updateFunc(request, stopList)
                if funcName == "pybake":
                    recent[RID] = database.views.pybake(request)
            activeList[pid] = ''
            threadDict.pop(RID, 0)
        sleep(1)


def datfuncCall(request):
    global activeList, stopList, stopDict, statDict, qList
    RID = request.POST['RID']
    request.POST['stopList'] = stopList
    funcName = request.POST['funcName']
    qDict = {'RID': RID, 'funcName': funcName, 'request': request}
    qList.append(qDict)
    q.put(qDict, True)
    statDict[RID] = int(q.qsize())
    while True:
        try:
            results = recent[RID]
            recent.pop(RID, 0)
            statDict.pop(RID, 0)
            stopDict.pop(RID, 0)
            return results
        except KeyError:
            if RID in stopList:
                statDict.pop(RID, 0)
        except Exception as e:
            print "Unexpected exception: "+str(e)
        sleep(1)


def datstat(RID):
    try:
        stopDict[RID] = stopDict[RID]
        return -1024
    except Exception:
        return statDict[RID] - stopped


def getStops():
    return stopList

