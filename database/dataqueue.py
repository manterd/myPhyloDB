from django.http import HttpResponse, HttpResponseNotFound
from Queue import Queue
import simplejson
from time import sleep
import threading

from views import uploadFunc, updateFunc, pybake
from parsers import reanalyze
from database.utils import threads


q = Queue(maxsize=0)
activeList = [0]
stopList = [0]  # change if adding more data threads
recent = {}
qList = []
threadDict = {}
statDict = {}
stopDict = {}
stopped = 0


def datstop(request):
    print "Told to stop a thing"
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
                stop = simplejson.dumps(myDict)
                return HttpResponse(stop, content_type='application/json')
    except:
        myDict = {'error': 'Analysis not running'}
        stop = simplejson.dumps(myDict)
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
                    recent[RID] = uploadFunc(request, stopList)
                if funcName == "reanalyze":
                    recent[RID] = reanalyze(request, stopList)
                if funcName == "updateFunc":
                    recent[RID] = updateFunc(request, stopList)
                if funcName == "pybake":
                    recent[RID] = pybake(request)
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
            # if results is None:  # if stopped as active process
            #     return HttpResponseNotFound()  # return rnf (better than None)
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
    except:
        return statDict[RID] - stopped


def getStops():
    return stopList

