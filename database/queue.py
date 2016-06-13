from django.http import HttpResponse, HttpResponseNotFound
from Queue import Queue
import simplejson
from time import sleep
import threading

from anova import anova_graphs
from diffabund import diffabund_graphs
from norm import norm_graphs
from pcoa import pcoa_graphs
from pca import pca_graphs
from gage import gage_graphs
from spls import spls_graphs
from wgcna import wgcna_graphs
from spac import spac_graphs
from soil_index import soil_index_graphs
from database.utils import threads


q = Queue(maxsize=0)
activeList = [0] * threads()
stopList = [0] * threads()
recent = {}
qList = []
threadDict = {}
statDict = {}
stopDict = {}
stopped = 0


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
                stop = simplejson.dumps(myDict)
                return HttpResponse(stop, content_type='application/json')
    except:
        myDict = {'error': 'Analysis not running'}
        stop = simplejson.dumps(myDict)
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
                    recent[RID] = norm_graphs.getNorm(request, RID, stopList, pid)
                if funcName == "getCatUnivData":
                    recent[RID] = anova_graphs.getCatUnivData(request, RID, stopList, pid)
                if funcName == "getQuantUnivData":
                    recent[RID] = anova_graphs.getQuantUnivData(request, RID, stopList, pid)
                if funcName == "getPCA":
                    recent[RID] = pca_graphs.getPCA(request, stopList, RID, pid)
                if funcName == "getPCoA":
                    recent[RID] = pcoa_graphs.getPCoA(request, stopList, RID, pid)
                if funcName == "getDiffAbund":
                    recent[RID] = diffabund_graphs.getDiffAbund(request, stopList, RID, pid)
                if funcName == "getGAGE":
                    recent[RID] = gage_graphs.getGAGE(request, stopList, RID, pid)
                if funcName == "getSPLS":
                    recent[RID] = spls_graphs.getSPLS(request, stopList, RID, pid)
                if funcName == "getWGCNA":
                    recent[RID] = wgcna_graphs.getWGCNA(request, stopList, RID, pid)
                if funcName == "getSpAC":
                    recent[RID] = spac_graphs.getSpAC(request, stopList, RID, pid)
                if funcName == "getsoil_index":
                    recent[RID] = soil_index_graphs.getsoil_index(request, stopList, RID, pid)
            activeList[pid] = ''
            threadDict.pop(RID, 0)
        sleep(1)


def funcCall(request):
    global activeList, stopList, stopDict, statDict, qList
    allJson = request.GET["all"]
    data = simplejson.loads(allJson)
    RID = data['RID']
    funcName = data['funcName']
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
            if results is None:  # if stopped as active process
                return HttpResponseNotFound()  # return rnf (better than None)
            return results
        except KeyError:
            if RID in stopList:
                try:
                    pid = stopList.index(RID)
                    # stopList[pid] = 0
                    # set stopList back to 0 after it is read in primary methods, here relies on timing to work
                except:
                    pass
                statDict.pop(RID, 0)
                # print "funcCall returning NF"
                return HttpResponseNotFound()
        sleep(1)


def stat(RID):
    return statDict[RID] - stopped

