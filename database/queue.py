from django.http import HttpResponse
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

"""
This still doesn't work quite right.

Check out the following pages.
http://eli.thegreenplace.net/2011/05/18/code-sample-socket-client-thread-in-python/
http://eli.thegreenplace.net/2011/12/27/python-threads-communication-and-stopping
"""


q = Queue()
num_threads = 3
activeList = [0] * num_threads
stopList = [0] * num_threads
recent = {}
qList = []
threadDict = {}


def stop(request):
    global activeList, stopList, threadDict
    RID = str(request.GET['all'])
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


def process(pid):
    global activeList, stopList, threadDict
    while True:
        data = q.get(block=True, timeout=None)
        RID = data['RID']
        funcName = data['funcName']
        request = data['request']
        activeList[pid] = RID
        thread = threading.current_thread()
        threadDict[RID] = thread.name
        if activeList[pid] == RID:
            if funcName == "getNorm":
                recent[RID] = norm_graphs.getNorm(request, RID, pid)
            if funcName == "getCatUnivData":
                recent[RID] = anova_graphs.getCatUnivData(request, RID, pid)
        sleep(1)
        activeList[pid] = ''
        threadDict.pop(RID, 0)


def funcCall(request):
    global activeList, stopList
    allJson = request.GET["all"]
    data = simplejson.loads(allJson)
    RID = data['RID']
    funcName = data['funcName']
    qDict = {'RID': RID, 'funcName': funcName, 'request': request}
    qList.append(qDict)
    q.put(qDict, True)

    while True:
        try:
            results = recent[RID]
            recent.pop(RID, 0)
            return results
        except KeyError:
            if RID in stopList:
                pid = stopList.index(RID)
                stopList[pid] = 0
                return HttpResponse('')
        sleep(1)
        print "active:", activeList
        print "stop:", stopList


