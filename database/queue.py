import Queue
from time import sleep
from anova import anova_graphs
from diffabund import diffabund_graphs
from norm import norm_graphs
from pcoa import pcoa_graphs
from pca import pca_graphs
from gage import gage_graphs
from spls import spls_graphs
from wgcna import wgcna_graphs
import simplejson
from django.http import HttpResponse


myQueue = Queue.Queue()  # so many vowels!
recent = {}
stopDict = {}
stopList = []
activeList = []
statDict = {}
stopped = 0


def funcCall(request):
    global stopDict
    allJson = request.GET["all"]
    data = simplejson.loads(allJson)
    RID = data['RID']
    funcName = data['funcName']
    stopDict[RID] = False
    queueDict = {'RID': RID, 'funcName': funcName, 'request': request}
    myQueue.put(queueDict, True)
    statDict[RID] = myQueue.qsize() - stopped
    while True:
        try:
            results = recent[RID]
            recent.pop(RID, 0)
            stopDict.pop(RID, 0)
            statDict.pop(RID, 0)
            return results
        except KeyError:
            if RID not in activeList and stopDict[RID]:
                stopDict.pop(RID, 0)
                myDict = {'error': 'none', 'message': 'Your analysis has been stopped!'}
                stop = simplejson.dumps(myDict)
                statDict.pop(RID, 0)
                return HttpResponse(stop, content_type='application/json')
        sleep(1)
        #print "I'm tired, RID: "+str(RID)+", funcName: "+str(funcName)


def stat(RID):
    return statDict[RID]


def stop(request):
    global stopped
    # find index of RID in activeList, set stopList at same index to True
    RID = str(request.GET["all"])
    try:
        if RID in activeList:
            stopList[activeList.index(RID)] = True
        stopDict[RID] = True
    except:
        stopped += 1

    myDict = {}
    myDict['error'] = 'none'
    myDict['message'] = 'Your analysis has been stopped!'
    stop = simplejson.dumps(myDict)

    return HttpResponse(stop, content_type='application/json')


def decQueue():
    for key in statDict:
        statDict[key] -= 1


def process(stop, PID):
    global stopped
    stopList.append(False)
    activeList.append(0)
    while stop[0] == 'True':
        try:
            curDict = myQueue.get(False)
            decQueue()
            myRID = curDict['RID']
            myFunc = curDict['funcName']
            myRequest = curDict['request']

            if not stopDict[myRID]:
                activeList[PID] = myRID
                # standard args are: data, stop, request ID, process ID
                if myFunc == "getNorm":
                    recent[curDict['RID']] = norm_graphs.getNorm(myRequest, stopList, myRID, PID)
                if myFunc == "getCatUnivData":
                    recent[curDict['RID']] = anova_graphs.getCatUnivData(myRequest, stopList, myRID, PID)
                if myFunc == "getQuantUnivData":
                    recent[curDict['RID']] = anova_graphs.getQuantUnivData(myRequest, stopList, myRID, PID)
                if myFunc == "getPCoA":
                    recent[curDict['RID']] = pcoa_graphs.getPCoA(myRequest, stopList, myRID, PID)
                if myFunc == "getPCA":
                    recent[curDict['RID']] = pca_graphs.getPCA(myRequest, stopList, myRID, PID)
                if myFunc == "getDiffAbund":
                    recent[curDict['RID']] = diffabund_graphs.getDiffAbund(myRequest, stopList, myRID, PID)
                if myFunc == "getGAGE":
                    recent[curDict['RID']] = gage_graphs.getGAGE(myRequest, stopList, myRID, PID)
                if myFunc == "getSPLS":
                    recent[curDict['RID']] = spls_graphs.getSPLS(myRequest, stopList, myRID, PID)
                if myFunc == "getWGCNA":
                    recent[curDict['RID']] = wgcna_graphs.getWGCNA(myRequest, stopList, myRID, PID)
                stopList[PID] = False
                activeList[PID] = 0
            else:
                stopped -= 1
        except:
            sleep(1)
    return

