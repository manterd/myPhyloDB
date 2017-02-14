from django.http import HttpResponse
from Queue import Queue
import json
from time import sleep, time
import threading

import anova.anova_graphs
import diffabund.diffabund_graphs
import norm.norm_graphs
import pcoa.pcoa_graphs
import pca.pca_graphs
import rf.rf_graphs
import gage.gage_graphs
import spls.spls_graphs
import wgcna.wgcna_graphs
import spac.spac_graphs
import soil_index.soil_index_graphs
import utils
from models import UserProfile


q = Queue(maxsize=0)
activeList = [0] * utils.threads()
stopList = [0] * utils.threads()
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
                    recent[RID] = norm.norm_graphs.getNorm(request, RID, stopList, pid)
                if funcName == "getCatUnivData":
                    recent[RID] = anova.anova_graphs.getCatUnivData(request, RID, stopList, pid)
                if funcName == "getQuantUnivData":
                    recent[RID] = anova.anova_graphs.getQuantUnivData(request, RID, stopList, pid)
                if funcName == "getPCA":
                    recent[RID] = pca.pca_graphs.getPCA(request, stopList, RID, pid)
                if funcName == "getPCoA":
                    recent[RID] = pcoa.pcoa_graphs.getPCoA(request, stopList, RID, pid)
                if funcName == "getRF":
                    recent[RID] = rf.rf_graphs.getRF(request, stopList, RID, pid)
                if funcName == "getDiffAbund":
                    recent[RID] = diffabund.diffabund_graphs.getDiffAbund(request, stopList, RID, pid)
                if funcName == "getGAGE":
                    recent[RID] = gage.gage_graphs.getGAGE(request, stopList, RID, pid)
                if funcName == "getSPLS":
                    recent[RID] = spls.spls_graphs.getSPLS(request, stopList, RID, pid)
                if funcName == "getWGCNA":
                    recent[RID] = wgcna.wgcna_graphs.getWGCNA(request, stopList, RID, pid)
                if funcName == "getSpAC":
                    recent[RID] = spac.spac_graphs.getSpAC(request, stopList, RID, pid)
                if funcName == "getsoil_index":
                    recent[RID] = soil_index.soil_index_graphs.getsoil_index(request, stopList, RID, pid)
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

