from django.http import HttpResponse
import json
from Queue import Queue
from time import sleep
import threading
import time

import functions
import database.views


datQ = Queue(maxsize=0)
datActiveList = [0]
datStopList = [0]
datRecent = {}
datQList = []
datThreadDict = {}
datStatDict = {}
datStopDict = {}
datStopped = 0


def datstop(request):   # this function needs reworked, some code needs to exit instead of being terminated, for cleanup
    global datActiveList, datStopList, datThreadDict, datStopDict, datStopped
    RID = str(request.GET['all'])
    datStopDict[RID] = True
    try:
        pid = datActiveList.index(RID)
        datStopList[pid] = RID
        datActiveList[pid] = 0
        threadName = datThreadDict[RID]
        functions.termP()   # only time this is called in entire code
        threads = threading.enumerate()
        for thread in threads:  # for each thread
            if thread.name == threadName:   # check if thread name matches
                #thread.terminate()  # terminate it
                datStopList[pid] = RID
                myDict = {'error': 'none', 'message': 'Your database changes have been cancelled!'}  # canceled is valid
                stop = json.dumps(myDict)
                # datStopped += 1   # might be redundant
                return HttpResponse(stop, content_type='application/json')
    except Exception:
        myDict = {'error': 'Analysis not running'}
        stop = json.dumps(myDict)
        return HttpResponse(stop, content_type='application/json')


def decremQ():
    for thing in datStatDict:
        datStatDict[thing] -= 1


def dataprocess(pid):
    global datActiveList, datThreadDict, datStopped, datRecent
    while True:
        data = datQ.get(block=True, timeout=None)
        decremQ()
        RID = data['RID']
        if RID in datStopDict:
            datStopDict.pop(RID, 0)
        else:
            funcName = data['funcName']
            request = data['request']
            datActiveList[pid] = RID
            thread = threading.current_thread()
            datThreadDict[RID] = thread.name
            if datActiveList[pid] == RID:
                if funcName == "uploadFunc":
                    # save that this is an upload in a dict somewhere, in stopdict, check if upload, removeproj if true
                    datRecent[RID] = database.views.uploadFunc(request, datStopList)
                if funcName == "reanalyze":
                    resp = functions.reanalyze(request, datStopList)
                    if resp is None:
                        resp = database.views.reprocess(request)
                    if resp == "Stopped":
                        resp = database.views.reprocess(request)
                        resp['error'] = "Reprocessing stopped"
                    datRecent[RID] = resp
                if funcName == "updateFunc":
                    datRecent[RID] = database.views.updateFunc(request, datStopList)
                if funcName == "pybake":
                    datRecent[RID] = database.views.pybake(request)
            datActiveList[pid] = ''
            datThreadDict.pop(RID, 0)
        sleep(1)


def datfuncCall(request):
    global datActiveList, datStopList, datStopDict, datStatDict, datQList
    RID = request.POST['RID']
    request.POST['stopList'] = datStopList
    funcName = request.POST['funcName']
    qDict = {'RID': RID, 'funcName': funcName, 'request': request}
    datQList.append(qDict)
    datQ.put(qDict, True)
    datStatDict[RID] = int(datQ.qsize())
    while True:
        try:
            results = datRecent[RID]
            datRecent.pop(RID, 0)
            datStatDict.pop(RID, 0)
            datStopDict.pop(RID, 0)
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
        return -1024
    except Exception:
        return datStatDict[RID] - datStopped


def getStops():
    return datStopList

