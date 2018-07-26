from django.http import HttpResponse
import json
from Queue import Queue
from time import sleep
import threading

import functions
import database.views

from database.models import Reference

# only needed for code mapping
#from pycallgraph import PyCallGraph
#from pycallgraph.output import GephiOutput
#from pycallgraph import Config


datQ = Queue(maxsize=0)
datActiveList = [0]
datQueueList = {}
datQueueFuncs = {}
datStopList = [0]
datRecent = {}
datQList = []
datThreadDict = {}
datStatDict = {}
datStopDict = {}
datStopped = 0


def datstop(request):
    global datActiveList, datQueueList, datQueueFuncs, datStopList, datThreadDict, datStopDict, datStopped, datRecent
    RID = str(request.GET['all'])
    datStopDict[RID] = RID
    try:
        pid = datActiveList.index(RID)
        datStopList[pid] = RID
        datActiveList[pid] = 0
        functions.termP()   # only time this is called in entire code
        datStopList[pid] = RID
        myDict = {'error': 'none', 'message': 'Your database changes have been cancelled!'}  # canceled is valid
        stop = json.dumps(myDict)
        return HttpResponse(stop, content_type='application/json')
    except Exception:   # not in activelist
        try:
            # check queuelist
            # wait for status, then cleanup
            # get function from funccall, return correct page rendering
            pid = datQueueList[RID]
            thisFunc = datQueueFuncs[pid]

            # uploadFunc reanalyze updateFunc pybake
            functions.log(request, "QSTOP", thisFunc)
            retStuff = None
            if thisFunc == "uploadFunc":
                retStuff = database.views.uploadFunc(request, datStopList)
            elif thisFunc == "reanalyze":
                retStuff = database.views.reprocess(request)
            elif thisFunc == "updateFunc":
                retStuff = database.views.updateFunc(request, datStopList)
            elif thisFunc == "pybake":
                retStuff = database.views.pybake(request)

            if retStuff is None:
                print "This is a problem"
            datRecent[pid] = retStuff

            myDict = {'error': 'none', 'message': 'Your database changes have been cancelled!'}
            stop = json.dumps(myDict)
            return HttpResponse(stop, content_type='application/json')

        except Exception as er:   # not in queuelist either. actual error
            print "Error with datStop:", er
            myDict = {'error': 'Analysis not running'}
            stop = json.dumps(myDict)
            return HttpResponse(stop, content_type='application/json')


def decremQ():
    for thing in datStatDict:
        datStatDict[thing] -= 1


def dataprocess(pid):
    global datActiveList, datQueueList, datThreadDict, datStopped, datRecent
    #count = 0   # counts number of functions cycled through
    #graphConfig = Config(max_depth=3)
    while True:
        #graphOut = GephiOutput(output_file='datFuncCall'+str(count)+'.gdf')
        #count += 1
        #with PyCallGraph(output=graphOut, config=graphConfig):  # when not debugging, tab back under here
        data = datQ.get(block=True, timeout=None)
        decremQ()
        RID = data['RID']
        if RID in datStopDict:
            datStopDict.pop(RID, 0)
        else:
            funcName = data['funcName']
            request = data['request']
            stopList = data['stop']
            datQueueList.pop(RID, 0)
            datActiveList[pid] = RID
            thread = threading.current_thread()
            datThreadDict[RID] = thread.name
            functions.log(request, "QSTART", funcName)
            if datActiveList[pid] == RID:
                if funcName == "uploadFunc":
                    datRecent[RID] = database.views.uploadFunc(request, stopList)
                if funcName == "fileUpFunc":
                    datRecent[RID] = database.views.fileUpFunc(request, stopList)
                if funcName == "reanalyze":
                    resp = functions.reanalyze(request, datStopList)
                    if resp is None:
                        resp = database.views.reprocess(request)
                    if resp == "Stopped":
                        resp = database.views.reprocess(request)
                        resp['error'] = "Reprocessing stopped"
                    datRecent[RID] = resp
                if funcName == "updateFunc":
                    datRecent[RID] = database.views.updateFunc(request, stopList)
                if funcName == "geneParse":
                    print 'starting'
                    datRecent[RID] = functions.geneParse(request)
                if funcName == "koParse":
                    datRecent[RID] = functions.koParse(request)
                if funcName == "nzParse":
                    datRecent[RID] = functions.nzParse(request)
            datActiveList[pid] = ''
            datThreadDict.pop(RID, 0)
            datStopDict.pop(RID, 0)  # clean this up when done or stopped, needs to be here since funcCall can end early
            functions.log(request, "QFINISH", funcName)
            sleep(1)    # TODO dos vuln? moved to after longer processes


def datfuncCall(request):
    global datActiveList, datQueueList, datQueueFuncs, datStopList, datStopDict, datStatDict, datQList
    RID = request.POST['RID']
    funcName = request.POST['funcName']
    datQueueList[RID] = RID  # add to queuelist, remove when processed or stopped
    datQueueFuncs[RID] = funcName
    qDict = {'RID': RID, 'funcName': funcName, 'request': request, 'stop': datStopList}
    datQList.append(qDict)
    datQ.put(qDict, True)
    datStatDict[RID] = int(datQ.qsize())

    # print log info, need to write this to a file somewhere
    functions.log(request, "QADD", funcName)
    while True:
        try:
            results = datRecent[RID]
            datRecent.pop(RID, 0)
            datStatDict.pop(RID, 0)
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
        if RID in datQueueList:
            return -512
        else:
            return -1024
    except Exception:
        # if not in stopDict (ie not called to stop) queuepos is current (minus number of pre stops)
        return datStatDict[RID] - datStopped


def getStops():
    return datStopList

