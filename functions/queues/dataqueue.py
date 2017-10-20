from django.http import HttpResponse
import json
from Queue import Queue
from time import sleep
import threading

import functions
import database.views

# imported for dataQueue fixes 2017/09/28
from django.shortcuts import render
from database.models import Reference
from database.forms import UploadForm1, UploadForm2, UploadForm4, UploadForm5, UploadForm6, UploadForm7, UploadForm8, UploadForm10

#from pycallgraph import PyCallGraph
#from pycallgraph.output import GephiOutput
#from pycallgraph import Config

import datetime

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
        try:    # check queuelist
            pid = datQueueList[RID]

            # wait for status, then cleanup
            # get function from funccall, return correct page rendering
            thisFunc = datQueueFuncs[pid]
            # can't populate these from here
            projects = Reference.objects.none()
            if request.user.is_superuser:
                projects = Reference.objects.all().order_by('projectid__project_name', 'path')
            elif request.user.is_authenticated():
                projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
            # uploadFunc reanalyze updateFunc pybake
            retStuff = None
            if thisFunc == "uploadFunc":
                print "Upload stopping early"
                retStuff = render(
                        request,
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': ""}
                    )
            elif thisFunc == "reanalyze":
                retStuff = render(
                        request,
                        'reprocess.html',
                        {'form4': UploadForm4,
                         'mform': UploadForm10},
                    )
            elif thisFunc == "updateFunc":
                state = ''
                retStuff = render(
                        request,
                        'update.html',
                        {'form5': UploadForm5,
                         'state': state}
                    )
            elif thisFunc == "pybake":
                form6 = UploadForm6(request.POST, request.FILES)
                form7 = UploadForm7(request.POST, request.FILES)
                form8 = UploadForm8(request.POST, request.FILES)

                if form6.is_valid():
                    file1 = request.FILES['taxonomy']
                    file2 = request.FILES['precalc_16S']
                    file3 = request.FILES['precalc_KEGG']
                    functions.geneParse(file1, file2, file3)

                if form7.is_valid():
                    file4 = request.FILES['ko_htext']
                    functions.koParse(file4)

                if form8.is_valid():
                    file5 = request.FILES['nz_htext']
                    functions.nzParse(file5)

                retStuff = render(
                    request,
                    'pybake.html',
                    {'form6': UploadForm6,
                     'form7': UploadForm7,
                     'form8': UploadForm8}
                )

            if retStuff == None:
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
            datQueueList.pop(RID, 0)
            datActiveList[pid] = RID
            thread = threading.current_thread()
            datThreadDict[RID] = thread.name
            msg = str(datetime.datetime.now())+": Starting on request from " + str(request.user.username) + " for " + str(funcName)
            functions.log(msg)
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
            datStopDict.pop(RID, 0)  # clean this up when done or stopped, needs to be here since funcCall can end early
            msg = str(datetime.datetime.now())+": Finished request from " + str(request.user.username) + " for " + str(funcName)
            functions.log(msg)
        sleep(1)


def datfuncCall(request):
    global datActiveList, datQueueList, datQueueFuncs, datStopList, datStopDict, datStatDict, datQList
    RID = request.POST['RID']
    request.POST['stopList'] = datStopList
    funcName = request.POST['funcName']
    datQueueList[RID] = RID  # add to queuelist, remove when processed or stopped
    datQueueFuncs[RID] = funcName
    qDict = {'RID': RID, 'funcName': funcName, 'request': request}
    datQList.append(qDict)
    datQ.put(qDict, True)
    datStatDict[RID] = int(datQ.qsize())

    # print log info, need to write this to a file somewhere
    msg = str(datetime.datetime.now())+": Received request from " + str(request.user.username) + " for " + str(funcName)
    functions.log(msg)

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
        return datStatDict[RID] - datStopped


def getStops():
    return datStopList

