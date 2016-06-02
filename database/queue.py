import Queue
from time import sleep
from anova import anova_graphs
import simplejson
from django.http import HttpResponse


myQueue = Queue.Queue()  # so many vowels!
recent = {}
stopDict = {}
stopList = [False]
activeList = []
statDict = {}


def funcCall(request):
    allJson = request.GET["all"]
    data = simplejson.loads(allJson)
    RID = data['RID']
    funcName = data['funcName']
    stopDict[RID] = False
    # print "RID: ", RID
    # print "funcName: ", funcName
    # print "request: ", request
    queueDict = {'RID': RID, 'funcName': funcName, 'request': request}
    myQueue.put(queueDict, True)
    statDict[RID] = myQueue.qsize()
    while True:
        if stopDict[RID]:
            myDict = {}
            myDict['error'] = 'none'
            myDict['message'] = 'Your analysis has been stopped!'
            stop = simplejson.dumps(myDict)
            return HttpResponse(stop, content_type='application/json')
        try:
            results = recent[RID]
            recent.pop(RID, 0)
            return results
        except Exception as e:
            pass
        sleep(1)


def stat(RID):
    return statDict[RID]


def stop(request):
    # find index of RID in activeList, set stopList at same index to True
    RID = str(request.GET["all"])
    # print "RID: ", RID
    try:
        if RID in activeList:
            stopList[activeList.index(RID)] = True
    except:
        pass

    stopDict[RID] = True
    # print "stop! ", stopDict
    myDict = {}
    myDict['error'] = 'none'
    myDict['message'] = 'Your analysis has been stopped!'
    stop = simplejson.dumps(myDict)
    return HttpResponse(stop, content_type='application/json')


def decQueue():
    for key in statDict:
        statDict[key] -= 1


def process(stop):
    while stop[0] == 'True':
        try:
            curDict = myQueue.get(False)
            # print "Found something to run"
            decQueue()
            # get main values
            myRID = curDict['RID']
            myFunc = curDict['funcName']
            myRequest = curDict['request']

            PID = 0  # which active process is about to be started?

            # check RID stops for skip
            # print "? ", stopDict[myRID]
            if stopDict[myRID]:
                stopDict.pop(myRID, 0)
            else:
                # run func
                # print "Not skipping"
                activeList.append(myRID)
                if myFunc == "getCatUnivData":
                    # run cat anova
                    # print "Starting Anova Cat"
                    recent[curDict['RID']] = anova_graphs.getCatUnivData(myRequest, stopList, myRID, PID)
                    # print "Finished AC"
                    stopList[0] = False
                    # standard args are: data, stop, request ID, process ID
                if myFunc == "getQuantUnivData":
                    # print "Starting Anova Quant"
                    recent[curDict['RID']] = anova_graphs.getQuantUnivData(myRequest, stopList, myRID, PID)
                    # print "Finished AQ"
                    stopList[0] = False
                activeList.remove(myRID)
            # print "Finished a loop"
        except:
            sleep(1)
    return
# I changed a thing

