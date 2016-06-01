import Queue
from time import sleep
import cherrypy
# four functions, enqueue, dequeue, process, stop

myQueue = Queue.Queue()  # so many vowels!
recent = {}
stops = []


def funcCall(RID, funcName, data):
    queueDict = {'RID': RID, 'funcName': funcName, 'data': data}
    enqueue(queueDict)
    while True:
        if RID in stops:
            return
        try:
            results = recent[RID]
            recent.pop(RID, 0)
            return results
        except Exception as e:
            pass  # not currently valid, check again
        sleep(1)


def enqueue(obj):
    # add dict to queue
    # obj assumed to be dictionary with keys funcName, RID, and data
    myQueue.put(obj, True)
    return


def dequeue():
    # remove first dict, return it
    return myQueue.get(False)


def stop(RID):
    print "Not yet implemented!"
    stops.append(RID)
    # add RID to stop list
    # when starting a new process, check if RID is in this list
    # if so, skip it
    # if RID matches an active process, kill it
    # kill process via stop tag passed down to each level of child process, with periodic/looping stop checks


def process(stop):
    while stop[0] == 'True':
        try:
            curDict = dequeue()
            sleep(5)  # simulate run time
            # get main values
            myRID = curDict['RID']
            myFunc = curDict['funcName']
            myData = curDict['data']
            # check RID stops for skip
            if myRID in stops:
                print "Skipping RID "+myRID
                continue
            else:
                # run func
                print "Actually running, calling "+str(myFunc)
                recent[curDict['RID']] = "results!"
        except:
            pass  # empty queue?
        sleep(1)
    return

