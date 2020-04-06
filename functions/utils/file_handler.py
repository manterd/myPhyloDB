import json
from django.http import HttpResponse
import hashlib
import os
import math
import datetime
import time
#from database.views import subDirList
from functions.utils.debug import debug
from functions.utils.utils_df import log, securityLog, errorLog
from threading import Lock

# k:v RID:(dict containing chunk data), for storing the raw data chunks sent for given RID
chunkStorageDict = {}
# k:v RID:(dict for saving chunk metrics), for storing duplicate send counts associated with an RID
chunkTallyDict = {}
# k:v RID:lastChunkNumSaved, for saving the index of the last (largest) chunk written to temp storage.
# Don't accept chunks below this number for given RID
bunchedChunkNumDict = {}
# dict for storing how many bunch files have been written for given RID
bunchCountDict = {}
# how many times can a chunk be resent before raising an error?
chunkRepeatLimit = 10
# how big is a chunk allowed to be? (can be smaller, for the most part this will match the data coming in)
chunkSize = 20000  # 20 KB chunk max, make sure this is synced with files.html value
# how much data will be stored in memory before being written to temp file?
maxMemory = 100000  # 100 KB max before offloading to temp file, this must be >= chunkSize
# derived value, how many chunks should be saved in memory before attempting to write to temp storage?
maximumMemorizedChunks = maxMemory / chunkSize

# lock for securing bunch-writing critical section of fileUploadChunk
bunchLock = Lock()

# need a way to keep track of how many times a duplicate chunk gets sent
mostRecentChunkTallyDict = {}
# metrics global variable, first value can be safely ignored (just making sure the type matches)
lastEndTime = time.time()

# TODO 1.3 need a process that cleans up requests which haven't finished after a period of time with no updating

import psutil
process = psutil.Process(os.getpid())
prevMem = 0.0
memTally = 0
enableMemDiff = False   # set this back to true for memory analysis during exploding_panda

# when adding new file types, just add name used in html (and for directory) to subDirList
subDirList = ['meta', 'shared', 'taxa', 'sequence', 'script',
              'sff', 'oligos', 'files', 'fna', 'qual', 'contig', 'fastq']


def memDiff():
    if enableMemDiff:
        global prevMem, memTally
        curMem = process.memory_info().rss / 1024.0
        print memTally, ":", curMem - prevMem, "KiB"
        prevMem = curMem
        memTally += 1
    # else just skip this function entirely


def fileUploadChunk(request):
    global lastEndTime
    # this function can cause significant memory leaks if chunk size is too large (http header size scales too fast)
    chunkStartTime = time.time()
    #print "Client:", chunkStartTime - lastEndTime
    global memTally
    memTally = 0
    memDiff()
    global chunkStorageDict
    retDict = {"error": "none"}
    # given an RID and a chunk of a file, write chunk to temp file under directory labelled by RID, name chunk by sequence
    # this function should also have some check for how many RID's are active from this user in this process
    # we want to limit upload files and how many uploads can occur at a time (ideally just one at a time per user)
    try:
        #myToken = request.POST['csrfmiddlewaretoken']

        if not request.user.is_authenticated():
            securityLog(request, "fileUploadChunk", "Authentication")
            raise Exception("Error: invalid user credentials")
        # TODO 1.3 user based performance usage cap

        myData = json.loads(request.body.rsplit('&', 1)[0])
        myRID = myData['RID']
        myChunkNum = myData['chunkNum']
        memDiff()   # the act of pulling myData from the request body is what takes up memory, myChunk is just a reference

        # do security checks before processing data
        if myRID in chunkStorageDict.keys():    # check that we've already started on this file
            # lets use a new dictionary for chunk tally, RID:dict where dict is chunkNum:chunkTally
            if myChunkNum in chunkTallyDict[myRID].keys():
                # this chunk was received before, check tally dict
                if chunkTallyDict[myRID][myChunkNum] >= chunkRepeatLimit:
                    securityLog(request, "fileUploadChunk", "TooManyDuplicates")
                    raise Exception("Error: Too many re-sent chunks")
        # need to limit how many uploads a user can have at a time (ideally one) TODO 1.3
        myChunk = myData['chunk']
        # need to hash the chunk and return the hash
        # at this point, myChunk is a Uint8 array made into a dictionary, convert back to array to make hash, then
        # send it back to binary entirely for the chunk saving / file writing later
        chunkData = []
        myKeys = myChunk.keys()
        keyList = []
        for k in myKeys:
            keyList.append(int(k))
        keyList = sorted(keyList)
        for index in keyList:
            chunkData.append(myChunk[str(index)])  # have to string case since its not actually an index
        myBuffer = bytearray(chunkData)
        # verify this chunk fits size specification
        if len(myBuffer) > chunkSize:
            raise Exception("Error: Received too much data")

        # now that we have the buffer, we can wipe myChunk and chunkData from memory
        memDiff()
        del chunkData
        del myChunk
        myData['chunk'].clear()
        memDiff()
        myHash = hashlib.md5(myBuffer).hexdigest()
        retDict['hash'] = myHash
        # save in memory up to a certain amount, then write that chunk down as a temp file. Cleanup during uploadComplete
        if myRID not in chunkStorageDict.keys():
            chunkStorageDict[myRID] = {}
            chunkTallyDict[myRID] = {}
            bunchCountDict[myRID] = 0
            bunchedChunkNumDict[myRID] = 0
            debug(request.user.username, "started an upload!")
            log(request, "fileUploadChunk", "NewFileUpload")
        myStorage = chunkStorageDict[myRID]
        myTallies = chunkTallyDict[myRID]

        # TODO 1.3 tweak this to support async calls so we can send many chunks at a time instead of drip-feeding 10KB for many GB files
        # TODO 1.3 atm more concerned with the file corruption seen on the archive test, missing exactly one sub file ?!
        # only broken on archive file, changing bunch and chunk sizes still makes correct files of non-archive types
        # TODO CURRENT this needs to ignore or report when chunkNum is less than lastBunched, since we already wrote that section down
        # bunch writing section should look for consecutive chunk writing available via sorted keys list
        # once we can write a section going from lastBunch of at least maxChunkCount, do such
        # this section should be put behind a lock (the if to get into it even?), as I'm not confident enough in the GIL

        # TODO 1.3 new memory leak issue (likely security risk as well) if we don't error out from missing chunks
        # like if the user skips chunk 1 and just keeps adding new ones, we'd never write a bunch file
        # need to detect and resolve this case, quite important
        #print "Tally!"
        if myChunkNum in myTallies.keys():
            # this chunk is a repeat
            myTallies[myChunkNum] += 1
        else:
            # first time we encountered this chunk
            myTallies[myChunkNum] = 0
        # either way we need to put myBuffer in this spot of myStorage
        myStorage[myChunkNum] = myBuffer

        # now how to we handle bunches? get sorted key list and see how far we can go before a gap
        # also these memory objects are shared, and this can get called multiple times, so we need to make sure only one
        # thread is allowed into the bunch-writing and dictionary-clearing section (critical section)
        #print "Hey wait we made it to bunch"
        # where do we put the lock for this section?
        ''' New bunch code '''
        # get number of chunks currently in memory
        curChunkMem = len(myStorage.keys())
        if curChunkMem >= maximumMemorizedChunks:
            # acquire lock to proceed
            #print "Capacity", myChunkNum
            if bunchLock.acquire(False):   # TODO 1.3 this lock should be unique to RID, atm its global
                # we have enough chunks to at least check for continuous chunks
                #print "LOCK", myChunkNum
                chunkKeys = myStorage.keys()
                newBunchStart = bunchedChunkNumDict[myRID]
                curChunkNum = newBunchStart
                while curChunkNum in chunkKeys:
                    # we have the data, this loop will exit when we run out of consecutive chunks
                    curChunkNum += 1
                #print "Got the lock, curChunkMem:", curChunkMem, ":", newBunchStart, "->", curChunkNum, "/", maximumMemorizedChunks
                # now newBunchStart -> curChunkNum (not including end) is the range of continuous chunks
                if curChunkNum-newBunchStart >= maximumMemorizedChunks:
                    # we have enough chunks to write to file, this is definitely a critical section
                    bunchPath = "temp/bunches/" + myRID
                    bunchNumber = bunchCountDict[myRID]
                    if not os.path.exists(bunchPath):
                        try:
                            os.makedirs(bunchPath)
                        except Exception as errr:
                            raise Exception("Error with making new bunch directory:" + str(errr))
                    # need to get correct chunk subset into sorted order
                    # could do writing directly from memory, we know how to iterate across the chunks correctly
                    # main thing is to delete each chunk after the bunch file is closed, so two loops
                    bunchFilePath = bunchPath + "/" + str(bunchNumber)
                    with open(bunchFilePath, "wb") as bunchFile:
                        # this loop iterator is the main change
                        iterVal = newBunchStart
                        while iterVal < curChunkNum:
                            bunchFile.write(myStorage[iterVal])
                            iterVal += 1
                        bunchFile.close()
                    # reset loop and cleanup memory chunks
                    iterVal = newBunchStart
                    while iterVal < curChunkNum:
                        myStorage.pop(iterVal)
                        iterVal += 1
                    bunchedChunkNumDict[myRID] = curChunkNum
                    bunchCountDict[myRID] += 1
                # release the lock at the end of the critical section (regardless of result)
                #print "Releasing", myChunkNum
                bunchLock.release()
            # else case from acquire means lock is busy, just return to client

        chunkStorageDict[myRID] = myStorage
        memDiff()

    except Exception as e:
        errorLog(request, "fileUploadChunk", e)
        retDict['error'] = str(e)
    res = json.dumps(retDict)

    lastEndTime = time.time()
    #print "Server:", lastEndTime-chunkStartTime

    return HttpResponse(res, content_type='json')


def writeBunchFile(chunks, path, number):
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except Exception as errr:
            print "Error with making new bunch directory:", errr
            return False
    bunchFilePath = path + "/" + str(number)
    with open(bunchFilePath, "wb") as bunchFile:
        for chunk in chunks:
            bunchFile.write(chunk)
        bunchFile.close()
        return True
    return False


def fileUploadComplete(request):

    retDict = {"error": "none"}

    try:

        if not request.user.is_authenticated():
            securityLog(request, "fileUploadComplete", "Authentication")
            raise Exception("Error: invalid user credentials")

        # myToken = request.POST['csrfmiddlewaretoken'] TODO 1.3 are we supposed to check these or does Django on its own?
        myData = json.loads(request.body.split('&')[0])
        myType = myData['type']
        myName = myData['name']
        myRID = myData['RID']
        log(request, "fileUploadComplete", myName)
        bunchDir = "temp/bunches/"+myRID

        uploadDest = str(os.getcwd()) + "/user_uploads/" + str(request.user.username) + "/" + datetime.date.today().isoformat()

        if not os.path.exists(uploadDest):
            os.makedirs(uploadDest)

        if myType in subDirList:
            subDest = os.path.join(uploadDest, myType)
            if not os.path.exists(subDest):
                os.makedirs(subDest)
            if len(myName.split('/')) != 1 or len(myName.split('\\')) != 1 or len(myName.split('..')) != 1:
                securityLog(request, "fileUploadComplete", "Traversal attempt with filename: "+str(myName))
                raise Exception("Invalid file name: "+str(myName))
            myFullPath = subDest + "/" + myName

            with open(myFullPath, "wb") as myFile:
                if os.path.exists(bunchDir):
                    # we have actually saved bunch data for this upload, let's handle that first
                    bunchFileNames = [f for f in os.listdir(bunchDir) if os.path.isfile(os.path.join(bunchDir, f))]
                    bunchNums = []
                    for bunchName in bunchFileNames:
                        bunchNums.append(int(bunchName))
                    bunchNums = sorted(bunchNums)
                    for num in bunchNums:
                        # read .bunch file named by num, write its data into myFile, then close and delete bunch
                        numPath = bunchDir + "/" + str(num)
                        # should double check the size of the bunchFile
                        #assert os.stat(numPath).st_size <= maxMemory, ".bunch file size error"  # TODO 1.4 better error msg
                        # this assertion breaks courtesy of thread-safe bunch writer (bunches are sometimes slightly bigger)

                        with open(numPath, "rb") as bunchFile:
                            myFile.write(bunchFile.read())  # assert says the file will fit in memory entirely
                            # so we can just write directly from the input, so long as maxMemory is configured correctly

                            bunchFile.close()   # more of a formality to call this but we are done, so close() it

                        # delete this file
                        os.remove(numPath)

                    # done with all bunches, delete directory
                    os.rmdir(bunchDir)

                # Get currently loaded chunks from memory
                myChunkData = chunkStorageDict[myRID]
                # need to update loop to go from lastBunchEnd onward, can assume all chunks are present TODO 1.3
                startChunk = bunchedChunkNumDict[myRID]
                # just run the rest of the list in sorted order, maybe throw an error if a sequence issue shows up
                if len(myChunkData.keys()) > 0:
                    myRemKeys = []
                    for key in myChunkData.keys():
                        myRemKeys.append(int(key))
                    prevKey = myRemKeys[0]
                    myRemKeys = sorted(myRemKeys)
                    for numKey in myRemKeys:
                        if numKey < startChunk:
                            # we have a problem, a chunk was never written to bunch, so this file will be corrupted
                            raise Exception("Chunk start error")
                        if prevKey != myRemKeys[0] and prevKey+1 != numKey:
                            # we have a different problem, non-sequential sequence of chunks (file will be corrupted)
                            raise Exception("Chunk sequencing error")
                        # safe to write
                        myFile.write(myChunkData[numKey])
                        prevKey = numKey
                myFile.close()
                '''
                for chunk in myChunkData:
                    myFile.write(chunk)
                myFile.close()
                '''
                debug(request.user.username, "finished a file upload:", myFullPath)

    except Exception as e:
        debug("fileUploadComplete:", request.user.username, "encountered error", e)
        errorLog(request, "fileUploadComplete", e)
        retDict['error'] = e
    res = json.dumps(retDict)
    return HttpResponse(res, content_type='application/json')
