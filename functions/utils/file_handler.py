import json
from django.http import HttpResponse
import hashlib
import os
import math
import datetime
from database.views import subDirList
from functions.utils.debug import debug

chunkStorageDict = {}
# chunkStorageDict uses key RID with value being a list, since we know which order we should be getting these in per RID
chunkSize = 10000  # 10 KB chunk max, make sure this is synced with files.html value
maxMemory = 10000000  # 10 MB max before offloading to temp file
maximumMemorizedChunks = maxMemory / chunkSize

# TODO 1.3 need a process that cleans up requests which haven't finished after a period of time with no updating


def fileUploadChunk(request):
    global chunkStorageDict
    retDict = {"error": "none"}
    # given an RID and a chunk of a file, write chunk to temp file under directory labelled by RID, name chunk by sequence
    # this function should also have some check for how many RID's are active from this user in this process
    # we want to limit upload files and how many uploads can occur at a time (ideally just one at a time per user)
    try:
        #myToken = request.POST['csrfmiddlewaretoken']

        if not request.user.is_authenticated():
            # TODO 1.3 security flag
            raise Exception("Error: invalid user credentials")

        myData = json.loads(request.body.rsplit('&', 1)[0])
        myRID = myData['RID']
        myChunkNum = myData['chunkNum']
        # TODO 1.3 Frequency control (security!)
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
        myHash = hashlib.md5(myBuffer).hexdigest()
        retDict['hash'] = myHash
        # now actually store the data somewhere, either write to temp file or hoard the whole thing in memory (might be huge)
        # save in memory up to a certain amount, then write that chunk down as a temp file. Cleanup during uploadComplete
        if myRID not in chunkStorageDict.keys():
            chunkStorageDict[myRID] = []
            debug(request.user.username, "started an upload!")
        myStorage = chunkStorageDict[myRID]

        # verify this chunk fits size specification  # TODO 1.3 implement
        if len(myBuffer) > chunkSize:
            raise Exception("Error: Received too much data")

        if myChunkNum % maximumMemorizedChunks < len(myStorage):  # by the time we've reached max, we can assume reset
            # we have a duplicate chunk send, assume the hash mismatched from before and overwrite
            # BUT we also should keep track of how many times a duplicate runs through, limit 10
            print "DUPLICATE CHUNK"  # TODO 1.3 implement tracker, dupes suggest hash mismatch or weird DoS attack
            myStorage[myChunkNum] = myBuffer
        else:
            # normal case where this chunk is new
            myStorage.append(myBuffer)

            # memory to temp file handler
            if len(myStorage) == maximumMemorizedChunks:
                # bunch starts sounding like a weird word the more I see it here, refactor? TODO 1.3 but really though
                # if above max, write bunch down as a temp file labelled RID-A.bunch or somesuch
                bunchPath = "temp/bunches/"+myRID
                bunchNumber = int(math.floor(myChunkNum / maximumMemorizedChunks))  # how many times have we done this before for this RID?
                if not os.path.exists(bunchPath):
                    try:
                        os.makedirs(bunchPath)
                    except Exception as errr:
                        raise Exception("Error with making new bunch directory:"+str(errr))
                if writeBunchFile(myStorage, bunchPath, bunchNumber):
                    # file wrote successfully
                    # empty old list from memory after writing file
                    del myStorage[:]
                    myStorage = []
                else:
                    raise Exception("Error with bunch writer")

        chunkStorageDict[myRID] = myStorage

    except Exception as e:
        print e  # TODO 1.3 log this
        retDict['error'] = e
    res = json.dumps(retDict)
    return HttpResponse(res, content_type='json')


def writeBunchFile(chunks, path, number):
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except Exception as errr:
            print "Error with making new bunch directory:", errr
            return False
    bunchFilePath = path + "/" + str(number)
    with open(bunchFilePath, "w+b") as bunchFile:
        for chunk in chunks:
            bunchFile.write(chunk)
        bunchFile.close()
        return True
    return False


def fileUploadComplete(request):
    # myToken = request.POST['csrfmiddlewaretoken'] TODO 1.3 are we supposed to check these or does Django on its own?
    myData = json.loads(request.body.split('&')[0])
    myType = myData['type']
    myName = myData['name']
    myRID = myData['RID']
    retDict = {"error": "none"}

    try:

        if not request.user.is_authenticated():
            # TODO 1.3 security flag
            raise Exception("Error: invalid user credentials")

        bunchDir = "temp/bunches/"+myRID

        uploadDest = str(os.getcwd()) + "/user_uploads/" + str(request.user.username) + "/" + datetime.date.today().isoformat()

        if not os.path.exists(uploadDest):
            os.makedirs(uploadDest)

        if myType in subDirList:
            subDest = os.path.join(uploadDest, myType)
            if not os.path.exists(subDest):
                os.makedirs(subDest)
            if len(myName.split('/')) != 1 or len(myName.split('\\')) != 1 or len(myName.split('..')) != 1:
                # TODO 1.3 security flag
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
                        assert os.stat(numPath).st_size <= maxMemory, ".bunch file size error"  # TODO 1.3 better error msg

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
                myFileData = ""
                for chunk in myChunkData:
                    myFileData += chunk
                    myFile.write(chunk)
                myFile.close()
                debug(request.user.username, "finished a file upload:", myFullPath)

    except Exception as e:
        debug("fileUploadComplete:", request.user.username, "encountered error", e)
        retDict['error'] = e
    res = json.dumps(retDict)
    return HttpResponse(res, content_type='application/json')
