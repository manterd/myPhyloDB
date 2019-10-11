from database.views import files
import json
from django.http import HttpResponse
import hashlib
import sys

chunkStorageDict = {}
# chunkStorageDict uses key RID with value being a list, since we know which order we should be getting these in per RID

def fileUploadChunk(request):
    #print "fileUploadChunk"
    global chunkStorageDict
    retDict = {"error": "none"}
    # given an RID and a chunk of a file, write chunk to temp file under directory labelled by RID, name chunk by sequence
    # this function should also have some check for how many RID's are active from this user in this process
    # we want to limit upload files and how many uploads can occur at a time (ideally just one at a time per user)
    try:
        #myToken = request.POST['csrfmiddlewaretoken']
        # TODO 1.3 put user check here, don't even start this if the user isn't logged in or the token is not valid
        #print "Trying"
        myData = json.loads(request.body.rsplit('&', 1)[0])
        myRID = myData['RID']
        myChunkNum = myData['chunkNum']
        # TODO 1.3 need to verify chunk is at or below expected size, also some sort of frequency control (security!)
        # TODO 1.3 file size control, write to file when enough chunks are in memory, reset memory. Rebuild final file from temp files + current memory
        myChunk = myData['chunk']   # in theory knowing if this is a string type file or not could fix reassembly
        #myChunk = request.POST['chunkBinaryData']
        # need to hash the chunk and return the hash
        #print "Encoding and hashing", myChunk
        # at this point, myChunk is an Int32 array made into a dictionary, how to convert back to binary (or hex should work too)
        chunkData = []
        myKeys = myChunk.keys()
        keyList = []
        for k in myKeys:
            keyList.append(int(k))
        keyList = sorted(keyList)
        for index in keyList:
            chunkData.append(myChunk[str(index)])  # have to string case since its not actuall an index


        myBuffer = bytearray(chunkData)
        #print "Buffer:", myBuffer
        myHash = hashlib.md5(myBuffer).hexdigest()
        #print "Hash:", myHash
        retDict['hash'] = myHash
        # now actually store the data somewhere, either write to temp file or hoard the whole thing in memory (might be huge)
        # save in memory up to a certain amount, then write that chunk down as a temp file. Cleanup during uploadComplete
        if myRID not in chunkStorageDict.keys():
            chunkStorageDict[myRID] = []
            print "New file being uploaded"
        myStorage = chunkStorageDict[myRID]

        # verify this chunk fits size specification  # TODO 1.3 implement

        #myChunk = myChunk.decode('utf-8')
        #myChunk = bytes(myChunk)
        #myChunk = myChunk.encode('utf-8')
        #myChunk = "".join(myChunk).decode('utf-8')

        #myChunk = myChunk.encode("ASCII")

        #print myChunk

        if myChunkNum < len(myStorage):
            # we have a duplicate chunk send, assume the hash mismatched from before and overwrite
            # BUT we also should keep track of how many times a duplicate runs through, limit 10
            print "DUPLICATE CHUNK"  # TODO 1.3 implement tracker
            myStorage[myChunkNum] = myBuffer
        else:
            # normal case where this chunk is new
            #print "Saving chunk"
            #print myChunk
            myStorage.append(myBuffer)
            #print "saved chunk"

        # TODO 1.3 check myStorage is below max bunch size (maximum number of chunks before writing to disk)
        # if above max, write bunch down as a temp file labelled RID-A.bunch or somesuch
        #print "placing into storage"
        chunkStorageDict[myRID] = myStorage

    except Exception as e:
        print e
        retDict['error'] = e
    res = json.dumps(retDict)
    #print "returning", HttpResponse(res, content_type='json')
    return HttpResponse(res, content_type='json')


def fileUploadComplete(request):
    print "Called fileUploadComplete"
    # myToken = request.POST['csrfmiddlewaretoken']
    myData = json.loads(request.body.split('&')[0])
    myType = myData['type']
    myName = myData['name']
    myRID = myData['RID']
    retDict = {"error": "none"}
    try:
        # This function needs to take an RID and find the temp folder holding chunks of a file, put them all back together
        # Use the old fileUpFunc code to get the right location to store the final file
        # TODO 1.3 Need to check for .bunch files with RID at start of name
        # Get currently loaded chunks from memory
        myChunkData = chunkStorageDict[myRID]
        myFileData = ""
        # TODO 1.3 holy cow we need to run the filename through some checks before doing this
        with open(myName, "w+b") as myFile:
            for chunk in myChunkData:
                myFileData += chunk
                myFile.write(chunk)
            myFile.close()
        print "Got file named", myName, ": it is", len(myFileData), "long (units?)"
        # so we're writing each byte as a character, we have the right data, just need to write it down AS BYTES, for some reason the string typing adds data
    except Exception as e:
        retDict['error'] = e
    res = json.dumps(retDict)
    return HttpResponse(res, content_type='application/json')

    '''
    errorText = ""
    RID = ''    # TODO 1.4 stopList and RID NYI, can go in with progress bar (same recoding basically)
    # There exists a django form plugin which can accomplish the progress tracking issue (as well as allow bigger file uploads overall)
    # TODO 1.3 Progress bar for uploading files. This part is hard to test locally, local uploads are somewhat fast
    # TODO 1.3 progress bar should come along with chunk based uploader
    # RID not particularly necessary for a one-way operation such as this. May be handy for logging of some description
    # Stoplist is the main reason to use RID, assuming stop function is implemented on the front end
    # Stopping a file upload is awkward since this is coded with only a single line technically performing the download
    # Need to interrupt the upload somehow if stop will be relevant
    # Further, progress bar AND stop both involve having a tighter control over this upload process, not sure if Django
    # supports accessing upload percentages and executing commands on operations I believe it assumes are done by now

    errorText = ""
    try:
        RID = request.POST['RID']
    except:
        pass
    date = datetime.date.today().isoformat()
    username = str(request.user.username)
    fileForm = BulkForm(request.POST, request.FILES)    # TODO 1.3 form overhaul for chunk uploader
    if fileForm.is_valid():

        try:

            uploadDest = str(os.getcwd()) + "/user_uploads/" + str(username) + "/" + date

            if not os.path.exists(uploadDest):
                os.makedirs(uploadDest)

            subFilesDict = {}
            for subDir in subDirList:
                subFilesDict[subDir] = request.FILES.getlist(str(subDir) + "files")
                subDest = os.path.join(uploadDest, subDir)
                if len(subFilesDict[subDir]) > 0:
                    if not os.path.exists(subDest):
                        os.makedirs(subDest)
                    try:
                        for each_file in subFilesDict[subDir]:
                            functions.handle_uploaded_file(each_file, subDest, each_file.name)
                    except Exception as upfileError:
                        errorText = "Error: " + str(upfileError)
        except Exception as genericUploadError:
            errorText = "Error: " + str(genericUploadError)
    else:
        errorText = "Error: invalid form"

    return files(request, errorText)
    '''

