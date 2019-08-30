from collections import defaultdict
import ctypes
import datetime
from django import forms
from django.core.exceptions import ValidationError
from django.http import HttpResponse
import inspect
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import re
import shutil
import json
import threading
import time
import zipfile
import math
from natsort import natsorted
from PyPDF2 import PdfFileReader, PdfFileMerger

from database.models import Project, Reference, Profile

import config.local_cfg
import functions
from functions.utils.debug import debug
from database import perms

from Queue import Queue
import psutil
import pickle

process = psutil.Process(os.getpid())


pd.set_option('display.max_colwidth', -1)

# logging code
loggerRunning = 0
preLogBackLog = Queue(maxsize=0)    # queue containing entries to write to log file
# second log for console page, should store most recent N entries, so make a queue and pop it when we reach N
consoleLog = []  # using array for full access without popping
# can't use queue because we want to keep the entries after reading them, will have to maintain order manually
oldestLogEntry = 0  # integer marking the index in consoleLog containing the first entry
# when we hit N entries, the next addition will overwrite oldestLogEntry and increment it by one
# when oldestLogEntry reaches N, reset it to 0
# this ensures we have a list of N log entries we can always read in correct order
maximumLogEntryCount = 50   # change this to influence logging history size (N)

server_start_time = -1


def log(request, reqType, name):
    global preLogBackLog
    text = ""
    text += str(datetime.datetime.now())
    typeLength = str(reqType).__len__()
    while typeLength < 8:
        reqType += " "
        typeLength += 1
    text += "  User: " + str(request.user.username) + "  \tType: " + str(reqType) + " \tName: " + str(name)
    # formatting subject to change
    preLogBackLog.put(text, True)
    # actual console output
    print text
    # add to secondary log for console page, this stores the N most recent log functions
    addToConsoleLog(text)


def addToConsoleLog(text):
    # check if list is full, if not just add. If full, overwrite oldest entry and move iterator for oldest by one
    global consoleLog, oldestLogEntry
    if len(consoleLog) >= maximumLogEntryCount:
        consoleLog[oldestLogEntry] = text
        oldestLogEntry += 1
        if oldestLogEntry >= maximumLogEntryCount:
            oldestLogEntry = 0
    else:
        consoleLog.append(text)


def getConsoleLog(request):
    if request.is_ajax():
        if not request.user.is_superuser or not request.user.is_authenticated:
            output = json.dumps("Invalid Permissions")
            return HttpResponse(output, content_type='application/json')
        # get console log for page, use oldestLogEntry to track the oldest value, loop if passing maximum count
        output = ""
        try:
            for iter in range(0, min(maximumLogEntryCount, len(consoleLog))):
                if oldestLogEntry + iter < maximumLogEntryCount:
                    output += str(consoleLog[oldestLogEntry+iter]) + "\n"
                else:
                    output += str(consoleLog[oldestLogEntry+iter-maximumLogEntryCount]) + "\n"
        except Exception as e:
            print "Exception during console log:", e
        res = json.dumps(output)
        return HttpResponse(res, content_type='application/json')


def setServerStartTime():
    global server_start_time
    server_start_time = time.time()


def getServerMetrics(request):
    # at present this only gets server run time, can add whatever is needed here and it'll get displayed on the page
    if request.is_ajax():
        if not request.user.is_superuser or not request.user.is_authenticated:
            output = json.dumps("Invalid Permissions")
            return HttpResponse(output, content_type='application/json')
        days = math.trunc((time.time() - server_start_time) / 86400)
        hours = math.trunc(((time.time() - server_start_time) % 86400) / 3600)
        minutes = math.trunc((((time.time() - server_start_time) % 86400) % 3600) / 60)
        seconds = math.trunc((((time.time() - server_start_time) % 86400) % 3600) % 60)
        output = "Server has been running for: " + str(days) + " days, " + str(hours) + " hours, " + str(minutes) + \
                 " minutes, " + str(seconds) + " seconds"
        res = json.dumps(output)
        return HttpResponse(res, content_type='application/json')


def startLogger():
    global loggerRunning, preLogBackLog
    if loggerRunning == 0:
        loggerRunning = 1
        while True:     # needs an exit in case server shutdown is called
            try:
                newText = preLogBackLog.get(block=True, timeout=None)
                logFile = open('server_log.txt', 'a')
                logFile.write(newText+"\n")
                logFile.close()
            except Exception as e:
                print "Error while logging: ", e
    else:
        print "Error: logger already active"


def ordered_set(seq, idfun=None):
    if idfun is None:
        def idfun(x):
            return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen:
            continue
        seen[marker] = 1
        result.append(item)
    return result


def handle_uploaded_file(f, path, name):  # move file from memory to disc
    try:
        if not os.path.exists(path):
            os.makedirs(path)
        dest = "/".join([str(path), str(name)])
        with open(str(dest), 'wb+') as destination:
            for chunk in f.chunks():
                destination.write(chunk)
    except Exception as e:
        print "Error with handling file:", e
        print "Was moving", f, "to", path, "with name", name


def mergePDF(path, finalFile):

    # Combining Pdf files from path into single finalFile
    pdf_files = [f for f in os.listdir(path) if f.endswith("pdf")]
    pdf_files = natsorted(pdf_files, key=lambda y: y.lower())
    merger = PdfFileMerger()
    for filename in pdf_files:
        merger.append(PdfFileReader(file(os.path.join(path, filename), 'rb')))
    merger.write(finalFile)


def remove_list(refList, user):
    #p_uuid = list(Reference.objects.filter(path__in=refList).values_list('projectid', flat=True))
    # Remove record from reference table
    for path in refList:
        debug("Remove_list:", path, user.username)
        ret = remove_proj(path, user)
        if ret == 1:
            debug("Encountered an issue while removing project")
            return 1
    return 0


def remove_proj(path, user):  # remove selected project path, after verifying user has correct permissions
    # so a permissions check is not entirely necessary here, though if its called from another context, maybe add perms
    # remove project, files and database
    # TODO 1.4 security alerts: some way to push security-related errors to current admins, no console babysitting required
    # TODO 1.3 for now, send security errors to a secondary log, which displays separately to the normal console
    pid = Reference.objects.get(path=path).projectid.projectid
    # Verify user has correct permissions for this step (either the owner or a superuser, no exceptions for deletion)
    myProj = Project.objects.get(projectid=pid)
    debug("RemProj:", myProj.owner.username, user.username, user.is_superuser)    # TODO 1.3 error:Exception page
    debug("ProjDetails:", pid, myProj.project_name, myProj)
    if myProj.owner == user or user.is_superuser:
        # remove current reference object, if its the last/only one, remove the project as a whole
        # actually remove project from database
        debug("Removing full project")
        # update perms before deleting project
        Reference.objects.get(path=path).delete()
        if myProj.status == "private":
            perms.remPriv(pid)
        else:
            perms.remPub(pid)
        myProj.delete()

        path = "/".join(["uploads", str(pid)])
        if os.path.exists(path):
            shutil.rmtree(path)

        return 0
    else:
        print "SECURITY: User lacks permission to remove this project!", Project.objects.get(projectid=pid).project_name, user
        # TODO 1.3 put this into security alerts feed (another console stat?)
        return 1

def multidict(ordered_pairs):
    d = defaultdict(list)
    for k, v in ordered_pairs:
        d[k].append(v)

    for k, v in d.items():
        if len(v) == 1:
            d[k] = v[0]
    return dict(d)


def taxaProfileDF(mySet):
    qs1 = Profile.objects.filter(sampleid__in=mySet).values('sampleid', 'kingdomid', 'phylaid', 'classid', 'orderid', 'familyid', 'genusid', 'speciesid', 'otuid', 'count')
    df = pd.DataFrame.from_records(qs1, columns=['sampleid', 'kingdomid', 'phylaid', 'classid', 'orderid', 'familyid', 'genusid', 'speciesid', 'otuid', 'count'])
    df = df.groupby(['sampleid', 'kingdomid', 'phylaid', 'classid', 'orderid', 'familyid', 'genusid', 'speciesid', 'otuid'])['count'].sum()
    df = df.to_frame(name='count')
    df2 = df.unstack(['sampleid']).fillna(0).stack(['sampleid'])
    df3 = df2.unstack(['sampleid'])
    taxaDF = df3['count']
    return taxaDF


def PCoA(dm):   # TODO 1.3 is this used anywhere? It appears not, can probably remove. Ask
    E_matrix = make_E_matrix(dm)
    F_matrix = make_F_matrix(E_matrix)  # pycharm complains about this function call because it cannot be 100% sure
    # about the variable type of E_matrix
    eigvals, eigvecs = np.linalg.eigh(F_matrix)
    negative_close_to_zero = np.isclose(eigvals, 0)  # another duck-typing related warning, in practice this is fine
    eigvals[negative_close_to_zero] = 0
    idxs_descending = eigvals.argsort()[::-1]
    eigvals = eigvals[idxs_descending]
    eigvecs = eigvecs[:, idxs_descending]
    eigvals, coordinates, proportion_explained = scores(eigvals, eigvecs)
    return eigvals, coordinates, proportion_explained


def make_E_matrix(dist_matrix):
    return dist_matrix * dist_matrix / -2.0


def make_F_matrix(E_matrix):
    row_means = E_matrix.mean(axis=1, keepdims=True, dtype=np.float64)
    col_means = E_matrix.mean(axis=0, keepdims=True, dtype=np.float64)
    matrix_mean = E_matrix.mean(dtype=np.float64)
    return E_matrix - row_means - col_means + matrix_mean


def scores(eigvals, eigvecs):
    num_positive = (eigvals >= 0).sum()  # TODO 1.3 fix
    eigvecs[:, num_positive:] = np.zeros(eigvecs[:, num_positive:].shape)
    eigvals[num_positive:] = np.zeros(eigvals[num_positive:].shape)
    coordinates = eigvecs * np.core.umath.sqrt(eigvals)
    proportion_explained = eigvals / eigvals.sum()
    return eigvals, coordinates, proportion_explained


def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            os.remove(os.path.join(dir, f))


class MultiFileInput(forms.FileInput):
    def render(self, name, value, attrs={}):    # TODO 1.3 fix
        attrs['multiple'] = 'multiple'
        return super(MultiFileInput, self).render(name, None, attrs=attrs)

    def value_from_datadict(self, data, files, name):
        if hasattr(files, 'getlist'):
            return files.getlist(name)
        else:
            return [files.get(name)]


class MultiFileField(forms.FileField):
    widget = MultiFileInput
    default_error_messages = {
        'min_num': u"Ensure at least %(min_num)s files are uploaded (received %(num_files)s).",
        'max_num': u"Ensure at most %(max_num)s files are uploaded (received %(num_files)s).",
        'file_size': u"File: %(uploaded_file_name)s, exceeded maximum upload size."
    }

    def __init__(self, *args, **kwargs):
        self.min_num = kwargs.pop('min_num', 0)
        self.max_num = kwargs.pop('max_num', None)
        self.maximum_file_size = kwargs.pop('maximum_file_size', None)
        super(MultiFileField, self).__init__(*args, **kwargs)

    def to_python(self, data):
        ret = []
        for item in data:
            ret.append(super(MultiFileField, self).to_python(item))
        return ret

    def validate(self, data):
        super(MultiFileField, self).validate(data)
        num_files = len(data)
        if len(data) and not data[0]:
            num_files = 0
        if num_files < self.min_num:
            raise ValidationError(self.error_messages['min_num'] % {'min_num': self.min_num, 'num_files': num_files})
        elif self.max_num and num_files > self.max_num:
            raise ValidationError(self.error_messages['max_num'] % {'max_num': self.max_num, 'num_files': num_files})
        for uploaded_file in data:
            if uploaded_file.size > self.maximum_file_size:
                raise ValidationError(self.error_messages['file_size'] % {'uploaded_file_name': uploaded_file.name})


def _async_raise(tid, exctype):
    if not inspect.isclass(exctype):
        raise TypeError("Only types can be raised (not instances)")
    res = ctypes.pythonapi.PyThreadState_SetAsyncExc(ctypes.c_long(tid), ctypes.py_object(exctype))

    if res == 0:
        raise ValueError("invalid thread id")
    elif res != 1:
        ctypes.pythonapi.PyThreadState_SetAsyncExc(ctypes.c_long(tid), 0)
        raise SystemError("PyThreadState_SetAsyncExc failed")


class stoppableThread(threading.Thread):
    def _get_my_tid(self):
        if not self.isAlive():
            raise threading.ThreadError("the thread is not active")

        if hasattr(self, "_thread_id"):
            return self._thread_id

        for tid, tobj in threading._active.items():
            if tobj is self:
                self._thread_id = tid
                return tid

        raise AssertionError("could not determine the thread's id")

    def raise_exc(self, exctype):
        _async_raise(self._get_my_tid(), exctype)

    def terminate(self):
        self.raise_exc(SystemExit)


def remove(path):
    if os.path.isdir(path):
        try:
            os.rmdir(path)
        except OSError:
            print "Unable to remove folder: %s" % path
    else:
        try:
            if os.path.exists(path):
                os.remove(path)
        except OSError:
            print "Unable to remove file: %s" % path


def wOdum(data, alpha):
    length = data.__len__()
    numrows, numcols = np.shape(data)
    dists = np.zeros((numrows, numrows))

    for i in xrange(length):
        for j in xrange(length):
            dist = 0.0
            if i == j:
                dists[i, j] = dist
                dists[j, i] = dist
            else:
                num, den = 0.0, 0.0
                otus = data[i].__len__()
                for l in xrange(otus):
                    u = data[i, l]
                    v = data[j, l]
                    if u + v > 0:
                        num += abs(u-v)/(u+v)*(u+v)**alpha
                        den += (u+v)**alpha
                dist = num / den
                dists[i, j] = dist
                dists[j, i] = dist
    return dists


def dictSum(*dicts):
    ret = defaultdict(int)
    for d in dicts:
        for k, v in d.items():
            ret[k] += v
    return dict(ret)


def cleanup(path):
    age = 3 * 60 * 60 * 24

    if os.path.exists(path):
        for file in os.listdir(path):
            now = time.time()
            filepath = os.path.join(path, file)
            modified = os.stat(filepath).st_mtime
            if modified < now - age:
                if os.path.isfile(filepath):
                    os.remove(filepath)


def analysisThreads():
    try:
        usr_threads = int(config.local_cfg.usr_num_threads)
    except Exception:
        usr_threads = 1

    if usr_threads <= 0:
        num_threads = 1
    elif usr_threads > mp.cpu_count():
        num_threads = int(mp.cpu_count())
    else:
        num_threads = usr_threads
    return int(num_threads)


def getRawDataTab(request):
    debug("getRawDataTab")
    if request.is_ajax():
        RID = request.GET["all"]
        func = request.GET["func"]
        treeType = int(request.GET["treeType"])
        myDir = 'myPhyloDB/media/temp/analysis/'
        path = str(myDir) + str(RID) + '.biom'

        debug("getRawDataTab: savedDF")
        savedDF = exploding_panda2(path, treeType)
        debug("getRawDataTab: panda")

        myDir = 'myPhyloDB/media/temp/' + str(func) + '/'
        fileName = str(myDir) + str(request.user.username) + '.' + str(func) + '.csv'

        debug("getRawDataTab: to_csv")
        savedDF.to_csv(fileName)

        myDir = 'myPhyloDB/media/temp/' + str(func) + '/'
        fileName2 = str(myDir) + str(request.user.username) + '.' + str(func) + '.gz'
        zf = zipfile.ZipFile(fileName2, "w", zipfile.ZIP_DEFLATED, allowZip64=True)

        debug("getRawDataTab: write csv")
        zf.write(fileName, 'usr_data.csv')
        zf.close()

        myDir = '../../myPhyloDB/media/temp/' + str(func) + '/'
        fileName2 = str(myDir) + str(request.user.username) + '.' + str(func) + '.gz'
        myDict = {'name': str(fileName2)}

        debug("getRawDataTab: return")
        res = json.dumps(myDict)
        return HttpResponse(res, content_type='application/json')


def getRawDataBiom(request):
    if request.is_ajax():
        RID = request.GET["all"]
        func = request.GET["func"]

        myDir = '../../myPhyloDB/media/temp/' + str(func) + '/'
        fileName = str(myDir) + str(RID) + '.biom'

        myDict = {'name': str(fileName)}
        res = json.dumps(myDict)
        return HttpResponse(res, content_type='application/json')


def getCoreBiom(request):
    if request.is_ajax():
        RID = request.GET["all"]

        myDir = '../../myPhyloDB/media/temp/core/Biom/' + str(RID) + "/"
        # all biom files for this RID are in biom/rid directory
        # rid.biom.zip is the archive containing the biom files
        fileName = str(myDir) + 'core.biom.zip'

        myDict = {'name': str(fileName)}
        res = json.dumps(myDict)
        return HttpResponse(res, content_type='application/json')


def removeFiles(request):   # DEPRECATED, see cleanup function in queue
    if request.is_ajax():
        RID = request.GET["all"]
        func = request.GET["func"]

        file = "myPhyloDB/media/temp/" + str(func) + "/Rplots/" + str(RID) + "." + str(func) + ".pdf"
        if os.path.exists(file):
            os.remove(file)

        file = "myPhyloDB/media/temp/" + str(func) + "/" + str(RID) + ".biom"
        if os.path.exists(file):
            os.remove(file)

        file = "myPhyloDB/media/temp/" + str(func) + "/" + str(RID) + ".csv"
        if os.path.exists(file):
            os.remove(file)

        return HttpResponse()


def excel_to_dict(wb, headerRow=1, nRows=1, sheet='Sheet1'):
    ws = wb[sheet]
    headerDict = dict()
    for col in xrange(1, ws.max_column+1):
        if ws.cell(row=headerRow, column=col).value:
            headerDict[col] = ws.cell(row=headerRow, column=col).value

    result_dict = []
    for row in range(headerRow+1, headerRow+1+nRows):
        line = dict()
        for key, value in headerDict.iteritems():
            cell_value = ws.cell(row=row, column=key).value
            if cell_value is None:
                cell_value = np.nan
            line[headerDict[key]] = cell_value
        result_dict.append(line)
    return result_dict


def getMetaDF(username, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, DepVar, levelDep=False):
    # we now use categorical data type to compress pandas (yes now we explode, implode, AND compress the poor fellas)
    # as well as filtering which data is necessary for exploding
    catFields = []
    catValues = []
    if metaValsCat:
        metaDictCat = json.JSONDecoder(object_pairs_hook=multidict).decode(metaValsCat)
        for key in sorted(metaDictCat):
            catFields.append(key)
            catValues.extend(metaDictCat[key])
    catSampleIDs = []
    if metaIDsCat:
        idDictCat = json.JSONDecoder(object_pairs_hook=multidict).decode(metaIDsCat)
        for key in sorted(idDictCat):
            myIDs = idDictCat[key]
            if type(myIDs) == list:
                print "List!"
                for ID in myIDs:
                    catSampleIDs.append(ID)
            else:
                print "Not list!"
                catSampleIDs.append(myIDs)
    quantFields = []
    quantValues = []
    if metaValsQuant:
        if isinstance(metaValsQuant, list):
            quantFields = metaValsQuant
        else:
            metaDictQuant = json.JSONDecoder(object_pairs_hook=multidict).decode(metaValsQuant)
            for key in sorted(metaDictQuant):
                quantFields.append(key)
                quantValues.extend(metaDictQuant[key])
    debug("getMetaDF: quantFields:", quantFields)
    quantSampleLists = []
    if metaIDsQuant:
        idDictQuant = json.JSONDecoder(object_pairs_hook=multidict).decode(metaIDsQuant)
        for key in sorted(idDictQuant):
            quantSampleLists.append(idDictQuant[key])
    quantSampleIDs = []
    if quantSampleLists:
        quantSampleIDs = list(set.intersection(*map(set, quantSampleLists)))
    if catSampleIDs and not quantSampleIDs:
        finalSampleIDs = list(set(catSampleIDs))
    elif quantSampleIDs and not catSampleIDs:
        finalSampleIDs = list(set(quantSampleIDs))
    else:
        finalSampleIDs = list(set(catSampleIDs) & set(quantSampleIDs))
    myDir = 'myPhyloDB/media/usr_temp/' + str(username) + '/'
    path = str(myDir) + 'myphylodb.biom'
    # if there are memory issues in analysis, blame the exploding panda
    debug("getMetaDF: Exploding Panda")
    savedDF, metaDF, remCatFields = exploding_panda(path, finalSampleIDs=finalSampleIDs, catFields=catFields, quantFields=quantFields, levelDep=levelDep, depVar=DepVar)
    return savedDF, metaDF, finalSampleIDs, catFields, remCatFields, quantFields, catValues, quantValues


def transformDF(transform, DepVar, finalDF):
    # replace zeros before transformation
    if transform != 0:
        if DepVar == 0:
            myList = finalDF.abund.tolist()
            nonZero = filter(lambda a: a != 0, myList)
            value = min(nonZero) / 2.0
            finalDF.abund.replace(to_replace=0, value=value, inplace=True)

        elif DepVar == 1:
            myList = finalDF.rel_abund.tolist()
            nonZero = filter(lambda a: a != 0, myList)
            value = min(nonZero) / 2.0
            finalDF.rel_abund.replace(to_replace=0, value=value, inplace=True)

        elif DepVar == 2:
            myList = finalDF.rich.tolist()
            nonZero = filter(lambda a: a != 0, myList)
            value = min(nonZero) / 2.0
            finalDF.rich.replace(to_replace=0, value=value, inplace=True)

        elif DepVar == 3:
            myList = finalDF.diversity.tolist()
            nonZero = filter(lambda a: a != 0, myList)
            value = min(nonZero) / 2.0
            finalDF.diversity.replace(to_replace=0, value=value, inplace=True)

        elif DepVar == 4:
            myList = finalDF.abund_16S.tolist()
            nonZero = filter(lambda a: a != 0, myList)
            value = min(nonZero) / 2.0
            finalDF.abund_16S.replace(to_replace=0, value=value, inplace=True)

    if transform == 1:
        if DepVar == 0:
            finalDF['abund'] = np.log(finalDF.abund)
        elif DepVar == 1:
            finalDF['rel_abund'] = np.log(finalDF.rel_abund)
        elif DepVar == 2:
            finalDF['rich'] = np.log(finalDF.rich)
        elif DepVar == 3:
            finalDF['diversity'] = np.log(finalDF.diversity)
        elif DepVar == 4:
            finalDF['abund_16S'] = np.log(finalDF.abund_16S)
    elif transform == 2:
        if DepVar == 0:
            finalDF['abund'] = np.log10(finalDF.abund)
        elif DepVar == 1:
            finalDF['rel_abund'] = np.log10(finalDF.rel_abund)
        elif DepVar == 2:
            finalDF['rich'] = np.log10(finalDF.rich)
        elif DepVar == 3:
            finalDF['diversity'] = np.log10(finalDF.diversity)
        elif DepVar == 4:
            finalDF['abund_16S'] = np.log10(finalDF.abund_16S)
    elif transform == 3:
        if DepVar == 0:
            finalDF['abund'] = np.sqrt(finalDF.abund)
        elif DepVar == 1:
            finalDF['rel_abund'] = np.sqrt(finalDF.rel_abund)
        elif DepVar == 2:
            finalDF['rich'] = np.sqrt(finalDF.rich)
        elif DepVar == 3:
            finalDF['diversity'] = np.sqrt(finalDF.diversity)
        elif DepVar == 4:
            finalDF['abund_16S'] = np.sqrt(finalDF.abund_16S)
    elif transform == 4:
        if DepVar == 0:
            finalDF['abund'] = np.log10(finalDF.abund/(1-finalDF.abund))
        elif DepVar == 1:
            finalDF['rel_abund'] = np.log10(finalDF.rel_abund/(1-finalDF.rel_abund))
        elif DepVar == 2:
            finalDF['rich'] = np.log10(finalDF.rich/(1-finalDF.rich))
        elif DepVar == 3:
            finalDF['diversity'] = np.log10(finalDF.diversity/(1-finalDF.diversity))
        elif DepVar == 4:
            finalDF['abund_16S'] = np.log10(finalDF.abund_16S/(1-finalDF.abund_16S))
    elif transform == 5:
        if DepVar == 0:
            finalDF['abund'] = np.arcsin(finalDF.abund)
        elif DepVar == 1:
            finalDF['rel_abund'] = np.arcsin(finalDF.rel_abund)
        elif DepVar == 2:
            finalDF['rich'] = np.arcsin(finalDF.rich)
        elif DepVar == 3:
            finalDF['diversity'] = np.arcsin(finalDF.diversity)
        elif DepVar == 4:
            finalDF['abund_16S'] = np.arcsin(finalDF.abund_16S)

    return finalDF


def categorize(dataFrame):
    colTotal = 0
    colCat = 0
    #print "Precat:", dataFrame.memory_usage().sum() / 1024.0 ** 2
    for col in dataFrame:
        try:
            if float(100.0 * len(dataFrame[col].unique()) / len(dataFrame[col])) < 50.0:
                dataFrame[col] = dataFrame[col].astype("category")
                if "NaN" not in dataFrame[col].cat.categories:
                    dataFrame[col] = dataFrame[col].cat.add_categories("NaN")
                colCat += 1
        except Exception as exc:
            print "Exception with categorize:", exc
            pass
        #else:
            #print ">=50% uniques:", col
        colTotal += 1
    #print "Postcat:", dataFrame.memory_usage().sum() / 1024.0 ** 2
    #print "Categorized:", colCat, "/", colTotal, "->", colCat*100.0/colTotal, "%"


prevMem = 0.0
memTally = 0
enableMemDiff = False   # set this back to true for memory analysis during exploding_panda

def memDiff():
    if enableMemDiff:
        global prevMem, memTally
        curMem = process.memory_info().rss / 1024.0 ** 2
        print memTally, ":", curMem - prevMem, "MB"
        prevMem = curMem
        memTally += 1
    # else just skip this function entirely


def exploding_panda(path, finalSampleIDs=[], catFields=[], quantFields=[], levelDep=False, depVar=10):
    # this function is the primary memory hog for analysis, and the largest chunk of normalize too
    # this uses the myphylodb.biom file normalization made
    global prevMem, memTally
    memTally = 0
    prevMem = process.memory_info().rss / 1024.0 ** 2
    # Load file
    file = open(path)
    data = json.load(file)      # assuming data is the main spike for this area, ideally be done with it ASAP (before peak), for cleanup sake
    file.close()
    debug("Exploding panda loaded json", len(data))
    memDiff()   # 139.8

    # Get metadata
    datCols = data['columns']
    datRows = data['rows']
    debug("Cols and Rows", len(datCols), len(datRows))
    # get metadata into dicionary keyed by id
    metaDict = {}
    for i in datCols:
        metaDict[str(i['id'])] = i['metadata']
    # get taxa data into dictionary
    taxaDict = {}
    for i in datRows:
        tempDict = {}
        taxon = i['metadata']['taxonomy']
        taxLength = len(taxon)
        if taxLength > 0:
            tempDict['kingdomid'] = taxon[0].split(': ')[0]
            tempDict['kingdomName'] = taxon[0].split(': ')[1]
        if taxLength > 1:
            tempDict['phylaid'] = taxon[1].split(': ')[0]
            tempDict['phylaName'] = taxon[1].split(': ')[1]
        if taxLength > 2:
            tempDict['classid'] = taxon[2].split(': ')[0]
            tempDict['className'] = taxon[2].split(': ')[1]
        if taxLength > 3:
            tempDict['orderid'] = taxon[3].split(': ')[0]
            tempDict['orderName'] = taxon[3].split(': ')[1]
        if taxLength > 4:
            tempDict['familyid'] = taxon[4].split(': ')[0]
            tempDict['familyName'] = taxon[4].split(': ')[1]
        if taxLength > 5:
            tempDict['genusid'] = taxon[5].split(': ')[0]
            tempDict['genusName'] = taxon[5].split(': ')[1]
        if taxLength > 6:
            tempDict['speciesid'] = taxon[6].split(': ')[0]
            tempDict['speciesName'] = taxon[6].split(': ')[1]
        if taxLength > 7:
            tempDict['otuid'] = taxon[7].split(': ')[0]
            tempDict['otuName'] = taxon[7].split(': ')[1]
        taxaDict[str(i['id'])] = tempDict
    # Get count data and calculate various dependent variables
    sampleids = [col['id'] for col in datCols]
    taxaids = [col['id'] for col in datRows]
    memDiff()   # 70.38
    mat = np.asarray(data['data']).T.tolist()
    memDiff()   # 72.13

    # cleanup memory
    del(datCols[:])
    del(datRows[:])
    data.clear()

    memDiff()   # -0.996

    metaDF = pd.DataFrame.from_dict(metaDict, orient='index').reset_index()
    metaDF.rename(columns={'index': 'sampleid'}, inplace=True)
    metaDF.set_index('sampleid', inplace=True)
    # this /\ is duplicating the index somehow, and verify_integrity doesn't fix the issue

    memDiff()   # 0.465
    categorize(metaDF)
    memDiff()   # 0.438

    projectIDs = metaDF['projectid'].unique()
    projects = Project.objects.filter(projectid__in=projectIDs)
    projectDict = {}
    for project in projects:
        projectDict[project.projectid] = project.project_name
    metaDF['project_name'] = metaDF['projectid'].map(projectDict)

    memDiff()   # 0.0

    if finalSampleIDs:
        debug("Exploding_panda: finalSampleIDs:", finalSampleIDs)
        metaDF = metaDF.loc[finalSampleIDs]     # TODO 1.3 error here when small sample count is run in soil health (?)

    metaDF.dropna(axis=1, how='all', inplace=True)
    metaDF.dropna(axis=0, how='all', inplace=True)

    memDiff()   # 0.574
    df = pd.DataFrame(mat, index=sampleids, columns=taxaids)
    memDiff()   # 20.56

    if finalSampleIDs:
        df = df.loc[finalSampleIDs]
        df.dropna(axis=1, how='all', inplace=True)
        df.dropna(axis=0, how='all', inplace=True)
    memDiff()   # 22.04
    # check depVar, 0 for abundDF, 1 for rel_abundDF, 2 for richDF, 3 for diversityDF, 4 for abund_16SDF, ALL for all
    abundDF = df.reset_index(drop=False)
    abundDF.rename(columns={'index': 'sampleid'}, inplace=True)
    abundDF = pd.melt(abundDF, id_vars='sampleid', value_vars=taxaids)
    abundDF.set_index('sampleid', inplace=True)
    categorize(abundDF)

    countDF = pd.DataFrame()
    rel_abundDF = pd.DataFrame()
    richDF = pd.DataFrame()
    diversityDF = pd.DataFrame()

    if depVar == 0:
        countDF = pd.DataFrame({
            'taxaid': abundDF['variable'],
            'abund': abundDF['value']
        }, index=abundDF.index)
    memDiff()   # 35.65

    if depVar == 1 or depVar == 10:
        rel_abundDF = df.div(df.sum(axis=1), axis=0)
        rel_abundDF.reset_index(drop=False, inplace=True)
        rel_abundDF.rename(columns={'index': 'sampleid'}, inplace=True)
        rel_abundDF = pd.melt(rel_abundDF, id_vars='sampleid', value_vars=taxaids)
        rel_abundDF.set_index('sampleid', inplace=True)

        if depVar == 1:
            countDF = pd.DataFrame({
                'taxaid': abundDF['variable'],
                'abund': abundDF['value'],
                'rel_abund': rel_abundDF['value']
            }, index=abundDF.index)
    memDiff()   # 123.5

    if depVar == 2 or depVar == 10:
        richDF = df / df
        richDF.fillna(0, inplace=True)
        richDF.reset_index(drop=False, inplace=True)
        richDF.rename(columns={'index': 'sampleid'}, inplace=True)
        richDF = pd.melt(richDF, id_vars='sampleid', value_vars=taxaids)
        richDF.set_index('sampleid', inplace=True)

        if depVar == 2:
            countDF = pd.DataFrame({
                'taxaid': abundDF['variable'],
                'abund': abundDF['value'],
                'rich': richDF['value']
            }, index=abundDF.index)
    memDiff()   # 70.69

    if depVar == 3 or depVar == 10:
        diversityDF = -df.div(df.sum(axis=1), axis=0) * np.log(df.div(df.sum(axis=1), axis=0))  # RuntimeWarning: divide by zero encountered
        # logically, is DF/DF.sum * log(DF/DF.sum)); so if DF.sum is 0 or equivalent, divide by zero occurs
        diversityDF[np.isinf(diversityDF)] = np.nan
        diversityDF.fillna(0.0, inplace=True)
        diversityDF.reset_index(drop=False, inplace=True)
        diversityDF.rename(columns={'index': 'sampleid'}, inplace=True)
        diversityDF = pd.melt(diversityDF, id_vars='sampleid', value_vars=taxaids)
        diversityDF.set_index('sampleid', inplace=True)

        if depVar == 3:
            countDF = pd.DataFrame({
                'taxaid': abundDF['variable'],
                'abund': abundDF['value'],
                'diversity': diversityDF['value']
            }, index=abundDF.index)
    memDiff()   # 141.5

    if depVar == 4 or depVar == 10:
        if 'rRNA_copies' in metaDF.columns:
            metaDF.replace(to_replace='nan', value=0.0, inplace=True)
            metaDF['rRNA_copies'] = metaDF['rRNA_copies'].astype(float)
            abund_16SDF = df.div(df.sum(axis=1), axis=0).multiply(metaDF['rRNA_copies'] / 1000.0, axis=0)
            abund_16SDF.fillna(0.0, inplace=True)
            abund_16SDF.reset_index(drop=False, inplace=True)
            abund_16SDF.rename(columns={'index': 'sampleid'}, inplace=True)
            abund_16SDF = pd.melt(abund_16SDF, id_vars='sampleid', value_vars=taxaids)
            abund_16SDF.set_index('sampleid', inplace=True)
        else:
            rows, cols = abundDF.shape
            abund_16SDF = np.zeros(rows)

        if depVar == 4:
            countDF = pd.DataFrame({
                'taxaid': abundDF['variable'],
                'abund': abundDF['value'],
                'abund_16S': abund_16SDF['value']
            }, index=abundDF.index)

        if depVar == 10:
            countDF = pd.DataFrame({
                'taxaid': abundDF['variable'],
                'abund': abundDF['value'],
                'rel_abund': rel_abundDF['value'],
                'rich': richDF['value'],
                'diversity': diversityDF['value'],
                'abund_16S': abund_16SDF['value']
            }, index=abundDF.index)
    memDiff()   # 0.0

    # Check if there is at least one categorical variable with multiple levels
    # Remove fields with only 1 level
    remCatFields = []
    if levelDep:
        if catFields:
            tempList = catFields[:]
            for i in tempList:
                noLevels = len(list(pd.unique(metaDF[i])))
                if noLevels < 2:
                    catFields.remove(i)
                    remCatFields.append(i)
    memDiff()   # 0.003

    # merging meta and counts into saved
    categorize(metaDF)
    categorize(countDF)
    savedDF = pd.merge(metaDF, countDF, left_index=True, right_index=True, how='inner')
    memDiff()
    categorize(savedDF)
    memDiff()   # 14.76

    savedDF.reset_index(drop=False, inplace=True)

    memDiff()   # 0.004

    taxaDF = pd.DataFrame.from_dict(taxaDict, orient='index').reset_index()
    taxaDF.set_index('index', inplace=True)
    taxaDF.reset_index(drop=False, inplace=True)
    taxaDF.rename(columns={'index': 'taxaid'}, inplace=True)
    memDiff()   # 0.0
    #categorize(taxaDF)     # categorizing taxaDF breaks kegg search, also only saves a small amount of memory (< 1%)
    # SPIKE HEREISH
    savedDF = pd.merge(savedDF, taxaDF, left_on='taxaid', right_on='taxaid', how='inner')   # this merge is the biggest memory spike by far, seems necessary though
    memDiff()   # 427.2

    #categorize(savedDF)  # Cannot categorize this either, same kegg issue later (cuts savedDF to 1/3 size though, rouhly 5% overall savings)

    memDiff()   # 0.0

    # make sure column types are correct
    metaDF[catFields] = metaDF[catFields].astype(str)
    try:
        metaDF[quantFields] = metaDF[quantFields].astype(float)
    except:
        pass
    savedDF.reset_index(drop=True, inplace=True)

    memDiff()   # 0.0

    return savedDF, metaDF, remCatFields


def exploding_panda2(path, treeType):
    try:
        # Load file     # TODO 1.3 revert changes to this function based on getRawTabData's needs
        debug("exploding_panda2")
        file = open(path)
        data = json.load(file)
        file.close()
        # Get metadata
        d = data['columns']
        metaDict = {}
        for i in d:
            metaDict[str(i['id'])] = i['metadata']
        debug("exploding_panda2: metaDict")
        metaDF = pd.DataFrame.from_dict(metaDict, orient='index').reset_index()
        metaDF.rename(columns={'index': 'sampleid'}, inplace=True)
        metaDF.set_index('sampleid', inplace=True)
        metaDF.dropna(axis=1, how='all', inplace=True)
        metaDF.dropna(axis=0, how='all', inplace=True)
        debug("exploding_panda2: metaDF")
        # Get count data and calculate various dependent variables
        sampleids = [col['id'] for col in data['columns']]
        taxaids = [col['id'] for col in data['rows']]

        mat = data['data']
        mat = np.asarray(mat).T.tolist()
        df = pd.DataFrame(mat, index=sampleids, columns=taxaids)
        debug("exploding_panda2: df")
        # sampleids are the row labels, for some reason turning to NaN in abundDF, but present here
        abundDF = df.reset_index(drop=False)    # changed drop back to false, sampleid stays in abund now
        abundDF.rename(columns={'index': 'sampleid'}, inplace=True)
        abundDF = pd.melt(abundDF, id_vars='sampleid', value_vars=taxaids)
        abundDF.set_index('sampleid', inplace=True)
        debug("exploding_panda2: abundDF")
        countDF = pd.DataFrame({
            'taxaid': abundDF['variable'],
            'DepVar': abundDF['value']
        }, index=abundDF.index)
        debug("exploding_panda2: countDF")
        savedDF = pd.merge(metaDF, countDF, left_index=True, right_index=True, how='inner')
        savedDF.reset_index(drop=False, inplace=True)
        #savedDF.drop(['0', '1', '2'], axis=1)
        # this is empty with three index columns, depvar, and taxaid (no actual index)
        debug("exploding_panda2: savedDF")
        d = data['rows']
        taxaDict = {}
        for i in d:
            tempDict = {}
            taxon = i['metadata']['taxonomy']
            length = len(taxon)
            if treeType == 1:
                if length > 0:
                    tempDict['Kingdom'] = taxon[0]
                if length > 1:
                    tempDict['Phylum'] = taxon[1]
                if length > 2:
                    tempDict['Class'] = taxon[2]
                if length > 3:
                    tempDict['Order'] = taxon[3]
                if length > 4:
                    tempDict['Family'] = taxon[4]
                if length > 5:
                    tempDict['Genus'] = taxon[5]
                if length > 6:
                    tempDict['Species'] = taxon[6]
                if length > 7:
                    tempDict['OTU99'] = taxon[7]
            elif treeType != 1:
                if length > 0:
                    tempDict['Level 1'] = taxon[0]
                if length > 1:
                    tempDict['Level 2'] = taxon[1]
                if length > 2:
                    tempDict['Level 3'] = taxon[2]
                if length > 3:
                    tempDict['Level 4'] = taxon[3]
                if length > 4:
                    tempDict['Level 5'] = taxon[4]

            taxaDict[str(i['id'])] = tempDict
        taxaDF = pd.DataFrame.from_dict(taxaDict, orient='index').reset_index()
        taxaDF.set_index('index', inplace=True)
        taxaDF.reset_index(drop=False, inplace=True)
        taxaDF.rename(columns={'index': 'taxaid'}, inplace=True)
        debug("exploding_panda2: taxaDF")
        try:
            savedDF = pd.merge(savedDF, taxaDF, left_on='taxaid', right_on='taxaid', how='inner')
            savedDF.reset_index(drop=True, inplace=True)
        except Exception as ACK:
            print "Exception with savedDF in explodingPanda2:", ACK
        debug("exploding_panda2: savedDF")
        # missing abundances? or does depvar cover that?
        return savedDF
    except Exception as exce:
        print "Exploding panda 2 exception:" + str(exce)


def imploding_panda(path, treeType, DepVar, myList, metaDF, finalDF):
    debug("IMPLOSION")
    myBiom = {}
    nameList = []
    tempDF = metaDF.set_index('sampleid', inplace=False)
    myList.sort()
    debug("SORT")
    for i in myList:
        try:
            nameList.append({"id": str(i), "metadata": tempDF.loc[i].to_dict()})
        except KeyError:
            #print "IMPLODING KEY ERROR:", i
            pass    # this exception occurs an alarming number of times (like 25% of myList)

    debug("NAMELIST")
    # get list of lists with abundances
    abundDF = pd.DataFrame()
    if DepVar == 0:
        abundDF = finalDF.pivot(index='rank_id', columns='sampleid', values='abund')
    if DepVar == 1:
        abundDF = finalDF.pivot(index='rank_id', columns='sampleid', values='rel_abund')
    elif DepVar == 2:
        abundDF = finalDF.pivot(index='rank_id', columns='sampleid', values='rich')
    elif DepVar == 3:
        abundDF = finalDF.pivot(index='rank_id', columns='sampleid', values='diversity')
    elif DepVar == 4:
        abundDF = finalDF.pivot(index='rank_id', columns='sampleid', values='abund_16S')

    abundDF.sort_index(axis=0, inplace=True)
    abundList = abundDF.values.tolist()

    debug("ABUND")
    # get list of taxa
    rank_id = abundDF.index.values.tolist()
    taxDict = {}
    if treeType == 1:
        taxDict = functions.getFullTaxonomy(rank_id)
    elif treeType == 2:
        taxDict = functions.getFullKO(rank_id)
    elif treeType == 3:
        taxDict = functions.getFullNZ(rank_id)

    for i in taxDict:
        if taxDict[i] == 'N|o|t| |f|o|u|n|d':
            taxDict[i] = 'N/A'

    debug("FULL TAXA")
    namesDF = pd.DataFrame.from_dict(taxDict, orient='index')
    namesDF.sort_index(axis=0, inplace=True)

    taxaList = []
    for index, row in namesDF.iterrows():
        rowList = row[0].split('|')
        metaDict = {'taxonomy':  rowList}
        taxaList.append({"id": index, "metadata": metaDict})

    shape = [len(taxaList), len(nameList)]
    myBiom['id'] = 'None'
    myBiom['format'] = 'myPhyloDB v.1.2.0'
    myBiom['format_url'] = 'http://biom-format.org'
    myBiom['generated_by'] = 'myPhyloDB'
    myBiom['type'] = 'OTU table'
    myBiom['date'] = str(datetime.datetime.now())
    myBiom['matrix_type'] = 'dense'
    myBiom['matrix_element_type'] = 'float'
    myBiom["shape"] = shape
    myBiom['rows'] = taxaList
    myBiom['columns'] = nameList
    myBiom['data'] = abundList

    debug("MYBIOM")
    with open(path, 'w') as outfile:
        json.dump(myBiom, outfile, ensure_ascii=True, indent=4)


def recLabels(lists, level):
    if lists.__len__() == 0:
        return {}

    first = lists
    splitset = []
    for i in range(0, level):
        children = []
        parents = []
        for set in first:
            children.append(set[set.__len__()-1])
            parents.append(set[0:set.__len__()-1])
        first = parents
        splitset.append(children)
    return makeLabels(" ", splitset)


# TODO 1.4 add documentation for old code like this
def makeLabels(name, list):
    retDict = {}
    if list.__len__() == 1:
        # final layer
        retDict['name'] = name
        retDict['categories'] = list[0]
        return retDict

    # change here
    children = []
    first = list[list.__len__()-1][0]
    iter = 0
    start = 0
    for stuff in list[list.__len__()-1]:
        if stuff != first:
            sublist = []
            for otherstuff in list[0:list.__len__()-1]:
                sublist.append(otherstuff[start:iter])
            children.append(makeLabels(first, sublist))
            first = stuff
            start = iter
        iter += 1

    # Repeat else condition at the end of the list
    sublist = []
    for otherstuff in list[0:list.__len__()-1]:
        sublist.append(otherstuff[start:iter])
    children.append(makeLabels(first, sublist))

    retDict['name'] = name
    retDict['categories'] = children
    return retDict


def remove_string_from_file_name(file_path, string_to_remove, dry_run=False):
    path, file_name = os.path.split(file_path)
    new_name = file_name.replace(string_to_remove, '')
    new_path = os.path.join(path, new_name)
    if dry_run:
        print "Would rename {} to {}".format(file_path, new_path)
    else:
        os.rename(file_path, new_path)


def rewrite_biom(samples, biomPath, selected):
    headers = []
    for thing in samples:
        headers.append(thing)
    with open(biomPath) as f:
        myJson = json.load(f)

    for part in myJson['rows']:
        # id is database id of lowest level of CORE taxa
        # metadata is a list of names of taxa, going from kingdom down to lowest level
        # need to take in ids and send them to function, can change function to run one at a time
        myID = part['id'][0]
        myParents = functions.getFullTaxaFromID(myID, selected)  # where is level from at this point
        myMeta = part['metadata'][0]
        newMeta = []
        newIter = 0
        # actually want k:v of metadata: {taxonomy:[list of strings]}
        for parent in myMeta:
            newMeta.append(str(myParents[newIter]) + ": " + str(parent))
            newIter += 1
        part['id'] = myID
        part['metadata'] = {'taxonomy': newMeta}
    for part in myJson['columns']:
        # part[metadata] is a list of values, each one needs its column
        iter = 0
        newPart = {}
        for subpart in part['metadata']:
            if headers[iter] != "sampleid":
                newPart[headers[iter]] = subpart
            iter += 1
        part['id'] = part['id'][0]
        part['metadata'] = newPart
    myJson['format'] = "myPhyloDB v.1.2.0"
    myJson['generated_by'] = "myPhyloDB"
    myJson['matrix_type'] = myJson['matrix_type'][0]
    myJson['format_url'] = myJson['format_url'][0]
    myJson['date'] = myJson['date'][0]
    myJson['type'] = myJson['type'][0]
    myJson['id'] = "None"
    myJson['matrix_element_type'] = myJson['matrix_element_type'][0]
    with open(biomPath, 'w') as outfile:
        json.dump(myJson, outfile, ensure_ascii=True, indent=4)


def write_taxa_summary(path, data):
    # write lists of unique taxa for each level based on given data
    # output file is exactly path (includes file name)
    # should check that path is valid first, directory exists etc

    # find each unique value among the levels of taxa
    # final output should be one row per level of taxa available, each row contains a list of unique ids
    # we are not necessarily given all levels, but each taxa will go to the SAME level for this call (not useful)
    # data should be a list of dictionaries
    # using dictionaries to store unique values we find as keys, since duplicate keys affect the same entry, and hashing
    # should have much better performance than searching each list for each entry
    uniqueKingdoms = {}
    uniquePhyla = {}
    uniqueClasses = {}
    uniqueOrders = {}
    uniqueFamilies = {}
    uniqueGenera = {}
    uniqueSpecies = {}
    uniqueOTUs = {}

    for taxonomy in data:
        # taxonomy should be a dictionary, id:lowestid (not relevant) metadata:{taxonomy:[taxaLevelList]}
        localIter = 0   # using this to track level
        for id in taxonomy['metadata']['taxonomy']:
            id = id.split(": ")[0]
            if localIter == 0:
                # kingdom
                uniqueKingdoms[id] = 0
            if localIter == 1:
                # phyla
                uniquePhyla[id] = 0
            if localIter == 2:
                # class
                uniqueClasses[id] = 0
            if localIter == 3:
                # order
                uniqueOrders[id] = 0
            if localIter == 4:
                # family
                uniqueFamilies[id] = 0
            if localIter == 5:
                # genus
                uniqueGenera[id] = 0
            if localIter == 6:
                # species
                uniqueSpecies[id] = 0
            if localIter == 7:
                # otu
                uniqueOTUs[id] = 0
            localIter += 1
            # TODO 1.3 PGPR?

    uniqueKingdoms = uniqueKingdoms.keys()
    uniquePhyla = uniquePhyla.keys()
    uniqueClasses = uniqueClasses.keys()
    uniqueOrders = uniqueOrders.keys()
    uniqueFamilies = uniqueFamilies.keys()
    uniqueGenera = uniqueGenera.keys()
    uniqueSpecies = uniqueSpecies.keys()
    uniqueOTUs = uniqueOTUs.keys()

    # now we have lists of unique ids for each level, write to file one line per level
    # alternatively, make a dictionary and just tell pickle to do it
    taxa_summary_dict = {'kingdom': uniqueKingdoms, 'phyla': uniquePhyla, 'class': uniqueClasses, 'order': uniqueOrders,
                         'family': uniqueFamilies, 'genus': uniqueGenera, 'species': uniqueSpecies, 'otu': uniqueOTUs}

    with open(path, 'wb') as f:
        pickle.dump(taxa_summary_dict, f)

