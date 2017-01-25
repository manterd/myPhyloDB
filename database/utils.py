from collections import defaultdict
import ctypes
from django import forms
from django.db.models import Q
from django.core.exceptions import ValidationError
from django.http import HttpResponse
import inspect
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import re
import shutil
import simplejson
import threading
import time

from models import Project, Reference, Profile
from config import local_cfg
from utils_kegg import getFullTaxonomy, getFullKO, getFullNZ


pd.set_option('display.max_colwidth', -1)


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
    if not os.path.exists(path):
        os.makedirs(path)
    dest = "/".join([str(path), str(name)])
    with open(str(dest), 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def remove_list(request):
    paths = request.POST.getlist('chkbx')
    p_uuid = list(Reference.objects.filter(path__in=paths).values_list('projectid', flat=True))

    # Remove record from reference table
    for path in paths:
        Reference.objects.filter(path=path).delete()
        if os.path.exists(path):
            shutil.rmtree(path, ignore_errors=True)

    # Remove record from project table and path (if necessary)
    for pid in p_uuid:
        qs = Project.objects.all().filter(projectid=pid)
        for item in qs:
            refLength = len(item.reference_set.values_list('refid', flat=True))
            if refLength == 0:
                Project.objects.filter(projectid=pid).delete()
                path = os.path.join("uploads", str(pid))
                if os.path.exists(path):
                    shutil.rmtree(path, ignore_errors=True)
            else:
                pass


def remove_proj(path):
    pid = Reference.objects.get(path=path).projectid.projectid

    # Remove paths (if necessary)
    qs = Project.objects.all().filter(projectid=pid)
    for item in qs:
        refLength = len(item.reference_set.values_list('refid', flat=True))
        if refLength == 0:
            path = os.path.join("uploads", str(pid))
            if os.path.exists(path):
                shutil.rmtree(path, ignore_errors=True)
        else:
            pass


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


def PCoA(dm):
    E_matrix = make_E_matrix(dm)
    F_matrix = make_F_matrix(E_matrix)
    eigvals, eigvecs = np.linalg.eigh(F_matrix)
    negative_close_to_zero = np.isclose(eigvals, 0)
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
    num_positive = (eigvals >= 0).sum()
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
    def render(self, name, value, attrs={}):
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


def threads():
    try:
        usr_threads = int(local_cfg.usr_num_threads)
    except:
        usr_threads = 1

    if usr_threads <= 0:
        num_threads = 1
    elif usr_threads > mp.cpu_count():
        num_threads = int(mp.cpu_count())
    else:
        num_threads = usr_threads
    return int(num_threads)


def getRawData(request):
    if request.is_ajax():
        RID = request.GET["all"]
        func = request.GET["func"]
        treeType = int(request.GET["treeType"])

        myDir = 'myPhyloDB/media/temp/' + str(func) + '/'
        fileName = str(myDir) + str(RID) + '.pkl'
        savedDF = pd.read_pickle(fileName)

        fRow, fCol = savedDF.shape
        if treeType == 1:
            idList = getFullTaxonomy(list(savedDF.rank_id.unique()))
            savedDF['Taxonomy'] = savedDF['rank_id'].map(idList)
        elif treeType == 2:
            idList = getFullKO(list(savedDF.rank_id.unique()))
            savedDF['Taxonomy'] = savedDF['rank_id'].map(idList)
        elif treeType == 3:
            idList = getFullNZ(list(savedDF.rank_id.unique()))
            savedDF['Taxonomy'] = savedDF['rank_id'].map(idList)

        savedDF.replace(to_replace='N/A', value=np.nan, inplace=True)
        savedDF.dropna(axis=1, how='all', inplace=True)
        savedDF.drop('rank_name', axis=1, inplace=True)

        myDir = 'myPhyloDB/media/temp/' + str(func) + '/'
        fileName = str(myDir) + str(RID) + '.csv'
        savedDF.to_csv(fileName)

        myDict = {}
        myDir = '/myPhyloDB/media/temp/' + str(func) + '/'
        fileName = str(myDir) + str(RID) + '.csv'
        myDict['name'] = str(fileName)
        res = simplejson.dumps(myDict)

        return HttpResponse(res, content_type='application/json')


def removeFiles(request):
    if request.is_ajax():
        RID = request.GET["all"]
        func = request.GET["func"]

        file = "myPhyloDB/media/temp/" + str(func) + "/Rplots/" + str(RID) + "." + str(func) + ".pdf"
        if os.path.exists(file):
            os.remove(file)

        file = "myPhyloDB/media/temp/" + str(func) + "/" + str(RID) + ".pkl"
        if os.path.exists(file):
            os.remove(file)

        file = "myPhyloDB/media/temp/" + str(func) + "/" + str(RID) + ".csv"
        if os.path.exists(file):
            os.remove(file)

        return HttpResponse()


def excel_to_dict(wb, headerRow=1, nRows=1, sheet='Sheet1'):
    ws = wb.get_sheet_by_name(sheet)
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


def getMetaDF(savedDF, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, DepVar):

    catFields = []
    catValues = []
    if metaValsCat:
        metaDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaValsCat)
        for key in sorted(metaDictCat):
            catFields.append(key)
            catValues.extend(metaDictCat[key])

    catSampleLists = []
    if metaIDsCat:
        idDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaIDsCat)
        for key in sorted(idDictCat):
            catSampleLists.append(idDictCat[key])

    catSampleIDs = []
    if catSampleLists:
        catSampleIDs = list(set.intersection(*map(set, catSampleLists)))

    quantFields = []
    quantValues = []
    if metaValsQuant:
        metaDictQuant = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaValsQuant)
        for key in sorted(metaDictQuant):
            quantFields.append(key)
            quantValues.extend(metaDictQuant[key])

    quantSampleLists = []
    if metaIDsQuant:
        idDictQuant = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaIDsQuant)
        for key in sorted(idDictQuant):
            quantSampleLists.append(idDictQuant[key])

    quantSampleIDs = []
    if quantSampleLists:
        quantSampleIDs = list(set.intersection(*map(set, quantSampleLists)))

    finalSampleIDs = list(set(catSampleIDs) | set(quantSampleIDs))

    # remove samples not selected
    if finalSampleIDs:
        savedDF = savedDF.loc[savedDF['sampleid'].isin(finalSampleIDs)]
    metaDF = savedDF.drop_duplicates(subset='sampleid', take_last=True)

    # Check if there is at least one categorical variable with multiple levels
    # Remove fields with only 1 level
    remCatFields = []
    if catFields:
        tempList = catFields[:]
        for i in tempList:
            noLevels = len(list(pd.unique(metaDF[i])))
            if noLevels < 2:
                catFields.remove(i)
                remCatFields.append(i)
        if remCatFields:
            allFields = catFields + quantFields
            wantedList = allFields + ['sampleid', 'sample_name']
            metaDF = metaDF[wantedList]

    # remove samples that do not have rRNA copy number data available
    if DepVar == 4:
        rnaDF = savedDF.loc[savedDF['abund_16S'] != 0]
        finalSampleIDs = list(pd.unique(rnaDF['sampleid']))
        savedDF = savedDF.loc[savedDF['sampleid'].isin(finalSampleIDs)]
        metaDF = metaDF.loc[metaDF['sampleid'].isin(finalSampleIDs)]

    # remove unnecessary fields
    wantedList = catFields + quantFields + ['sampleid', 'kingdomid', 'kingdomName', 'phylaid', 'phylaName', 'classid', 'className', 'orderid', 'orderName', 'familyid', 'familyName', 'genusid', 'genusName', 'speciesid', 'speciesName', 'otuid', 'otuName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']
    savedDF = savedDF[wantedList]

    wantedList = catFields + quantFields + ['sampleid', 'sample_name']
    wantedList = list(set(wantedList))
    metaDF = metaDF[wantedList]
    metaDF.set_index('sampleid', drop=True, inplace=True)

    # make sure column types are correct
    metaDF[catFields] = metaDF[catFields].astype(str)
    metaDF[quantFields] = metaDF[quantFields].astype(float)

    savedDF.dropna(axis=0, inplace=True)
    metaDF.dropna(axis=0, inplace=True)
    finalSampleIDs = metaDF.index.tolist()

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


def getViewProjects(request):  # JUMP

    projects = Project.objects.none()
    if request.user.is_superuser:
        projects = Project.objects.all().order_by('project_name')
    elif request.user.is_authenticated():
        # run through list of projects, when valid project is found, append filterIDS with ID
        # projects will be a queryset set to all projects, then filtered by ids in filterIDS
        filterIDS = []
        for proj in Project.objects.all():
            good = False  # good to add to list
            if proj.owner == request.user:
                good = True
            if proj.status == 'public':
                good = True
            checkList = proj.whitelist_view.split(';')
            for name in checkList:
                if name == request.user.username:
                    good = True
            if good:
                filterIDS.append(proj.projectid)
        projects = Project.objects.all().filter(projectid__in=filterIDS)
    if not request.user.is_superuser and not request.user.is_authenticated():
        # impossible to have guest user be on whitelist (hopefully), so public only check
        projects = Project.objects.all().filter( Q(status='public') ).order_by('project_name')
    return projects


def getEditProjects(request):

    projects = Project.objects.none()
    if request.user.is_superuser:
        projects = Project.objects.all().order_by('project_name')

    elif request.user.is_authenticated():
        # run through list of projects, when valid project is found, append filterIDS with ID
        # projects will be a queryset set to all projects, then filtered by ids in filterIDS
        filterIDS = []
        for proj in Project.objects.all():
            good = False  # good to add to list
            if proj.owner == request.user:
                good = True
            checkList = proj.whitelist_edit.split(';')
            for name in checkList:
                if name == request.user.username:
                    good = True
            if good:
                filterIDS.append(proj.projectid)
        projects = Project.objects.all().filter(projectid__in=filterIDS)

    return projects
