from collections import defaultdict
import ctypes
from django import forms
from django.core.exceptions import ValidationError
from django.db import transaction
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


def handle_uploaded_file(f, path, name):
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
    qs1 = Profile.objects.filter(sampleid__in=mySet).values('sampleid', 'kingdomid', 'phylaid', 'classid', 'orderid', 'familyid', 'genusid', 'speciesid', 'count')
    df = pd.DataFrame.from_records(qs1, columns=['sampleid', 'kingdomid', 'phylaid', 'classid', 'orderid', 'familyid', 'genusid', 'speciesid', 'count'])
    df = df.groupby(['sampleid', 'kingdomid', 'phylaid', 'classid', 'orderid', 'familyid', 'genusid', 'speciesid'])['count'].sum()
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

        myDir = 'myPhyloDB/media/temp/' + str(func) + '/'
        fileName = str(myDir) + str(RID) + '.pkl'
        savedDF = pd.read_pickle(fileName)

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
    for col in xrange(1, ws.max_column):
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
