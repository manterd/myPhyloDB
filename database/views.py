import ast
import datetime
from django.contrib.auth.models import User
from django.contrib.auth.hashers import check_password
from django.contrib.auth.decorators import login_required, user_passes_test
from django.contrib.auth.signals import user_logged_in
from django.db import transaction
from django.db.models import Q
from django.http import *
from django_pandas.io import read_frame
from django.shortcuts import render_to_response
from django.template import RequestContext
import fileinput
import logging
import multiprocessing as mp
import os
import pandas as pd
import pickle
import shutil
import simplejson
import tarfile
import time
import zipfile
from uuid import uuid4

from forms import UploadForm1, UploadForm2, UploadForm4, UploadForm5, \
    UploadForm6, UploadForm7, UploadForm8, UploadForm9, UserRegForm
from models import Project, Reference, Sample, Air, Microbial, Human_Associated, Soil, Water, UserDefined, OTU_99
from models import ko_lvl1, ko_lvl2, ko_lvl3, ko_entry
from models import nz_lvl1, nz_lvl2, nz_lvl3, nz_lvl4, nz_entry
from models import PICRUSt, UserProfile
from parsers import mothur, projectid, parse_project, parse_sample, parse_taxonomy, parse_profile
from utils import handle_uploaded_file, remove_list, remove_proj
from models import addQueue, getQueue, subQueue
from database.pybake.pybake import koParse, nzParse
from database.pybake.pybake import geneParse
from database.utils import cleanup


rep_project = ''
pd.set_option('display.max_colwidth', -1)
LOG_FILENAME = 'error_log.txt'


def home(request):
    return render_to_response(
        'home.html',
        context_instance=RequestContext(request)
    )


def map(request):
    return render_to_response(
        'map.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/myPhyloDB/accounts/login/')
def upload(request):
    projects = Reference.objects.none()
    if request.user.is_superuser:
        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
    elif request.user.is_authenticated():
        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)

    return render_to_response(
        'upload.html',
        {'projects': projects,
         'form1': UploadForm1,
         'form2': UploadForm2,
         'error': ""},
        context_instance=RequestContext(request)
    )

@login_required(login_url='/myPhyloDB/accounts/login/')
def download(request):

    projects = Reference.objects.none()
    if request.user.is_superuser:
        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
    elif request.user.is_authenticated():
        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)

    return render_to_response(
        'download.html',
        {'projects': projects,
         'type': "GET",
         'paths': ''},
        context_instance=RequestContext(request)
    )


@login_required(login_url='/myPhyloDB/accounts/login/')
def profile(request):

    projects = Project.objects.none()
    if request.user.is_superuser:
        projects = Project.objects.all().order_by('project_name')
    elif request.user.is_authenticated():
        path_list = Reference.objects.filter(Q(author=request.user)).values_list('projectid_id')
        projects = Project.objects.all().filter( Q(projectid__in=path_list) | Q(status='public') ).order_by('project_name')
    if not request.user.is_superuser and not request.user.is_authenticated():
        projects = Project.objects.all().filter( Q(status='public') ).order_by('project_name')

    return render_to_response(
        'profile.html',
        {'projects': projects},
        context_instance=RequestContext(request)
    )


def getProjectFiles(request):
    all = request.GET['all']
    selKeys = simplejson.loads(all)

    paths = Reference.objects.filter(refid__in=selKeys).values_list('path', flat=True)
    zip_file = os.path.join(os.getcwd(), 'myPhyloDB', 'media', 'usr_temp', request.user.username, 'final_data.zip')
    zf = zipfile.ZipFile(zip_file, "w")
    foundCount = 0
    errorMsg = "none"
    for path in paths:
        if os.path.exists(path):
            foundCount += 1
            for root, dirs, files in os.walk(path):
                for file in files:
                    if file.endswith('.shared') or file.endswith('.taxonomy') or file.endswith('.fasta') or file.endswith('.names') or file.endswith('.groups') or file.endswith('.logfile') or file.endswith('.xlsx') or file.endswith('.xls'):
                        zf.write(os.path.join(root, file), zipfile.ZIP_DEFLATED)
    zf.close()

    if foundCount != len(paths):
        errorMsg = "Failed to locate project files"

    results = {'files': zip_file, 'error': errorMsg}
    myJson = simplejson.dumps(results, ensure_ascii=False)
    return HttpResponse(myJson)


def upStop(request):
    # cleanup mid upload project!
    projects = Reference.objects.none()
    if request.user.is_superuser:
        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
    elif request.user.is_authenticated():
        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
    return render_to_response(
            'upload.html',
            {'projects': projects,
             'form1': UploadForm1,
             'form2': UploadForm2,
             'error': "Upload stopped"
             },
            context_instance=RequestContext(request)
    )


def upErr(msg, request, dest, sid):
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
    logging.exception(myDate)
    remove_proj(dest)
    transaction.savepoint_rollback(sid)

    if request.user.is_superuser:
        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
    elif request.user.is_authenticated():
        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
    else:
        projects = None
    return render_to_response(
        'upload.html',
        {'projects': projects,
         'form1': UploadForm1,
         'form2': UploadForm2,
         'error': msg
         },
        context_instance=RequestContext(request)
    )


@transaction.atomic
def uploadFunc(request, stopList):
    ### create a savepoint
    sid = transaction.savepoint()
    curUser = User.objects.get(username=request.user.username)

    ### start of main upload function
    projects = Reference.objects.none()
    if request.method == 'POST' and 'Upload' in request.POST:
        start = datetime.datetime.now()
        form1 = UploadForm1(request.POST, request.FILES)
        source = str(request.POST['source'])
        userID = str(request.user.id)
        processors = int(request.POST['processors'])
        RID = request.POST['RID']
        PID = 0  # change if adding additional data threads

        if stopList[PID] == RID:
            return upStop(request)

        if form1.is_valid():
            file1 = request.FILES['docfile1']
            try:
                p_uuid, pType, num_samp = projectid(file1)
                if not num_samp:
                    return render_to_response(
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "Error: Please open and save your meta file in Excel/OpenOffice in order to perform the necessary calculations..."
                         },
                        context_instance=RequestContext(request)
                    )

            except:
                logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                logging.exception(myDate)

                if request.user.is_superuser:
                    projects = Reference.objects.all().order_by('projectid__project_name', 'path')
                elif request.user.is_authenticated():
                    projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
                return render_to_response(
                    'upload.html',
                    {'projects': projects,
                     'form1': UploadForm1,
                     'form2': UploadForm2,
                     'error': "There was an error parsing your meta file:" + str(file1.name)
                     },
                    context_instance=RequestContext(request)
                )

            date = datetime.date.today().isoformat()
            hour = datetime.datetime.now().hour
            minute = datetime.datetime.now().minute
            second = datetime.datetime.now().second
            timestamp = ".".join([str(hour), str(minute), str(second)])
            datetimestamp = "_".join([str(date), str(timestamp)])
            dest = "/".join(["uploads", str(p_uuid), str(datetimestamp)])
            metaName = 'final_meta.xlsx'
            metaFile = '/'.join([dest, metaName])

            try:
                handle_uploaded_file(file1, dest, metaName)
                res = parse_project(metaFile, p_uuid, curUser)
                if res == "none":
                    print "Parsed project with no errors"
                else:
                    print "Encountered error while parsing meta file: "+res
            except:
                return upErr("There was an error parsing your meta file:" + str(file1.name), request, dest, sid)

            if source == 'mothur':
                raw = False
                batch = 'blank'
            elif source == '454_sff':
                raw = True
                batch = request.FILES['docfile7']
            elif source == '454_fastq':
                raw = True
                batch = request.FILES['docfile7']
            elif source == 'miseq':
                raw = True
                batch = request.FILES['docfile15']
            else:
                raw = False
                batch = ''

            if stopList[PID] == RID:
                remove_proj(dest)
                transaction.savepoint_rollback(sid)
                return upStop(request)

            try:
                refDict = parse_sample(metaFile, p_uuid, pType, num_samp, dest, batch, raw, source, userID)
            except:
                return upErr("There was an error parsing your meta file:" + str(file1.name), request, dest, sid)

            if stopList[PID] == RID:
                remove_proj(dest)
                transaction.savepoint_rollback(sid)
                return upStop(request)

            if source == 'mothur':
                try:
                    file3 = request.FILES['docfile3']
                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)

                    if request.user.is_superuser:
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
                    elif request.user.is_authenticated():
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
                    return render_to_response(
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "Cannot open taxonomy file"
                         },
                        context_instance=RequestContext(request)
                    )
                handle_uploaded_file(file3, dest, file3.name)

                try:
                    parse_taxonomy(file3)
                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)

                    if request.user.is_superuser:
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
                    elif request.user.is_authenticated():
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
                    return render_to_response(
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "There was an error parsing your taxonomy file:" + str(file3.name)
                         },
                        context_instance=RequestContext(request)
                    )

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                try:
                    file4 = request.FILES['docfile4']
                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)

                    if request.user.is_superuser:
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
                    elif request.user.is_authenticated():
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
                    return render_to_response(
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "Cannot open shared file"
                         },
                        context_instance=RequestContext(request)
                    )
                handle_uploaded_file(file4, dest, file4.name)

                try:
                    parse_profile(file3, file4, p_uuid, refDict)
                    end = datetime.datetime.now()
                    print 'Total time for upload:', end-start
                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)

                    if request.user.is_superuser:
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
                    elif request.user.is_authenticated():
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
                    return render_to_response(
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "There was an error parsing your shared file:" + str(file4.name)
                         },
                        context_instance=RequestContext(request)
                    )

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

            elif source == '454_sff':
                mothurdest = 'mothur/temp'
                if not os.path.exists(mothurdest):
                    os.makedirs(mothurdest)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                file_list = request.FILES.getlist('sff_files')
                for file in file_list:  # .tar.gz not included in this somehow?
                    if file.name.endswith(('.gz', '.tgz', '.tar')):
                        handle_uploaded_file(file, dest, file.name)
                        tar = tarfile.open(os.path.join(dest, file.name))
                        tar.extractall(path=mothurdest)
                        tar.close()
                    elif file.name.endswith('.zip'):
                        handle_uploaded_file(file, dest, file.name)
                        zip = zipfile.ZipFile(os.path.join(dest, file.name))
                        zip.extractall(mothurdest)
                        zip.close()
                    else:
                        handle_uploaded_file(file, mothurdest, file.name)
                    handle_uploaded_file(file, dest, file.name)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                file_list = request.FILES.getlist('oligo_files')
                for file in file_list:
                    if file.name.endswith(('.gz', '.tgz', '.tar')):
                        handle_uploaded_file(file, dest, file.name)
                        tar = tarfile.open(os.path.join(dest, file.name))
                        tar.extractall(path=mothurdest)
                        tar.close()
                    elif file.name.endswith('.zip'):
                        handle_uploaded_file(file, dest, file.name)
                        zip = zipfile.ZipFile(os.path.join(dest, file.name))
                        zip.extractall(mothurdest)
                        zip.close()
                    else:
                        handle_uploaded_file(file, mothurdest, file.name)
                    handle_uploaded_file(file, dest, file.name)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                file5 = request.FILES['docfile5']
                handle_uploaded_file(file5, mothurdest, 'temp.txt')
                handle_uploaded_file(file5, dest, 'temp.txt')

                batch = 'mothur.batch'
                file7 = request.FILES['docfile7']

                avail_proc = mp.cpu_count()
                use_proc = min(avail_proc, processors)
                actual_proc = 'processors=' + str(use_proc)

                handle_uploaded_file(file7, mothurdest, batch)

                for line in fileinput.input('mothur/temp/mothur.batch', inplace=1):
                    print line.replace("processors=X", actual_proc),

                handle_uploaded_file(file7, dest, batch)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                # check queue for mothur availability, then run or wait
                # if mothur is available, run it
                # else sleep, loop back to check again
                addQueue()
                while getQueue() > 1:
                    time.sleep(5)
                try:
                    mothur(dest, source)
                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)

                    if request.user.is_superuser:
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
                    elif request.user.is_authenticated():
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
                    return render_to_response(
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "There was an error with your mothur batch file:" + str(file7.name)
                         },
                        context_instance=RequestContext(request)
                    )

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                subQueue()

                try:
                    with open('% s/final.taxonomy' % dest, 'rb') as file3:
                        parse_taxonomy(file3)
                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)

                    if request.user.is_superuser:
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
                    elif request.user.is_authenticated():
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
                    return render_to_response(
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "There was an error parsing your taxonomy file: final.taxonomy"
                         },
                        context_instance=RequestContext(request)
                    )

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                try:
                    with open('% s/final.taxonomy' % dest, 'rb') as file3:
                        with open('% s/final.shared' % dest, 'rb') as file4:
                            parse_profile(file3, file4, p_uuid, refDict)
                            end = datetime.datetime.now()
                    print 'Total time for upload:', end-start

                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)

                    if request.user.is_superuser:
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
                    elif request.user.is_authenticated():
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
                    return render_to_response(
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "There was an error parsing your shared file: final.shared"
                         },
                        context_instance=RequestContext(request)
                    )

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

            elif source == '454_fastq':
                mothurdest = 'mothur/temp'

                if not os.path.exists(mothurdest):
                    os.makedirs(mothurdest)

                file_list = request.FILES.getlist('fna_files')
                tempList = []
                if len(file_list) > 1:
                    for file in file_list:
                        if file.name.endswith(('.gz', '.tgz', '.tar')):
                            handle_uploaded_file(file, dest, file.name)
                            tar = tarfile.open(os.path.join(dest, file.name))
                            tar.extractall(path=mothurdest)
                            tar.close()
                        elif file.name.endswith('.zip'):
                            handle_uploaded_file(file, dest, file.name)
                            zip = zipfile.ZipFile(os.path.join(dest, file.name))
                            zip.extractall(mothurdest)
                            zip.close()
                        else:
                            handle_uploaded_file(file, mothurdest, file.name)
                        handle_uploaded_file(file, dest, file.name)
                        if os.name == 'nt':
                            myStr = "mothur\\temp\\" + str(file.name)
                        else:
                            myStr = "mothur/temp/" + str(file.name)
                        tempList.append(myStr)
                    inputList = "-".join(tempList)
                    if os.name == 'nt':
                        os.system('"mothur\\mothur-win\\mothur.exe \"#merge.files(input=%s, output=mothur\\temp\\temp.fasta)\""' % inputList)
                    else:
                        os.system("mothur/mothur-linux/mothur \"#merge.files(input=%s, output=mothur/temp/temp.fasta)\"" % inputList)
                else:
                    for file in file_list:
                        fasta = 'temp.fasta'  # potentially change file.name for fasta on zip saves
                        if file.name.endswith(('.gz', '.tgz')):
                            handle_uploaded_file(file, dest, file.name)
                            tar = tarfile.open(os.path.join(dest, file.name))
                            tar.extractall(path=mothurdest)
                            tar.close()
                        elif file.name.endswith('.zip'):
                            handle_uploaded_file(file, dest, file.name)
                            zip = zipfile.ZipFile(os.path.join(dest, file.name))
                            zip.extractall(mothurdest)
                            zip.close()
                        else:
                            handle_uploaded_file(file, mothurdest, fasta)
                        handle_uploaded_file(file, dest, file.name)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                file_list = request.FILES.getlist('qual_files')
                tempList = []
                if len(file_list) > 1:
                    for file in file_list:
                        if file.name.endswith(('.gz', '.tgz', '.tar')):
                            handle_uploaded_file(file, dest, file.name)
                            tar = tarfile.open(os.path.join(dest, file.name))
                            tar.extractall(path=mothurdest)
                            tar.close()
                        elif file.name.endswith('.zip'):
                            handle_uploaded_file(file, dest, file.name)
                            zip = zipfile.ZipFile(os.path.join(dest, file.name))
                            zip.extractall(mothurdest)
                            zip.close()
                        else:
                            handle_uploaded_file(file, mothurdest, file.name)
                        if os.name == 'nt':
                            myStr = "mothur\\temp\\" + str(file.name)
                        else:
                            myStr = "mothur/temp/" + str(file.name)
                        tempList.append(myStr)
                        handle_uploaded_file(file, dest, file.name)
                    inputList = "-".join(tempList)
                    if os.name == 'nt':
                        os.system('"mothur\\mothur-win\\mothur.exe \"#merge.files(input=%s, output=mothur\\temp\\temp.qual)\""' % inputList)
                    else:
                        os.system("mothur/mothur-linux/mothur \"#merge.files(input=%s, output=mothur/temp/temp.qual)\"" % inputList)
                else:
                    for file in file_list:
                        qual = 'temp.qual'
                        if file.name.endswith(('.gz', '.tgz')):  # change file.name to qual potentially
                            handle_uploaded_file(file, dest, file.name)
                            tar = tarfile.open(os.path.join(dest, file.name))
                            tar.extractall(path=mothurdest)
                            tar.close()
                        elif file.name.endswith('.zip'):
                            handle_uploaded_file(file, dest, file.name)
                            zip = zipfile.ZipFile(os.path.join(dest, file.name))
                            zip.extractall(mothurdest)
                            zip.close()
                        else:
                            handle_uploaded_file(file, mothurdest, qual)
                        handle_uploaded_file(file, dest, file)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                oligo = 'temp.oligos'
                file6 = request.FILES['docfile6']
                handle_uploaded_file(file6, mothurdest, oligo)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                handle_uploaded_file(file6, dest, file6.name)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                batch = 'mothur.batch'
                file7 = request.FILES['docfile7']

                avail_proc = mp.cpu_count()
                use_proc = min(avail_proc, processors)
                actual_proc = 'processors=' + str(use_proc)

                handle_uploaded_file(file7, mothurdest, batch)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                for line in fileinput.input('mothur/temp/mothur.batch', inplace=1):
                    print line.replace("processors=X", actual_proc),

                handle_uploaded_file(file7, dest, batch)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                try:
                    mothur(dest, source)

                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)

                    if request.user.is_superuser:
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
                    elif request.user.is_authenticated():
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
                    return render_to_response(
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "There was an error with your mothur batch file: " + str(file7.name)
                         },
                        context_instance=RequestContext(request)
                    )

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                try:
                    with open('% s/final.taxonomy' % dest, 'rb') as file3:
                        parse_taxonomy(file3)
                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)

                    if request.user.is_superuser:
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
                    elif request.user.is_authenticated():
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
                    return render_to_response(
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "There was an error parsing your taxonomy file: final.taxonomy"
                         },
                        context_instance=RequestContext(request)
                    )

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                try:
                    with open('% s/final.taxonomy' % dest, 'rb') as file3:
                        with open('% s/final.shared' % dest, 'rb') as file4:
                            parse_profile(file3, file4, p_uuid, refDict)
                            end = datetime.datetime.now()
                    print 'Total time for upload:', end-start
                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)

                    if request.user.is_superuser:
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
                    elif request.user.is_authenticated():
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
                    return render_to_response(
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "There was an error parsing your shared file: final.shared"
                         },
                        context_instance=RequestContext(request)
                    )

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

            elif source == 'miseq':
                mothurdest = 'mothur/temp'
                if not os.path.exists(mothurdest):
                    os.makedirs(mothurdest)

                fastq = 'temp.files'
                file13 = request.FILES['docfile13']
                handle_uploaded_file(file13, mothurdest, fastq)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                handle_uploaded_file(file13, dest, fastq)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                file_list = request.FILES.getlist('fastq_files')

                for file in file_list:
                    if file.name.endswith(('.gz', '.tgz', '.tar')):
                        handle_uploaded_file(file, dest, file.name)
                        tar = tarfile.open(os.path.join(dest, file.name))
                        tar.extractall(path=mothurdest)
                        tar.close()
                    elif file.name.endswith('.zip'):
                        handle_uploaded_file(file, dest, file.name)
                        zip = zipfile.ZipFile(os.path.join(dest, file.name))
                        zip.extractall(mothurdest)
                        zip.close()
                    else:
                        handle_uploaded_file(file, mothurdest, file.name)
                    handle_uploaded_file(file, dest, file.name)

                    if stopList[PID] == RID:
                        remove_proj(dest)
                        transaction.savepoint_rollback(sid)
                        return upStop(request)

                batch = 'mothur.batch'
                file15 = request.FILES['docfile15']

                avail_proc = mp.cpu_count()
                use_proc = min(avail_proc, processors)
                actual_proc = 'processors=' + str(use_proc)

                handle_uploaded_file(file15, mothurdest, batch)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                for line in fileinput.input('mothur/temp/mothur.batch', inplace=1):
                    print line.replace("processors=X", actual_proc),

                handle_uploaded_file(file15, dest, batch)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                try:
                    mothur(dest, source)
                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)

                    if request.user.is_superuser:
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
                    elif request.user.is_authenticated():
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
                    return render_to_response(
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "There was an error with your mothur batch file: " + str(file15.name)
                         },
                        context_instance=RequestContext(request)
                    )

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                try:
                    with open('% s/final.taxonomy' % dest, 'rb') as file3:
                        parse_taxonomy(file3)
                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)

                    if request.user.is_superuser:
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
                    elif request.user.is_authenticated():
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
                    return render_to_response(
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "There was an error parsing taxonomy file: final.taxonomy"},
                        context_instance=RequestContext(request)
                    )

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

                try:
                    with open('% s/final.taxonomy' % dest, 'rb') as file3:
                        with open('% s/final.shared' % dest, 'rb') as file4:
                            parse_profile(file3, file4, p_uuid, refDict)
                            end = datetime.datetime.now()
                    print 'Total time for upload:', end-start
                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)

                    if request.user.is_superuser:
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
                    elif request.user.is_authenticated():
                        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
                    return render_to_response(
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "There was an error parsing your shared file: final.shared"},
                        context_instance=RequestContext(request)
                    )

                if stopList[PID] == RID:
                    remove_proj(dest)
                    transaction.savepoint_rollback(sid)
                    return upStop(request)

            else:
                print ('Please check that all necessary files have been selected.')

    elif request.method == 'POST' and 'clickMe' in request.POST:
        remove_list(request)

    if request.user.is_superuser:
        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
    elif request.user.is_authenticated():
        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)

    return render_to_response(
        'upload.html',
        {'projects': projects,
         'form1': UploadForm1,
         'form2': UploadForm2,
         'error': ""},
        context_instance=RequestContext(request)
    )


def projectTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        if selSamples:
            if request.user.is_superuser:
                qs = Sample.objects.filter(sampleid__in=selSamples)
            elif request.user.is_authenticated():
                path_list = Sample.objects.filter(refid__author=request.user).values_list('sampleid', flat=True)
                qs = Sample.objects.filter( Q(sampleid__in=path_list) | Q(projectid__status='public') ).filter(sampleid__in=selSamples)
            else:
                qs = Sample.objects.none()
        else:
            qs = Sample.objects.none()

        results = {}
        qs1 = qs.values_list(
            "projectid",
            "projectid__status",
            "projectid__project_name",
            "sample_name",
            "projectid__project_desc",
            "projectid__start_date",
            "projectid__end_date",
            "projectid__pi_last",
            "projectid__pi_first",
            "projectid__pi_affiliation",
            "projectid__pi_email",
            "projectid__pi_phone"
        )

        qs1 = [[u'nan' if x is None else x for x in c] for c in qs1]
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


def sampleTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        if selSamples:
            if request.user.is_superuser:
                qs = Sample.objects.filter(sampleid__in=selSamples)
            elif request.user.is_authenticated():
                path_list = Sample.objects.filter(refid__author=request.user).values_list('sampleid', flat=True)
                qs = Sample.objects.filter( Q(sampleid__in=path_list) | Q(projectid__status='public') ).filter(sampleid__in=selSamples)
            else:
                qs = Sample.objects.none()
        else:
            qs = Sample.objects.none()

        results = {}
        qs1 = qs.values_list(
            "projectid__project_name",
            "refid",
            "sampleid",
            "sample_name",
            "organism",
            "collection_date",
            "depth",
            "elev",
            "seq_platform",
            "seq_gene",
            "seq_gene_region",
            "seq_barcode",
            "seq_for_primer",
            "seq_rev_primer",
            "env_biome",
            "env_feature",
            "env_material",
            "geo_loc_country",
            "geo_loc_state",
            "geo_loc_city",
            "geo_loc_farm",
            "geo_loc_plot",
            "latitude",
            "longitude",
            "annual_season_precpt",
            "annual_season_temp"
        )

        qs1 = [[u'nan' if x is None else x for x in c] for c in qs1]
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


def referenceTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        if selSamples:
            if request.user.is_superuser:
                qs = Sample.objects.filter(sampleid__in=selSamples)
            elif request.user.is_authenticated():
                path_list = Sample.objects.filter(refid__author=request.user).values_list('sampleid', flat=True)
                qs = Sample.objects.filter( Q(sampleid__in=path_list) | Q(projectid__status='public') ).filter(sampleid__in=selSamples)
            else:
                qs = Sample.objects.none()
        else:
            qs = Sample.objects.none()

        results = {}
        qs1 = qs.values_list(
            "projectid__project_name",
            "projectid",
            "refid",
            "sample_name",
            "sampleid",
            "refid__raw",
            "refid__path",
            "refid__source",
            "refid__alignDB",
            "refid__templateDB",
            "refid__taxonomyDB",
            "refid__author_id__username"
        )

        qs1 = [[u'nan' if x is None else x for x in c] for c in qs1]
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


def airTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        if selSamples:
            if request.user.is_superuser:
                qs = Air.objects.filter(sampleid__in=selSamples)
            elif request.user.is_authenticated():
                path_list = Air.objects.filter(refid__author=request.user).values_list('sampleid', flat=True)
                qs = Air.objects.filter( Q(sampleid__in=path_list) | Q(projectid__status='public') ).filter(sampleid__in=selSamples)
            else:
                qs = Air.objects.none()
        else:
            qs = Air.objects.none()

        results = {}
        qs1 = qs.values_list(
            "projectid__project_name",
            "sampleid__sample_name",
            "barometric_press",
            "carb_dioxide",
            "carb_monoxide",
            "chem_admin_term",
            "chem_admin_time",
            "elev",
            "humidity",
            "methane",
            "organism_type",
            "organism_count",
            "oxy_stat_samp",
            "oxygen",
            "perturbation_type",
            "perturbation_interval",
            "pollutants_type",
            "pollutants_concentration",
            "rel_to_oxygen",
            "resp_part_matter_substance",
            "resp_part_matter_concentration",
            "samp_collect_device",
            "samp_mat_process",
            "samp_salinity",
            "samp_size",
            "samp_store_dur",
            "samp_store_loc",
            "samp_store_temp",
            "solar_irradiance",
            "temp",
            "ventilation_rate",
            "ventilation_type",
            "volatile_org_comp_name",
            "volatile_org_comp_concentration",
            "wind_direction",
            "wind_speed",
        )

        qs1 = [[u'nan' if x is None else x for x in c] for c in qs1]
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


def associatedTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        if selSamples:
            if request.user.is_superuser:
                qs = Human_Associated.objects.filter(sampleid__in=selSamples)
            elif request.user.is_authenticated():
                path_list = Human_Associated.objects.filter(refid__author=request.user).values_list('sampleid', flat=True)
                qs = Human_Associated.objects.filter( Q(sampleid__in=path_list) | Q(projectid__status='public') ).filter(sampleid__in=selSamples)
            else:
                qs = Human_Associated.objects.none()
        else:
            qs = Human_Associated.objects.none()

        results = {}
        qs1 = qs.values_list(
            "projectid__project_name",
            "sampleid__sample_name",
            "samp_collect_device",
            "samp_mat_process",
            "samp_size",
            "samp_store_temp",
            "samp_store_dur",
            "samp_type",
            "samp_location",
            "samp_temp",
            "samp_ph",
            "samp_oxy_stat",
            "samp_salinity",
            "host_subject_id",
            "host_age",
            "host_pulse",
            "host_gender",
            "host_ethnicity",
            "host_height",
            "host_weight",
            "host_bmi",
            "host_weight_loss_3_month",
            "host_body_temp",
            "host_occupation",
            "pet_farm_animal",
            "smoker",
            "diet_type",
            "diet_duration",
            "diet_frequency",
            "diet_last_six_month",
            "last_meal",
            "medic_hist_perform",
            "disease_type",
            "disease_location",
            "disease_duration",
            "organism_count",
            "tumor_location",
            "tumor_mass",
            "tumor_stage",
            "drug_usage",
            "drug_type",
            "drug_duration",
            "drug_frequency",
            "perturbation",
            "pert_type",
            "pert_duration",
            "pert_frequency",
            "fetal_health_stat",
            "amniotic_fluid_color",
            "gestation_stat",
            "maternal_health_stat"
        )

        qs1 = [[u'nan' if x is None else x for x in c] for c in qs1]
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


def microbialTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        if selSamples:
            if request.user.is_superuser:
                qs = Microbial.objects.filter(sampleid__in=selSamples)
            elif request.user.is_authenticated():
                path_list = Microbial.objects.filter(refid__author=request.user).values_list('sampleid', flat=True)
                qs = Microbial.objects.filter( Q(sampleid__in=path_list) | Q(projectid__status='public') ).filter(sampleid__in=selSamples)
            else:
                qs = Microbial.objects.none()
        else:
            qs = Microbial.objects.none()

        results = {}
        qs1 = qs.values_list(
            "projectid__project_name",
            "sampleid__sample_name",
            "alkalinity",
            "alkyl_diethers",
            "altitude",
            "aminopept_act",
            "ammonium",
            "bacteria_carb_prod",
            "biomass_part_name",
            "biomass_amount",
            "bishomohopanol",
            "bromide",
            "calcium",
            "carb_nitro_ratio",
            "chem_administration_term",
            "chem_administration_time",
            "chloride",
            "chlorophyll",
            "diether_lipids_name",
            "diether_lipids_concentration",
            "diss_carb_dioxide",
            "diss_hydrogen",
            "diss_inorg_carb",
            "diss_org_carb",
            "diss_org_nitro",
            "diss_oxygen",
            "glucosidase_act",
            "magnesium",
            "mean_frict_vel",
            "mean_peak_frict_vel",
            "methane",
            "n_alkanes_name",
            "n_alkanes_concentration",
            "nitrate",
            "nitrite",
            "nitro",
            "org_carb",
            "org_matter",
            "org_nitro",
            "organism_name",
            "organism_count",
            "oxy_stat_samp",
            "part_org_carb",
            "perturbation_type",
            "perturbation_interval",
            "petroleum_hydrocarb",
            "ph",
            "phaeopigments_type",
            "phaeopigments_concentration",
            "phosphate",
            "phosplipid_fatt_acid_name",
            "phosplipid_fatt_acid_concentration",
            "potassium",
            "pressure",
            "redox_potential",
            "rel_to_oxygen",
            "salinity",
            "samp_collect_device",
            "samp_mat_process",
            "samp_size",
            "samp_store_dur",
            "samp_store_loc",
            "samp_store_temp",
            "silicate",
            "sodium",
            "sulfate",
            "sulfide",
            "temp",
            "tot_carb",
            "tot_nitro",
            "tot_org_carb",
            "turbidity",
            "water_content"
        )

        qs1 = [[u'nan' if x is None else x for x in c] for c in qs1]
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


def soilTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        if selSamples:
            if request.user.is_superuser:
                qs = Soil.objects.filter(sampleid__in=selSamples)
            elif request.user.is_authenticated():
                path_list = Soil.objects.filter(refid__author=request.user).values_list('sampleid', flat=True)
                qs = Soil.objects.filter( Q(sampleid__in=path_list) | Q(projectid__status='public') ).filter(sampleid__in=selSamples)
            else:
                qs = Soil.objects.none()
        else:
            qs = Soil.objects.none()

        results = {}
        qs1 = qs.values_list(
            "projectid__project_name",
            "sampleid__sample_name",
            "samp_collection_device",
            "samp_size",
            "samp_depth",
            "samp_prep",
            "samp_sieve_size",
            "samp_store_dur",
            "samp_store_loc",
            "samp_store_temp",
            "samp_weight_dna_ext",
            "pool_dna_extracts",
            "fao_class",
            "local_class",
            "texture_class",
            "porosity",
            "profile_position",
            "slope_aspect",
            "slope_gradient",
            "bulk_density",
            "drainage_class",
            "water_content_soil",
            "cur_land_use",
            "cur_vegetation",
            "cur_crop",
            "cur_cultivar",
            "crop_rotation",
            "cover_crop",
            "fert_amendment_class",
            "fert_placement",
            "fert_type",
            "fert_tot_amount",
            "fert_N_tot_amount",
            "fert_P_tot_amount",
            "fert_K_tot_amount",
            "irrigation_type",
            "irrigation_tot_amount",
            "residue_removal",
            "residue_growth_stage",
            "residue_removal_percent",
            "tillage_event",
            "tillage_event_depth",
            "amend1_class",
            "amend1_active_ingredient",
            "amend1_tot_amount",
            "amend2_class",
            "amend2_active_ingredient",
            "amend2_tot_amount",
            "amend3_class",
            "amend3_active_ingredient",
            "amend3_tot_amount",
            "rRNA_copies",
            "microbial_biomass_C",
            "microbial_biomass_N",
            "microbial_respiration",
            "soil_pH",
            "soil_EC",
            "soil_C",
            "soil_OM",
            "soil_N",
            "soil_NO3_N",
            "soil_NH4_N",
            "soil_P",
            "soil_K",
            "soil_S",
            "soil_Zn",
            "soil_Fe",
            "soil_Cu",
            "soil_Mn",
            "soil_Ca",
            "soil_Mg",
            "soil_Na",
            "soil_B",
            "plant_C",
            "plant_N",
            "plant_P",
            "plant_K",
            "plant_Ca",
            "plant_Mg",
            "plant_S",
            "plant_Na",
            "plant_Cl",
            "plant_Al",
            "plant_B",
            "plant_Cu",
            "plant_Fe",
            "plant_Mn",
            "plant_Zn",
            "crop_tot_biomass_fw",
            "crop_tot_biomass_dw",
            "crop_tot_above_biomass_fw",
            "crop_tot_above_biomass_dw",
            "crop_tot_below_biomass_fw",
            "crop_tot_below_biomass_dw",
            "harv_fraction",
            "harv_fresh_weight",
            "harv_dry_weight",
            "ghg_chamber_placement",
            "ghg_N2O",
            "ghg_CO2",
            "ghg_NH4",
            "soil_water_cap",
            "soil_surf_hard",
            "soil_subsurf_hard",
            "soil_agg_stability",
            "soil_ACE_protein",
            "soil_active_C"
        )

        qs1 = [[u'nan' if x is None else x for x in c] for c in qs1]
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


def waterTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        if selSamples:
            if request.user.is_superuser:
                qs = Water.objects.filter(sampleid__in=selSamples)
            elif request.user.is_authenticated():
                path_list = Water.objects.filter(refid__author=request.user).values_list('sampleid', flat=True)
                qs = Water.objects.filter( Q(sampleid__in=path_list) | Q(projectid__status='public') ).filter(sampleid__in=selSamples)
            else:
                qs = Water.objects.none()
        else:
            qs = Water.objects.none()

        results = {}
        qs1 = qs.values_list(
            "projectid__project_name",
            "sampleid__sample_name",
            "alkalinity",
            "alkyl_diethers",
            "altitude",
            "aminopept_act",
            "ammonium",
            "atmospheric_data",
            "bac_prod",
            "bac_resp",
            "bacteria_carb_prod",
            "biomass_part_name",
            "biomass_amount",
            "bishomohopanol",
            "bromide",
            "calcium",
            "carb_nitro_ratio",
            "chem_administration_name",
            "chem_administration_time",
            "chloride",
            "chlorophyll",
            "conduc",
            "density",
            "diether_lipids",
            "diss_carb_dioxide",
            "diss_hydrogen",
            "diss_inorg_carb",
            "diss_inorg_nitro",
            "diss_inorg_phosp",
            "diss_org_carb",
            "diss_org_nitro",
            "diss_oxygen",
            "down_par",
            "elev",
            "fluor",
            "glucosidase_act",
            "light_intensity",
            "magnesium",
            "mean_frict_vel",
            "mean_peak_frict_vel",
            "n_alkanes",
            "nitrate",
            "nitrite",
            "nitro",
            "org_carb",
            "org_matter",
            "org_nitro",
            "organism_name",
            "organism_count",
            "oxy_stat_samp",
            "part_org_carb",
            "part_org_nitro",
            "perturbation_type",
            "perturbation_interval",
            "pretroleum_hydrocarb",
            "ph",
            "phaeopigments",
            "phosphate",
            "phosplipid_fatt_acid",
            "photon_flux",
            "potassium",
            "pressure",
            "primary_prod",
            "redox_potential",
            "rel_to_oxygen",
            "samp_mat_process",
            "samp_salinity",
            "samp_size",
            "samp_store_dur",
            "samp_store_loc",
            "samp_store_temp",
            "samp_vol_we_dna_ext",
            "silicate",
            "sodium",
            "soluble_react_phosp",
            "source_material_id",
            "sulfate",
            "sulfide",
            "suspen_part_matter",
            "temp",
            "tidal_stage",
            "tot_depth_water_col",
            "tot_diss_nitro",
            "tot_inorg_nitro",
            "tot_nitro",
            "tot_part_carb",
            "tot_phosp",
            "water_current_direction",
            "water_current_magnitude"
        )

        qs1 = [[u'nan' if x is None else x for x in c] for c in qs1]
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


def userTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        if selSamples:
            if request.user.is_superuser:
                qs = UserDefined.objects.filter(sampleid__in=selSamples)
            elif request.user.is_authenticated():
                path_list = UserDefined.objects.filter(refid__author=request.user).values_list('sampleid', flat=True)
                qs = UserDefined.objects.filter( Q(sampleid__in=path_list) | Q(projectid__status='public') ).filter(sampleid__in=selSamples)
            else:
                qs = UserDefined.objects.none()
        else:
            qs = UserDefined.objects.none()

        results = {}
        qs1 = qs.values_list(
            "projectid__project_name",
            "sampleid__sample_name",
            "usr_cat1",
            "usr_cat2",
            "usr_cat3",
            "usr_cat4",
            "usr_cat5",
            "usr_cat6",
            "usr_quant1",
            "usr_quant2",
            "usr_quant3",
            "usr_quant4",
            "usr_quant5",
            "usr_quant6"
        )

        qs1 = [[u'nan' if x is None else x for x in c] for c in qs1]
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


@login_required(login_url='/myPhyloDB/accounts/login/')
def select(request):
    if request.method == 'POST':
        uploaded = request.FILES["normFile"]
        savedDF = pd.read_csv(uploaded, index_col=0, sep='\t')
        selList = list(set(savedDF['sampleid']))

        # save selected samples file to users temp/ folder
        myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
        path = str(myDir) + 'usr_sel_samples.pkl'

        if not os.path.exists(myDir):
            os.makedirs(myDir)

        with open(path, 'wb') as f:
            pickle.dump(selList, f)

        # save normalized file to users temp/ folder
        myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
        path = str(myDir) + 'usr_norm_data.csv'

        if not os.path.exists(myDir):
            os.makedirs(myDir)
        savedDF.to_csv(path, sep='\t')

        projectList = list(set(savedDF['projectid']))

        if not Project.objects.filter(projectid__in=projectList).exists():
            return render_to_response(
                'select.html',
                {'form9': UploadForm9,
                 'selList': '',
                 'normpost': 'error',
                 'method': 'POST'},
                 context_instance=RequestContext(request)
            )

        if not Reference.objects.filter(projectid__in=projectList).exists():
            return render_to_response(
                'select.html',
                {'form9': UploadForm9,
                 'selList': '',
                 'normpost': 'error',
                 'method': 'POST'},
                context_instance=RequestContext(request)
            )

        if not Sample.objects.filter(sampleid__in=selList).exists():
            return render_to_response(
                'select.html',
                {'form9': UploadForm9,
                 'selList': '',
                 'normpost': 'error',
                 'method': 'POST'},
                context_instance=RequestContext(request)
            )

        biome = {}
        tempDF = savedDF.drop_duplicates(subset='sampleid', take_last=True)
        tempDF.set_index('sampleid', drop=True, inplace=True)

        metaDF = tempDF.drop(['kingdomName', 'phylaName', 'className', 'orderName', 'familyName', 'genusName', 'speciesName', 'otuName', 'otuid', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity'], axis=1)
        myList = list(tempDF.index.values)

        nameList = []
        for i in myList:
            nameList.append({"id": str(i), "metadata": metaDF.loc[i].to_dict()})

        # get list of lists with abundances
        taxaOnlyDF = savedDF.loc[:, ['sampleid', 'kingdomName', 'phylaName', 'className', 'orderName', 'familyName', 'genusName', 'speciesName', 'otuName', 'otuid', 'abund']]
        taxaOnlyDF = taxaOnlyDF.pivot(index='otuid', columns='sampleid', values='abund')
        dataList = taxaOnlyDF.values.tolist()

        # get list of taxa
        namesDF = savedDF.loc[:, ['sampleid', 'otuid']]
        namesDF['taxa'] = savedDF.loc[:, ['kingdomName', 'phylaName', 'className', 'orderName', 'familyName', 'genusName', 'speciesName', 'otuName']].values.tolist()
        namesDF = namesDF.pivot(index='otuid', columns='sampleid', values='taxa')

        taxaList = []
        for index, row in namesDF.iterrows():
            metaDict = {}
            metaDict['taxonomy'] = row[0]
            taxaList.append({"id": index, "metadata": metaDict})

        biome['format'] = 'Biological Observation Matrix 0.9.1-dev'
        biome['format_url'] = 'http://biom-format.org/documentation/format_versions/biom-1.0.html'
        biome['type'] = 'OTU table'
        biome['generated_by'] = 'myPhyloDB'
        biome['date'] = str(datetime.datetime.now())
        biome['matrix_type'] = 'dense'
        biome['matrix_element_type'] = 'float'
        biome['rows'] = taxaList
        biome['columns'] = nameList
        biome['data'] = dataList

        myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
        path = str(myDir) + 'usr_norm_data.biom'
        with open(path, 'w') as outfile:
            simplejson.dump(biome, outfile, ensure_ascii=True, indent=4, sort_keys=True)

        return render_to_response(
            'select.html',
            {'form9': UploadForm9,
             'selList': selList,
             'normpost': 'Success',
             'method': 'POST'},
            context_instance=RequestContext(request)
        )

    else:
        return render_to_response(
            'select.html',
            {'form9': UploadForm9,
             'selList': '',
             'normpost': '',
             'method': 'GET'},
            context_instance=RequestContext(request))


def taxaJSON(request):
    results = {}
    qs1 = OTU_99.objects.values_list('kingdomid', 'kingdomid__kingdomName', 'phylaid', 'phylaid__phylaName', 'classid', 'classid__className', 'orderid', 'orderid__orderName', 'familyid', 'familyid__familyName', 'genusid', 'genusid__genusName', 'speciesid', 'speciesid__speciesName', 'otuid', 'otuName')
    results['data'] = list(qs1)
    myJson = simplejson.dumps(results, ensure_ascii=False)
    return HttpResponse(myJson)


def taxa(request):
    return render_to_response(
        'taxa.html',
        context_instance=RequestContext(request)
    )


def pathJSON(request):
    results = {}
    qs1 = ko_entry.objects.using('picrust').values_list('ko_lvl1_id', 'ko_lvl1_id__ko_lvl1_name', 'ko_lvl2_id', 'ko_lvl2_id__ko_lvl2_name', 'ko_lvl3_id', 'ko_lvl3_id__ko_lvl3_name', 'ko_lvl4_id', 'ko_orthology', 'ko_name', 'ko_desc')
    results['data'] = list(qs1)
    myJson = simplejson.dumps(results, ensure_ascii=False)
    return HttpResponse(myJson)


def pathTaxaJSON(request):
    if request.is_ajax():
        wanted = request.GET['key']
        wanted = wanted.replace('"', '')
        koList = []

        if ko_lvl1.objects.using('picrust').filter(ko_lvl1_id=wanted).exists():
            record = ko_lvl1.objects.using('picrust').get(ko_lvl1_id=wanted)
            koList = record.ko_entry_set.values_list('ko_orthology', flat=True)
        elif ko_lvl2.objects.using('picrust').filter(ko_lvl2_id=wanted).exists():
            record = ko_lvl2.objects.using('picrust').get(ko_lvl2_id=wanted)
            koList = record.ko_entry_set.values_list('ko_orthology', flat=True)
        elif ko_lvl3.objects.using('picrust').filter(ko_lvl3_id=wanted).exists():
            record = ko_lvl3.objects.using('picrust').get(ko_lvl3_id=wanted)
            koList = record.ko_entry_set.values_list('ko_orthology', flat=True)
        elif ko_entry.objects.using('picrust').filter(ko_lvl4_id=wanted).exists():
            koList = ko_entry.objects.using('picrust').filter(ko_lvl4_id=wanted).values_list('ko_orthology', flat=True)

        finalotuList = []
        if koList:
            qs = PICRUSt.objects.using('picrust')
            picrustDF = read_frame(qs, fieldnames=['otuid__otuid', 'geneCount'], verbose=False)
            picrustDF.rename(columns={'otuid__otuid': 'otuid'}, inplace=True)
            otuList = picrustDF['otuid'].tolist()
            picrustDF.set_index('otuid', drop=True, inplace=True)
            picrustDF['gene'] = 'No'

            finalotuList = []
            for otu in otuList:
                cell = picrustDF.at[otu, 'geneCount']
                d = ast.literal_eval(cell)
                present = 'No'
                for key, value in d.items():
                    if key in koList and value > 0:
                        present = 'Yes'
                        break
                if present == 'Yes':
                    finalotuList.append(otu)

        results = {}
        qs1 = OTU_99.objects.filter(otuid__in=finalotuList).values_list('kingdomid', 'kingdomid__kingdomName', 'phylaid', 'phylaid__phylaName', 'classid', 'classid__className', 'orderid', 'orderid__orderName', 'familyid', 'familyid__familyName', 'genusid', 'genusid__genusName', 'speciesid', 'speciesid__speciesName', 'otuName', 'otuid')
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


def kegg_path(request):
    return render_to_response(
        'kegg_path.html',
        context_instance=RequestContext(request)
    )


def nzJSON(request):
    qs1 = nz_entry.objects.using('picrust').values_list('nz_lvl1_id','nz_lvl1_id__nz_lvl1_name', 'nz_lvl2_id', 'nz_lvl2_id__nz_lvl2_name', 'nz_lvl3_id', 'nz_lvl3_id__nz_lvl3_name', 'nz_lvl4_id', 'nz_lvl4_id__nz_lvl4_name', 'nz_lvl5_id', 'nz_orthology', 'nz_name', 'nz_desc')
    results = {}
    results['data'] = list(qs1)
    myJson = simplejson.dumps(results, ensure_ascii=False)
    return HttpResponse(myJson)


def nzTaxaJSON(request):
    if request.is_ajax():
        wanted = request.GET['key']
        wanted = wanted.replace('"', '')
        koList = []

        if nz_lvl1.objects.using('picrust').filter(nz_lvl1_id=wanted).exists():
            record = nz_lvl1.objects.using('picrust').get(nz_lvl1_id=wanted)
            koList = record.nz_entry_set.values_list('nz_orthology', flat=True)
        elif nz_lvl2.objects.using('picrust').filter(nz_lvl2_id=wanted).exists():
            record = nz_lvl2.objects.using('picrust').get(nz_lvl2_id=wanted)
            koList = record.nz_entry_set.values_list('nz_orthology', flat=True)
        elif nz_lvl3.objects.using('picrust').filter(nz_lvl3_id=wanted).exists():
            record = nz_lvl3.objects.using('picrust').get(nz_lvl3_id=wanted)
            koList = record.nz_entry_set.values_list('nz_orthology', flat=True)
        elif nz_lvl4.objects.using('picrust').filter(nz_lvl4_id=wanted).exists():
            record = nz_lvl4.objects.using('picrust').get(nz_lvl4_id=wanted)
            koList = record.nz_entry_set.values_list('nz_orthology', flat=True)
        elif nz_entry.objects.using('picrust').filter(nz_lvl5_id=wanted).exists():
            koList = nz_entry.objects.using('picrust').filter(nz_lvl5_id=wanted).values_list('nz_orthology', flat=True)

        finalotuList = []
        if koList:
            qs = PICRUSt.objects.using('picrust')
            picrustDF = read_frame(qs, fieldnames=['otuid__otuid', 'geneCount'], verbose=False)
            picrustDF.rename(columns={'otuid__otuid': 'otuid'}, inplace=True)
            otuList = picrustDF['otuid'].tolist()
            picrustDF.set_index('otuid', drop=True, inplace=True)
            picrustDF['gene'] = 'No'

            finalotuList = []
            for otu in otuList:
                cell = picrustDF.at[otu, 'geneCount']
                d = ast.literal_eval(cell)
                present = 'No'
                for key, value in d.items():
                    if key in koList and value > 0:
                        present = 'Yes'
                        break
                if present == 'Yes':
                    finalotuList.append(otu)

        results = {}
        qs1 = OTU_99.objects.filter(otuid__in=finalotuList).values_list('kingdomid', 'kingdomid__kingdomName', 'phylaid', 'phylaid__phylaName', 'classid', 'classid__className', 'orderid', 'orderid__orderName', 'familyid', 'familyid__familyName', 'genusid', 'genusid__genusName', 'speciesid', 'speciesid__speciesName', 'otuid', 'otuName')
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


def kegg_enzyme(request):
    return render_to_response(
        'kegg_enzyme.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/myPhyloDB/accounts/login/')
def norm(request):
    cleanup('myPhyloDB/media/temp/norm')

    return render_to_response(
        'norm.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/myPhyloDB/accounts/login/')
def ANOVA(request):
    cleanup('myPhyloDB/media/temp/anova')

    return render_to_response(
        'anova.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/myPhyloDB/accounts/login/')
def rich(request):
    cleanup('myPhyloDB/media/temp/spac')

    return render_to_response(
        'SpAC.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/myPhyloDB/accounts/login/')
def soil_index(request):
    cleanup('myPhyloDB/media/temp/soil_index')

    return render_to_response(
        'soil_index.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/myPhyloDB/accounts/login/')
def DiffAbund(request):
    cleanup('myPhyloDB/media/temp/diffabund')

    return render_to_response(
        'diff_abund.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/myPhyloDB/accounts/login/')
def GAGE(request):
    cleanup('myPhyloDB/media/temp/gage')

    return render_to_response(
        'gage.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/myPhyloDB/accounts/login/')
def PCA(request):
    cleanup('myPhyloDB/media/temp/pca')

    return render_to_response(
        'pca.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/myPhyloDB/accounts/login/')
def PCoA(request):
    cleanup('myPhyloDB/media/temp/pcoa')

    return render_to_response(
        'pcoa.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/myPhyloDB/accounts/login/')
def RF(request):
    cleanup('myPhyloDB/media/temp/rf')

    return render_to_response(
        'rf.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/myPhyloDB/accounts/login/')
def SPLS(request):
    cleanup('myPhyloDB/media/temp/spls')

    return render_to_response(
        'spls.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/myPhyloDB/accounts/login/')
def WGCNA(request):
    cleanup('myPhyloDB/media/temp/wgcna')

    return render_to_response(
        'wgcna.html',
        context_instance=RequestContext(request)
    )


def saveSampleList(request):
    if request.is_ajax():
        allJson = request.GET["all"]
        selList = simplejson.loads(allJson)

        # remove old files
        try:
            dirpath = 'myPhyloDB/media/usr_temp/' + str(request.user)
            shutil.rmtree(dirpath, ignore_errors=True)
        except Exception as e:
            print e

        # save file
        myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
        path = str(myDir) + 'usr_sel_samples.pkl'

        if not os.path.exists(myDir):
            os.makedirs(myDir)

        with open(path, 'wb') as f:
            pickle.dump(selList, f)

        thisUser = User.objects.filter(username=request.user.username)
        if UserProfile.objects.filter(user=thisUser).exists():
            myData = UserProfile.objects.get(user=thisUser)
            myData.dataID = uuid4().hex
            myData.save()
        else:
            UserProfile.objects.create(user=request.user, dataID=uuid4().hex)

        text = 'Selected sample(s) have been recorded!'
        res = simplejson.dumps(text, encoding="Latin-1")
        return HttpResponse(res, content_type='application/json')


@login_required(login_url='/myPhyloDB/accounts/login/')
def reprocess(request):
    try:
        alignFile = request.FILES['docfile8']
        alignDB = request.FILES['docfile8'].name
        handle_uploaded_file(alignFile, 'mothur/reference/align', alignDB)
    except:
        pass

    try:
        templateFile = request.FILES['docfile9']
        templateDB = request.FILES['docfile9'].name
        handle_uploaded_file(templateFile, 'mothur/reference/template', templateDB)
    except:
        pass

    try:
        taxonomyFile = request.FILES['docfile10']
        taxonomyDB = request.FILES['docfile10'].name
        handle_uploaded_file(taxonomyFile, 'mothur/reference/taxonomy', taxonomyDB)
    except:
        pass

    alignDB = sorted(os.listdir('mothur/reference/align/'))
    templateDB = sorted(os.listdir('mothur/reference/template/'))
    taxonomyDB = sorted(os.listdir('mothur/reference/taxonomy/'))

    return render_to_response(
        'reprocess.html',
        {'form4': UploadForm4,
         'alignDB': alignDB,
         'templateDB': templateDB,
         'taxonomyDB': taxonomyDB},
        context_instance=RequestContext(request)
    )


@login_required(login_url='/myPhyloDB/accounts/login/')
def update(request):
    state = ''
    return render_to_response(
        'update.html',
        {'form5': UploadForm5,
         'state': state},
        context_instance=RequestContext(request)
    )


def updaStop(request):
    state = "Update stopped"
    return render_to_response(
        'update.html',
        {'form5': UploadForm5,
         'state': state},
        context_instance=RequestContext(request)
    )


def updateFunc(request, stopList):
    form5 = UploadForm5(request.POST, request.FILES)
    state = ''
    curUser = User.objects.get(username=request.user.username)

    PID = 0
    RID = request.POST['RID']

    if stopList[PID] == RID:
        return updaStop(request)

    if form5.is_valid():
        refid = request.POST['refid']
        file1 = request.FILES['docfile11']

        if stopList[PID] == RID:
            return updaStop(request)

        ref = Reference.objects.get(refid=refid)
        dest = ref.path
        p_uuid, pType, num_samp = projectid(file1)
        if not num_samp:
            return render_to_response(
                'update.html',
                {'form5': UploadForm5,
                 'state': "Error: Please open and save your meta file in Excel/OpenOffice in order to perform the necessary calculations..."
                 },
                context_instance=RequestContext(request)
            )

        if stopList[PID] == RID:
            return updaStop(request)

        try:
            metaFile = '/'.join([dest, file1.name])
            handle_uploaded_file(file1, dest, file1.name)
            parse_project(metaFile, p_uuid, curUser)
        except Exception as e:
            state = "There was an error parsing your metafile: " + str(file1.name) + "\nError info: "+str(e)
            return render_to_response(
                'update.html',
                {'form5': UploadForm5,
                 'state': state},
                context_instance=RequestContext(request)
            )

        if stopList[PID] == RID:
            return updaStop(request)

        try:
            bat = 'mothur.batch'
            batPath = '/'.join([dest, bat])
            reference = Reference.objects.get(refid=refid)
            raw = reference.raw
            source = reference.source
            userID = str(request.user.id)

            if os.path.exists(batPath):
                with open(batPath, 'rb') as batFile:
                    parse_sample(metaFile, p_uuid, pType, num_samp, dest, batFile, raw, source, userID)
            else:
                batFile = 'you do not really need me'
                parse_sample(metaFile, p_uuid, pType, num_samp, dest, batFile, raw, source, userID)

        except Exception as e:
            state = "There was an error parsing your metafile: " + str(file1.name) + "\nError info: "+str(e)
            return render_to_response(
                'update.html',
                {'form5': UploadForm5,
                 'state': state},
                context_instance=RequestContext(request)
            )

        if stopList[PID] == RID:
            return updaStop(request)

        state = 'Path: ' + str(dest) + ' is finished parsing!'

    if stopList[PID] == RID:
        return updaStop(request)

    return render_to_response(
        'update.html',
        {'form5': UploadForm5,
         'state': state},
        context_instance=RequestContext(request)
    )


@login_required(login_url='/myPhyloDB/accounts/login/')
@user_passes_test(lambda u: u.is_superuser)
def pybake(request):
    # split up for queue
    form6 = UploadForm6(request.POST, request.FILES)
    form7 = UploadForm7(request.POST, request.FILES)
    form8 = UploadForm8(request.POST, request.FILES)

    if form6.is_valid():
        file1 = request.FILES['taxonomy']
        file2 = request.FILES['precalc_16S']
        file3 = request.FILES['precalc_KEGG']
        geneParse(file1, file2, file3)

    if form7.is_valid():
        file4 = request.FILES['ko_htext']
        koParse(file4)

    if form8.is_valid():
        file5 = request.FILES['nz_htext']
        nzParse(file5)

    return render_to_response(
        'pybake.html',
        {'form6': UploadForm6,
         'form7': UploadForm7,
         'form8': UploadForm8},
        context_instance=RequestContext(request)
    )


def uploadNorm(request):
    uploaded = request.FILES["normFile"]
    savedDF = pd.read_csv(uploaded, index_col=0, sep='\t')
    selList = list(set(savedDF['sampleid']))

    # save selected samples file to users temp/ folder
    myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
    path = str(myDir) + 'usr_sel_samples.pkl'

    if not os.path.exists(myDir):
        os.makedirs(myDir)

    with open(path, 'wb') as f:
        pickle.dump(selList, f)

    # save normalized file to users temp/ folder
    myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
    path = str(myDir) + 'usr_norm_data.csv'

    if not os.path.exists(myDir):
        os.makedirs(myDir)
    savedDF.to_csv(path, sep='\t')

    projectList = list(set(savedDF['projectid']))

    if Project.objects.filter(projectid__in=projectList).exists():
        projects = Project.objects.filter(projectid__in=projectList)

    else:
        return render_to_response(
            'select.html',
            {'normpost': 'One or more projects were not found in your database!'},
            context_instance=RequestContext(request)
        )

    if Reference.objects.filter(projectid__in=projectList).exists():
        refs = Reference.objects.filter(projectid__in=projectList)
    else:
        return render_to_response(
            'select.html',
            {'normpost': 'One or more projects were not found in your database!'},
            context_instance=RequestContext(request)
        )

    if Sample.objects.filter(sampleid__in=selList).exists():
        samples = Sample.objects.filter(sampleid__in=selList)
    else:
        return render_to_response(
            'select.html',
            {'normpost': 'One or more samples were not found in your database!'},
            context_instance=RequestContext(request)
        )

    biome = {}
    tempDF = savedDF.drop_duplicates(subset='sampleid', take_last=True)
    tempDF.set_index('sampleid', drop=True, inplace=True)

    metaDF = tempDF.drop(['kingdomName', 'phylaName', 'className', 'orderName', 'familyName', 'genusName', 'speciesName', 'otuName', 'otuid', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity'], axis=1)
    myList = list(tempDF.index.values)

    nameList = []
    for i in myList:
        nameList.append({"id": str(i), "metadata": metaDF.loc[i].to_dict()})

    # get list of lists with abundances
    taxaOnlyDF = savedDF.loc[:, ['sampleid', 'kingdomName', 'phylaName', 'className', 'orderName', 'familyName', 'genusName', 'speciesName', 'otuName', 'otuid', 'abund']]
    taxaOnlyDF = taxaOnlyDF.pivot(index='otuid', columns='sampleid', values='abund')
    dataList = taxaOnlyDF.values.tolist()

    # get list of taxa
    namesDF = savedDF.loc[:, ['sampleid', 'otuid']]
    namesDF['taxa'] = savedDF.loc[:, ['kingdomName', 'phylaName', 'className', 'orderName', 'familyName', 'genusName', 'speciesName', 'otuName']].values.tolist()
    namesDF = namesDF.pivot(index='otuid', columns='sampleid', values='taxa')

    taxaList = []
    for index, row in namesDF.iterrows():
        metaDict = {}
        metaDict['taxonomy'] = row[0]
        taxaList.append({"id": index, "metadata": metaDict})

    biome['format'] = 'Biological Observation Matrix 0.9.1-dev'
    biome['format_url'] = 'http://biom-format.org/documentation/format_versions/biom-1.0.html'
    biome['type'] = 'OTU table'
    biome['generated_by'] = 'myPhyloDB'
    biome['date'] = str(datetime.datetime.now())
    biome['matrix_type'] = 'dense'
    biome['matrix_element_type'] = 'float'
    biome['rows'] = taxaList
    biome['columns'] = nameList
    biome['data'] = dataList

    myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
    path = str(myDir) + 'usr_norm_data.biom'
    with open(path, 'w') as outfile:
        simplejson.dump(biome, outfile, ensure_ascii=True, indent=4, sort_keys=True)

    return render_to_response(
        'select.html',
        {'form9': UploadForm9,
         'projects': projects,
         'refs': refs,
         'samples': samples,
         'selList': selList,
         'normpost': 'Success'},
        context_instance=RequestContext(request)
    )


# Function to create a user folder at login
def login_usr_callback(sender, user, request, **kwargs):
    user = request.user
    if not os.path.exists('myPhyloDB/media/usr_temp/'+str(user)):
        os.makedirs('myPhyloDB/media/usr_temp/'+str(user))
    cleanup('myPhyloDB/media/usr_temp/'+str(user))

user_logged_in.connect(login_usr_callback)  # JUMP does this belong in a function?


# Function has been added to context processor in setting file
def usrFiles(request):
    selFiles = os.path.exists('myPhyloDB/media/usr_temp/' + str(request.user) + '/usr_sel_samples.pkl')
    normFiles = os.path.exists('myPhyloDB/media/usr_temp/' + str(request.user) + '/usr_norm_data.h5')

    try:
        dataID = UserProfile.objects.get(user=request.user).dataID  # UUID from last selection
    except:
        dataID = "nodataidfound"

    return {
        'selFiles': selFiles,
        'normFiles': normFiles,
        'dataID': dataID
    }


def changeuser(request):
    return render_to_response(
        'changeuser.html',
        {"form": UserRegForm,
            "error": "none"},
        context_instance=RequestContext(request)
    )


def updateInfo(request):
    stuff = request.POST
    pword = stuff['password1']
    error = "none"
    verified = check_password(pword, request.user.password)
    if verified:
        thisUser = User.objects.get(username=request.user.username)
        myData = UserProfile.objects.get(user=thisUser)
        for key in stuff:
            thisUser.__setattr__(key, stuff[key])
            myData.__setattr__(key, stuff[key])
        thisUser.save()
        myData.save()
    else:
        error = "Incorrect password"


    return render_to_response(
        'changeuser.html',
        {"form": UserRegForm,
            "error": error},
        context_instance=RequestContext(request)
    )


def addPerms(request):
    if request.is_ajax():
        allJson = simplejson.loads(request.GET["all"])
        # get selected projects list (files? not samples as subsets though)
        selList = allJson['keys']
        # get list of names to add
        nameList = allJson['names']
        nameList = nameList.split(';')
        # get permission mode (view vs edit)  * view is redundant if project is public
        permLevel = allJson['mode']
        # permLevel 0 means view only
        # permLevel 1 means editing as well

        thisUser = User.objects.get(username=request.user.username)

        # loop through nameList, verify each exists
        errorList = []
        finalNameList = nameList[:]
        for name in nameList:
            if not User.objects.filter(username=name).exists():
                errorList.append(name)
                finalNameList.remove(name)
                # remove name from nameList

        # get project by id from selList  # JUMP
        for pid in selList:
            curProj = Project.objects.get(projectid=pid)
            # check if thisUser is owner of project (or super, check old edit perms basically)
            if thisUser == curProj.owner:
                # loop through nameList, check if name is present on permList already
                for curName in finalNameList:
                    viewList = curProj.whitelist_view.split(';')
                    found = False
                    for name in viewList:
                        if name == curName:
                            found = True
                    if not found:
                        # add current name to current project's viewPerms
                        curProj.whitelist_view += curName+";"
                    # check if mode grants edit perm as well
                    if permLevel:
                        editList = curProj.whitelist_edit.split(';')
                        found = False
                        for name in editList:
                            if name == curName:
                                found = True
                        if not found:
                            # add current name to current project's editPerms
                            curProj.whitelist_edit += curName+";"
            # save changes
            curProj.save()
        # return new profile page (boxes cleared, perhaps success alert message)

        # check if errorList has entries
        if len(errorList) > 0:
            # set return text to be an error reporting failed names
            text = "Name(s) not found: "
            iter = 1
            for name in errorList:
                if iter == len(errorList):
                    text += name
                else:
                    text += name + ", "
                iter += 1
        else:
            text = 'Users have been added to selected project(s)'
        retDict = {"error": text}
        res = simplejson.dumps(retDict, encoding="Latin-1")
        return HttpResponse(res, content_type='application/json')


def remPerms(request):
    print "NYI"
    return
