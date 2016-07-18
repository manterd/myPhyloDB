import datetime
from django.contrib.auth.decorators import login_required
from django.contrib.auth.signals import user_logged_in
from django.db.models import Q
from django.http import *
from django.shortcuts import render_to_response
from django.template import RequestContext
from django.views.decorators.csrf import csrf_exempt
import fileinput
import logging
import multiprocessing as mp
import os
import pandas as pd
import pickle
import simplejson
import xlrd
import time

from forms import UploadForm1, UploadForm2, UploadForm4, UploadForm5, \
    UploadForm6, UploadForm7, UploadForm8, UploadForm9
from models import Project, Reference, Sample, Species
from models import ko_entry, nz_entry
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


def upStop(request):
    # cleanup mid upload project!
    projects = Reference.objects.none()
    if request.user.is_superuser:
        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
    elif request.user.is_authenticated():
        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
    thing =render_to_response(
            'upload.html',
            {'projects': projects,
             'form1': UploadForm1,
             'form2': UploadForm2,
             'error': "Upload stopped"
             },
            context_instance=RequestContext(request)
        )
    return thing


def uploadFunc(request, stopList):
    projects = Reference.objects.none()
    if request.method == 'POST' and 'Upload' in request.POST:
        start = datetime.datetime.now()
        form1 = UploadForm1(request.POST, request.FILES)
        source = str(request.POST['source'])
        userID = str(request.User__id)
        processors = int(request.POST['processors'])
        RID = request.POST['RID']
        PID = 0  # change if adding additional data threads

        if stopList[PID] == RID:
            return upStop(request)

        if form1.is_valid():
            file1 = request.FILES['docfile1']
            try:
                p_uuid, pType, num_samp = projectid(file1)
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
            metaName = 'final_meta.xls'
            metaFile = '/'.join([dest, metaName])

            try:
                handle_uploaded_file(file1, dest, metaName)
                parse_project(metaFile, p_uuid)
            except:
                logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                logging.exception(myDate)

                remove_proj(dest)

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
                return upStop(request)

            try:
                refDict = parse_sample(metaFile, p_uuid, pType, num_samp, dest, batch, raw, source, userID)
            except:
                logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                logging.exception(myDate)

                remove_proj(dest)

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

            if stopList[PID] == RID:
                remove_proj(dest)
                return upStop(request)

            if source == 'mothur':
                file3 = request.FILES['docfile3']
                handle_uploaded_file(file3, dest, file3.name)

                try:
                    parse_taxonomy(file3)
                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)

                    remove_proj(dest)

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
                    return upStop(request)

                file4 = request.FILES['docfile4']
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
                    return upStop(request)

            elif source == '454_sff':
                mothurdest = 'mothur/temp'
                if not os.path.exists(mothurdest):
                    os.makedirs(mothurdest)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    return upStop(request)

                file_list = request.FILES.getlist('sff_files')
                for each in file_list:
                    file = each
                    handle_uploaded_file(file, mothurdest, each)
                    handle_uploaded_file(file, dest, each)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    return upStop(request)

                file_list = request.FILES.getlist('oligo_files')
                for each in file_list:
                    file = each
                    handle_uploaded_file(file, mothurdest, each)
                    handle_uploaded_file(file, dest, each)

                if stopList[PID] == RID:
                    remove_proj(dest)
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
                    return upStop(request)

            elif source == '454_fastq':
                mothurdest = 'mothur/temp'

                if not os.path.exists(mothurdest):
                    os.makedirs(mothurdest)

                file_list = request.FILES.getlist('fna_files')
                tempList = []
                if len(file_list) > 1:
                    for each in file_list:
                        file = each
                        handle_uploaded_file(file, mothurdest, each)
                        handle_uploaded_file(file, dest, each)
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
                    for each in file_list:
                        file = each
                        fasta = 'temp.fasta'
                        handle_uploaded_file(file, mothurdest, fasta)
                        handle_uploaded_file(file, dest, each)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    return upStop(request)

                file_list = request.FILES.getlist('qual_files')
                tempList = []
                if len(file_list) > 1:
                    for each in file_list:
                        file = each
                        handle_uploaded_file(file, mothurdest, each)
                        handle_uploaded_file(file, dest, each)
                        if os.name == 'nt':
                            myStr = "mothur\\temp\\" + str(file.name)
                        else:
                            myStr = "mothur/temp/" + str(file.name)
                        tempList.append(myStr)
                    inputList = "-".join(tempList)
                    if os.name == 'nt':
                        os.system('"mothur\\mothur-win\\mothur.exe \"#merge.files(input=%s, output=mothur\\temp\\temp.qual)\""' % inputList)
                    else:
                        os.system("mothur/mothur-linux/mothur \"#merge.files(input=%s, output=mothur/temp/temp.qual)\"" % inputList)
                else:
                    for each in file_list:
                        file = each
                        qual = 'temp.qual'
                        handle_uploaded_file(file, mothurdest, qual)
                        handle_uploaded_file(file, dest, each)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    return upStop(request)

                oligo = 'temp.oligos'
                file6 = request.FILES['docfile6']
                handle_uploaded_file(file6, mothurdest, oligo)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    return upStop(request)

                handle_uploaded_file(file6, dest, file6.name)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    return upStop(request)

                batch = 'mothur.batch'
                file7 = request.FILES['docfile7']

                avail_proc = mp.cpu_count()
                use_proc = min(avail_proc, processors)
                actual_proc = 'processors=' + str(use_proc)

                handle_uploaded_file(file7, mothurdest, batch)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    return upStop(request)

                for line in fileinput.input('mothur/temp/mothur.batch', inplace=1):
                    print line.replace("processors=X", actual_proc),

                handle_uploaded_file(file7, dest, batch)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    return upStop(request)

                try:
                    mothur(dest, source)

                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)

                    remove_proj(dest)

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
                    return upStop(request)

                try:
                    with open('% s/final.taxonomy' % dest, 'rb') as file3:
                        parse_taxonomy(file3)
                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)

                    remove_proj(dest)

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
                    return upStop(request)

                handle_uploaded_file(file13, dest, fastq)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    return upStop(request)

                file_list = request.FILES.getlist('fastq_files')
                for each in file_list:
                    file = each
                    handle_uploaded_file(file, mothurdest, each)
                    if stopList[PID] == RID:
                        remove_proj(dest)
                        return upStop(request)
                    handle_uploaded_file(file, dest, each)
                    if stopList[PID] == RID:
                        remove_proj(dest)
                        return upStop(request)

                batch = 'mothur.batch'
                file15 = request.FILES['docfile15']

                avail_proc = mp.cpu_count()
                use_proc = min(avail_proc, processors)
                actual_proc = 'processors=' + str(use_proc)

                handle_uploaded_file(file15, mothurdest, batch)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    return upStop(request)

                for line in fileinput.input('mothur/temp/mothur.batch', inplace=1):
                    print line.replace("processors=X", actual_proc),

                handle_uploaded_file(file15, dest, batch)

                if stopList[PID] == RID:
                    remove_proj(dest)
                    return upStop(request)

                try:
                    mothur(dest, source)
                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)

                    remove_proj(dest)

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
                    return upStop(request)

                try:
                    with open('% s/final.taxonomy' % dest, 'rb') as file3:
                        parse_taxonomy(file3)
                except:
                    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
                    myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
                    logging.exception(myDate)

                    remove_proj(dest)

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

        qs = Sample.objects.none()
        if request.user.is_superuser:
            qs = Sample.objects.filter(sampleid__in=selSamples)
        elif request.user.is_authenticated():
            path_list = Reference.objects.filter(Q(author=request.user)).values_list('sampleid_id')
            qs = Sample.objects.all().filter( Q(sampleid__in=path_list) | Q(status='public') ).filter(sampleid__in=selSamples)

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
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        print "projects: ", myJson
        return HttpResponse(myJson)

# need to finish the rest of the tables for the select page...
def sampleTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        qs = Sample.objects.none()
        if request.user.is_superuser:
            qs = Sample.objects.filter(sampleid__in=selSamples)
        elif request.user.is_authenticated():
            path_list = Reference.objects.filter(Q(author=request.user)).values_list('sampleid_id')
            qs = Sample.objects.all().filter( Q(sampleid__in=path_list) | Q(status='public') ).filter(sampleid__in=selSamples)

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
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


def referenceTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        qs = Sample.objects.none()
        if request.user.is_superuser:
            qs = Sample.objects.filter(sampleid__in=selSamples)
        elif request.user.is_authenticated():
            path_list = Reference.objects.filter(Q(author=request.user)).values_list('sampleid_id')
            qs = Sample.objects.all().filter( Q(sampleid__in=path_list) | Q(status='public') ).filter(sampleid__in=selSamples)

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
            "refid__author"
        )
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


def airTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        qs = Sample.objects.none()
        if request.user.is_superuser:
            qs = Sample.objects.filter(sampleid__in=selSamples)
        elif request.user.is_authenticated():
            path_list = Reference.objects.filter(Q(author=request.user)).values_list('sampleid_id')
            qs = Sample.objects.all().filter( Q(sampleid__in=path_list) | Q(status='public') ).filter(sampleid__in=selSamples)

        results = {}
        qs1 = qs.values_list(
            "air__projectid__project_name",
            "air__sampleid__sample_name",
            "air__barometric_press",
            "air__carb_dioxide",
            "air__carb_monoxide",
            "air__chem_admin_term",
            "air__chem_admin_time",
            "air__elev",
            "air__humidity",
            "air__methane",
            "air__organism_type",
            "air__organism_count",
            "air__oxy_stat_samp",
            "air__oxygen",
            "air__perturbation_type",
            "air__perturbation_interval",
            "air__pollutants_type",
            "air__pollutants_concentration",
            "air__rel_to_oxygen",
            "air__resp_part_matter_substance",
            "air__resp_part_matter_concentration",
            "air__samp_collect_device",
            "air__samp_mat_process",
            "air__samp_salinity",
            "air__samp_size",
            "air__samp_store_dur",
            "air__samp_store_loc",
            "air__samp_store_temp",
            "air__solar_irradiance",
            "air__temp",
            "air__ventilation_rate",
            "air__ventilation_type",
            "air__volatile_org_comp_name",
            "air__volatile_org_comp_concentration",
            "air__wind_direction",
            "air__wind_speed"
        )
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        print "air: ", myJson
        return HttpResponse(myJson)


def associatedTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        qs = Sample.objects.none()
        if request.user.is_superuser:
            qs = Sample.objects.filter(sampleid__in=selSamples)
        elif request.user.is_authenticated():
            path_list = Reference.objects.filter(Q(author=request.user)).values_list('sampleid_id')
            qs = Sample.objects.all().filter( Q(sampleid__in=path_list) | Q(status='public') ).filter(sampleid__in=selSamples)

        results = {}
        qs1 = qs.values_list(
            "human_associated__projectid__project_name",
            "human_associated__sampleid__sample_name",
            "human_associated__samp_collect_device",
            "human_associated__samp_mat_process",
            "human_associated__samp_size",
            "human_associated__samp_store_temp",
            "human_associated__samp_store_dur",
            "human_associated__samp_type",
            "human_associated__samp_location",
            "human_associated__samp_temp",
            "human_associated__samp_ph",
            "human_associated__samp_oxy_stat",
            "human_associated__samp_salinity",
            "human_associated__host_subject_id",
            "human_associated__host_age",
            "human_associated__host_pulse",
            "human_associated__host_gender",
            "human_associated__host_ethnicity",
            "human_associated__host_height",
            "human_associated__host_weight",
            "human_associated__host_bmi",
            "human_associated__host_weight_loss_3_month",
            "human_associated__host_body_temp",
            "human_associated__host_occupation",
            "human_associated__pet_farm_animal",
            "human_associated__smoker",
            "human_associated__diet_type",
            "human_associated__diet_duration",
            "human_associated__diet_frequency",
            "human_associated__diet_last_six_month",
            "human_associated__last_meal",
            "human_associated__medic_hist_perform",
            "human_associated__disease_type",
            "human_associated__disease_location",
            "human_associated__disease_duration",
            "human_associated__organism_count",
            "human_associated__tumor_location",
            "human_associated__tumor_mass",
            "human_associated__tumor_stage",
            "human_associated__drug_usage",
            "human_associated__drug_type",
            "human_associated__drug_duration",
            "human_associated__drug_frequency",
            "human_associated__perturbation",
            "human_associated__pert_type",
            "human_associated__pert_duration",
            "human_associated__pert_frequency",
            "human_associated__fetal_health_stat",
            "human_associated__amniotic_fluid_color",
            "human_associated__gestation_stat",
            "human_associated__maternal_health_stat"
        )
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


def microbialTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        qs = Sample.objects.none()
        if request.user.is_superuser:
            qs = Sample.objects.filter(sampleid__in=selSamples)
        elif request.user.is_authenticated():
            path_list = Reference.objects.filter(Q(author=request.user)).values_list('sampleid_id')
            qs = Sample.objects.all().filter( Q(sampleid__in=path_list) | Q(status='public') ).filter(sampleid__in=selSamples)

        results = {}
        qs1 = qs.values_list(
            "microbial__projectid__project_name",
            "microbial__sampleid__sample_name",
            "microbial__alkalinity",
            "microbial__alkyl_diethers",
            "microbial__altitude",
            "microbial__aminopept_act",
            "microbial__ammonium",
            "microbial__bacteria_carb_prod",
            "microbial__biomass_part_name",
            "microbial__biomass_amount",
            "microbial__bishomohopanol",
            "microbial__bromide",
            "microbial__calcium",
            "microbial__carb_nitro_ratio",
            "microbial__chem_administration_term",
            "microbial__chem_administration_time",
            "microbial__chloride",
            "microbial__chlorophyll",
            "microbial__diether_lipids_name",
            "microbial__diether_lipids_concentration",
            "microbial__diss_carb_dioxide",
            "microbial__diss_hydrogen",
            "microbial__diss_inorg_carb",
            "microbial__diss_org_carb",
            "microbial__diss_org_nitro",
            "microbial__diss_oxygen",
            "microbial__glucosidase_act",
            "microbial__magnesium",
            "microbial__mean_frict_vel",
            "microbial__mean_peak_frict_vel",
            "microbial__methane",
            "microbial__n_alkanes_name",
            "microbial__n_alkanes_concentration",
            "microbial__nitrate",
            "microbial__nitrite",
            "microbial__nitro",
            "microbial__org_carb",
            "microbial__org_matter",
            "microbial__org_nitro",
            "microbial__organism_name",
            "microbial__organism_count",
            "microbial__oxy_stat_samp",
            "microbial__part_org_carb",
            "microbial__perturbation_type",
            "microbial__perturbation_interval",
            "microbial__petroleum_hydrocarb",
            "microbial__ph",
            "microbial__phaeopigments_type",
            "microbial__phaeopigments_concentration",
            "microbial__phosphate",
            "microbial__phosplipid_fatt_acid_name",
            "microbial__phosplipid_fatt_acid_concentration",
            "microbial__potassium",
            "microbial__pressure",
            "microbial__redox_potential",
            "microbial__rel_to_oxygen",
            "microbial__salinity",
            "microbial__samp_collect_device",
            "microbial__samp_mat_process",
            "microbial__samp_size",
            "microbial__samp_store_dur",
            "microbial__samp_store_loc",
            "microbial__samp_store_temp",
            "microbial__silicate",
            "microbial__sodium",
            "microbial__sulfate",
            "microbial__sulfide",
            "microbial__temp",
            "microbial__tot_carb",
            "microbial__tot_nitro",
            "microbial__tot_org_carb",
            "microbial__turbidity",
            "microbial__water_content"
        )
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


def soilTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        qs = Sample.objects.none()
        if request.user.is_superuser:
            qs = Sample.objects.filter(sampleid__in=selSamples)
        elif request.user.is_authenticated():
            path_list = Reference.objects.filter(Q(author=request.user)).values_list('sampleid_id')
            qs = Sample.objects.all().filter( Q(sampleid__in=path_list) | Q(status='public') ).filter(sampleid__in=selSamples)

        results = {}
        qs1 = qs.values_list(
            "soil__projectid__project_name",
            "soil__sampleid__sample_name",
            "soil__samp_collection_device",
            "soil__samp_size",
            "soil__samp_depth",
            "soil__samp_prep",
            "soil__samp_sieve_size",
            "soil__samp_store_dur",
            "soil__samp_store_loc",
            "soil__samp_store_temp",
            "soil__samp_weight_dna_ext",
            "soil__pool_dna_extracts",
            "soil__fao_class",
            "soil__local_class",
            "soil__texture_class",
            "soil__porosity",
            "soil__profile_position",
            "soil__slope_aspect",
            "soil__slope_gradient",
            "soil__bulk_density",
            "soil__drainage_class",
            "soil__water_content_soil",
            "soil__cur_land_use",
            "soil__cur_vegetation",
            "soil__cur_crop",
            "soil__cur_cultivar",
            "soil__crop_rotation",
            "soil__cover_crop",
            "soil__fert_amendment_class",
            "soil__fert_placement",
            "soil__fert_type",
            "soil__fert_tot_amount",
            "soil__fert_N_tot_amount",
            "soil__fert_P_tot_amount",
            "soil__fert_K_tot_amount",
            "soil__irrigation_type",
            "soil__irrigation_tot_amount",
            "soil__residue_removal",
            "soil__residue_growth_stage",
            "soil__residue_removal_percent",
            "soil__tillage_event",
            "soil__tillage_event_depth",
            "soil__amend1_class",
            "soil__amend1_active_ingredient",
            "soil__amend1_tot_amount",
            "soil__amend2_class",
            "soil__amend2_active_ingredient",
            "soil__amend2_tot_amount",
            "soil__amend3_class",
            "soil__amend3_active_ingredient",
            "soil__amend3_tot_amount",
            "soil__rRNA_copies",
            "soil__microbial_biomass_C",
            "soil__microbial_biomass_N",
            "soil__microbial_respiration",
            "soil__soil_pH",
            "soil__soil_EC",
            "soil__soil_C",
            "soil__soil_OM",
            "soil__soil_N",
            "soil__soil_NO3_N",
            "soil__soil_NH4_N",
            "soil__soil_P",
            "soil__soil_K",
            "soil__soil_S",
            "soil__soil_Zn",
            "soil__soil_Fe",
            "soil__soil_Cu",
            "soil__soil_Mn",
            "soil__soil_Ca",
            "soil__soil_Mg",
            "soil__soil_Na",
            "soil__soil_B",
            "soil__plant_C",
            "soil__plant_N",
            "soil__plant_P",
            "soil__plant_K",
            "soil__plant_Ca",
            "soil__plant_Mg",
            "soil__plant_S",
            "soil__plant_Na",
            "soil__plant_Cl",
            "soil__plant_Al",
            "soil__plant_B",
            "soil__plant_Cu",
            "soil__plant_Fe",
            "soil__plant_Mn",
            "soil__plant_Zn",
            "soil__crop_tot_biomass_fw",
            "soil__crop_tot_biomass_dw",
            "soil__crop_tot_above_biomass_fw",
            "soil__crop_tot_above_biomass_dw",
            "soil__crop_tot_below_biomass_fw",
            "soil__crop_tot_below_biomass_dw",
            "soil__harv_fraction",
            "soil__harv_fresh_weight",
            "soil__harv_dry_weight",
            "soil__ghg_chamber_placement",
            "soil__ghg_N2O",
            "soil__ghg_CO2",
            "soil__ghg_NH4",
            "soil__soil_water_cap",
            "soil__soil_surf_hard",
            "soil__soil_subsurf_hard",
            "soil__soil_agg_stability",
            "soil__soil_ACE_protein",
            "soil__soil_active_C"
        )
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


def waterTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        qs = Sample.objects.none()
        if request.user.is_superuser:
            qs = Sample.objects.filter(sampleid__in=selSamples)
        elif request.user.is_authenticated():
            path_list = Reference.objects.filter(Q(author=request.user)).values_list('sampleid_id')
            qs = Sample.objects.all().filter( Q(sampleid__in=path_list) | Q(status='public') ).filter(sampleid__in=selSamples)

        results = {}
        qs1 = qs.values_list(
            "water__projectid__project_name",
            "water__sampleid__sample_name",
            "water__alkalinity",
            "water__alkyl_diethers",
            "water__altitude",
            "water__aminopept_act",
            "water__ammonium",
            "water__atmospheric_data",
            "water__bac_prod",
            "water__bac_resp",
            "water__bacteria_carb_prod",
            "water__biomass_part_name",
            "water__biomass_amount",
            "water__bishomohopanol",
            "water__bromide",
            "water__calcium",
            "water__carb_nitro_ratio",
            "water__chem_administration_name",
            "water__chem_administration_time",
            "water__chloride",
            "water__chlorophyll",
            "water__conduc",
            "water__density",
            "water__diether_lipids",
            "water__diss_carb_dioxide",
            "water__diss_hydrogen",
            "water__diss_inorg_carb",
            "water__diss_inorg_nitro",
            "water__diss_inorg_phosp",
            "water__diss_org_carb",
            "water__diss_org_nitro",
            "water__diss_oxygen",
            "water__down_par",
            "water__elev",
            "water__fluor",
            "water__glucosidase_act",
            "water__light_intensity",
            "water__magnesium",
            "water__mean_frict_vel",
            "water__mean_peak_frict_vel",
            "water__n_alkanes",
            "water__nitrate",
            "water__nitrite",
            "water__nitro",
            "water__org_carb",
            "water__org_matter",
            "water__org_nitro",
            "water__organism_name",
            "water__organism_count",
            "water__oxy_stat_samp",
            "water__part_org_carb",
            "water__part_org_nitro",
            "water__perturbation_type",
            "water__perturbation_interval",
            "water__pretroleum_hydrocarb",
            "water__ph",
            "water__phaeopigments",
            "water__phosphate",
            "water__phosplipid_fatt_acid",
            "water__photon_flux",
            "water__potassium",
            "water__pressure",
            "water__primary_prod",
            "water__redox_potential",
            "water__rel_to_oxygen",
            "water__samp_mat_process",
            "water__samp_salinity",
            "water__samp_size",
            "water__samp_store_dur",
            "water__samp_store_loc",
            "water__samp_store_temp",
            "water__samp_vol_we_dna_ext",
            "water__silicate",
            "water__sodium",
            "water__soluble_react_phosp",
            "water__source_material_id",
            "water__sulfate",
            "water__sulfide",
            "water__suspen_part_matter",
            "water__temp",
            "water__tidal_stage",
            "water__tot_depth_water_col",
            "water__tot_diss_nitro",
            "water__tot_inorg_nitro",
            "water__tot_nitro",
            "water__tot_part_carb",
            "water__tot_phosp",
            "water__water_current_direction",
            "water__water_current_magnitude"
        )
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


def userTableJSON(request):
    if request.is_ajax():
        jsonSamples = request.GET['key']
        selSamples = simplejson.loads(jsonSamples)

        qs = Sample.objects.none()
        if request.user.is_superuser:
            qs = Sample.objects.filter(sampleid__in=selSamples)
        elif request.user.is_authenticated():
            path_list = Reference.objects.filter(Q(author=request.user)).values_list('sampleid_id')
            qs = Sample.objects.all().filter( Q(sampleid__in=path_list) | Q(status='public') ).filter(sampleid__in=selSamples)

        results = {}
        qs1 = qs.values_list(
            "userdefined__projectid__project_name",
            "userdefined__sampleid__sample_name",
            "userdefined__usr_cat1",
            "userdefined__usr_cat2",
            "userdefined__usr_cat3",
            "userdefined__usr_cat4",
            "userdefined__usr_cat5",
            "userdefined__usr_cat6",
            "userdefined__usr_quant1",
            "userdefined__usr_quant2",
            "userdefined__usr_quant3",
            "userdefined__usr_quant4",
            "userdefined__usr_quant5",
            "userdefined__usr_quant6"
        )
        results['data'] = list(qs1)
        myJson = simplejson.dumps(results, ensure_ascii=False)
        return HttpResponse(myJson)


@login_required(login_url='/myPhyloDB/accounts/login/')
def select(request):
    return render_to_response(
        'select.html',
        {'form9': UploadForm9},
        context_instance=RequestContext(request))


def taxaJSON(request):
    results = {}
    qs1 = Species.objects.values_list('kingdomid', 'kingdomid__kingdomName', 'phylaid', 'phylaid__phylaName', 'classid', 'classid__className', 'orderid', 'orderid__orderName', 'familyid', 'familyid__familyName', 'genusid', 'genusid__genusName', 'speciesid', 'speciesName')
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
    qs1 = ko_entry.objects.using('picrust').values_list('ko_lvl1_id__ko_lvl1_name', 'ko_lvl1_id', 'ko_lvl2_id__ko_lvl2_name', 'ko_lvl2_id', 'ko_lvl3_id__ko_lvl3_name', 'ko_lvl3_id', 'ko_lvl4_id', 'ko_orthology', 'ko_name', 'ko_desc')
    results['data'] = list(qs1)
    myJson = simplejson.dumps(results, ensure_ascii=False)
    return HttpResponse(myJson)


def kegg_path(request):
    return render_to_response(
        'kegg_path.html',
        context_instance=RequestContext(request)
    )


def nzJSON(request):
    qs1 = nz_entry.objects.using('picrust').values_list('nz_lvl1_id__nz_lvl1_name', 'nz_lvl1_id', 'nz_lvl2_id__nz_lvl2_name', 'nz_lvl2_id', 'nz_lvl3_id__nz_lvl3_name', 'nz_lvl3_id', 'nz_lvl4_id__nz_lvl4_name', 'nz_lvl4_id', 'nz_lvl5_id', 'nz_orthology', 'nz_name', 'nz_desc')
    results = {}
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

        # save file to users temp/ folder
        myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
        path = str(myDir) + 'usr_sel_samples.pkl'

        if not os.path.exists(myDir):
            os.makedirs(myDir)

        with open(path, 'wb') as f:
            pickle.dump(selList, f)

        # remove normalized file
        try:
            os.remove('myPhyloDB/media/usr_temp/' + str(request.user) + '/usr_norm_data.csv')
        except OSError:
            pass

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
        placeholder = ''

    try:
        templateFile = request.FILES['docfile9']
        templateDB = request.FILES['docfile9'].name
        handle_uploaded_file(templateFile, 'mothur/reference/template', templateDB)
    except:
        placeholder = ''

    try:
        taxonomyFile = request.FILES['docfile10']
        taxonomyDB = request.FILES['docfile10'].name
        handle_uploaded_file(taxonomyFile, 'mothur/reference/taxonomy', taxonomyDB)
    except:
        placeholder = ''

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

        if stopList[PID] == RID:
            return updaStop(request)

        try:
            metaFile = '/'.join([dest, file1.name])
            handle_uploaded_file(file1, dest, file1.name)
            parse_project(metaFile, p_uuid)
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
            f = xlrd.open_workbook(metaFile)
            sheet = f.sheet_by_name('Project')
            num_samp = int(sheet.cell_value(rowx=5, colx=0))

            bat = 'mothur.batch'
            batPath = '/'.join([dest, bat])
            reference = Reference.objects.get(refid=refid)
            raw = reference.raw
            source = reference.source
            userID = str(request.User__id)

            # add stops to parse functions?
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

    metaDF = tempDF.drop(['kingdomName', 'phylaName', 'className', 'orderName', 'familyName', 'genusName', 'speciesName', 'speciesid', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity'], axis=1)
    myList = list(tempDF.index.values)

    nameList = []
    for i in myList:
        nameList.append({"id": str(i), "metadata": metaDF.loc[i].to_dict()})

    # get list of lists with abundances
    taxaOnlyDF = savedDF.loc[:, ['sampleid', 'kingdomName', 'phylaName', 'className', 'orderName', 'familyName', 'genusName', 'speciesName', 'speciesid', 'abund']]
    taxaOnlyDF = taxaOnlyDF.pivot(index='speciesid', columns='sampleid', values='abund')
    dataList = taxaOnlyDF.values.tolist()

    # get list of taxa
    namesDF = savedDF.loc[:, ['sampleid', 'speciesid']]
    namesDF['taxa'] = savedDF.loc[:, ['kingdomName', 'phylaName', 'className', 'orderName', 'familyName', 'genusName', 'speciesName']].values.tolist()
    namesDF = namesDF.pivot(index='speciesid', columns='sampleid', values='taxa')

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

user_logged_in.connect(login_usr_callback)


# Function has been added to context processor in setting file
def usrFiles(request):

    selFiles = os.path.exists('myPhyloDB/media/usr_temp/' + str(request.user) + '/usr_sel_samples.pkl')
    normFiles = os.path.exists('myPhyloDB/media/usr_temp/' + str(request.user) + '/usr_norm_data.csv')

    return {
        'selFiles': selFiles,
        'normFiles': normFiles
    }
