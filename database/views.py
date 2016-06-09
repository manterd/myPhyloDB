import datetime
from django.contrib.auth.decorators import login_required
from django.contrib.auth.signals import user_logged_in
from django.db.models import Q
from django.http import *
from django.shortcuts import render_to_response, HttpResponseRedirect
from django.template import RequestContext
import fileinput
import logging
import multiprocessing as mp
import os
import pandas as pd
import pickle
import simplejson
import xlrd
import time
from uuid import uuid4

from forms import UploadForm1, UploadForm2, UploadForm4, UploadForm5, \
    UploadForm6, UploadForm7, UploadForm8, UploadForm9
from models import Project, Reference, Sample, Species
from models import ko_entry, nz_entry, nz_lvl1
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


@login_required(login_url='/accounts/login/')
def upload(request):
    projects = Reference.objects.none()
    if request.method == 'POST' and 'Upload' in request.POST:
        start = datetime.datetime.now()
        form1 = UploadForm1(request.POST, request.FILES)
        source = str(request.POST['source'])
        userID = str(request.user.id)
        processors = int(request.POST['processors'])

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

            elif source == '454_sff':
                mothurdest = 'mothur/temp'

                if not os.path.exists(mothurdest):
                    os.makedirs(mothurdest)

                file_list = request.FILES.getlist('sff_files')
                for each in file_list:
                    file = each
                    handle_uploaded_file(file, mothurdest, each)
                    handle_uploaded_file(file, dest, each)

                file_list = request.FILES.getlist('oligo_files')
                for each in file_list:
                    file = each
                    handle_uploaded_file(file, mothurdest, each)
                    handle_uploaded_file(file, dest, each)

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

                oligo = 'temp.oligos'
                file6 = request.FILES['docfile6']
                handle_uploaded_file(file6, mothurdest, oligo)
                handle_uploaded_file(file6, dest, file6.name)

                batch = 'mothur.batch'
                file7 = request.FILES['docfile7']

                avail_proc = mp.cpu_count()
                use_proc = min(avail_proc, processors)
                actual_proc = 'processors=' + str(use_proc)

                handle_uploaded_file(file7, mothurdest, batch)

                for line in fileinput.input('mothur/temp/mothur.batch', inplace=1):
                    print line.replace("processors=X", actual_proc),

                handle_uploaded_file(file7, dest, batch)

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

            elif source == 'miseq':
                mothurdest = 'mothur/temp'
                if not os.path.exists(mothurdest):
                    os.makedirs(mothurdest)

                fastq = 'temp.files'
                file13 = request.FILES['docfile13']
                handle_uploaded_file(file13, mothurdest, fastq)
                handle_uploaded_file(file13, dest, fastq)

                file_list = request.FILES.getlist('fastq_files')
                for each in file_list:
                    file = each
                    handle_uploaded_file(file, mothurdest, each)
                    handle_uploaded_file(file, dest, each)

                batch = 'mothur.batch'
                file15 = request.FILES['docfile15']

                avail_proc = mp.cpu_count()
                use_proc = min(avail_proc, processors)
                actual_proc = 'processors=' + str(use_proc)

                handle_uploaded_file(file15, mothurdest, batch)

                for line in fileinput.input('mothur/temp/mothur.batch', inplace=1):
                    print line.replace("processors=X", actual_proc),

                handle_uploaded_file(file15, dest, batch)

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


@login_required(login_url='/accounts/login/')
def select(request):
    projects = Project.objects.none()
    if request.user.is_superuser:
        projects = Project.objects.all()
    elif request.user.is_authenticated():
        path_list = Reference.objects.filter(Q(author=request.user)).values_list('projectid_id')
        projects = Project.objects.all().filter( Q(projectid__in=path_list) | Q(status='public') )
    if not request.user.is_superuser and not request.user.is_authenticated():
        projects = Project.objects.all().filter( Q(status='public') )

    refs = Reference.objects.filter(projectid__in=projects)
    samples = Sample.objects.filter(projectid__in=projects)

    return render_to_response(
        'select.html',
        {'form9': UploadForm9,
         'projects': projects,
         'refs': refs,
         'samples': samples},
        context_instance=RequestContext(request)
    )


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


@login_required(login_url='/accounts/login/')
def norm(request):
    cleanup('media/temp/norm')

    return render_to_response(
        'norm.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/accounts/login/')
def ANOVA(request):
    cleanup('media/temp/anova')

    return render_to_response(
        'anova.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/accounts/login/')
def rich(request):
    cleanup('media/temp/spac')

    return render_to_response(
        'SpAC.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/accounts/login/')
def soil_index(request):
    cleanup('media/temp/soil_index')

    return render_to_response(
        'soil_index.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/accounts/login/')
def DiffAbund(request):
    cleanup('media/temp/diffabund')

    return render_to_response(
        'diff_abund.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/accounts/login/')
def GAGE(request):
    cleanup('media/temp/gage')

    return render_to_response(
        'gage.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/accounts/login/')
def PCA(request):
    cleanup('media/temp/pca')

    return render_to_response(
        'pca.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/accounts/login/')
def PCoA(request):
    cleanup('media/temp/pcoa')

    return render_to_response(
        'pcoa.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/accounts/login/')
def SPLS(request):
    cleanup('media/temp/spls')

    return render_to_response(
        'spls.html',
        context_instance=RequestContext(request)
    )


@login_required(login_url='/accounts/login/')
def WGCNA(request):
    cleanup('media/temp/wgcna')

    return render_to_response(
        'wgcna.html',
        context_instance=RequestContext(request)
    )


def saveSampleList(request):
    if request.is_ajax():
        allJson = request.GET["all"]
        selList = simplejson.loads(allJson)

        # save file to users temp/ folder
        myDir = 'media/usr_temp/' + str(request.user) + '/'
        path = str(myDir) + 'usr_sel_samples.pkl'

        if not os.path.exists(myDir):
            os.makedirs(myDir)

        with open(path, 'wb') as f:
            pickle.dump(selList, f)

        # remove normalized file
        try:
            os.remove('media/usr_temp/' + str(request.user) + '/usr_norm_data.csv')
        except OSError:
            pass

        text = 'Selected sample(s) have been recorded!'
        res = simplejson.dumps(text, encoding="Latin-1")
        return HttpResponse(res, content_type='application/json')


@login_required(login_url='/accounts/login/')
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


@login_required(login_url='/accounts/login/')
def update(request):
    form5 = UploadForm5(request.POST, request.FILES)
    state = ''

    if form5.is_valid():
        refid = request.POST['refid']
        file1 = request.FILES['docfile11']

        ref = Reference.objects.get(refid=refid)
        dest = ref.path
        p_uuid, pType, num_samp = projectid(file1)

        try:
            metaFile = '/'.join([dest, file1.name])
            handle_uploaded_file(file1, dest, file1.name)
            parse_project(metaFile, p_uuid)
        except:
            state = "There was an error parsing your metafile: " + str(file1.name)
            return render_to_response(
                'update.html',
                {'form5': UploadForm5,
                 'state': state},
                context_instance=RequestContext(request)
            )

        try:
            f = xlrd.open_workbook(metaFile)
            sheet = f.sheet_by_name('Project')
            num_samp = int(sheet.cell_value(rowx=5, colx=0))

            bat = 'mothur.batch'
            batPath = '/'.join([dest, bat])
            reference = Reference.objects.get(refid=refid)
            raw = reference.raw
            source = reference.source
            userID = str(request.user.id)

            with open(batPath, 'rb') as batFile:
                refDict = parse_sample(metaFile, p_uuid, pType, num_samp, dest, batFile, raw, source, userID)

        except Exception as e:
            state = "There was an error parsing your metafile: " + str(file1.name)
            return render_to_response(
                'update.html',
                {'form5': UploadForm5,
                 'state': state},
                context_instance=RequestContext(request)
            )

        state = 'Path: ' + str(dest) + ' is finished parsing!'

    return render_to_response(
        'update.html',
        {'form5': UploadForm5,
         'state': state},
        context_instance=RequestContext(request)
    )


@login_required(login_url='/accounts/login/')
def pybake(request):
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
    myDir = 'media/usr_temp/' + str(request.user) + '/'
    path = str(myDir) + 'usr_sel_samples.pkl'

    if not os.path.exists(myDir):
        os.makedirs(myDir)

    with open(path, 'wb') as f:
        pickle.dump(selList, f)

    # save normalized file to users temp/ folder
    myDir = 'media/usr_temp/' + str(request.user) + '/'
    path = str(myDir) + 'usr_norm_data.csv'

    if not os.path.exists(myDir):
        os.makedirs(myDir)
    savedDF.to_csv(path, sep='\t')

    projects = list(set(savedDF['projectid']))

    if Reference.objects.filter(projectid__in=projects).exists():
        refs = Reference.objects.filter(projectid__in=projects).values_list('refid', flat=True)
    else:
        refs = uuid4().hex

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

    myDir = 'media/usr_temp/' + str(request.user) + '/'
    path = str(myDir) + 'usr_norm_data.biom'
    with open(path, 'w') as outfile:
        simplejson.dump(biome, outfile, ensure_ascii=True, indent=4, sort_keys=True)

    return render_to_response(
        'select.html',
        {'form9': UploadForm9,
         'projects': projects,
         'refs': refs,
         'samples': selList,
         'normpost': 'Success'},
        context_instance=RequestContext(request)
    )


# Function to create a user folder at login
def login_usr_callback(sender, user, request, **kwargs):
    user = request.user
    if not os.path.exists('media/usr_temp/'+str(user)):
        os.makedirs('media/usr_temp/'+str(user))
    cleanup('media/usr_temp/'+str(user))

user_logged_in.connect(login_usr_callback)


# Add info on user files to base
# Function needs to be added to settings file
def usrFiles(request):

    selFiles = os.path.exists('media/usr_temp/' + str(request.user) + '/usr_sel_samples.pkl')
    normFiles = os.path.exists('media/usr_temp/' + str(request.user) + '/usr_norm_data.csv')

    return {
        'selFiles': selFiles,
        'normFiles': normFiles
    }