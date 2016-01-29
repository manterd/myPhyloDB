import datetime
from django.contrib.auth.decorators import login_required
from django.db.models import Q
from django.http import *
from django.shortcuts import render
import fileinput
import logging
import multiprocessing as mp
import os
import pandas as pd
import pickle
import simplejson
import xlrd
import time

from forms import UploadForm1, UploadForm2, UploadForm4, UploadForm5
from models import Project, Reference, Sample, Species
from models import addQueue, getQueue, subQueue
from parsers import mothur, projectid, parse_project, parse_sample, parse_taxonomy, parse_profile
from utils import handle_uploaded_file, remove_list, remove_proj


LOG_FILENAME = 'error_log.txt'
rep_project = ''


def home(request):
    return render(
        request,
        'home.html',
        {}
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
                return render(
                    request,
                    'upload.html',
                    {'projects': projects,
                     'form1': UploadForm1,
                     'form2': UploadForm2,
                     'error': "TmyPhyloDB could not parse your Excel meta file\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
                     }
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
                return render(
                    request,
                    'upload.html',
                    {'projects': projects,
                     'form1': UploadForm1,
                     'form2': UploadForm2,
                     'error': "TmyPhyloDB could not parse your Excel meta file\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
                     }
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
                return render(
                    request,
                    'upload.html',
                    {'projects': projects,
                     'form1': UploadForm1,
                     'form2': UploadForm2,
                     'error': "TmyPhyloDB could not parse your Excel meta file\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
                     }
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
                    return render(
                        request,
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "TmyPhyloDB could not parse your taxonomy file\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
                         }
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
                    return render(
                        request,
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "TmyPhyloDB could not parse your taxonomy file\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
                         }
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
                    return render(
                        request,
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "TmyPhyloDB found an error with your Mothur batch file\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
                         }
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
                    return render(
                        request,
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "TmyPhyloDB could not parse your taxonomy file\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
                         }
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
                    return render(
                        request,
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "TmyPhyloDB could not parse your shared file\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
                         }
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
                    return render(
                        request,
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "TmyPhyloDB found an error with your Mothur batch file\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
                         }
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
                    return render(
                        request,
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "TmyPhyloDB could not parse your taxonomy file\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
                         }
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
                    return render(
                        request,
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "TmyPhyloDB could not parse your shared file\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
                         }
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
                    return render(
                        request,
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "TmyPhyloDB found an error with your Mothur batch file\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
                         }
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
                    return render(
                        request,
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "TmyPhyloDB could not parse your taxonomy file\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
                         }
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
                    return render(
                        request,
                        'upload.html',
                        {'projects': projects,
                         'form1': UploadForm1,
                         'form2': UploadForm2,
                         'error': "TmyPhyloDB could not parse your shared file\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
                         }
                    )

            else:
                print ('Please check that all necessary files have been selected.')

    elif request.method == 'POST' and 'clickMe' in request.POST:
        remove_list(request)

    if request.user.is_superuser:
        projects = Reference.objects.all().order_by('projectid__project_name', 'path')
    elif request.user.is_authenticated():
        projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)

    return render(
        request,
        'upload.html',
        {'projects': projects,
         'form1': UploadForm1,
         'form2': UploadForm2,
         'error': ""
         }
    )


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

    return render(
        request,
        'select.html',
        {'projects': projects,
         'refs': refs,
         'samples': samples
         }
    )


def taxa(request):
    qs1 = Species.objects.values('kingdomid__kingdomName', 'kingdomid', 'phylaid__phylaName', 'phylaid', 'classid__className', 'classid', 'orderid__orderName', 'orderid', 'familyid__familyName', 'familyid', 'genusid__genusName', 'genusid', 'speciesName', 'speciesid')
    df = pd.DataFrame.from_records(qs1)
    df.rename(columns={'kingdomid__kingdomName': 'Kingdom Name', 'kingdomid': 'Kingdom ID'}, inplace=True)
    df.rename(columns={'phylaid__phylaName': 'Phylum Name', 'phylaid': 'Phylum ID'}, inplace=True)
    df.rename(columns={'classid__className': 'Class Name', 'classid': 'Class ID'}, inplace=True)
    df.rename(columns={'orderid__orderName': 'Order Name', 'orderid': 'Order ID'}, inplace=True)
    df.rename(columns={'familyid__familyName': 'Family Name', 'familyid': 'Family ID'}, inplace=True)
    df.rename(columns={'genusid__genusName': 'Genus Name', 'genusid': 'Genus ID'}, inplace=True)
    df.rename(columns={'speciesName': 'Species Name', 'speciesid': 'Species ID'}, inplace=True)
    table = df.to_html(classes="table display", columns=['Kingdom Name', 'Kingdom ID', 'Phylum Name', 'Phylum ID', 'Class Name', 'Class ID', 'Order Name', 'Order ID', 'Family Name', 'Family ID', 'Genus Name', 'Genus ID', 'Species Name', 'Species ID'])
    table = table.replace('border="1"', 'border="0"')

    return render(
        request,
        'taxa.html',
        {'table': table}
    )


def norm(request):
    platform = os.name

    return render(
        request,
        'norm.html',
        {'platform': platform}
    )


def ANOVA(request):
    return render(
        request,
        'anova.html',
        {}
    )


def DiffAbund(request):
    return render(
        request,
        'diff_abund.html',
        {}
    )


def PCoA(request):
    name = request.user
    ip = request.META.get('REMOTE_ADDR')
    fileStr = str(name) + "." + str(ip)

    return render(
        request,
        'pcoa.html',
        {'fileStr': fileStr}
    )


def SPLS(request):
    name = request.user
    ip = request.META.get('REMOTE_ADDR')
    fileStr = str(name) + "." + str(ip)

    return render(
        request,
        'spls.html',
        {'fileStr': fileStr}
    )


def saveSampleCookie(request):
    if request.is_ajax():
        allJson = request.GET["all"]
        selList = simplejson.loads(allJson)
        qs = Sample.objects.all().filter(sampleid__in=selList).values_list('sampleid')
        request.session['selected_samples'] = pickle.dumps(qs.query)

        text = 'Selected sample(s) have been recorded!'
        res = simplejson.dumps(text, encoding="Latin-1")
        return HttpResponse(res, content_type='application/json')


def getSampleCookie(request):
    myDict = {}
    if "selected_samples" in request.session:
        myDict['select'] = 'yes'
    else:
        myDict['select'] = 'no'
    if "savedDF" in request.session:
        myDict['norm'] = 'yes'
    else:
        myDict['norm'] = 'no'
    if "NormMeth" in request.session:
        myDict['NormMeth'] = request.session['NormMeth']
    else:
        myDict['NormMeth'] = 'no'
    res = simplejson.dumps(myDict, encoding="Latin-1")
    return HttpResponse(res, content_type='application/json')


def clearNormCookie(request):
    request.session.pop('NormMeth', None)
    request.session.pop('savedDF', None)
    return HttpResponse()


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

    return render(
        request,
        'reprocess.html',
        {'form4': UploadForm4,
         'alignDB': alignDB,
         'templateDB': templateDB,
         'taxonomyDB': taxonomyDB
         }
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
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)

            state = "TmyPhyloDB could not parse your Excel meta file\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            return render(
                request,
                'update.html',
                {'form5': UploadForm5,
                 'state': state}
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
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)

            state = "TmyPhyloDB could not parse your Excel meta file\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            return render(
                request,
                'update.html',
                {'form5': UploadForm5,
                 'state': state
                 }
            )

        state = 'Path: ' + str(dest) + ' is finished parsing!'

    return render(
        request,
        'update.html',
        {'form5': UploadForm5,
         'state': state
         }
    )

