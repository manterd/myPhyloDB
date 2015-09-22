import datetime
import os
import pandas as pd
import pickle
import simplejson
from django.http import *
from django.shortcuts import render_to_response
from django.template import RequestContext
from forms import UploadForm1, UploadForm2, UploadForm4, UploadForm5
from models import Project, Reference, Sample, Species
from parsers import mothur, projectid, parse_project, parse_reference, parse_sample, parse_taxonomy, parse_profile
from utils import handle_uploaded_file, remove_list, remove_proj
from django.contrib.auth.decorators import login_required
from django.contrib.auth import authenticate, login, logout
from uuid import uuid4
from django.db.models import Q
from django.contrib.auth.models import User as Users

rep_project = ''


def home(request):
    return render_to_response('home.html')


@login_required(login_url='/myPhyloDB/login/')
def upload(request):
    projects = Reference.objects.none()

    if request.method == 'POST' and 'Upload' in request.POST:
        form1 = UploadForm1(request.POST, request.FILES)
        source = str(request.POST['source'])
        userID = str(request.user.id)

        if form1.is_valid():
            try:
                file1 = request.FILES['docfile1']
                p_uuid, pType = projectid(file1)
            except Exception as e:
                print("Error with project file: " + str(e))

                if request.user.is_superuser:
                    projects = Reference.objects.all().order_by('projectid__project_name', 'path')
                elif request.user.is_authenticated():
                    projects = Reference.objects.all().order_by('projectid__project_name', 'path').filter(author=request.user)
                return render_to_response(
                    'upload.html',
                    {'projects': projects,
                     'form1': UploadForm1,
                     'form2': UploadForm2,
                     'error': "There was an error parsing your Project file"},
                    context_instance=RequestContext(request)
                )

            date = datetime.date.today().isoformat()
            hour = datetime.datetime.now().hour
            minute = datetime.datetime.now().minute
            second = datetime.datetime.now().second
            timestamp = ".".join([str(hour), str(minute), str(second)])
            datetimestamp = "_".join([str(date), str(timestamp)])
            dest = "/".join(["uploads", str(p_uuid), str(datetimestamp)])
            refid = uuid4().hex

            try:
                parse_project(file1, p_uuid)

            except Exception as e:
                print("Error with project file: " + str(e))
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
                     'error': "There was an error parsing your Project file"},
                    context_instance=RequestContext(request)
                )

            try:
                if source == 'mothur':
                    file7 = 'blank'
                    raw = False
                    parse_reference(p_uuid, refid, dest, file7, raw, source, userID)
                elif source == '454':
                    file7 = request.FILES['docfile7']
                    raw = True
                    parse_reference(p_uuid, refid, dest, file7, raw, source, userID)
                elif source == 'miseq':
                    file15 = request.FILES['docfile15']
                    raw = True
                    parse_reference(p_uuid, refid, dest, file15, raw, source, userID)

            except Exception as e:
                print("Error with project file: " + str(e))
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
                     'error': "There was an error parsing your Project file"},
                    context_instance=RequestContext(request)
                )

            try:
                parse_sample(file1, p_uuid, refid, dest, pType)
            except Exception as e:
                print("Error with sample file: " + str(e))
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
                     'error': "There was an error parsing your Sample file"},
                    context_instance=RequestContext(request)
                )

            if source == 'mothur':
                taxonomy = ".".join(["mothur", "taxonomy"])
                file3 = request.FILES['docfile3']
                handle_uploaded_file(file3, dest, taxonomy)

                try:
                    with open('% s/mothur.taxonomy' % dest, 'rb') as file3:
                        parse_taxonomy(file3)
                except Exception as e:
                    print("Error with taxonomy file: " + str(e))
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
                         'error': "There was an error parsing your Taxonomy file"},
                        context_instance=RequestContext(request)
                    )

                shared = ".".join(["mothur", "shared"])
                file4 = request.FILES['docfile4']
                handle_uploaded_file(file4, dest, shared)

                try:
                    with open('% s/mothur.taxonomy' % dest, 'rb') as file3:
                        with open('% s/mothur.shared' % dest, 'rb') as file4:
                            parse_profile(file3, file4, p_uuid, refid)
                except Exception as e:
                    print("Error with shared file: " + str(e))
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
                         'error': "There was an error parsing your Shared file"},
                        context_instance=RequestContext(request)
                    )

            elif source == '454':
                mothurdest = 'mothur/temp'

                if not os.path.exists(mothurdest):
                    os.makedirs(mothurdest)

                sff = 'temp.sff'
                file5 = request.FILES['docfile5']
                handle_uploaded_file(file5, mothurdest, sff)

                oligo = 'temp.oligos'
                file6 = request.FILES['docfile6']
                handle_uploaded_file(file6, mothurdest, oligo)

                batch = 'mothur.batch'
                file7 = request.FILES['docfile7']
                handle_uploaded_file(file7, mothurdest, batch)

                try:
                    mothur(dest, source)

                except Exception as e:
                    print("Encountered error with Mothur: " + str(e))
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
                         'error': "There was an error with Mothur"},
                        context_instance=RequestContext(request)
                    )

                try:
                    with open('% s/mothur.taxonomy' % dest, 'rb') as file3:
                        parse_taxonomy(file3)
                except Exception as e:
                    print("Error with post-mothur taxonomy file: " + str(e))
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
                         'error': "There was an error parsing your Taxonomy file (post-Mothur)"},
                        context_instance=RequestContext(request)
                    )

                try:
                    with open('% s/mothur.taxonomy' % dest, 'rb') as file3:
                        with open('% s/mothur.shared' % dest, 'rb') as file4:
                            parse_profile(file3, file4, p_uuid, refid)
                except Exception as e:
                    print("Error with parsing post-mothur profile: " + str(e))
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
                         'error': "There was an error parsing your Profile (post-Mothur)"},
                        context_instance=RequestContext(request)
                    )

            elif source == 'miseq':
                mothurdest = 'mothur/temp'
                if not os.path.exists(mothurdest):
                    os.makedirs(mothurdest)

                fastq = 'temp.files'
                file13 = request.FILES['docfile13']
                handle_uploaded_file(file13, mothurdest, fastq)

                file_list = request.FILES.getlist('files')
                for each in file_list:
                    file = each
                    handle_uploaded_file(file, mothurdest, each)

                batch = 'mothur.batch'
                file15 = request.FILES['docfile15']
                handle_uploaded_file(file15, mothurdest, batch)

                try:
                    mothur(dest, source)
                except Exception as e:
                    print("Encountered error with Mothur: " + str(e))
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
                         'error': "There was an error with Mothur"},
                        context_instance=RequestContext(request)
                    )

                try:
                    with open('% s/mothur.taxonomy' % dest, 'rb') as file3:
                        parse_taxonomy(file3)
                except Exception as e:
                    print("Error with post-mothur taxonomy file: " + str(e))
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
                         'error': "There was an error parsing your Taxonomy file (post-Mothur)"},
                        context_instance=RequestContext(request)
                    )

                try:
                    with open('% s/mothur.taxonomy' % dest, 'rb') as file3:
                        with open('% s/mothur.shared' % dest, 'rb') as file4:
                            parse_profile(file3, file4, p_uuid, refid)
                except Exception as e:
                    print("Error with parsing post-mothur profile: " + str(e))
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
                         'error': "There was an error parsing your Profile (post-Mothur)"},
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


def select(request):
    projects = Project.objects.none()
    if request.user.is_superuser:
        projects = Project.objects.all()
    if request.user.is_authenticated():
        path_list = Reference.objects.filter(Q(author=request.user)).values_list('projectid_id')
        projects = Project.objects.all().filter( Q(projectid__in=path_list) | Q(status='public') )
    if not request.user.is_superuser and not request.user.is_authenticated():
        projects = Project.objects.all().filter( Q(status='public') )

    refs = Reference.objects.filter(projectid__in=projects)
    samples = Sample.objects.filter(projectid__in=projects)

    return render_to_response(
        'select.html',
        {'projects': projects,
         'refs': refs,
         'samples': samples},
        context_instance=RequestContext(request)
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

    return render_to_response(
        'taxa.html',
        {'table': table},
        context_instance=RequestContext(request)
    )


def ANOVA(request):
    return render_to_response(
        'anova.html',
        context_instance=RequestContext(request)
    )


def DiffAbund(request):
    return render_to_response(
        'diff_abund.html',
        context_instance=RequestContext(request)
    )


def PCoA(request):
    return render_to_response(
        'pcoa.html',
        context_instance=RequestContext(request)
    )


def saveCookie(request):
    if request.is_ajax():
        allJson = request.GET["all"]
        selList = simplejson.loads(allJson)
        qs = Sample.objects.all().filter(sampleid__in=selList).values_list('sampleid')
        request.session['selected_samples'] = pickle.dumps(qs.query)

        text = 'Selected sample(s) have been recorded!'
        res = simplejson.dumps(text, encoding="Latin-1")
        return HttpResponse(res, content_type='application/json')


def getCookie(request):
    samples = Sample.objects.all()
    try:
        samples.query = pickle.loads(request.session['selected_samples'])
        return HttpResponse('yes', content_type='application/text')
    except:
        return HttpResponse('no', content_type='application/text')


@login_required(login_url='/myPhyloDB/login/')
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


@login_required(login_url='/myPhyloDB/login/')
def update(request):
    form5 = UploadForm5(request.POST, request.FILES)
    state = ''

    if form5.is_valid():
        refid = request.POST['refid']
        file1 = request.FILES['docfile11']
        file2 = request.FILES['docfile12']

        project = Reference.objects.get(refid=refid)
        p_uuid = project.projectid.projectid
        pType = project.projectid.projectType
        dest = project.path

        try:
            parse_project(file1, dest, p_uuid, pType)
        except Exception as e:
            state = "Error with project file: " + str(e)
            return render_to_response(
                'update.html',
                {'form5': UploadForm5,
                 'state': state},
                context_instance=RequestContext(request)
            )

        try:
            parse_sample(file2, p_uuid, refid, dest, pType)
        except Exception as e:
            state = "Error with sample file: " + str(e)
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
