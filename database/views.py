import datetime
import os
import pandas as pd
import pickle
import simplejson
from django.core.servers.basehttp import FileWrapper
from django.http import HttpResponse
from django.shortcuts import render_to_response
from django.template import RequestContext
from forms import UploadForm1, UploadForm2, UploadForm3
from models import Project, Sample, Species
from parsers import mothur, projectid, parse_project, parse_sample, parse_taxonomy, parse_profile
from utils import handle_uploaded_file, remove_list
from django.contrib.auth.decorators import login_required
from django.contrib.auth import logout


def home(request):
    return render_to_response('home.html')


def logout_view(request):
    logout(request)


@login_required(login_url='/myPhyloDB/login/')
def users(request):
    return render_to_response('users.html')


@login_required(login_url='/myPhyloDB/login/')
def upload(request):
    if request.method == 'POST' and 'Upload' in request.POST:
        form1 = UploadForm1(request.POST, request.FILES)
        form2 = UploadForm2(request.POST, request.FILES)
        form3 = UploadForm3(request.POST, request.FILES)

        if form1.is_valid():
            project = ".".join(["project", "csv"])
            file1 = request.FILES['docfile1']
            p_uuid = projectid(file1)

            date = datetime.date.today().isoformat()
            hour = datetime.datetime.now().hour
            minute = datetime.datetime.now().minute
            second = datetime.datetime.now().second
            timestamp = ".".join([str(hour), str(minute), str(second)])
            datetimestamp = "_".join([str(date), str(timestamp)])

            dest = "/".join(["uploads", str(p_uuid), str(datetimestamp)])
            handle_uploaded_file(file1, dest, project)
            parse_project(file1, dest, p_uuid)

            sample = ".".join(["sample", "csv"])
            file2 = request.FILES['docfile2']
            handle_uploaded_file(file2, dest, sample)
            parse_sample(file2, p_uuid)

            if form2.is_valid():
                taxonomy = ".".join(["mothur", "taxonomy"])
                file3 = request.FILES['docfile3']
                handle_uploaded_file(file3, dest, taxonomy)
                parse_taxonomy(file3, "1")

                shared = ".".join(["mothur", "shared"])
                file4 = request.FILES['docfile4']
                handle_uploaded_file(file4, dest, shared)
                parse_profile(file3, file4, p_uuid)

            if form3.is_valid():
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

                mothur(dest)

                with open('% s/mothur.taxonomy' % dest, 'rb') as file3:
                    parse_taxonomy(file3, "no")

                with open('% s/mothur.taxonomy' % dest, 'rb') as file3:
                    with open('% s/mothur.shared' % dest, 'rb') as file4:
                        parse_profile(file3, file4, p_uuid)

    elif request.method == 'POST' and 'clickMe' in request.POST:
        remove_list(request)

    projects = Project.objects.all().order_by('project_name')
    return render_to_response(
        'upload.html',
        {'projects': projects,
         'form1': UploadForm1,
         'form2': UploadForm2,
         'form3': UploadForm3},
        context_instance=RequestContext(request)
    )


def select(request):
    projects = Project.objects.all()
    samples = Sample.objects.all()

    return render_to_response(
        'select.html',
        {'projects': projects,
         'samples': samples},
        context_instance=RequestContext(request)
    )


def taxa(request):
    qs1 = Species.objects.values('kingdomid__kingdomName', 'kingdomid', 'phylaid__phylaName', 'phylaid', 'classid__className', 'classid', 'orderid__orderName', 'orderid', 'familyid__familyName', 'familyid', 'genusid__genusName', 'genusid', 'speciesName', 'speciesid')
    speciesDF = pd.DataFrame.from_records(qs1)
    speciesDF.rename(columns={'kingdomid__kingdomName': 'Kingdom Name', 'kingdomid': 'Kingdom ID'}, inplace=True)
    speciesDF.rename(columns={'phylaid__phylaName': 'Phylum Name', 'phylaid': 'Phylum ID'}, inplace=True)
    speciesDF.rename(columns={'classid__className': 'Class Name', 'classid': 'Class ID'}, inplace=True)
    speciesDF.rename(columns={'orderid__orderName': 'Order Name', 'orderid': 'Order ID'}, inplace=True)
    speciesDF.rename(columns={'familyid__familyName': 'Family Name', 'familyid': 'Family ID'}, inplace=True)
    speciesDF.rename(columns={'genusid__genusName': 'Genus Name', 'genusid': 'Genus ID'}, inplace=True)
    speciesDF.rename(columns={'speciesName': 'Species Name', 'speciesid': 'Species ID'}, inplace=True)
    table = speciesDF.to_html(classes="table display", columns=['Kingdom Name', 'Kingdom ID', 'Phylum Name', 'Phylum ID', 'Class Name', 'Class ID', 'Order Name', 'Order ID', 'Family Name', 'Family ID', 'Genus Name', 'Genus ID', 'Species Name', 'Species ID'])
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


def instructions(request):
    filename = "instructions/Manual.pdf"
    wrapper = FileWrapper(file(filename, 'rb'))
    response = HttpResponse(wrapper, content_type='application/pdf')
    response['Content-Disposition'] = 'attachment; filename="Manual.pdf"'
    response['Content-Length'] = os.path.getsize(filename)
    return response


def project_file1(request):
    filename = "sample_files/Example1.project.csv"
    wrapper = FileWrapper(file(filename))
    response = HttpResponse(wrapper, content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="Example1.project.csv"'
    response['Content-Length'] = os.path.getsize(filename)
    return response


def sample_file1(request):
    filename = "sample_files/Example1.sample.csv"
    wrapper = FileWrapper(file(filename))
    response = HttpResponse(wrapper, content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="Example1.sample.csv"'
    response['Content-Length'] = os.path.getsize(filename)
    return response


def shared_file(request):
    filename = "sample_files/final.pds.wang.tx.shared"
    wrapper = FileWrapper(file(filename))
    response = HttpResponse(wrapper, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="final.shared"'
    response['Content-Length'] = os.path.getsize(filename)
    return response


def taxonomy_file(request):
    filename = "sample_files/final.pds.wang.tx.1.cons.taxonomy"
    wrapper = FileWrapper(file(filename))
    response = HttpResponse(wrapper, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="final.taxonomy"'
    response['Content-Length'] = os.path.getsize(filename)
    return response


def project_file2(request):
    filename = "sample_files/Example2.project.csv"
    wrapper = FileWrapper(file(filename))
    response = HttpResponse(wrapper, content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="Example2.project.csv"'
    response['Content-Length'] = os.path.getsize(filename)
    return response


def sample_file2(request):
    filename = "sample_files/Example2.sample.csv"
    wrapper = FileWrapper(file(filename))
    response = HttpResponse(wrapper, content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="Example2.sample.csv"'
    response['Content-Length'] = os.path.getsize(filename)
    return response


def sff(request):
    filename = "sample_files/Example2.sff"
    wrapper = FileWrapper(file(filename))
    response = HttpResponse(wrapper, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="Example2.sff"'
    response['Content-Length'] = os.path.getsize(filename)
    return response


def oligos(request):
    filename = "sample_files/Example2.oligos"
    wrapper = FileWrapper(file(filename))
    response = HttpResponse(wrapper, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="Example2.oligos"'
    response['Content-Length'] = os.path.getsize(filename)
    return response


def batch1(request):
    filename = "sample_files/Example2.mothur_win.batch"
    wrapper = FileWrapper(file(filename))
    response = HttpResponse(wrapper, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="Example2.mothur_win.batch"'
    response['Content-Length'] = os.path.getsize(filename)
    return response


def batch2(request):
    filename = "sample_files/Example2.mothur_linux.batch"
    wrapper = FileWrapper(file(filename))
    response = HttpResponse(wrapper, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="Example2.mothur_linux.batch"'
    response['Content-Length'] = os.path.getsize(filename)
    return response