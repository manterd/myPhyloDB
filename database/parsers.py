import csv
import pandas as pd
import re
import simplejson
from django.http import HttpResponse
from models import Project, Reference, Sample, Soil, Human_Associated, UserDefined
from models import Kingdom, Phyla, Class, Order, Family, Genus, Species, Profile
from utils import remove_proj, purge
from uuid import uuid4
import numpy as np
from numpy import *
import glob
import os
import shutil
from django.contrib.auth.models import User
import xlrd
from xlutils.copy import copy
import xlwt
from xlwt import Style


stage = ''
perc = 0
rep_project = ''


def mothur(dest, source):
    global stage, perc
    stage = "Step 3 of 5: Running mothur...please check your host terminal for progress!"
    perc = 0

    if os.name == 'nt':
        try:
            os.system('"mothur\\mothur-win\\mothur.exe mothur\\temp\\mothur.batch"')
        except Exception as e:
            print "Mothur failed: " + str(e)
    else:
        try:
            os.system("mothur\/mothur-linux\/mothur mothur\/temp\/mothur.batch")
        except Exception as e:
            print "Mothur failed: " + str(e)

    if source == '454':
        shutil.copy('mothur/temp/temp.sff', '% s/mothur.sff' % dest)
        shutil.copy('mothur/temp/temp.oligos', '% s/mothur.oligos' % dest)
        shutil.copy('mothur/temp/mothur.batch', '% s/mothur.batch' % dest)
        shutil.copy('mothur/temp/final.fasta', '% s/final.fasta' % dest)
        shutil.copy('mothur/temp/final.names', '% s/final.names' % dest)
        shutil.copy('mothur/temp/final.groups', '% s/final.groups' % dest)
        shutil.copy('mothur/temp/final.taxonomy', '% s/mothur.taxonomy' % dest)
        shutil.copy('mothur/temp/final.shared', '% s/mothur.shared' % dest)

        shutil.rmtree('mothur/temp')

    if source == 'miseq':
        shutil.copy('mothur/temp/temp.files', '% s/final.files' % dest)
        shutil.copy('mothur/temp/mothur.batch', '% s/mothur.batch' % dest)
        shutil.copy('mothur/temp/final.fasta', '% s/final.fasta' % dest)
        shutil.copy('mothur/temp/final.names', '% s/final.names' % dest)
        shutil.copy('mothur/temp/final.groups', '% s/final.groups' % dest)
        shutil.copy('mothur/temp/final.taxonomy', '% s/mothur.taxonomy' % dest)
        shutil.copy('mothur/temp/final.shared', '% s/mothur.shared' % dest)

        for afile in glob.glob(r'mothur/temp/*.fastq'):
            shutil.copy(afile, dest)

        shutil.rmtree('mothur/temp')

    purge('mothur/reference/align', '.8mer')
    purge('mothur/reference/taxonomy', 'numNonZero')
    purge('mothur/reference/taxonomy', '.8mer.prob')
    purge('mothur/reference/taxonomy', '.tree.sum')
    purge('mothur/reference/taxonomy', '.tree.train')
    purge('mothur/reference/template', '.8mer')
    purge('mothur/reference/template', '.summary')

    for afile in glob.glob(r'*.logfile'):
        shutil.move(afile, dest)


def status(request):
    if request.is_ajax():
        dict = {}
        dict['stage'] = stage
        dict['perc'] = perc
        dict['project'] = rep_project
        json_data = simplejson.dumps(dict, encoding="Latin-1")
        return HttpResponse(json_data, content_type='application/json')


def projectid(Document):
    f = xlrd.open_workbook(file_contents=Document.read())
    sheet = f.sheet_by_name('Project')
    pType = sheet.cell_value(rowx=5, colx=2)

    projectid = sheet.cell_value(rowx=5, colx=3)
    num_samp = int(sheet.cell_value(rowx=5, colx=0))

    if projectid == '':
        return uuid4().hex, pType, num_samp
    else:
        return projectid, pType, num_samp


def parse_project(Document, p_uuid):
    global stage, perc
    stage = "Step 1 of 5: Parsing project file..."
    perc = 0

    df = pd.read_excel(Document, skiprows=4, sheetname='Project')
    rowDict = df.to_dict(outtype='records')[0]
    rowDict.pop('num_samp')
    for key in rowDict.keys():
        if key == 'projectid':
            rowDict[key] = p_uuid

    m = Project(**rowDict)
    m.save()


def parse_reference(p_uuid, refid, path, file7, raw, source, userid):
    project = Project.objects.get(projectid=p_uuid)
    author = User.objects.get(id=userid)

    if raw:
        align_ref = ''
        template_ref = ''
        taxonomy_ref = ''
        for row in file7:
            if "align.seqs" in row:
                for item in row.split(','):
                    if "reference=" in item:
                        string = item.split('=')
                        align_ref = string[1].replace('mothur/reference/align/', '')
            if "classify.seqs" in row:
                for item in row.split(','):
                    if "template=" in item:
                        string = item.split('=')
                        template_ref = string[1].replace('mothur/reference/template/', '')
                    if "taxonomy=" in item:
                        string = item.split('=')
                        taxonomy_ref = string[1].replace('mothur/reference/taxonomy/', '')

        m = Reference(refid=refid, projectid=project, path=path, source=source, raw=True, alignDB=align_ref, templateDB=template_ref, taxonomyDB=taxonomy_ref, author=author)
        m.save()
    else:
        m = Reference(refid=refid, projectid=project, path=path, source=source, raw=False, alignDB='null', templateDB='null', taxonomyDB='null', author=author)
        m.save()


def parse_sample(Document, p_uuid, refid, pType, total):
    global stage, perc
    stage = "Step 2 of 5: Parsing sample file..."
    perc = 0
    step = 0

    project = Project.objects.get(projectid=p_uuid)
    ref = Reference.objects.get(refid=refid)

    df1 = pd.read_excel(Document, skiprows=5, sheetname='MIMARKs')
    if pType == 'human associated':
        df2 = pd.read_excel(Document, skiprows=5, sheetname='Human Associated')
    elif pType == 'soil':
        df2 = pd.read_excel(Document, skiprows=5, sheetname='Soil')
    else:
        df2 = pd.DataFrame()

    df3 = pd.read_excel(Document, skiprows=5, sheetname='User')

    idList = []
    for i in xrange(total):
        step += 1.0
        perc = int((step / total/2) * 100)

        row = df1.iloc[[i]].to_dict(outtype='records')[0]
        s_uuid = row['sampleid']
        row.pop('seq_method')
        row.pop('geo_loc_name')
        row.pop('lat_lon')

        if not Sample.objects.filter(sampleid=s_uuid).exists():
            s_uuid = uuid4().hex
            row['sampleid'] = s_uuid
            idList.append(s_uuid)
            m = Sample(projectid=project, refid=ref, **row)
            m.save()
        else:
            idList.append(s_uuid)
            m = Sample(projectid=project, refid=ref, **row)
            m.save()

        sample = Sample.objects.get(sampleid=s_uuid)

        if pType == "human associated":
            row = df2.iloc[[i]].to_dict(outtype='records')[0]
            row.pop('sampleid')
            row.pop('sample_name')
            m = Human_Associated(projectid=project, refid=ref, sampleid=sample, **row)
            m.save()
        elif pType == "soil":
            row = df2.iloc[[i]].to_dict(outtype='records')[0]
            row.pop('sampleid')
            row.pop('sample_name')
            m = Soil(projectid=project, refid=ref, sampleid=sample, **row)
            m.save()
        else:
            placeholder = ''

        row = df3.iloc[[i]].to_dict(outtype='records')[0]
        row.pop('sampleid')
        row.pop('sample_name')
        m = UserDefined(projectid=project, refid=ref, sampleid=sample, **row)
        m.save()

    rb = xlrd.open_workbook(Document, formatting_info=True)
    nSheets = rb.nsheets
    wb = copy(rb)

    style = xlwt.XFStyle()

    # border
    borders = xlwt.Borders()
    borders.bottom = xlwt.Borders.THIN
    borders.top = xlwt.Borders.THIN
    borders.left = xlwt.Borders.THIN
    borders.right = xlwt.Borders.THIN
    style.borders = borders

    # font
    font = xlwt.Font()
    font.name = 'Calibri'
    font.height = 11 * 20
    style.font = font

    # background color
    pattern = xlwt.Pattern()
    pattern.pattern = xlwt.Pattern.SOLID_PATTERN
    pattern.pattern_fore_colour = xlwt.Style.colour_map['gray25']
    style.pattern = pattern

    for each in xrange(nSheets):
        ws = wb.get_sheet(each)
        if ws.name == 'Project':
            ws.write(5, 2, p_uuid, style)

        if ws.name == 'MIMARKs':
            for i in xrange(total):
                j = i + 6
                ws.write(j, 0, idList[i], style)

        if pType == 'human associated':
            if ws.name == 'Human Associated':
                for i in xrange(total):
                    j = i + 6
                    ws.write(j, 0, idList[i], style)
        elif pType == 'soil':
            if ws.name == 'Soil':
                for i in xrange(total):
                    j = i + 6
                    ws.write(j, 0, idList[i], style)
        else:
            placeholder = ''

        if ws.name == 'User':
            for i in xrange(total):
                j = i + 6
                ws.write(j, 0, idList[i], style)

        wb.save(Document)
        perc += 50/nSheets


def parse_taxonomy(Document):
    global stage, perc
    stage = "Step 4 of 5: Parsing taxonomy file..."
    perc = 0

    f = csv.reader(Document, delimiter='\t')
    f.next()
    total = 0.0
    for row in f:
        if row:
            total += 1.0
    Document.seek(0)
    f.next()
    step = 0.0
    for row in f:
        if row:
            step += 1.0
            perc = int(step / total * 100)
            subbed = re.sub(r'(\(.*?\)|k__|p__|c__|o__|f__|g__|s__)', '', row[2])
            subbed = subbed[:-1]
            taxon = subbed.split(';')

            if not Kingdom.objects.filter(kingdomName=taxon[0]).exists():
                kid = uuid4().hex
                record = Kingdom(kingdomid=kid, kingdomName=taxon[0])
                record.save()

            k = Kingdom.objects.get(kingdomName=taxon[0]).kingdomid
            if not Phyla.objects.filter(kingdomid_id=k, phylaName=taxon[1]).exists():
                pid = uuid4().hex
                record = Phyla(kingdomid_id=k, phylaid=pid, phylaName=taxon[1])
                record.save()

            p = Phyla.objects.get(kingdomid_id=k, phylaName=taxon[1]).phylaid
            if not Class.objects.filter(kingdomid_id=k, phylaid_id=p, className=taxon[2]).exists():
                cid = uuid4().hex
                record = Class(kingdomid_id=k, phylaid_id=p, classid=cid, className=taxon[2])
                record.save()

            c = Class.objects.get(kingdomid_id=k, phylaid_id=p, className=taxon[2]).classid
            if not Order.objects.filter(kingdomid_id=k, phylaid_id=p, classid_id=c, orderName=taxon[3]).exists():
                oid = uuid4().hex
                record = Order(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid=oid, orderName=taxon[3])
                record.save()

            o = Order.objects.get(kingdomid_id=k, phylaid_id=p, classid_id=c, orderName=taxon[3]).orderid
            if not Family.objects.filter(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyName=taxon[4]).exists():
                fid = uuid4().hex
                record = Family(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid=fid, familyName=taxon[4])
                record.save()

            f = Family.objects.get(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyName=taxon[4]).familyid
            if not Genus.objects.filter(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusName=taxon[5]).exists():
                gid = uuid4().hex
                record = Genus(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusid=gid, genusName=taxon[5])
                record.save()

            try:
                g = Genus.objects.get(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusName=taxon[5]).genusid
                if not Species.objects.filter(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusid_id=g, speciesName=taxon[6]).exists():
                    sid = uuid4().hex
                    record = Species(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusid_id=g, speciesid=sid, speciesName=taxon[6])
                    record.save()
            except:
                g = Genus.objects.get(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusName=taxon[5]).genusid
                if not Species.objects.filter(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusid_id=g, speciesName='unclassified').exists():
                    sid = uuid4().hex
                    record = Species(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusid_id=g, speciesid=sid, speciesName='unclassified')
                    record.save()


def parse_profile(file3, file4, p_uuid, refid):
    global stage, perc
    stage = "Step 5 of 5: Parsing shared file..."
    perc = 0

    data1 = genfromtxt(file3, delimiter='\t', dtype=None, autostrip=True)
    arr1 = np.delete(data1, 1, axis=1)
    df1 = pd.DataFrame(arr1[1:, 1:], index=arr1[1:, 0], columns=arr1[0, 1:])
    df1 = df1[df1.index != 'False']
    a = df1.columns.values.tolist()
    a = [x for x in a if x != 'False']
    df1 = df1[a]
    file3.close()

    data2 = genfromtxt(file4, delimiter='\t', dtype=None, autostrip=True)
    arr2 = data2.T
    arr2 = np.delete(arr2, 0, axis=0)
    arr2 = np.delete(arr2, 1, axis=0)
    df2 = pd.DataFrame(arr2[1:, 1:], index=arr2[1:, 0], columns=arr2[0, 1:])
    df2 = df2[df2.index != 'False']
    b = df2.columns.values.tolist()
    b = [x for x in b if x != 'False']
    df2 = df2[b]
    file4.close()

    df3 = df1.join(df2, how='outer')
    df3['Taxonomy'].replace(to_replace='(\(.*?\)|k__|p__|c__|o__|f__|g__|s__)', value='', regex=True, inplace=True)
    df3.reset_index(drop=True, inplace=True)
    del df1, df2

    indexList = df3.index.values.tolist()
    total = len(indexList)
    sampleList = df3.columns.values.tolist()
    sampleList.remove('Taxonomy')

    step = 0.0
    for index, row in df3.iterrows():
        step += 1.0
        perc = int(step / total * 100)
        taxon = str(row['Taxonomy'])
        taxon = taxon[:-1]
        taxaList = taxon.split(';')
        k = taxaList[0]
        p = taxaList[1]
        c = taxaList[2]
        o = taxaList[3]
        f = taxaList[4]
        g = taxaList[5]
        try:
            s = taxaList[6]
        except:
            s = 'unclassified'

        t_kingdom = Kingdom.objects.get(kingdomName=k)
        t_phyla = Phyla.objects.get(kingdomid_id=t_kingdom, phylaName=p)
        t_class = Class.objects.get(kingdomid_id=t_kingdom, phylaid_id=t_phyla, className=c)
        t_order = Order.objects.get(kingdomid_id=t_kingdom, phylaid_id=t_phyla, classid_id=t_class, orderName=o)
        t_family = Family.objects.get(kingdomid_id=t_kingdom, phylaid_id=t_phyla, classid_id=t_class, orderid_id=t_order, familyName=f)
        t_genus = Genus.objects.get(kingdomid_id=t_kingdom, phylaid_id=t_phyla, classid_id=t_class, orderid_id=t_order, familyid_id=t_family, genusName=g)
        t_species = Species.objects.get(kingdomid_id=t_kingdom, phylaid_id=t_phyla, classid_id=t_class, orderid_id=t_order, familyid_id=t_family, genusid_id=t_genus, speciesName=s)

        for name in sampleList:
            count = int(row[str(name)])
            if count > 0:
                project = Project.objects.get(projectid=p_uuid)
                ref = Reference.objects.get(refid=refid)
                sample = Sample.objects.filter(projectid=p_uuid).get(sample_name=name)

                if Profile.objects.filter(sampleid_id=sample, kingdomid_id=t_kingdom, phylaid_id=t_phyla, classid_id=t_class, orderid_id=t_order, familyid_id=t_family, genusid_id=t_genus, speciesid_id=t_species).exists():
                    t = Profile.objects.get(sampleid_id=sample, kingdomid_id=t_kingdom, phylaid_id=t_phyla, classid_id=t_class, orderid_id=t_order, familyid_id=t_family, genusid_id=t_genus, speciesid_id=t_species)
                    old = t.count
                    new = old + int(count)
                    t.count = new
                    t.save()
                else:
                    record = Profile(projectid=project, refid=ref, sampleid=sample, kingdomid=t_kingdom, phylaid=t_phyla, classid=t_class, orderid=t_order, familyid=t_family, genusid=t_genus, speciesid=t_species, count=count)
                    record.save()


def reanalyze(request):
    global rep_project

    if request.is_ajax():
        mothurdest = 'mothur/temp'
        if not os.path.exists(mothurdest):
            os.makedirs(mothurdest)

        allJson = request.GET["all"]
        all = simplejson.loads(allJson)
        ids = all["ids"]
        new_align = 'reference=mothur/reference/align/' + str(all['alignDB'])
        new_taxonomy = 'taxonomy=mothur/reference/taxonomy/' + str(all['taxonomyDB'])
        new_tax = str(all['taxonomyDB'])
        new_tax_tag = str(new_tax.split('.')[-2:-1][0])
        new_template = 'template=mothur/reference/template/' + str(all['templateDB'])

        projects = Reference.objects.all().filter(refid__in=ids)

        for project in projects:
            rep_project = 'myPhyloDB is currently reprocessing project: ' + str(project.projectid.project_name)
            dest = project.path
            source = project.source

            if source == '454':
                shutil.copy('% s/mothur.sff' % dest, '% s/temp.sff' % mothurdest)
                shutil.copy('% s/mothur.oligos' % dest, '% s/temp.oligos' % mothurdest)
            if source == 'miseq':
                shutil.copy('% s/final.files' % dest, '% s/temp.files' % mothurdest)
                for afile in glob.glob(r'% s/*.fastq' % dest):
                    shutil.copy(afile, mothurdest)

            orig_align = 'reference=mothur/reference/align/' + str(project.alignDB)
            orig_taxonomy = 'taxonomy=mothur/reference/taxonomy/' + str(project.taxonomyDB)
            orig_tax = str(project.taxonomyDB)
            orig_tax_tag = str(orig_tax.split('.')[-2:-1][0])
            orig_template = 'template=mothur/reference/template/' + str(project.templateDB)

            try:
                with open("% s/mothur.batch" % dest, 'r+') as bat:
                    with open("% s/mothur.batch" % mothurdest, 'wb+') as destination:
                        method = 'wang'
                        foundClassify = False
                        orig_tag_meth = ''
                        new_tag_meth = ''
                        for line in bat:
                            if "align.seqs" in line:
                                line = line.replace(str(orig_align), str(new_align))
                            if "classify.seqs" in line:
                                line = line.replace(str(orig_template), str(new_template))
                                line = line.replace(str(orig_taxonomy), str(new_taxonomy))
                                cmds = line.split(',')
                                for item in cmds:
                                    if "method=" in item:
                                        method = item.split('=')[1]
                                orig_tag_meth = str(orig_tax_tag) + "." + str(method)
                                new_tag_meth = str(new_tax_tag) + "." + str(method)
                                foundClassify = True
                            if (str(orig_tag_meth) in line) and foundClassify:
                                line = line.replace(str(orig_tag_meth), str(new_tag_meth))
                            destination.write(line)
            except Exception as e:
                print("Error with batch file: ", e)

            p_uuid = project.projectid.projectid
            refid = project.refid
            dest = project.path
            pType = project.projectid.projectType

            shutil.copy('% s/final_meta.csv' % dest, '% s/final_meta.csv' % mothurdest)

            remove_proj(dest)

            if not os.path.exists(dest):
                os.makedirs(dest)

            shutil.copy('% s/final_meta.csv' % mothurdest, '% s/final_meta.csv' % dest)

            try:
                mothur(dest, source)
            except Exception as e:
                print("Error with mothur: " + str(e))
                return HttpResponse(
                    simplejson.dumps({"error": "yes"}),
                    content_type="application/json"
                )

            metaName = 'final_meta.xls'
            metaFile = '/'.join([dest, metaName])
            parse_project(metaFile, p_uuid)

            with open('% s/mothur.batch' % dest, 'rb') as file7:
                raw = True
                userID = str(request.user.id)
                parse_reference(p_uuid, refid, dest, file7, raw, source, userID)

            f = xlrd.open_workbook(file_contents=metaFile)
            sheet = f.sheet_by_name('Project')
            num_samp = int(sheet.cell_value(rowx=5, colx=0))
            parse_sample(metaFile, p_uuid, refid, pType, num_samp)

            with open('% s/mothur.taxonomy' % dest, 'rb') as file3:
                parse_taxonomy(file3)

            with open('% s/mothur.taxonomy' % dest, 'rb') as file3:
                with open('% s/mothur.shared' % dest, 'rb') as file4:
                    parse_profile(file3, file4, p_uuid, refid)

        return HttpResponse(
            simplejson.dumps({"error": "no"}),
            content_type="application/json"
        )
