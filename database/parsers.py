import csv
import pandas as pd
import re
import simplejson
from django.http import HttpResponse
from models import Project, Sample, Collect, Climate, Soil_class, Soil_nutrient, Management, Microbial, User
from models import Kingdom, Phyla, Class, Order, Family, Genus, Species, Profile
from pandas import Series
from uuid import uuid4


stage = ''
perc = 0


def status(request):
    global stage, perc
    if request.is_ajax():
        dict = {}
        dict['stage'] = stage
        dict['perc'] = perc
        json_data = simplejson.dumps(dict, encoding="Latin-1")
        return HttpResponse(json_data, content_type='application/json')


def projectid(Document):
    f = csv.DictReader(Document, delimiter=',')
    for row in f:
        if row:
            index = row['project_id']
            if Project.objects.filter(projectid=index).exists():
                return index
            else:
                return uuid4().hex


def parse_project(Document, path, p_uuid):
    global stage, perc
    stage = "Step 1 of 4: Parsing project file..."
    perc = 0

    f = csv.DictReader(Document, delimiter=',')
    for row in f:
        if row:
            perc = 100
            row_dict = dict((k, v) for k, v in row.iteritems() if v != '')
            row_dict.pop('project_id')
            m = Project(projectid=p_uuid, path=path, **row_dict)
            m.save()
    Document.close()


def parse_sample(Document, p_uuid):
    global stage, perc
    stage = "Step 2 of 4: Parsing sample file..."
    perc = 0
    step = 0.0
    total = 0.0
    f = csv.reader(Document)
    f.next()
    for row in f:
        if row:
            total += 1.0

    f = csv.DictReader(Document, delimiter=',')
    for row in f:
        if row:
            row_dict = dict((k, v) for k, v in row.iteritems() if v != '')
            step += 1.0
            perc = int(step / total * 100)
            index = row_dict['sample_id']
            if not Sample.objects.filter(sampleid=index).exists():
                s_uuid = uuid4().hex
                row_dict.pop('sample_id')

                project = Project.objects.get(projectid=p_uuid)
                wanted_keys = ['sample_name', 'organism', 'title', 'seq_method', 'collection_date', 'biome', 'feature', 'geo_loc_country', 'geo_loc_state', 'geo_loc_city', 'geo_loc_farm', 'geo_loc_plot', 'latitude', 'longitude', 'material', 'elevation']
                sampleDict = {x: row_dict[x] for x in wanted_keys if x in row_dict}
                m = Sample(projectid=project, sampleid=s_uuid, **sampleDict)
                m.save()

                sample = Sample.objects.get(sampleid=s_uuid)

                wanted_keys = ['depth', 'pool_dna_extracts', 'samp_size', 'samp_collection_device', 'samp_weight_dna_ext', 'sieving', 'storage_cond']
                collectDict = {x: row_dict[x] for x in wanted_keys if x in row_dict}
                m = Collect(projectid=project, sampleid=sample, **collectDict)
                m.save()

                wanted_keys = ['annual_season_precpt', 'annual_season_temp']
                climateDict = {x: row_dict[x] for x in wanted_keys if x in row_dict}
                m = Climate(projectid=project, sampleid=sample, **climateDict)
                m.save()

                wanted_keys = ['bulk_density', 'drainage_class', 'fao_class', 'horizon', 'local_class', 'porosity', 'profile_position', 'slope_aspect', 'slope_gradient', 'soil_type', 'texture_class', 'water_content_soil']
                soil_classDict = {x: row_dict[x] for x in wanted_keys if x in row_dict}
                m = Soil_class(projectid=project, sampleid=sample, **soil_classDict)
                m.save()

                wanted_keys = ['pH', 'EC', 'tot_C', 'tot_OM', 'tot_N', 'NO3_N', 'NH4_N', 'P', 'K', 'S', 'Zn', 'Fe', 'Cu', 'Mn', 'Ca', 'Mg', 'Na', 'B']
                soil_nutrDict = {x: row_dict[x] for x in wanted_keys if x in row_dict}
                m = Soil_nutrient(projectid=project, sampleid=sample, **soil_nutrDict)
                m.save()

                wanted_keys = ['agrochem_amendments', 'agrochem_amendments_desc', 'biological_amendments', 'biological_amendments_desc', 'cover_crop', 'crop_rotation', 'cur_land_use', 'cur_vegetation', 'cur_crop', 'cur_cultivar', 'organic', 'previous_land_use', 'soil_amendments', 'soil_amendments_desc', 'tillage']
                mgtDict = {x: row_dict[x] for x in wanted_keys if x in row_dict}
                m = Management(projectid=project, sampleid=sample, **mgtDict)
                m.save()

                wanted_keys = ['rRNA_copies', 'microbial_biomass_C', 'microbial_biomass_N', 'microbial_respiration']
                microbeDict = {x: row_dict[x] for x in wanted_keys if x in row_dict}
                m = Microbial(projectid=project, sampleid=sample, **microbeDict)
                m.save()

                wanted_keys = ['usr_cat1', 'usr_cat2', 'usr_cat3', 'usr_cat4', 'usr_cat5', 'usr_cat6', 'usr_quant1', 'usr_quant2', 'usr_quant3', 'usr_quant4', 'usr_quant5', 'usr_quant6']
                userDict = {x: row_dict[x] for x in wanted_keys if x in row_dict}
                m = User(projectid=project, sampleid=sample, **userDict)
                m.save()


def parse_taxonomy(Document):
    global stage, perc
    stage = "Step 3 of 4: Parsing taxonomy file..."
    perc = 0
    step = 0.0
    total = 0.0
    f = csv.reader(Document)
    f.next()
    for row in f:
        if row:
            total += 1.0

    f = csv.reader(Document, delimiter='\t')
    f.next()
    for row in f:
        if row:
            step += 1.0
            perc = int(step / total * 100)
            subbed = re.sub(r'(\(.*?\)|k__|p__|c__|o__|f__|g__|s__)', '', row[2])
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
            g = Genus.objects.get(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusName=taxon[5]).genusid

            if not Species.objects.filter(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusid_id=g, speciesName=taxon[6]).exists():
                sid = uuid4().hex
                record = Species(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusid_id=g, speciesid=sid, speciesName=taxon[6])
                record.save()


def parse_profile(file3, file4, p_uuid):
    global stage, perc
    stage = "Step 4 of 4: Parsing shared file..."
    perc = 0
    step = 0.0
    total = 0.0

    f1 = csv.DictReader(file3, delimiter='\t')
    rows_list1 = []
    for row in f1:
        if row:
            row_dict = dict((k, v) for k, v in row.iteritems() if v)
            rows_list1.append(row_dict)
    df1 = pd.DataFrame(rows_list1)
    df1.set_index('OTU', inplace=True)
    file3.close()
    df1.drop('Size', axis=1, inplace=True)

    f2 = csv.DictReader(file4, delimiter='\t')
    rows_list2 = []
    for row in f2:
        if row:
            row_dict = dict((k, v) for k, v in row.iteritems() if v != '')
            rows_list2.append(row_dict)
    df2 = pd.DataFrame(rows_list2)
    df2.set_index('Group', inplace=True)
    file4.close()
    df2.drop(['label', 'numOtus'], axis=1, inplace=True)
    df3 = df2.T

    df4 = df1.merge(df3, left_index=True, right_index=True, how='outer')
    df4['Taxonomy'].replace(to_replace='(\(.*?\)|k__|p__|c__|o__|f__|g__|s__)', value='', regex=True, inplace=True)
    df4.reset_index(drop=True, inplace=True)

    df4.set_index(['Taxonomy'], inplace=True)
    df5 = df4.unstack().reset_index(name='count')
    df5.rename(columns={'level_0': 'sample'}, inplace=True)

    df6 = df5.groupby('sample')['count'].sum().reset_index()
    df6.rename(columns={'count': 'total'}, inplace=True)
    df7 = df5.merge(df6, on='sample', how='outer')

    df8 = df7['Taxonomy'].str.split(';').apply(Series, 1)
    df8.rename(columns={0: 'kingdom', 1: 'phyla', 2: 'class', 3: 'order', 4: 'family', 5: 'genus', 6: 'species'}, inplace=True)
    df9 = df7.join(df8, how='outer')

    for row in df9.iterrows():
        total += 1.0

    for index, row in df9.iterrows():
        step += 1.0
        perc = int(step / total * 100)

        if row['count'] > 0:
            name = row['sample']
            k = row['kingdom']
            p = row['phyla']
            c = row['class']
            o = row['order']
            f = row['family']
            g = row['genus']
            s = row['species']
            count = row['count']
            project = Project.objects.get(projectid=p_uuid)
            sample = Sample.objects.filter(projectid=p_uuid).get(sample_name=name)
            t_kingdom = Kingdom.objects.get(kingdomName=k)
            t_phyla = Phyla.objects.get(kingdomid_id=t_kingdom, phylaName=p)
            t_class = Class.objects.get(kingdomid_id=t_kingdom, phylaid_id=t_phyla, className=c)
            t_order = Order.objects.get(kingdomid_id=t_kingdom, phylaid_id=t_phyla, classid_id=t_class, orderName=o)
            t_family = Family.objects.get(kingdomid_id=t_kingdom, phylaid_id=t_phyla, classid_id=t_class, orderid_id=t_order, familyName=f)
            t_genus = Genus.objects.get(kingdomid_id=t_kingdom, phylaid_id=t_phyla, classid_id=t_class, orderid_id=t_order, familyid_id=t_family, genusName=g)
            t_species = Species.objects.get(kingdomid_id=t_kingdom, phylaid_id=t_phyla, classid_id=t_class, orderid_id=t_order, familyid_id=t_family, genusid_id=t_genus, speciesName=s)

            if Profile.objects.filter(sampleid_id=sample, kingdomid_id=t_kingdom, phylaid_id=t_phyla, classid_id=t_class, orderid_id=t_order, familyid_id=t_family, genusid_id=t_genus, speciesid_id=t_species).exists():
                t = Profile.objects.get(sampleid_id=sample, kingdomid_id=t_kingdom, phylaid_id=t_phyla, classid_id=t_class, orderid_id=t_order, familyid_id=t_family, genusid_id=t_genus, speciesid_id=t_species)
                old = t.count
                new = old + int(count)
                t.count = new
                t.save()
            else:
                record = Profile(projectid=project, sampleid=sample, kingdomid=t_kingdom, phylaid=t_phyla, classid=t_class, orderid=t_order, familyid=t_family, genusid=t_genus, speciesid=t_species, count=count)
                record.save()