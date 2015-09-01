import csv
import pandas as pd
import re
import simplejson
from django.http import HttpResponse
from models import Project, Reference, Sample, Soil, Human_Gut, Microbial, User, Human_Associated, Air, Water
from models import Kingdom, Phyla, Class, Order, Family, Genus, Species, Profile
from utils import remove_proj
from uuid import uuid4
import numpy as np
from numpy import *
import subprocess
import glob
import os
import shutil


stage = ''
perc = 0
rep_project = ''


def mothur(dest):
    global stage, perc
    stage = "Step 3 of 5: Running mothur...please check your host terminal for progress!"
    perc = 0

    if os.name == 'nt':
        subprocess.call('mothur/mothur-win/mothur.exe mothur/temp/mothur.batch')
    else:
        subprocess.call('mothur/mothur-linux/mothur mothur/temp/mothur.batch', shell=True)

    shutil.copy('mothur/temp/temp.sff', '% s/mothur.sff' % dest)
    shutil.copy('mothur/temp/temp.oligos', '% s/mothur.oligos' % dest)
    shutil.copy('mothur/temp/mothur.batch', '% s/mothur.batch' % dest)
    shutil.copy('mothur/temp/final.fasta', '% s/final.fasta' % dest)
    shutil.copy('mothur/temp/final.names', '% s/final.names' % dest)
    shutil.copy('mothur/temp/final.groups', '% s/final.groups' % dest)
    shutil.copy('mothur/temp/final.taxonomy', '% s/mothur.taxonomy' % dest)
    shutil.copy('mothur/temp/final.shared', '% s/mothur.shared' % dest)

    shutil.rmtree('mothur/temp')

    for fl in glob.glob("*.logfile"):
        os.remove(fl)


def status(request):
    if request.is_ajax():
        dict = {}
        dict['stage'] = stage
        dict['perc'] = perc
        dict['project'] = rep_project
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


def parse_project(Document, path, p_uuid, pType):
    global stage, perc
    stage = "Step 1 of 5: Parsing project file..."
    perc = 0

    f = csv.DictReader(Document, delimiter=',')
    for row in f:
        if row:
            perc = 100
            row_dict = dict((k, v) for k, v in row.iteritems() if v != '')
            row_dict.pop('project_id')
            m = Project(projectType=pType, projectid=p_uuid, path=path, **row_dict)
            m.save()

    try:
        Document.seek(0)
        f = csv.reader(Document, delimiter=',')
        if not os.path.exists(path):
            os.makedirs(path)
        with open('% s/project.csv' % path, 'wb+') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',')
            for stuff in f:
                if str(stuff[0]) == "null":
                    stuff[0] = str(p_uuid)
                spamwriter.writerow(stuff)
    except Exception as e:
        print(e)


def parse_reference(p_uuid, path, file7, raw):
    refid = uuid4().hex
    project = Project.objects.get(projectid=p_uuid)

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

        m = Reference(refid=refid, projectid=project, path=path, alignDB=align_ref, templateDB=template_ref, taxonomyDB=taxonomy_ref)
        m.save()


def parse_sample(Document, p_uuid, path, pType):
    global stage, perc
    stage = "Step 2 of 5: Parsing sample file..."
    perc = 0

    sampleidlist = []

    f = csv.reader(Document, delimiter=',')
    f.next()
    total = 0.0
    for row in f:
        if row:
            total += 1.0

    Document.seek(0)
    f = csv.DictReader(Document, delimiter=',')
    step = 0.0
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

                sampleidlist.append(s_uuid)
                sample = Sample.objects.get(sampleid=s_uuid)

                if pType == "soil":
                    wanted_keys = ['depth', 'pool_dna_extracts', 'samp_size', 'samp_collection_device', 'samp_weight_dna_ext', 'sieving', 'storage_cond', 'annual_season_precpt', 'annual_season_temp', 'bulk_density', 'drainage_class', 'fao_class', 'horizon', 'local_class', 'porosity', 'profile_position', 'slope_aspect', 'slope_gradient', 'soil_type', 'texture_class', 'water_content_soil', 'pH', 'EC', 'tot_C', 'tot_OM', 'tot_N', 'NO3_N', 'NH4_N', 'P', 'K', 'S', 'Zn', 'Fe', 'Cu', 'Mn', 'Ca', 'Mg', 'Na', 'B', 'agrochem_amendments', 'agrochem_amendments_desc', 'biological_amendments', 'biological_amendments_desc', 'cover_crop', 'crop_rotation', 'cur_land_use', 'cur_vegetation', 'cur_crop', 'cur_cultivar', 'organic', 'previous_land_use', 'soil_amendments', 'soil_amendments_desc', 'tillage', 'rRNA_copies', 'microbial_biomass_C', 'microbial_biomass_N', 'microbial_respiration']
                    soilDict = {x: row_dict[x] for x in wanted_keys if x in row_dict}
                    m = Soil(projectid=project, sampleid=sample, **soilDict)
                    m.save()

                if pType == "human_gut":
                    wanted_keys = ['age', 'body_mass_index', 'body_product', 'chem_administration', 'diet', 'disease', 'ethnicity', 'family_relationship', 'gastrointest_disord', 'genotype', 'height', 'host_body_temp', 'host_subject_id', 'ihmc_medication_code', 'last_meal', 'liver_disord', 'medic_hist_perform', 'nose_throat_disord', 'occupation', 'organism_count', 'oxy_stat_samp', 'perturbation', 'phenotype', 'pulse', 'rel_to_oxygen', 'samp_collect_device', 'samp_mat_process', 'samp_salinity', 'samp_size', 'samp_store_loc', 'samp_store_temp', 'sex', 'special_diet', 'temp', 'tissue', 'tot_mass', 'user_defined']
                    gutDict = {x: row_dict[x] for x in wanted_keys if x in row_dict}
                    m = Human_Gut(projectid=project, sampleid=sample, **gutDict)
                    m.save()

                if pType == "microbial":
                    wanted_keys = ['alkalinity', 'alkyl_diethers', 'altitude', 'aminopept_act', 'ammonium', 'bacteria_carb_prod', 'biomass', 'bishomohopanol', 'bromide', 'calcium', 'carb_nitro_ratio', 'chem_administration', 'chloride', 'chlorophyll', 'diether_lipids', 'diss_carb_dioxide', 'diss_hydrogen', 'diss_inorg_carb', 'diss_org_carb', 'diss_org_nitro', 'diss_oxygen', 'glucosidase_act', 'magnesium', 'mean_frict_vel', 'mean_peak_frict_vel', 'methane', 'n_alkanes', 'nitrate', 'nitrite', 'nitro', 'org_carb', 'org_matter', 'org_nitro', 'organism_count', 'oxy_stat_samp', 'part_org_carb', 'perturbation', 'petroleum_hydrocarb', 'ph', 'phaeopigments', 'phosphate', 'phosplipid_fatt_acid', 'potassium', 'pressure', 'redox_potential', 'rel_to_oxygen', 'salinity', 'samp_collect_device', 'samp_mat_process', 'samp_size', 'samp_store_dur', 'samp_store_loc', 'samp_store_temp', 'silicate', 'sodium', 'sulfate', 'sulfide', 'temp', 'tot_carb', 'tot_nitro', 'tot_org_carb', 'turbidity', 'water_content', 'user_defined']
                    gutDict = {x: row_dict[x] for x in wanted_keys if x in row_dict}
                    m = Microbial(projectid=project, sampleid=sample, **gutDict)
                    m.save()

                if pType == "human_associated":
                    wanted_keys = ['age', 'amniotic_fluid_color', 'blood_blood_disord', 'body_mass_index', 'body_product', 'chem_administration', 'diet', 'diet_last_six_month', 'disease', 'drug_usage', 'ethnicity', 'family_relationship', 'fetal_health_stat', 'genotype', 'gestation_state', 'height', 'hiv_stat', 'host_body_temp', 'host_subject_id', 'ihmc_medication_code', 'kidney_disord', 'last_meal', 'maternal_health_stat', 'medic_hist_perform', 'nose_throat_disord', 'occupation', 'perturbation', 'pet_farm_animal', 'phenotype', 'pulmonary_disord', 'pulse', 'rel_to_oxygen', 'samp_collect_device', 'samp_mat_process', 'samp_salinity', 'samp_size', 'samp_store_dur', 'samp_store_loc', 'samp_store_temp', 'sex', 'smoker', 'study_complt_stat', 'temp', 'tissue', 'tot_mass', 'travel_out_six_month', 'twin_sibling', 'urine_collect_meth', 'urogenit_tract_disor', 'weight_loss_3_month', 'user_defined']
                    gutDict = {x: row_dict[x] for x in wanted_keys if x in row_dict}
                    m = Human_Associated(projectid=project, sampleid=sample, **gutDict)
                    m.save()

                if pType == "air":
                    wanted_keys = ['barometric_press', 'carb_dioxide', 'carb_monoxide', 'chem_administration', 'elev', 'humidity', 'methane', 'organism_count', 'oxy_stat_samp', 'oxygen', 'perturbation', 'pollutants', 'rel_to_oxygen', 'resp_part_matter', 'samp_collect_device', 'samp_mat_process', 'samp_salinity', 'samp_size', 'samp_store_dur', 'samp_store_loc', 'samp_store_temp', 'solar_irradiance', 'temp', 'ventilation_rate', 'ventilation_type', 'volatile_org_comp', 'wind_direction', 'wind_speed', 'user_defined']
                    gutDict = {x: row_dict[x] for x in wanted_keys if x in row_dict}
                    m = Air(projectid=project, sampleid=sample, **gutDict)
                    m.save()

                if pType == "water":
                    wanted_keys = ['alkalinity', 'alkyl_diethers', 'altitude', 'aminopept_act', 'ammonium', 'atmospheric_data', 'bac_prod', 'bac_resp', 'bacteria_carb_prod', 'biomass', 'bishomohopanol', 'bromide', 'calcium', 'carb_nitro_ratio', 'chem_administration', 'chloride', 'chlorophyll', 'conduc', 'density', 'diether_lipids', 'diss_carb_dioxide', 'diss_hydrogen', 'diss_inorg_carb', 'diss_inorg_nitro', 'diss_inorg_phosp', 'diss_org_carb', 'diss_org_nitro', 'diss_oxygen', 'down_par', 'elev', 'fluor', 'glucosidase_act', 'light_intensity', 'magnesium', 'mean_frict_vel', 'mean_peak_frict_vel', 'n_alkanes', 'nitrate', 'nitrite', 'nitro', 'org_carb', 'org_matter', 'org_nitro', 'organism_count', 'oxy_stat_samp', 'part_org_carb', 'part_org_nitro', 'perturbation', 'petroleum_hydrocarb', 'ph', 'phaeopigments', 'phosphate', 'phosplipid_fatt_acid', 'photon_flux', 'potassium', 'pressure', 'primary_prod', 'redox_potential', 'rel_to_oxygen', 'samp_mat_process', 'samp_salinity', 'samp_size', 'samp_store_dur', 'samp_store_loc', 'samp_store_temp', 'samp_vol_we_dna_ext', 'silicate', 'sodium', 'soluble_react_phosp', 'source_material_id', 'sulfate', 'sulfide', 'suspend_part_matter', 'temp', 'tidal_stage', 'tot_depth_water_col', 'tot_diss_nitro', 'tot_inorg_nitro', 'tot_nitro', 'tot_part_carb', 'tot_phosp', 'water_current', 'user_defined']
                    gutDict = {x: row_dict[x] for x in wanted_keys if x in row_dict}
                    m = Water(projectid=project, sampleid=sample, **gutDict)
                    m.save()

                wanted_keys = ['usr_cat1', 'usr_cat2', 'usr_cat3', 'usr_cat4', 'usr_cat5', 'usr_cat6', 'usr_quant1', 'usr_quant2', 'usr_quant3', 'usr_quant4', 'usr_quant5', 'usr_quant6']
                userDict = {x: row_dict[x] for x in wanted_keys if x in row_dict}
                m = User(projectid=project, sampleid=sample, **userDict)
                m.save()

    try:
        Document.seek(0)
        f = csv.reader(Document, delimiter=',')
        if not os.path.exists(path):
            os.makedirs(path)
        with open('% s/sample.csv' % path, 'wb+') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',')
            itera = 0
            for stuff in f:
                if str(stuff[0]) == "null":
                    stuff[0] = str(sampleidlist[itera])
                    itera += 1
                spamwriter.writerow(stuff)
    except Exception as e:
        print(e)


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


def parse_profile(file3, file4, p_uuid):
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
                sample = Sample.objects.filter(projectid=p_uuid).get(sample_name=name)

                if Profile.objects.filter(sampleid_id=sample, kingdomid_id=t_kingdom, phylaid_id=t_phyla, classid_id=t_class, orderid_id=t_order, familyid_id=t_family, genusid_id=t_genus, speciesid_id=t_species).exists():
                    t = Profile.objects.get(sampleid_id=sample, kingdomid_id=t_kingdom, phylaid_id=t_phyla, classid_id=t_class, orderid_id=t_order, familyid_id=t_family, genusid_id=t_genus, speciesid_id=t_species)
                    old = t.count
                    new = old + int(count)
                    t.count = new
                    t.save()
                else:
                    record = Profile(projectid=project, sampleid=sample, kingdomid=t_kingdom, phylaid=t_phyla, classid=t_class, orderid=t_order, familyid=t_family, genusid=t_genus, speciesid=t_species, count=count)
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

        projects = Reference.objects.all().filter(projectid_id__in=ids)
        for project in projects:
            rep_project = 'myPhyloDB is currently reprocessing project: ' + str(project.projectid.project_name)
            dest = project.projectid.path
            shutil.copy('% s/mothur.sff' % dest, '% s/temp.sff' % mothurdest)
            shutil.copy('% s/mothur.oligos' % dest, '% s/temp.oligos' % mothurdest)

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
            dest = project.path
            pType = project.projectid.projectType

            shutil.copy('% s/project.csv' % dest, '% s/project.csv' % mothurdest)
            shutil.copy('% s/sample.csv' % dest, '% s/sample.csv' % mothurdest)

            if not os.path.exists(dest):
                os.makedirs(dest)

            try:
                mothur(dest)
            except Exception as e:
                print("Error with mothur: " + str(e))
                return HttpResponse(
                    simplejson.dumps({"error": "yes"}),
                    content_type="application/json"
                )

            remove_proj(p_uuid)

            with open('mothur/temp/project.csv', 'rb') as file1:
                parse_project(file1, dest, p_uuid, pType)

            with open('mothur/temp/sample.csv', 'rb') as file2:
                parse_sample(file2, p_uuid, dest, pType)

            with open('mothur/temp/mothur.batch', 'rb') as file7:
                raw = True
                parse_reference(p_uuid, dest, file7, raw)

            with open('% s/mothur.taxonomy' % dest, 'rb') as file3:
                parse_taxonomy(file3)

            with open('% s/mothur.taxonomy' % dest, 'rb') as file3:
                with open('% s/mothur.shared' % dest, 'rb') as file4:
                    parse_profile(file3, file4, p_uuid)

        return HttpResponse(
            simplejson.dumps({"error": "no"}),
            content_type="application/json"
        )
