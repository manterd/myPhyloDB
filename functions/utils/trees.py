from django.http import HttpResponse
from django.db.models import Q
from django.contrib.auth.models import User
import operator
import pandas as pd
import pickle
from pyper import *
import json

from database.models import Project, Reference, Sample, Air, Human_Associated, Microbial, Soil, Water, UserDefined, \
    Kingdom, Phyla, Class, Order, Family, Genus, Species, OTU_99, Profile, \
    ko_lvl1, ko_entry, \
    nz_lvl1, nz_entry, \
    UserProfile, DaymetData

import functions

from database import perms

import traceback


pd.set_option('display.max_colwidth', -1)


def getProjectTree(request):    # get all projects this user has permission to view (run analysis on)
    myTree = {'title': 'All Projects', 'isFolder': True, 'expand': True, 'hideCheckbox': True, 'children': []}

    projects = perms.getViewProjects(request)
    for project in projects:
        myNode = {
            'title': project.project_name,
            'tooltip': "Project type: " + project.projectType + "\nDescription: " + project.project_desc + "\nID: " + project.projectid + "\nPI: " + project.pi_first + " " + project.pi_last + "\nAffiliation: " + project.pi_affiliation + "\nOwner: " + project.owner.username,
            'id': project.projectid,
            'isFolder': True,
            'isLazy': True,
            'wip': project.wip
        }
        myTree['children'].append(myNode)
    # Convert result list to a JSON string
    res = json.dumps(myTree)

    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')


def getProjectTreeChildren(request):    # get the samples belonging to projects the user has permission for
    if request.is_ajax():
        projectid = request.GET["id"]
        samples = Sample.objects.filter(projectid=projectid).values_list('sampleid', flat=True)

        nodes = []
        reads = Sample.objects.filter(sampleid__in=samples).order_by('sample_name').values('sampleid', 'sample_name', 'reads')
        for i in reads:
            if int(i['reads'] > 0):
                myNode = {
                    'title': 'Name: ' + str(i['sample_name']) + '; Reads: ' + str(i['reads']),
                    'tooltip': 'ID: ' + str(i['sampleid']),
                    'id': str(i['sampleid']),
                    'isFolder': False
                }
                nodes.append(myNode)
        res = json.dumps(nodes)
        return HttpResponse(res, content_type='application/json')


def getSampleCatTree(request):  # get categorical headers for projects user selected
    myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
    path = str(myDir) + 'usr_norm_samples.pkl'
    with open(path, 'rb') as f:
        samples = pickle.load(f)

    # There are lists for all fields associated with a project type, split into groups, here.
    # When adding new attributes to the database (models), one must also add them here somewhere (and/or in quant equivalent later)
    # The reason for the hard-coding of categories comes down to grouping: we don't want ALL soil fields to show up in the same folder
    # so we split the full list into functional groups, which must be specified manually in order to exist (no way to measure group)
    # The only way around this might be some sort of structure in the model itself for each group, like soilHealth as a model contained by soil

    # convert samples list into django queryset based on id's in current list
    projectList = Sample.objects.filter(sampleid__in=samples).values_list('projectid').distinct()
    typeList = Project.objects.filter(projectid__in=projectList).values_list('projectType', flat=True)

    myTree = {'title': 'Meta Data: Categorical', 'id': 'root', 'isFolder': False,  'hideCheckbox': True, 'expand': True, 'children': []}
    project = {'title': 'Projects', 'id': 'project', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    mimark = {'title': 'MIMARKs', 'id': 'mimark', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    air = {'title': 'Air', 'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    human_associated = {'title': 'Human Associated', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    microbial = {'title': 'Microbial', 'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    soil = {'title': 'Soil', 'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    water = {'title': 'Water', 'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    user = {'title': 'User-defined', 'id': 'user', 'isFolder': True,  'hideCheckbox': True, 'children': []}

    list = ['project_name']
    for i in range(len(list)):
        myNode = {'title': list[i], 'id': 'project', 'isFolder': True, 'pType': 'project', 'isLazy': True, 'children': []}
        project['children'].append(myNode)

    list = ['sample_name', 'organism', 'collection_date', 'depth', 'elev', 'seq_platform', 'seq_gene', 'seq_gene_region', 'seq_barcode', 'seq_for_primer', 'seq_rev_primer', 'env_biome', 'env_feature', 'env_material', 'geo_loc_country', 'geo_loc_state', 'geo_loc_city', 'geo_loc_farm', 'geo_loc_plot']
    for i in range(len(list)):
        myNode = {'title': list[i], 'id': 'mimark', 'isFolder': True, 'pType': 'mimark', 'isLazy': True, 'children': []}
        mimark['children'].append(myNode)

    if 'air' in typeList:
        list = ['chem_admin_term', 'chem_admin_time', 'organism_type', 'oxy_stat_samp', 'perturbation_type', 'perturbation_interval', 'pollutants_type', 'rel_to_oxygen', 'resp_part_matter_substance', 'samp_collect_device', 'samp_mat_process', 'samp_store_loc', 'ventilation_type', 'volatile_org_comp_name', 'wind_direction']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'air', 'isFolder': True, 'pType': 'air', 'isLazy': True, 'children': []}
            air['children'].append(myNode)

    if 'microbial' in typeList:
        list = ['biomass_part_name', 'chem_administration_term', 'chem_administration_time', 'diether_lipids_name', 'n_alkanes_name', 'organism_name', 'oxy_stat_samp', 'perturbation_type', 'perturbation_interval', 'phaeopigments_type', 'phosplipid_fatt_acid_name', 'rel_to_oxygen', 'samp_collect_device', 'samp_mat_process', 'samp_store_loc']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'microbial', 'isFolder': True, 'pType': 'microbial', 'isLazy': True, 'children': []}
            microbial['children'].append(myNode)

    if 'water' in typeList:
        list = ['biomass_part_name', 'chem_administration_name', 'chem_administration_time', 'organism_name', 'oxy_stat_samp', 'perturbation_type', 'perturbation_interal', 'rel_to_oxygen', 'samp_mat_process', 'samp_store_loc', 'source_material_id', 'tidal_stage', 'water_current_direction']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'water', 'isFolder': True, 'pType': 'water', 'isLazy': True, 'children': []}
            water['children'].append(myNode)

    if 'human associated' in typeList:
        samp_collect = {'title': 'Sample Collection', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['samp_collect_device', 'samp_mat_process', 'samp_store_loc']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human associated', 'isLazy': True, 'children': []}
            samp_collect['children'].append(myNode)
        human_associated['children'].append(samp_collect)

        samp_class = {'title': 'Sample Classification', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['samp_type', 'samp_location', 'samp_oxy_stat']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human associated', 'isLazy': True, 'children': []}
            samp_class['children'].append(myNode)
        human_associated['children'].append(samp_class)

        host = {'title': 'Host', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['host_subject_id', 'host_gender', 'host_ethnicity', 'host_occupation', 'pet_farm_animal', 'obesity', 'smoker']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human associated', 'isLazy': True, 'children': []}
            host['children'].append(myNode)
        human_associated['children'].append(host)

        diet = {'title': 'Diet', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['diet_type', 'diet_frequency', 'diet_last_six_month', 'last_meal']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human associated', 'isLazy': True, 'children': []}
            diet['children'].append(myNode)
        human_associated['children'].append(diet)

        disease = {'title': 'Disease', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['medic_hist_perform', 'disease_type', 'disease_location', 'tumor_location', 'tumor_stage']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human associated', 'isLazy': True, 'children': []}
            disease['children'].append(myNode)
        human_associated['children'].append(disease)

        drug_use = {'title': 'Drug Usage', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['drug_usage', 'drug_type', 'drug_frequency']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human associated', 'isLazy': True, 'children': []}
            drug_use['children'].append(myNode)
        human_associated['children'].append(drug_use)

        interven = {'title': 'Intervention', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['perturbation', 'pert_type', 'pert_frequency']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human associated', 'isLazy': True, 'children': []}
            interven['children'].append(myNode)
        human_associated['children'].append(interven)

        fetal = {'title': 'Fetal', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['fetal_health_stat', 'amniotic_fluid_color', 'gestation_stat', 'maternal_health_stat']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human_associated', 'isLazy': True, 'children': []}
            fetal['children'].append(myNode)
        human_associated['children'].append(fetal)

    if 'soil' in typeList:
        samp_collect = {'title': 'Sample Collection', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['samp_collection_device', 'samp_depth', 'samp_prep', 'samp_store_loc']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            samp_collect['children'].append(myNode)
        soil['children'].append(samp_collect)

        soil_class = {'title': 'Soil Classification', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['fao_class', 'local_class', 'texture_class', 'profile_position', 'slope_aspect', 'drainage_class']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            soil_class['children'].append(myNode)
        soil['children'].append(soil_class)

        crop_info = {'title': 'Crop Information', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['cur_land_use', 'cur_vegetation', 'cur_crop', 'cur_cultivar', 'crop_rotation', 'cover_crop']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            crop_info['children'].append(myNode)
        soil['children'].append(crop_info)

        fert_mgt = {'title': 'Fertilizer Mgt', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['fert_amendment_class', 'fert_placement', 'fert_type']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            fert_mgt['children'].append(myNode)
        soil['children'].append(fert_mgt)

        irrigation = {'title': 'Irrigation', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['irrigation_type']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            irrigation['children'].append(myNode)
        soil['children'].append(irrigation)

        residue = {'title': 'Residue Mgt', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['residue_removal', 'residue_growth_stage']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            residue['children'].append(myNode)
        soil['children'].append(residue)

        tillage = {'title': 'Tillage Mgt', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['tillage_event']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            tillage['children'].append(myNode)
        soil['children'].append(tillage)

        amend = {'title': 'Amendment Mgt', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['amend1_class', 'amend1_active_ingredient', 'amend2_class', 'amend2_active_ingredient', 'amend3_class', 'amend3_active_ingredient']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            amend['children'].append(myNode)
        soil['children'].append(amend)

        biomass = {'title': 'Crop Biomass', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['harv_fraction']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            biomass['children'].append(myNode)
        soil['children'].append(biomass)

        ghg = {'title': 'GHG Flux', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['ghg_chamber_placement']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            ghg['children'].append(myNode)
        soil['children'].append(ghg)

        health = {'title': 'Soil Health', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['soil_texture_sand', 'soil_texture_silt', 'soil_texture_clay', 'soil_water_cap', 'soil_water_cap_rating'
            , 'soil_surf_hardness', 'soil_surf_hardness_rating', 'soil_subsurf_hardness', 'soil_subsurf_hardness_rating'
            , 'soil_agg_stability', 'soil_agg_stability_rating', 'soil_organic_matter', 'soil_organic_matter_rating'
            , 'soil_ACE_protein_index', 'soil_ACE_protein_index_rating', 'soil_root_pathogen_pressure'
            , 'soil_root_pathogen_pressure_rating', 'soil_respiration_four_day', 'soil_soil_respiration_four_day_rating'
            , 'soil_active_C', 'soil_active_C_rating', 'soil_pH_rating', 'soil_p_rating', 'soil_k_rating'
            , 'soil_minor_elements_rating', 'CASH_SHI_rating']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            health['children'].append(myNode)
        soil['children'].append(health)

    list = ['usr_cat1', 'usr_cat2', 'usr_cat3', 'usr_cat4', 'usr_cat5', 'usr_cat6']
    for i in range(len(list)):
        myNode = {'title': list[i], 'id': 'user', 'isFolder': True, 'pType': 'user', 'isLazy': True, 'children': []}
        user['children'].append(myNode)

    myTree['children'].append(project)
    myTree['children'].append(mimark)
    if 'air' in typeList:
        myTree['children'].append(air)
    if 'microbial' in typeList:
        myTree['children'].append(microbial)
    if 'water' in typeList:
        myTree['children'].append(water)
    if 'human associated' in typeList:
        myTree['children'].append(human_associated)
    if 'soil' in typeList:
        myTree['children'].append(soil)
    myTree['children'].append(user)

    # Convert result list to a JSON string
    res = json.dumps(myTree)

    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')


def getSampleCatTreeChildren(request):  # populate selected nodes of categorical data with values and sample info
    if request.is_ajax():
        myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
        path = str(myDir) + 'usr_norm_samples.pkl'
        with open(path, 'rb') as f:
            selected = pickle.load(f)

        filtered = []
        reads = Sample.objects.filter(sampleid__in=selected).order_by('sample_name').values('sampleid', 'sample_name', 'reads')
        for i in reads:
            if int(i['reads']) > 0:
                filtered.append(i['sampleid'])

        field = request.GET["field"]
        pType = request.GET["pType"]

        your_fields = Sample._meta.local_fields
        mimark = [f.name for f in your_fields]

        your_fields = UserDefined._meta.local_fields
        user = [f.name for f in your_fields]

        your_fields = Air._meta.local_fields
        air = [f.name for f in your_fields]

        your_fields = Human_Associated._meta.local_fields
        human_associated = [f.name for f in your_fields]

        your_fields = Microbial._meta.local_fields
        microbial = [f.name for f in your_fields]

        your_fields = Soil._meta.local_fields
        soil = [f.name for f in your_fields]

        your_fields = Water._meta.local_fields
        water = [f.name for f in your_fields]

        myNode = []
        if field == 'project_name':
            table_field = 'projectid__' + field
            values = Sample.objects.values_list('projectid__project_name', flat='True').filter(sampleid__in=filtered).exclude(projectid__wip=True).distinct()
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        'isFolder': True,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        myNode2 = {
                            'title': str(item.sample_name) + ' (ID: ' + str(item.sampleid) + '; Reads: ' + str(item.reads) + ')',
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'project',
                            'tooltip': 'Project: ' + item.projectid.project_name + ' (ID: ' + item.projectid.projectid + ')',
                            'hideCheckbox': False,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        elif field in mimark:
            values = Sample.objects.values_list(field, flat='True').filter(sampleid__in=filtered).distinct().order_by(field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        'isFolder': True,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        myNode2 = {
                            'title': str(item.sample_name) + ' (ID: ' + str(item.sampleid) + '; Reads: ' + str(item.reads) + ')',
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'mimark',
                            'tooltip': 'Project: ' + item.projectid.project_name + ' (ID: ' + item.projectid.projectid + ')',
                            'hideCheckbox': False,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        elif field in air and pType == 'air':
            table_field = 'air__' + field
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=filtered).distinct().order_by(table_field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        'isFolder': True,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        myNode2 = {
                            'title': str(item.sample_name) + ' (ID: ' + str(item.sampleid) + '; Reads: ' + str(item.reads) + ')',
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'air',
                            'tooltip': 'Project: ' + item.projectid.project_name + ' (ID: ' + item.projectid.projectid + ')',
                            'hideCheckbox': False,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        elif field in human_associated and pType == 'human associated':
            table_field = 'human_associated__' + field
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=filtered).distinct().order_by(table_field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        'isFolder': True,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        myNode2 = {
                            'title': str(item.sample_name) + ' (ID: ' + str(item.sampleid) + '; Reads: ' + str(item.reads) + ')',
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'human_associated',
                            'tooltip': 'Project: ' + item.projectid.project_name + ' (ID: ' + item.projectid.projectid + ')',
                            'hideCheckbox': False,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        elif field in microbial and pType == 'microbial':
            table_field = 'microbial__' + field
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=filtered).distinct().order_by(table_field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        'isFolder': True,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        myNode2 = {
                            'title': str(item.sample_name) + ' (ID: ' + str(item.sampleid) + '; Reads: ' + str(item.reads) + ')',
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'microbial',
                            'tooltip': 'Project: ' + item.projectid.project_name + ' (ID: ' + item.projectid.projectid + ')',
                            'hideCheckbox': False,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        elif field in soil and pType == 'soil':
            table_field = 'soil__' + field
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=filtered).distinct().order_by(table_field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        'isFolder': True,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        myNode2 = {
                            'title': str(item.sample_name) + ' (ID: ' + str(item.sampleid) + '; Reads: ' + str(item.reads) + ')',
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'soil',
                            'tooltip': 'Project: ' + item.projectid.project_name + ' (ID: ' + item.projectid.projectid + ')',
                            'hideCheckbox': False,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        elif field in water and pType == 'water':
            table_field = 'water__' + field
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=filtered).distinct().order_by(table_field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        'isFolder': True,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        myNode2 = {
                            'title': str(item.sample_name) + ' (ID: ' + str(item.sampleid) + '; Reads: ' + str(item.reads) + ')',
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'water',
                            'tooltip': 'Project: ' + item.projectid.project_name + ' (ID: ' + item.projectid.projectid + ')',
                            'hideCheckbox': False,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        elif field in user:
            table_field = 'userdefined__' + field
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=filtered).distinct().order_by(table_field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        'isFolder': True,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        myNode2 = {
                            'title': str(item.sample_name) + ' (ID: ' + str(item.sampleid) + '; Reads: ' + str(item.reads) + ')',
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'user',
                            'tooltip': 'Project: ' + item.projectid.project_name + ' (ID: ' + item.projectid.projectid + ')',
                            'hideCheckbox': False,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        res = json.dumps(myNode)
        return HttpResponse(res, content_type='application/json')


def getSampleQuantTree(request):    # get variable names for quantitative data, grouped by type
    myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
    path = str(myDir) + 'usr_norm_samples.pkl'
    with open(path, 'rb') as f:
        samples = pickle.load(f)

    projectList = Sample.objects.filter(sampleid__in=samples).values_list('projectid').distinct()
    typeList = Project.objects.filter(projectid__in=projectList).values_list('projectType', flat=True)

    myTree = {'title': 'Meta Data: Quantitative', 'id': 'root', 'isFolder': False,  'hideCheckbox': True, 'expand': True, 'children': []}
    mimark = {'title': 'MIMARKs', 'id': 'mimark', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    air = {'title': 'Air', 'id': 'air', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    human_associated = {'title': 'Human Associated', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    microbial = {'title': 'Microbial', 'id': 'microbial', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    soil = {'title': 'Soil', 'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    water = {'title': 'Water', 'id': 'water', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    user = {'title': 'User-defined', 'id': 'user', 'isFolder': True,  'hideCheckbox': True, 'children': []}

    list = ['latitude', 'longitude', 'annual_season_temp', 'annual_season_precpt']
    for i in range(len(list)):
        myNode = {'title': list[i], 'id': 'mimark', 'isFolder': True, 'pType': 'mimark', 'isLazy': True, 'children': []}
        mimark['children'].append(myNode)

    if 'air' in typeList:
        list = ['barometric_press', 'carb_dioxide', 'carb_monoxide', 'elev', 'humidity', 'methane', 'organism_count', 'oxygen', 'pollutants_concentration', 'resp_part_matter_concentration', 'samp_salinity', 'samp_size', 'samp_store_dur', 'samp_store_temp', 'solar_irradiance', 'temp', 'ventilation_rate', 'volatile_org_comp_concentration', 'wind_speed']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'air', 'isFolder': True, 'pType': 'air', 'isLazy': True, 'children': []}
            air['children'].append(myNode)

    if 'microbial' in typeList:
        list = ['alkalinity', 'alkyl_diethers', 'altitude', 'aminopept_act', 'ammonium', 'bacteria_carb_prod', 'biomass_amount', 'bishomohopanol', 'bromide', 'calcium', 'carb_nitro_ratio', 'chloride', 'chlorophyll', 'diether_lipids_concentration', 'diss_carb_dioxide', 'diss_hydrogen', 'diss_inorg_carb', 'diss_org_carb', 'diss_org_nitro', 'diss_oxygen', 'glucosidase_act', 'magnesium', 'mean_frict_vel', 'mean_peak_frict_vel', 'methane', 'n_alkanes_concentration', 'nitrate', 'nitrite', 'nitro', 'org_carb', 'org_matter', 'org_nitro', 'organism_count', 'part_org_carb', 'petroleum_hydrocarb', 'ph', 'phaeopigments_concentration', 'phosphate', 'phosplipid_fatt_acid_concentration', 'potassium', 'pressure', 'redox_potential', 'salinity', 'samp_size', 'samp_store_dur', 'samp_store_temp', 'silicate', 'sodium', 'sulfate', 'sulfide', 'temp', 'tot_carb', 'tot_nitro', 'tot_org_carb', 'turbidity', 'water_content']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'microbial', 'isFolder': True, 'pType': 'microbial', 'isLazy': True, 'children': []}
            microbial['children'].append(myNode)

    if 'water' in typeList:
        list = ['alkalinity', 'alkyl_diethers', 'altitude', 'aminopept_act', 'ammonium', 'atmospheric_data', 'bac_prod', 'bac_resp', 'bacteria_carb_prod', 'biomass_amount', 'bishomohopanol', 'bromide', 'calcium', 'carb_nitro_ratio', 'chloride', 'chlorophyll', 'conduc', 'density', 'diether_lipids', 'diss_carb_dioxide', 'diss_hydrogen', 'diss_inorg_carb', 'diss_inorg_nitro', 'diss_inorg_phosp', 'diss_org_carb', 'diss_org_nitro', 'diss_oxygen', 'down_par', 'elev', 'fluor', 'glucosidase_act', 'light_intensity', 'magnesium', 'mean_frict_vel', 'mean_peak_frict_vel', 'n_alkanes', 'nitrate', 'nitrite', 'nitro', 'org_carb', 'org_matter', 'org_nitro', 'organism_count', 'part_org_carb', 'part_org_nitro', 'pretroleum_hydrocarb', 'ph', 'phaeopigments', 'phosphate', 'phosplipid_fatt_acid', 'photon_flux', 'potassium', 'pressure', 'primary_prod', 'redox_potential', 'samp_salinity', 'samp_size', 'samp_store_dur', 'samp_store_temp', 'samp_vol_we_dna_ext', 'silicate', 'sodium', 'soluble_react_phosp', 'sulfate', 'sulfide', 'suspen_part_matter', 'temp', 'tot_depth_water_col', 'tot_diss_nitro', 'tot_inorg_nitro', 'tot_nitro', 'tot_part_carb', 'tot_phosp', 'water_current_magnitude']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'water', 'isFolder': True, 'pType': 'water', 'isLazy': True, 'children': []}
            water['children'].append(myNode)

    if 'human associated' in typeList:
        samp_collect = {'title': 'Sample Collection', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['samp_size', 'samp_store_temp', 'samp_store_dur']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human associated', 'isLazy': True, 'children': []}
            samp_collect['children'].append(myNode)
        human_associated['children'].append(samp_collect)

        samp_class = {'title': 'Sample Classification', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['samp_temp', 'samp_ph', 'samp_salinity']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human associated', 'isLazy': True, 'children': []}
            samp_class['children'].append(myNode)
        human_associated['children'].append(samp_class)

        host = {'title': 'Host', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['host_age', 'host_pulse', 'host_height', 'host_weight', 'host_bmi', 'host_weight_loss_3_month', 'host_body_temp']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human associated', 'isLazy': True, 'children': []}
            host['children'].append(myNode)
        human_associated['children'].append(host)

        diet = {'title': 'Diet', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['diet_duration']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human associated', 'isLazy': True, 'children': []}
            diet['children'].append(myNode)
        human_associated['children'].append(diet)

        disease = {'title': 'Disease', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['disease_duration', 'organism_count', 'tumor_mass']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human associated', 'isLazy': True, 'children': []}
            disease['children'].append(myNode)
        human_associated['children'].append(disease)

        drug_use = {'title': 'Drug Usage', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['drug_duration']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human associated', 'isLazy': True, 'children': []}
            drug_use['children'].append(myNode)
        human_associated['children'].append(drug_use)

        interven = {'title': 'Intervention', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['pert_duration']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human associated', 'isLazy': True, 'children': []}
            interven['children'].append(myNode)
        human_associated['children'].append(interven)

    if 'soil' in typeList:
        samp_collect = {'title': 'Sample Collection', 'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['samp_size', 'samp_sieve_size', 'samp_store_dur', 'samp_store_temp', 'samp_weight_dna_ext', 'pool_dna_extracts']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            samp_collect['children'].append(myNode)
        soil['children'].append(samp_collect)

        soil_class = {'title': 'Soil Classification',  'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['porosity', 'slope_gradient', 'bulk_density', 'water_content_soil']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            soil_class['children'].append(myNode)
        soil['children'].append(soil_class)

        fertilizer = {'title': 'Fertilizer Mgt',  'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['fert_tot_amount', 'fert_N_tot_amount', 'fert_P_tot_amount', 'fert_K_tot_amount']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            fertilizer['children'].append(myNode)
        soil['children'].append(fertilizer)

        irrigation = {'title': 'Irrigation',  'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['irrigation_tot_amount']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            irrigation['children'].append(myNode)
        soil['children'].append(irrigation)

        residue = {'title': 'Residue Mgt',  'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['residue_removal_percent']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            residue['children'].append(myNode)
        soil['children'].append(residue)

        tillage = {'title': 'Tillage Mgt',  'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['tillage_event_depth']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            tillage['children'].append(myNode)
        soil['children'].append(tillage)

        amend = {'title': 'Amendment Mgt',  'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['amend1_tot_amount', 'amend2_tot_amount', 'amend3_tot_amount']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            amend['children'].append(myNode)
        soil['children'].append(amend)

        microbe = {'title': 'Microbial Biomass',  'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['rRNA_copies', 'microbial_biomass_C', 'microbial_biomass_N', 'microbial_respiration']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            microbe['children'].append(myNode)
        soil['children'].append(microbe)

        soil_nutrient = {'title': 'Soil Nutrients',  'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['soil_pH', 'soil_EC', 'soil_C', 'soil_OM', 'soil_N', 'soil_NO3_N', 'soil_NH4_N', 'soil_P', 'soil_K', 'soil_S', 'soil_Zn', 'soil_Fe', 'soil_Cu', 'soil_Mn', 'soil_Ca', 'soil_Mg', 'soil_Na', 'soil_B']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            soil_nutrient['children'].append(myNode)
        soil['children'].append(soil_nutrient)

        plant_nutrient = {'title': 'Plant Nutrients',  'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['plant_C', 'plant_N', 'plant_P', 'plant_K', 'plant_Ca', 'plant_Mg', 'plant_S', 'plant_Na', 'plant_Cl', 'plant_Al', 'plant_B', 'plant_Cu', 'plant_Fe', 'plant_Mn', 'plant_Zn']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            plant_nutrient['children'].append(myNode)
        soil['children'].append(plant_nutrient)

        biomass = {'title': 'Crop Biomass',  'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['crop_tot_biomass_fw', 'crop_tot_biomass_dw', 'crop_tot_above_biomass_fw', 'crop_tot_above_biomass_dw', 'crop_tot_below_biomass_fw', 'crop_tot_below_biomass_dw', 'harv_fresh_weight', 'harv_dry_weight']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            biomass['children'].append(myNode)
        soil['children'].append(biomass)

        ghg = {'title': 'GHG Flux',  'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['ghg_N2O', 'ghg_CO2', 'ghg_NH4']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            ghg['children'].append(myNode)
        soil['children'].append(ghg)

        health = {'title': 'Soil Health', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['soil_water_cap', 'soil_surf_hard', 'soil_subsurf_hard', 'soil_agg_stability', 'soil_ACE_protein', 'soil_active_C']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'soil', 'isFolder': True, 'pType': 'soil', 'isLazy': True, 'children': []}
            health['children'].append(myNode)
        soil['children'].append(health)

    list = ['usr_quant1', 'usr_quant2', 'usr_quant3', 'usr_quant4', 'usr_quant5', 'usr_quant6']
    for i in range(len(list)):
        myNode = {'title': list[i], 'id': 'user', 'isFolder': True, 'pType': 'user', 'isLazy': True, 'children': []}
        user['children'].append(myNode)

    myTree['children'].append(mimark)
    if 'air' in typeList:
        myTree['children'].append(air)
    if 'microbial' in typeList:
        myTree['children'].append(microbial)
    if 'water' in typeList:
        myTree['children'].append(water)
    if 'human associated' in typeList:
        myTree['children'].append(human_associated)
    if 'soil' in typeList:
        myTree['children'].append(soil)
    myTree['children'].append(user)

    # daymetData, should only be visible if user has a DaymetData object (atm keyed on username)
    try:
        # myDaymet line is to check if object exists
        if DaymetData.objects.filter(user=request.user).exists():
            daymet = {'title': 'Daymet Data', 'id': 'user', 'isFolder': True, 'hideCheckbox': True, 'children': []}
            # "year" "yday" "dayl" "prcp" "srad"  "swe"  "tmax"  "tmin"  "vp"
            list = ['year', 'yday', 'dayl', 'prcp', 'srad', 'swe', 'tmax', 'tmin', 'vp']
            for i in range(len(list)):
                myNode = {'title': list[i], 'id': 'user', 'isFolder': True, 'pType': 'user', 'isLazy': True, 'children': []}
                daymet['children'].append(myNode)
            myTree['children'].append(daymet)
    except Exception as ofcoursetheresanexceptionhere:
        print "Hit the exception on daymetdata:", ofcoursetheresanexceptionhere
        pass

    # Convert result list to a JSON string
    res = json.dumps(myTree)

    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')


def getSampleQuantTreeChildren(request):    # get actual values to populate parent tree with (for lazy loaded children)
    if request.is_ajax():
        myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
        path = str(myDir) + 'usr_norm_samples.pkl'
        with open(path, 'rb') as f:
            selected = pickle.load(f)

        filtered = []
        reads = Sample.objects.filter(sampleid__in=selected).order_by('sample_name').values('sampleid', 'sample_name', 'reads')
        for i in reads:
            if int(i['reads']) > 0:
                filtered.append(i['sampleid'])

        field = request.GET["field"]
        pType = request.GET["pType"]

        your_fields = Sample._meta.local_fields
        mimark = [f.name for f in your_fields]

        your_fields = UserDefined._meta.local_fields
        user = [f.name for f in your_fields]

        your_fields = Air._meta.local_fields
        air = [f.name for f in your_fields]

        your_fields = Human_Associated._meta.local_fields
        human_associated = [f.name for f in your_fields]

        your_fields = Microbial._meta.local_fields
        microbial = [f.name for f in your_fields]

        your_fields = Soil._meta.local_fields
        soil = [f.name for f in your_fields]

        your_fields = Water._meta.local_fields
        water = [f.name for f in your_fields]

        your_fields = DaymetData._meta.local_fields
        daymet = [f.name for f in your_fields]

        myNode = []
        if field in mimark:
            values = Sample.objects.values_list(field, flat='True').filter(sampleid__in=filtered).distinct().order_by(field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        'isFolder': True,
                        'hideCheckbox': False,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        myNode2 = {
                            'title': str(item.sample_name) + ' (ID: ' + str(item.sampleid) + '; Reads: ' + str(item.reads) + ')',
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'mimark',
                            'tooltip': 'Project: ' + item.projectid.project_name + ' (ID: ' + item.projectid.projectid + ')',
                            'hideCheckbox': False,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        elif field in air and pType == 'air':
            table_field = 'air__' + field
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=filtered).distinct().order_by(table_field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        'isFolder': True,
                        'hideCheckbox': False,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        myNode2 = {
                            'title': str(item.sample_name) + ' (ID: ' + str(item.sampleid) + '; Reads: ' + str(item.reads) + ')',
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'air',
                            'tooltip': 'Project: ' + item.projectid.project_name + ' (ID: ' + item.projectid.projectid + ')',
                            'hideCheckbox': False,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        elif field in human_associated and pType == 'human associated':
            table_field = 'human_associated__' + field
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=filtered).distinct().order_by(table_field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        'isFolder': True,
                        'hideCheckbox': False,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        myNode2 = {
                            'title': str(item.sample_name) + ' (ID: ' + str(item.sampleid) + '; Reads: ' + str(item.reads) + ')',
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'human_associated',
                            'tooltip': 'Project: ' + item.projectid.project_name + ' (ID: ' + item.projectid.projectid + ')',
                            'hideCheckbox': False,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        elif field in microbial and pType == 'microbial':
            table_field = 'microbial__' + field
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=filtered).distinct().order_by(table_field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        'isFolder': True,
                        'hideCheckbox': False,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        myNode2 = {
                            'title': str(item.sample_name) + ' (ID: ' + str(item.sampleid) + '; Reads: ' + str(item.reads) + ')',
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'microbial',
                            'tooltip': 'Project: ' + item.projectid.project_name + ' (ID: ' + item.projectid.projectid + ')',
                            'hideCheckbox': False,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        elif field in soil and pType == 'soil':
            table_field = 'soil__' + field
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=filtered).distinct().order_by(table_field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        'isFolder': True,
                        'hideCheckbox': False,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        myNode2 = {
                            'title': str(item.sample_name) + ' (ID: ' + str(item.sampleid) + '; Reads: ' + str(item.reads) + ')',
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'soil',
                            'tooltip': 'Project: ' + item.projectid.project_name + ' (ID: ' + item.projectid.projectid + ')',
                            'hideCheckbox': False,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        elif field in water and pType == 'water':
            table_field = 'water__' + field
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=filtered).distinct().order_by(table_field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        'isFolder': True,
                        'hideCheckbox': False,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        myNode2 = {
                            'title': str(item.sample_name) + ' (ID: ' + str(item.sampleid) + '; Reads: ' + str(item.reads) + ')',
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'water',
                            'tooltip': 'Project: ' + item.projectid.project_name + ' (ID: ' + item.projectid.projectid + ')',
                            'hideCheckbox': False,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        elif field in user:
            table_field = 'userdefined__' + field
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=filtered).distinct().order_by(table_field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        'isFolder': True,
                        'hideCheckbox': False,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        myNode2 = {
                            'title': str(item.sample_name) + ' (ID: ' + str(item.sampleid) + '; Reads: ' + str(item.reads) + ')',
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'user',
                            'tooltip': 'Project: ' + item.projectid.project_name + ' (ID: ' + item.projectid.projectid + ')',
                            'hideCheckbox': False,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        elif field in daymet:
            # the trick with daymet tree children is that a single daymet object contains all info needed, buried in strings
            # need to get the correct strings (sampleid and whatever node is being expanded upon) and pair values together
            # in theory just need to query object by user, then run through correct strings, format for tree, done
            #values = DaymetData.objects.values_list(field, flat='True').filter(sampleid__in=filtered).distinct().order_by(field)
            # populate values list from paired strings
            myDaymetData = DaymetData.objects.get(user=request.user)
            daySamps = myDaymetData.sampleIDs.split(";")[:-1]
            # "year" "yday" "dayl" "prcp" "srad"  "swe"  "tmax"  "tmin"  "vp"
            # for clarity and security, values list is selected based on a series of if/elifs
            # check which daymet data is being queried for, pull up corresponding list
            # then iterate through in order to make sampID key with value as queried val, use dictionary to make nodes
            valueString = ""
            #print "daymet tree children, looking for:", field
            if field == "year":
                valueString = myDaymetData.year
            elif field == "yday":
                valueString = myDaymetData.yday
            elif field == "dayl":
                valueString = myDaymetData.dayl
                field = "dayl..s."
            elif field == "prcp":
                valueString = myDaymetData.prcp
                field = "prcp..mm.day."
            elif field == "srad":
                valueString = myDaymetData.srad
                field = "srad..W.m.2."
            elif field == "swe":
                valueString = myDaymetData.swe
                field = "swe..kg.m.2."
            elif field == "tmax":
                valueString = myDaymetData.tmax
                field = "tmax..deg.c."
            elif field == "tmin":
                valueString = myDaymetData.tmin
                field = "tmin..deg.c."
            elif field == "vp":
                valueString = myDaymetData.vp
                field = "vp..Pa."
            else:
                print "User requested a non existent daymet field!!"

            # "year" "yday" "dayl..s." "prcp..mm.day." "srad..W.m.2."  "swe..kg.m.2."  "tmax..deg.c."  "tmin..deg.c."  "vp..Pa."
            values = valueString.split(';')[:-1]
            #print "Daymet for tree:", len(daySamps), ":", len(values)
            #print "Values:", values
            # group identical values (and their corresponding sampleids)
            valDict = {}
            sampIter = 0
            for val in values:
                if val not in valDict.keys():
                    valDict[val] = []

                valDict[val].append(Sample.objects.get(sampleid=daySamps[sampIter]))
                sampIter += 1
            # now just sync iterating through daySamps and values
            for keyVal in valDict.keys():
                if pd.notnull(val) and not val == 'nan':
                    myNode1 = {
                        'title': keyVal,
                        'id': field,
                        'isFolder': True,
                        'hideCheckbox': False,
                        'children': []
                    }
                    # folder is value, put all sampleids paired with this value in tree (no duplicates)

                    #args_list = []
                    #args_list.append(Q(**{field: val}))
                    #items = DaymetData.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    # item object needs samplename (from id) and sampleid (easy) plus projectid and projectname (from sample)
                    # value is still from list
                    for curSamp in valDict[keyVal]:
                        myNode2 = {
                            'title': str(curSamp.sample_name) + ' (ID: ' + str(curSamp.sampleid) + '; Reads: ' + str(
                                curSamp.reads) + ')',
                            'id': curSamp.sampleid,
                            'field': field,
                            'value': val,
                            'table': 'mimark',
                            'tooltip': 'Project: ' + curSamp.projectid.project_name + ' (ID: ' + curSamp.projectid.projectid + ')',
                            'hideCheckbox': False,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        res = json.dumps(myNode)
        return HttpResponse(res, content_type='application/json')


def getTaxaTree(request):   # Get main taxanomic datapoints (just kingdoms as folders really)
    myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
    path = str(myDir) + 'usr_norm_samples.pkl'
    with open(path, 'rb') as f:
        selected = pickle.load(f)

    filtered = []
    reads = Sample.objects.filter(sampleid__in=selected).values('sampleid', 'sample_name', 'reads')
    for i in reads:
        if int(i['reads']) > 0:
            filtered.append(i['sampleid'])

    selected_taxa = Profile.objects.filter(sampleid__in=filtered).values_list('kingdomid', flat=True)

    myTree = {'title': 'Taxa Name', 'tooltip': 'root', 'isFolder': False, 'hideCheckbox': True, 'expand': True, 'children': []}
    kingdoms = Kingdom.objects.filter(kingdomid__in=selected_taxa).order_by('kingdomName')
    for kingdom in kingdoms:
        myNode = {
            'title': kingdom.kingdomName,
            'id': kingdom.kingdomid,
            'tooltip': "Kingdom",
            'isFolder': True,
            'expand': False,
            'isLazy': True,
            'children': []
        }
        myTree['children'].append(myNode)

    # Convert result list to a JSON string
    res = json.dumps(myTree)

    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')


def getTaxaTreeChildren(request):   # get deeper layers of the taxa tree upon request (based on selected parent node)
    if request.is_ajax():
        myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
        path = str(myDir) + 'usr_norm_samples.pkl'
        with open(path, 'rb') as f:
            selected = pickle.load(f)

        filtered = []
        reads = Sample.objects.filter(sampleid__in=selected).values('sampleid', 'sample_name', 'reads')
        for i in reads:
            if i['reads'] > 0:
                filtered.append(i['sampleid'])

        selected_taxa = Profile.objects.filter(sampleid__in=filtered)

        taxa = request.GET["tooltip"]
        id = request.GET["id"]

        nodes = []
        if taxa == 'Kingdom':
            qs = Phyla.objects.filter(phylaid__in=selected_taxa.values_list('phylaid')).filter(**{'kingdomid': id}).order_by('phylaName')
            for item in qs:
                myNode = {
                    'title': item.phylaName,
                    'id': item.phylaid,
                    'tooltip': "Phyla",
                    'isFolder': True,
                    'isLazy': True
                }
                nodes.append(myNode)

        if taxa == 'Phyla':
            qs = Class.objects.filter(classid__in=selected_taxa.values_list('classid')).filter(**{'phylaid': id}).order_by('className')
            for item in qs:
                myNode = {
                    'title': item.className,
                    'id': item.classid,
                    'tooltip': "Class",
                    'isFolder': True,
                    'isLazy': True
                }
                nodes.append(myNode)

        elif taxa == 'Class':
            qs = Order.objects.filter(orderid__in=selected_taxa.values_list('orderid')).filter(**{'classid': id}).order_by('orderName')
            for item in qs:
                myNode = {
                    'title': item.orderName,
                    'id': item.orderid,
                    'tooltip': "Order",
                    'isFolder': True,
                    'isLazy': True
                }
                nodes.append(myNode)

        elif taxa == 'Order':
            qs = Family.objects.filter(familyid__in=selected_taxa.values_list('familyid')).filter(**{'orderid': id}).order_by('familyName')
            for item in qs:
                myNode = {
                    'title': item.familyName,
                    'id': item.familyid,
                    'tooltip': "Family",
                    'isFolder': True,
                    'isLazy': True
                }
                nodes.append(myNode)

        elif taxa == 'Family':
            qs = Genus.objects.filter(genusid__in=selected_taxa.values_list('genusid')).filter(**{'familyid': id}).order_by('genusName')
            for item in qs:
                myNode = {
                    'title': item.genusName,
                    'id': item.genusid,
                    'tooltip': "Genus",
                    'isFolder': True,
                    'isLazy': True
                }
                nodes.append(myNode)

        elif taxa == 'Genus':
            qs = Species.objects.filter(speciesid__in=selected_taxa.values_list('speciesid')).filter(**{'genusid': id}).order_by('speciesName')
            for item in qs:
                myNode = {
                    'title': item.speciesName,
                    'id': item.speciesid,
                    'tooltip': "Species",
                    'isFolder': True,
                    'isLazy': True
                }
                nodes.append(myNode)

        elif taxa == 'Species':
            qs = OTU_99.objects.filter(otuid__in=selected_taxa.values_list('otuid')).filter(**{'speciesid': id}).order_by('otuName')
            for item in qs:
                myNode = {
                    'title': item.otuName,
                    'id': item.otuid,
                    'tooltip': "OTU_99",
                    'isLazy': False
                }
                nodes.append(myNode)

        elif taxa == 'OTU_99':
            pass

        res = json.dumps(nodes)
        return HttpResponse(res, content_type='application/json')


def getKEGGTree(request):   # Create the first three layers of the kegg orthology tree (the static ones)
    myTree = {'title': 'KEGG Pathway', 'isFolder': False, 'hideCheckbox': True, 'expand': True, 'children': []}

    level1 = ko_lvl1.objects.using('picrust').all().order_by('ko_lvl1_name').distinct()
    for item in level1:
        myNode = {
            'title': item.ko_lvl1_name,
            'id': item.ko_lvl1_id,
            'tooltip': "Level1",
            'kegg': 'Level1',
            'isFolder': True,
            'expand': False,
            'children': []
        }
        level2 = item.ko_lvl2_set.all().exclude(ko_lvl2_name='Enzyme families').order_by('ko_lvl2_name').distinct()
        for thing1 in level2:
            myNode1 = {
                'title': thing1.ko_lvl2_name,
                'id': thing1.ko_lvl2_id,
                'tooltip': "Level2",
                'kegg': 'Level2',
                'isFolder': True,
                'expand': False,
                'children': []
            }
            level3 = thing1.ko_lvl3_set.all().filter(ko_lvl3_name__contains='[PATH:ko').order_by('ko_lvl3_name').distinct()
            for thing2 in level3:
                myNode2 = {
                    'title': thing2.ko_lvl3_name,
                    'id': thing2.ko_lvl3_id,
                    'tooltip': "Level3",
                    'kegg': 'Level3',
                    'isFolder': True,
                    'expand': False,
                    'isLazy': True
                }
                myNode1['children'].append(myNode2)
            myNode['children'].append(myNode1)
        myTree['children'].append(myNode)

    # Convert result list to a JSON string
    res = json.dumps(myTree)

    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')


def getKEGGTreeChildren(request):   # get lazy-loaded data for parent node in kegg orthology tree !UNUSED! see getKegg2
    if request.is_ajax():
        level = request.GET["tooltip"]
        id = request.GET["id"]

        nodes = []
        if level == 'Level3':
            qs = ko_entry.objects.using('picrust').filter(**{'ko_lvl3_id': id}).order_by('ko_name')
            for item in qs:
                myNode = {
                    'title': str(item.ko_orthology) + ": " + str(item.ko_name),
                    'id': item.ko_lvl4_id,
                    'tooltip': item.ko_desc,
                    'kegg': 'Level4',
                    'isFolder': False,
                    'isLazy': True
                }
                nodes.append(myNode)

        elif level == 'Level4':
            pass

        res = json.dumps(nodes)
        return HttpResponse(res, content_type='application/json')


def getKEGGTree2(request):  # get lazy-loaded data for parent node in kegg orthology tree
    myTree = {'title': 'KEGG Pathway', 'isFolder': False, 'hideCheckbox': True, 'expand': True, 'children': []}

    if os.name == 'nt':
        r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
    else:
        r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

    r("load('myPhyloDB/media/kegg/kegg.sets.ko.RData')")
    pathArray = r.get("names(kegg.sets.ko)")

    pathList = []
    for i in pathArray:
        pathList.append(i.split()[0])

    level1 = ko_lvl1.objects.using('picrust').all().order_by('ko_lvl1_name').distinct()
    for item in level1:
        myNode = {
            'title': item.ko_lvl1_name,
            'id': item.ko_lvl1_id,
            'hideCheckbox': True,
            'tooltip': "Level1",
            'kegg': 'Level1',
            'isFolder': True,
            'expand': False,
            'children': []
        }
        level2 = item.ko_lvl2_set.all().order_by('ko_lvl2_name').distinct()
        for thing1 in level2:
            myNode1 = {
                'title': thing1.ko_lvl2_name,
                'id': thing1.ko_lvl2_id,
                'hideCheckbox': True,
                'tooltip': "Level2",
                'kegg': 'Level2',
                'isFolder': True,
                'expand': False,
                'children': []
            }
            level3 = thing1.ko_lvl3_set.all().filter(ko_lvl3_name__contains='[PATH:ko').order_by('ko_lvl3_name').distinct()
            for thing2 in level3:
                fullName = thing2.ko_lvl3_name
                path = fullName.split('[PATH:')[1].split(']')[0]
                if path in pathList:
                    myNode2 = {
                        'title': thing2.ko_lvl3_name,
                        'id': thing2.ko_lvl3_id,
                        'tooltip': "Level3",
                        'kegg': 'Level3',
                        'isFolder': False,
                        'expand': False
                    }
                    myNode1['children'].append(myNode2)
            myNode['children'].append(myNode1)
        myTree['children'].append(myNode)

    # Convert result list to a JSON string
    res = json.dumps(myTree)

    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')


def getNZTree(request):  # populate first three layers of kegg enzyme tree with lazy-loaded folders
    myTree = {'title': 'KEGG Enzyme', 'isFolder': False, 'hideCheckbox': True, 'expand': True, 'children': []}

    level1 = nz_lvl1.objects.using('picrust').all().order_by('nz_lvl1_name').distinct()

    for item in level1:
        myNode = {
            'title': item.nz_lvl1_name,
            'id': item.nz_lvl1_id,
            'tooltip': "Level1",
            'kegg': 'Level1',
            'isFolder': True,
            'expand': False,
            'children': []
        }
        level2 = item.nz_lvl2_set.all().order_by('nz_lvl2_name').distinct()
        for thing1 in level2:
            myNode1 = {
                'title': thing1.nz_lvl2_name,
                'id': thing1.nz_lvl2_id,
                'tooltip': "Level2",
                'kegg': 'Level2',
                'isFolder': True,
                'expand': False,
                'children': []
            }
            level3 = thing1.nz_lvl3_set.all().order_by('nz_lvl3_name').distinct()
            for thing2 in level3:
                myNode2 = {
                    'title': thing2.nz_lvl3_name,
                    'id': thing2.nz_lvl3_id,
                    'tooltip': "Level3",
                    'kegg': 'Level3',
                    'isFolder': True,
                    'expand': False,
                    'children': []
                }
                level4 = thing2.nz_lvl4_set.all().order_by('nz_lvl4_name').distinct()
                for thing3 in level4:
                    myNode3 = {
                        'title': thing3.nz_lvl4_name,
                        'id': thing3.nz_lvl4_id,
                        'tooltip': "Level4",
                        'kegg': 'Level4',
                        'isFolder': True,
                        'expand': False,
                        'isLazy': True
                    }
                    myNode2['children'].append(myNode3)
                myNode1['children'].append(myNode2)
            myNode['children'].append(myNode1)
        myTree['children'].append(myNode)

    # Convert result list to a JSON string
    res = json.dumps(myTree)

    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')


def getNZTreeChildren(request):  # get child nodes from lazy-loaded kegg enzyme tree
    if request.is_ajax():
        level = request.GET["tooltip"]
        id = request.GET["id"]

        nodes = []
        if level == 'Level4':
            qs = nz_entry.objects.using('picrust').filter(**{'nz_lvl4_id': id}).order_by('nz_name')
            for item in qs:
                myNode = {
                    'title': str(item.nz_orthology) + ": " + str(item.nz_name),
                    'id': item.nz_lvl5_id,
                    'tooltip': item.nz_desc,
                    'kegg': 'Level5',
                    'isFolder': False,
                    'isLazy': True
                }
                nodes.append(myNode)

        elif level == 'Level5':
            pass

        res = json.dumps(nodes)
        return HttpResponse(res, content_type='application/json')


def makeUpdateTree(request):
    myTree = {'title': 'All Projects', 'isFolder': True, 'expand': True, 'hideCheckbox': True, 'children': []}

    projects = perms.getEditProjects(request)
    for project in projects:
        myNode = {
            'title': project.project_name,
            'tooltip': project.project_desc,
            'isFolder': True,
            'hideCheckbox': True,
            'children': [],
            'wip': project.wip
        }

        refids = project.reference_set.all().order_by('path')
        for ref in refids:
            myNode1={
                'title': "Path: " + str(ref.path),
                'tooltip': ref.projectid.project_desc,
                'id': ref.refid,
                'isFolder': True,
                'children': []
            }

            myNode['children'].append(myNode1)

        myTree['children'].append(myNode)

    # Convert result list to a JSON string
    res = json.dumps(myTree)

    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')


def makeFilesTree(request):
    # print "Populating file tree"
    userid = str(request.user.id)
    name = request.GET['name']
    # filter is true when we want specifically files with parent folder 'name', and don't want to display that in output
    filter = request.GET['filter']
    username = request.user.username
    if name == 'Files':
        myName = name
    else:
        myName = str(name) + ' Files'
    myTree = {'title': myName, 'isFolder': True, 'expand': False, 'hideCheckbox': True, 'children': []}
    cwd = os.getcwd()

    userList = []

    if request.user.is_authenticated:
        if request.user.is_superuser:
            # get all usernames
            users = User.objects.all()
            for user in users:
                userList.append(user.username)
        else:
            # unfiltered mode is for file deletion tree, should only be this user unless super? disregard perms?
            # get usernames which this user has permissions for, start by adding this user
            userList.append(username)
            # now we add usernames based on hasPermsFrom attribute  (assuming list is kept up to date properly)
            myPerms = UserProfile.objects.get(user=request.user).hasPermsFrom.split(';')
            for permname in myPerms:
                if permname != "":
                    userList.append(permname)

        if filter == "true":
            if name.lower() == "script":
                userList = []
                userList.append('admin')

        for user in userList:
            folderPath = cwd + "/user_uploads/" + user
            myTree['children'].append(
                fileTreeChildren(folderPath, user, filter, name))
    else:
        pass



    # given a userid folder, should be able to populate a node with said folder and its subdirs, with actual files there
    # print "MyTree:", myTree

    # Convert result list to a JSON string
    res = json.dumps(myTree)

    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')


def fileTreeChildren(path, name, filter, filterName):

    '''

    :param path: location to use as root
    :param name: display name of folder (and id)
    :return: a singular node, which contains its own children and grandchildren
    '''

    myNode = {
        'title': name,
        'id': name,
        'isFolder': True,
        'children': [],
        'hideCheckbox': True
    }
    try:
        for root, dirs, files in os.walk(path):
            for curdir in dirs:
                shouldAdd = False
                dirNode = {
                    'title': curdir,
                    'id': curdir,
                    'isFolder': True,
                    'children': []
                }
                if filter == "false":   # javascript false vs python False, vs 0 vs '0'
                    nextpath = os.path.join(path, curdir)
                    for temproot, tempdirs, tempfiles in os.walk(nextpath):
                        for tempdir in tempdirs:
                            # add folder node here
                            otherDirNode = {
                                'title': tempdir,
                                'id': tempdir,
                                'isFolder': True,
                                'children': []
                            }
                            for moreroot, moredirs, morefiles in os.walk(os.path.join(nextpath, tempdir)):
                                for morefile in morefiles:
                                    #print "Adding file", morefile
                                    myNode2 = {
                                        'title': morefile,
                                        'id': morefile,
                                        'isFolder': False
                                    }
                                    otherDirNode['children'].append(myNode2)
                                    shouldAdd = True

                            dirNode['children'].append(otherDirNode)
                else:
                    # if sff or oligos don't hide checkbox?, filter out ones with too few levels on python end
                    # python end filtering because control, JS can get messed with easily

                    if filterName.lower() not in "sff oligos fna qual":  # add all multis here
                        dirNode['hideCheckbox'] = True

                    for temproot, tempdirs, tempfiles in os.walk(os.path.join(path, curdir+"/"+filterName.lower())):
                        for curfile in tempfiles:
                            myNode2 = {
                                'title': curfile,
                                'id': curfile,
                                'isFolder': False,
                            }
                            dirNode['children'].append(myNode2)
                            shouldAdd = True
                if shouldAdd:
                    myNode['children'].append(dirNode)

    except Exception as exc:
        print "Error making file tree:", exc
    return myNode


def makeReproTree(request):
    myTree = {'title': 'All Uploads', 'isFolder': True, 'expand': True, 'hideCheckbox': True, 'children': []}

    projects = perms.getEditProjects(request)
    for project in projects:
        myNode = {
            'title': "Project: " + str(project.project_name),
            'tooltip': project.project_desc,
            'isFolder': True,
            'children': [],
            'wip': project.wip
        }

        refids = project.reference_set.all().order_by('path')
        for ref in refids:
            myNode1={
                'title': "Path: " + str(ref.path),
                'tooltip': ref.projectid.project_desc,
                'id': ref.refid,
                'isFolder': True,
                'children': []
            }
            myNode['children'].append(myNode1)

        myTree['children'].append(myNode)

    # Convert result list to a JSON string
    res = json.dumps(myTree)

    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')


def getDownloadTree(request):   # get list of available downloads for this user (project files they have perms for? TODO PERMS CHECK HERE
    myTree = {'title': 'All Projects', 'isFolder': True, 'expand': True, 'hideCheckbox': True, 'children': []}

    projects = Project.objects.none()
    if request.user.is_superuser:
        projects = Project.objects.all().order_by('project_name')
    elif request.user.is_authenticated():
        path_list = Reference.objects.filter(Q(author=request.user)).values_list('projectid_id')
        projects = Project.objects.all().filter( Q(projectid__in=path_list) | Q(status='public') ).order_by('project_name')
    if not request.user.is_superuser and not request.user.is_authenticated():
        projects = Project.objects.all().filter( Q(status='public') ).order_by('project_name')

    for project in projects:
        if Sample.objects.filter(projectid=project.projectid).exists():
            myNode = {
                'title': project.project_name,
                'tooltip': "Project type: " + project.projectType + "\nDescription: " + project.project_desc + "\nID: " + project.projectid + "\nPI: " + project.pi_first + " " + project.pi_last + "\nAffiliation: " + project.pi_affiliation,
                'id': project.projectid,
                'isFolder': True,
                'hideCheckbox': False,
                'isLazy': True,
                'wip': project.wip
            }
            myTree['children'].append(myNode)
    # Convert result list to a JSON string
    res = json.dumps(myTree)

    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')


def getDownloadTreeChildren(request):
    # get project children which are visible to current user (check
    if request.is_ajax():
        projectid = request.GET["id"]
        samples = Reference.objects.filter(projectid=projectid)

        nodes = []
        for sample in samples:
            myNode = {
                'title': sample.path,
                'id': sample.refid,
                'isFolder': False,
                'children': []
            }
            for root, dirs, files in os.walk(sample.path):
                for file in files:
                    myNode2 = {
                        'title': file,
                        'id': file,
                        'isFolder': False
                    }
                    myNode['children'].append(myNode2)
            nodes.append(myNode)

        res = json.dumps(nodes)
        return HttpResponse(res, content_type='application/json')


def getPermissionTree(request):
    myTree = {'title': 'My Projects', 'isFolder': True, 'expand': True, 'hideCheckbox': True, 'children': []}
    publicTree = {'title': 'Public', 'isFolder': True, 'expand': True, 'hideCheckbox': True, 'children': []}
    privateTree = {'title': 'Private', 'isFolder': True, 'expand': True, 'hideCheckbox': True, 'children': []}
    # split tree into private and public projects via status flag

    projects = Project.objects.none()

    if request.user.is_superuser:
        projects = Project.objects.all()
    elif request.user.is_authenticated():
        # get projects this user owns, we want to display them all (plus their whitelists)
        projects = Project.objects.all().filter(owner=request.user)
    # getPermissionTree is for viewing projects owned by this user and permissions which were granted for said projects
    for project in projects:
        viewList = project.whitelist_view.split(";")
        projPerms = []
        for uname in viewList:
            if uname != "":
                myPerm = {
                    'title': uname,
                    'tooltip': "Select to revoke viewing permission for this project",
                    'id': uname+";"+project.projectid,
                    'isFolder': False,
                    'wip': project.wip
                }
                projPerms.append(myPerm)
        myNode = {
            'title': project.project_name,
            'tooltip': "Project type: " + project.projectType + "\nDescription: " + project.project_desc + "\nID: " + project.projectid + "\nPI: " + project.pi_first + " " + project.pi_last + "\nAffiliation: " + project.pi_affiliation,
            'id': project.projectid,
            'isFolder': True,
            'isLazy': False,
            'wip': project.wip,
            'children': projPerms
        }
        if project.status == "public":  # should only be private and public, defaulting to private just in case
            publicTree['children'].append(myNode)
        else:
            privateTree['children'].append(myNode)

    myTree['children'].append(publicTree)
    myTree['children'].append(privateTree)
    # Convert result list to a JSON string
    res = json.dumps(myTree)

    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')


def getFilePermTree(request):   # given == True for gavePermsTo, else hasPermsFrom
    # this tree should be a single level, just a list of usernames (case sensitive)
    # user should be able to select usernames to revoke permissions for
    # security goals should be to not display more information than necessary (uuid, actual case of username?)

    given = request.GET['given']

    if given == "true":
        myTree = {'title': "Users you've given permission to", 'isFolder': True, 'expand': True, 'hideCheckbox': True, 'children': []}
        myPermList = UserProfile.objects.get(user=request.user).gavePermsTo.split(";")
    else:
        myTree = {'title': "", 'isFolder': True, 'expand': True, 'hideCheckbox': True, 'children': []}
        myPermList = UserProfile.objects.get(user=request.user).hasPermsFrom.split(";")

    for username in myPermList:
        if username != "":
            myNode = {
                'title': username,
                'id': username,
                'isFolder': False,
                'isLazy': False
            }
            if given != "true":
                myNode['hideCheckbox'] = True
            myTree['children'].append(myNode)
    # Convert result list to a JSON string
    res = json.dumps(myTree)

    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')


def getLocationSamplesTree(request):
    # request should contain marker data (aka lat lon pair, sampleids, projectnames)
    # use sampleids to build tree of sample names, parent node is projectname (selectable but not sent if selected)
    # tree should link sampleids of itself to ids in main select tree, mirror select/unselect
    # actual select code should be unaffected, just want grouped sample data

    # permissions currently only allow for check via project, so samples must be from projects with perms

    sampString = request.GET['sampleIDs']
    myTree = {'title': 'Samples at Location', 'isFolder': True, 'expand': True, 'hideCheckbox': True, 'children': []}
    sampIDs = sampString.split(":")
    projectDict = {}
    for sampID in sampIDs:
        if sampID != "":
            try:
                thisSample = Sample.objects.get(sampleid=sampID)
                thisProject = thisSample.projectid  # this attribute is misnamed, as its a full project object, not ID
                projID = str(thisProject.projectid)
                if projID not in projectDict:
                    projectDict[projID] = []
                projectDict[projID].append(thisSample)

            except Exception as e:
                print "Error during query for location tree:", e

    myProjects = perms.getViewProjects(request)

    for projID in projectDict:
        thisProject = Project.objects.get(projectid=projID)
        if thisProject in myProjects:   # simple enough security check, since samples are looped by project already
            myProjNode = {
                'title': thisProject.project_name,
                'tooltip': "Project type: " + thisProject.projectType + "\nDescription: " + thisProject.project_desc + "\nID: " + thisProject.projectid + "\nPI: " + thisProject.pi_first + " " + thisProject.pi_last + "\nAffiliation: " + thisProject.pi_affiliation + "\nOwner: " + thisProject.owner.username,
                'id': thisProject.projectid,
                'isFolder': True,
                'wip': thisProject.wip,
                'children': [],
                # 'expand': True # use this for auto expanded projects, currently inactive for compactness sake
            }
            for samp in projectDict[projID]:
                try:
                    mySampNode = {
                        'title': 'Name: ' + str(samp.sample_name) + '; Reads: ' + str(samp.reads),
                        'tooltip': 'ID: ' + str(samp.sampleid),
                        'id': str(samp.sampleid),
                        'isFolder': False
                    }
                    myProjNode['children'].append(mySampNode)
                except Exception as e:
                    print "Error during sample selection for location tree:", e

            myTree['children'].append(myProjNode)


    res = json.dumps(myTree)
    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')


def getFilterSamplesTree(request):
    # request should contain marker data (aka lat lon pair, sampleids, projectnames)
    # use sampleids to build tree of sample names, parent node is projectname (selectable but not sent if selected)
    # tree should link sampleids of itself to ids in main select tree, mirror select/unselect
    # actual select code should be unaffected, just want grouped sample data

    # permissions currently only allow for check via project, so samples must be from projects with perms

    pType = request.GET['ptype']
    fName = request.GET['fname']
    fVal = request.GET['fval']
    myTree = {'title': 'Filtered Data', 'isFolder': True, 'expand': True, 'hideCheckbox': True, 'children': []}
    # print "PTYPE:", pType
    # print "FNAME:", fName
    # print "FVAL:", fVal
    if pType != "" and pType is not None and fName != "" and fName is not None:
        projectSet = []
        myProjects = perms.getViewProjects(request)
        sampleDict = {}
        for proj in myProjects:
            if proj.projectType == pType or pType == "mimarks" or pType == "user_defined":
                projectSet.append(proj)
                sampleDict[proj.projectid] = []
        # print "Found", len(projectSet), "projects of correct type"
        # writing new queries here, should grab all sample objects with 1-2 queries (no loops necessary)
        #newTime = time.time()
        try:

            workSamps = None
            if pType.lower() == "mimarks":
                workSamps = Sample.objects.filter(projectid__in=projectSet)
            if pType.lower() == "soil":
                workSamps = Soil.objects.filter(projectid__in=projectSet)
            if pType.lower() == "human_associated":
                workSamps = Human_Associated.objects.filter(projectid__in=projectSet)
            if pType.lower() == "air":
                workSamps = Air.objects.filter(projectid__in=projectSet)
            if pType.lower() == "water":
                workSamps = Water.objects.filter(projectid__in=projectSet)
            if pType.lower() == "microbial":
                workSamps = Microbial.objects.filter(projectid__in=projectSet)
            if pType.lower() == "user_defined":
                workSamps = UserDefined.objects.filter(projectid__in=projectSet)
            # now check workSamps before checking values
            if workSamps is not None:
                # workSamps actually has data, time to filter by value
                if fName != "nofilterselected":
                    # a filter has actually been set
                    for workSamp in workSamps:
                        workVal = str(getattr(workSamp, fName))
                        thisValPasses = False
                        if fVal == "":
                            if workVal is not None and workVal != "" and workVal != "nan" and workVal != "None":
                                thisValPasses = True
                        elif workVal == fVal:
                            thisValPasses = True
                        if thisValPasses:
                            sampleDict[workSamp.projectid.projectid].append(workSamp)  # projectid still foreign project
                else:
                    #nofilter
                    # just return mySamples
                    for workSamp in workSamps:
                        sampleDict[workSamp.projectid.projectid].append(workSamp)

        except Exception as e:
            print "Error with filter subquerying", e
            traceback.print_exc()

        projectSet = []
        for projid in sampleDict.keys():
            if len(sampleDict[projid]) > 0:
                projectSet.append(Project.objects.get(projectid=projid))

        projectSet.sort(key=lambda x: x.project_name)   # this is actually quite fast

        for proj in projectSet:
            myProjNode = {
                'title': proj.project_name,
                'tooltip': "Project type: " + proj.projectType + "\nDescription: " + proj.project_desc + "\nID: " + proj.projectid + "\nPI: " + proj.pi_first + " " + proj.pi_last + "\nAffiliation: " + proj.pi_affiliation + "\nOwner: " + proj.owner.username,
                'id': proj.projectid,
                'isFolder': True,
                'wip': proj.wip,
                'children': [],
                # 'expand': True # use this for auto expanded projects, currently inactive for compactness sake
            }
            for samp in sampleDict[proj.projectid]:
                try:
                    if fName == "nofilterselected":
                        if pType == "mimarks":
                            mySampNode = {
                                'title': 'Name: ' + str(samp.sample_name) + '; Reads: ' + str(
                                    samp.reads),
                                'tooltip': 'ID: ' + str(samp.sampleid),
                                'id': str(samp.sampleid),
                                'isFolder': False
                            }
                        else:
                            mySampNode = {
                                'title': 'Name: ' + str(
                                    samp.sampleid.sample_name) + '; Reads: ' + str(samp.sampleid.reads),
                                'tooltip': 'ID: ' + str(samp.sampleid.sampleid),
                                'id': str(samp.sampleid.sampleid),
                                'isFolder': False
                            }
                    else:
                        myVal = str(getattr(samp, fName))
                        if pType == "mimarks":
                            mySampNode = {
                                'title': "Value: " + myVal + '; Name: ' + str(samp.sample_name) + '; Reads: ' + str(samp.reads),
                                'tooltip': 'ID: ' + str(samp.sampleid),
                                'id': str(samp.sampleid),
                                'isFolder': False
                            }
                        else:
                            mySampNode = {
                                'title': "Value: " + myVal + '; Name: ' + str(samp.sampleid.sample_name) + '; Reads: ' + str(samp.sampleid.reads),
                                'tooltip': 'ID: ' + str(samp.sampleid.sampleid),
                                'id': str(samp.sampleid.sampleid),
                                'isFolder': False
                            }
                    myProjNode['children'].append(mySampNode)
                except Exception as e:
                    print "Error during sample selection for location tree:", e

            myTree['children'].append(myProjNode)

    res = json.dumps(myTree)
    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')
