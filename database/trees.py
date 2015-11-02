import pickle
import operator
import simplejson
from django.http import HttpResponse
from django.db.models import Q, Sum
from models import Project, Sample, Reference
from models import Kingdom, Class, Order, Family, Genus, Species, Profile
from models import Air, Human_Associated, Microbial, Soil, Water, UserDefined
import numpy as np
import pandas as pd


def getProjectTree(request):
    myTree = {'title': 'All Projects', 'isFolder': True, 'expand': True, 'hideCheckbox': True, 'children': []}

    projects = Project.objects.none()
    if request.user.is_superuser:
        projects = Project.objects.all().order_by('project_name')
    if request.user.is_authenticated():
        path_list = Reference.objects.filter(Q(author=request.user)).values_list('projectid_id')
        projects = Project.objects.all().filter( Q(projectid__in=path_list) | Q(status='public') ).order_by('project_name')
    if not request.user.is_superuser and not request.user.is_authenticated():
        projects = Project.objects.all().filter( Q(status='public') ).order_by('project_name')

    for project in projects:
        myNode = {
            'title': project.project_name,
            'tooltip': project.project_desc,
            'id': project.projectid,
            'isFolder': True,
            'isLazy': True
        }
        myTree['children'].append(myNode)

    # Convert result list to a JSON string
    res = simplejson.dumps(myTree, encoding="Latin-1")

    # Support for the JSONP protocol.
    response_dict = {}
    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
    return HttpResponse(response_dict, content_type='application/json')

    response_dict = {}
    response_dict.update({'children': myTree})
    return HttpResponse(response_dict, content_type='application/javascript')


def getProjectTreeChildren(request):
    if request.is_ajax():
        projectid = request.GET["id"]
        samples = Sample.objects.filter(projectid=projectid)

        nodes = []
        for sample in samples:
            reads = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
            myNode = {
                'title': 'Sample: ' + sample.sample_name + '; Reads: ' + str(reads['count__sum']),
                'tooltip': sample.sample_name,
                'id': sample.sampleid,
                'isFolder': False
            }
            nodes.append(myNode)

        res = simplejson.dumps(nodes, encoding="Latin-1")
        return HttpResponse(res, content_type='application/json')


def getSampleCatTree(request):
    samples = Sample.objects.all()
    samples.query = pickle.loads(request.session['selected_samples'])

    projectList = samples.values_list('projectid').distinct()

    projectType = Project.objects.all().filter(projectid__in=projectList)
    typeList = []
    for p in projectType:
        typeList.append(p.projectType)

    myTree = {'title': 'Meta Data: Categorical', 'id': 'root', 'isFolder': False,  'hideCheckbox': True, 'expand': True, 'children': []}
    mimark = {'title': 'MIMARKs', 'id': 'mimark', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    air = {'title': 'Air', 'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    human_associated = {'title': 'Human Associated', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    microbial = {'title': 'Microbial', 'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    soil = {'title': 'Soil', 'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    water = {'title': 'Water', 'id': 'soil', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    user = {'title': 'User-defined', 'id': 'user', 'isFolder': True,  'hideCheckbox': True, 'children': []}

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

    if 'human_associated' in typeList:
        samp_collect = {'title': 'Sample Collection', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['samp_collect_device', 'samp_mat_process', 'samp_store_loc']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human_associated', 'isLazy': True, 'children': []}
            samp_collect['children'].append(myNode)
        human_associated['children'].append(samp_collect)

        samp_class = {'title': 'Sample Classification', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['samp_type', 'samp_location', 'samp_oxy_stat']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human_associated', 'isLazy': True, 'children': []}
            samp_class['children'].append(myNode)
        human_associated['children'].append(samp_class)

        host = {'title': 'Host', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['host_subject_id', 'host_gender', 'host_ethnicity', 'host_occupation', 'pet_farm_animal', 'obesity', 'smoker']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human_associated', 'isLazy': True, 'children': []}
            host['children'].append(myNode)
        human_associated['children'].append(host)

        diet = {'title': 'Host', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['diet_type', 'diet_frequency', 'diet_last_six_month', 'last_meal']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human_associated', 'isLazy': True, 'children': []}
            diet['children'].append(myNode)
        human_associated['children'].append(diet)

        disease = {'title': 'Disease', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['medic_hist_perform', 'disease_type', 'disease_location', 'tumor_location', 'tumor_stage']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human_associated', 'isLazy': True, 'children': []}
            disease['children'].append(myNode)
        human_associated['children'].append(disease)

        drug_use = {'title': 'Drug Usage', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['drug_usage', 'drug_type', 'drug_frequency']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human_associated', 'isLazy': True, 'children': []}
            drug_use['children'].append(myNode)
        human_associated['children'].append(drug_use)

        interven = {'title': 'Intervention', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['perturbation', 'pert_type', 'pert_frequency']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human_associated', 'isLazy': True, 'children': []}
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

    list = ['usr_cat1', 'usr_cat2', 'usr_cat3', 'usr_cat4', 'usr_cat5', 'usr_cat6']
    for i in range(len(list)):
        myNode = {'title': list[i], 'id': 'user', 'isFolder': True, 'pType': 'user', 'isLazy': True, 'children': []}
        user['children'].append(myNode)

    myTree['children'].append(mimark)
    if 'human_associated' in typeList:
        myTree['children'].append(human_associated)
    if 'soil' in typeList:
        myTree['children'].append(soil)

    myTree['children'].append(user)

    # Convert result list to a JSON string
    res = simplejson.dumps(myTree, encoding="Latin-1")

    # Support for the JSONP protocol.
    response_dict = {}
    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')

    response_dict = {}
    response_dict.update({'children': myTree})
    return HttpResponse(response_dict, content_type='application/javascript')


def getSampleCatTreeChildren(request):
    samples = Sample.objects.all()
    samples.query = pickle.loads(request.session['selected_samples'])
    selected = samples.values_list('sampleid')

    filterList = []
    for sample in selected:
        reads = Profile.objects.filter(sampleid=sample).aggregate(Sum('count'))
        if reads['count__sum'] is not None:
            filterList.append(sample[0])
    filtered = Sample.objects.all().filter(sampleid__in=filterList).values_list('sampleid')

    if request.is_ajax():
        field = request.GET["field"]
        pType = request.GET["pType"]
        mimark = Sample._meta.get_all_field_names()
        user = UserDefined._meta.get_all_field_names()
        air = Air._meta.get_all_field_names()
        human_associated = Human_Associated._meta.get_all_field_names()
        microbial = Microbial._meta.get_all_field_names()
        soil = Soil._meta.get_all_field_names()
        water = Water._meta.get_all_field_names()

        myNode = []
        if field in mimark:
            values = Sample.objects.values_list(field, flat='True').filter(sampleid__in=filtered).distinct().order_by(field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        #'tooltip': 'Value',
                        'isFolder': True,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        reads = Profile.objects.filter(sampleid=item.sampleid).aggregate(Sum('count'))
                        myNode2 = {
                            'title': 'Sample: ' + item.sample_name + '; Reads: ' + str(reads['count__sum']),
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'mimark',
                            'tooltip': 'Project: ' + item.projectid.project_name,
                            'hideCheckbox': True,
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
                        #'tooltip': 'Value',
                        'isFolder': True,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        reads = Profile.objects.filter(sampleid=item.sampleid).aggregate(Sum('count'))
                        myNode2 = {
                            'title': 'Sample: ' + item.sample_name + '; Reads: ' + str(reads['count__sum']),
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'air',
                            'tooltip': 'Project: ' + item.projectid.project_name,
                            'hideCheckbox': True,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        elif field in human_associated and pType == 'human_associated':
            table_field = 'human_associated__' + field
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=filtered).distinct().order_by(table_field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        #'tooltip': 'Value',
                        'isFolder': True,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        reads = Profile.objects.filter(sampleid=item.sampleid).aggregate(Sum('count'))
                        myNode2 = {
                            'title': 'Sample: ' + item.sample_name + '; Reads: ' + str(reads['count__sum']),
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'human_associated',
                            'tooltip': 'Project: ' + item.projectid.project_name,
                            'hideCheckbox': True,
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
                        #'tooltip': 'Value',
                        'isFolder': True,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        reads = Profile.objects.filter(sampleid=item.sampleid).aggregate(Sum('count'))
                        myNode2 = {
                            'title': 'Sample: ' + item.sample_name + '; Reads: ' + str(reads['count__sum']),
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'microbial',
                            'tooltip': 'Project: ' + item.projectid.project_name,
                            'hideCheckbox': True,
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
                        #'tooltip': 'Value',
                        'isFolder': True,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        reads = Profile.objects.filter(sampleid=item.sampleid).aggregate(Sum('count'))
                        myNode2 = {
                            'title': 'Sample: ' + item.sample_name + '; Reads: ' + str(reads['count__sum']),
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'soil',
                            'tooltip': 'Project: ' + item.projectid.project_name,
                            'hideCheckbox': True,
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
                        #'tooltip': 'Value',
                        'isFolder': True,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        reads = Profile.objects.filter(sampleid=item.sampleid).aggregate(Sum('count'))
                        myNode2 = {
                            'title': 'Sample: ' + item.sample_name + '; Reads: ' + str(reads['count__sum']),
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'water',
                            'tooltip': 'Project: ' + item.projectid.project_name,
                            'hideCheckbox': True,
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
                        #'tooltip': 'Value',
                        'isFolder': True,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        reads = Profile.objects.filter(sampleid=item.sampleid).aggregate(Sum('count'))
                        myNode2 = {
                            'title': 'Sample: ' + item.sample_name + '; Reads: ' + str(reads['count__sum']),
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'user',
                            'tooltip': 'Project: ' + item.projectid.project_name,
                            'hideCheckbox': True,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        res = simplejson.dumps(myNode, encoding="Latin-1")
        return HttpResponse(res, content_type='application/json')


def getSampleQuantTree(request):
    samples = Sample.objects.all()
    samples.query = pickle.loads(request.session['selected_samples'])

    projectList = samples.values_list('projectid').distinct()

    projectType = Project.objects.all().filter(projectid__in=projectList)
    typeList = []
    for p in projectType:
        typeList.append(p.projectType)

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

    if 'human_associated' in typeList:
        samp_collect = {'title': 'Sample Collection', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['samp_size', 'samp_store_temp', 'samp_store_dur']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human_associated', 'isLazy': True, 'children': []}
            samp_collect['children'].append(myNode)
        human_associated['children'].append(samp_collect)

        samp_class = {'title': 'Sample Classification', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['samp_temp', 'samp_ph', 'samp_salinity']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human_associated', 'isLazy': True, 'children': []}
            samp_class['children'].append(myNode)
        human_associated['children'].append(samp_class)

        host = {'title': 'Host', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['host_age', 'host_pulse', 'host_height', 'host_weight', 'host_bmi', 'host_weight_loss_3_month', 'host_body_temp']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human_associated', 'isLazy': True, 'children': []}
            host['children'].append(myNode)
        human_associated['children'].append(host)

        diet = {'title': 'Diet', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['diet_duration']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human_associated', 'isLazy': True, 'children': []}
            diet['children'].append(myNode)
        human_associated['children'].append(diet)

        disease = {'title': 'Disease', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['disease_duration', 'organism_count', 'tumor_mass']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human_associated', 'isLazy': True, 'children': []}
            disease['children'].append(myNode)
        human_associated['children'].append(disease)

        drug_use = {'title': 'Drug Usage', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['drug_duration']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human_associated', 'isLazy': True, 'children': []}
            drug_use['children'].append(myNode)
        human_associated['children'].append(drug_use)

        interven = {'title': 'Intervention', 'id': 'human_associated', 'isFolder': True,  'hideCheckbox': True, 'children': []}
        list = ['pert_duration']
        for i in range(len(list)):
            myNode = {'title': list[i], 'id': 'human_associated', 'isFolder': True, 'pType': 'human_associated', 'isLazy': True, 'children': []}
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

    list = ['usr_quant1', 'usr_quant2', 'usr_quant3', 'usr_quant4', 'usr_quant5', 'usr_quant6']
    for i in range(len(list)):
        myNode = {'title': list[i], 'id': 'user', 'isFolder': True, 'pType': 'user', 'isLazy': True, 'children': []}
        user['children'].append(myNode)

    myTree['children'].append(mimark)
    if 'human_associated' in typeList:
        myTree['children'].append(human_associated)
    if 'soil' in typeList:
        myTree['children'].append(soil)

    myTree['children'].append(user)

    # Convert result list to a JSON string
    res = simplejson.dumps(myTree, encoding="Latin-1")

    # Support for the JSONP protocol.
    response_dict = {}
    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')

    response_dict = {}
    response_dict.update({'children': myTree})
    return HttpResponse(response_dict, content_type='application/javascript')


def getSampleQuantTreeChildren(request):
    samples = Sample.objects.all()
    samples.query = pickle.loads(request.session['selected_samples'])
    selected = samples.values_list('sampleid')

    filterList = []
    for sample in selected:
        reads = Profile.objects.filter(sampleid=sample).aggregate(Sum('count'))
        if reads['count__sum'] is not None:
            filterList.append(sample[0])
    filtered = Sample.objects.all().filter(sampleid__in=filterList).values_list('sampleid')

    if request.is_ajax():
        field = request.GET["field"]
        pType = request.GET["pType"]
        mimark = Sample._meta.get_all_field_names()
        user = UserDefined._meta.get_all_field_names()
        air = Air._meta.get_all_field_names()
        human_associated = Human_Associated._meta.get_all_field_names()
        microbial = Microbial._meta.get_all_field_names()
        soil = Soil._meta.get_all_field_names()
        water = Water._meta.get_all_field_names()

        myNode = []
        if field in mimark:
            values = Sample.objects.values_list(field, flat='True').filter(sampleid__in=filtered).distinct().order_by(field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        #'tooltip': 'Value',
                        'isFolder': True,
                        'hideCheckbox': False,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        reads = Profile.objects.filter(sampleid=item.sampleid).aggregate(Sum('count'))
                        myNode2 = {
                            'title': 'Sample: ' + item.sample_name + '; Reads: ' + str(reads['count__sum']),
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'mimark',
                            'tooltip': 'Project: ' + item.projectid.project_name,
                            'hideCheckbox': True,
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
                        #'tooltip': 'Value',
                        'isFolder': True,
                        'hideCheckbox': False,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        reads = Profile.objects.filter(sampleid=item.sampleid).aggregate(Sum('count'))
                        myNode2 = {
                            'title': 'Sample: ' + item.sample_name + '; Reads: ' + str(reads['count__sum']),
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'air',
                            'tooltip': 'Project: ' + item.projectid.project_name,
                            'hideCheckbox': True,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        elif field in human_associated and pType == 'human_associated':
            table_field = 'human_associated__' + field
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=filtered).distinct().order_by(table_field)
            for j in range(len(values)):
                if pd.notnull(values[j]) and not values[j] == 'nan':
                    myNode1 = {
                        'title': values[j],
                        'id': field,
                        #'tooltip': 'Value',
                        'isFolder': True,
                        'hideCheckbox': False,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        reads = Profile.objects.filter(sampleid=item.sampleid).aggregate(Sum('count'))
                        myNode2 = {
                            'title': 'Sample: ' + item.sample_name + '; Reads: ' + str(reads['count__sum']),
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'human_associated',
                            'tooltip': 'Project: ' + item.projectid.project_name,
                            'hideCheckbox': True,
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
                        #'tooltip': 'Value',
                        'isFolder': True,
                        'hideCheckbox': False,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        reads = Profile.objects.filter(sampleid=item.sampleid).aggregate(Sum('count'))
                        myNode2 = {
                            'title': 'Sample: ' + item.sample_name + '; Reads: ' + str(reads['count__sum']),
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'microbial',
                            'tooltip': 'Project: ' + item.projectid.project_name,
                            'hideCheckbox': True,
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
                        #'tooltip': 'Value',
                        'isFolder': True,
                        'hideCheckbox': False,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        reads = Profile.objects.filter(sampleid=item.sampleid).aggregate(Sum('count'))
                        myNode2 = {
                            'title': 'Sample: ' + item.sample_name + '; Reads: ' + str(reads['count__sum']),
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'soil',
                            'tooltip': 'Project: ' + item.projectid.project_name,
                            'hideCheckbox': True,
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
                        #'tooltip': 'Value',
                        'isFolder': True,
                        'hideCheckbox': False,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        reads = Profile.objects.filter(sampleid=item.sampleid).aggregate(Sum('count'))
                        myNode2 = {
                            'title': 'Sample: ' + item.sample_name + '; Reads: ' + str(reads['count__sum']),
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'water',
                            'tooltip': 'Project: ' + item.projectid.project_name,
                            'hideCheckbox': True,
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
                        #'tooltip': 'Value',
                        'isFolder': True,
                        'hideCheckbox': False,
                        'children': []
                    }
                    args_list = []
                    args_list.append(Q(**{table_field: values[j]}))
                    items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=filtered).order_by('sample_name')
                    for item in items:
                        reads = Profile.objects.filter(sampleid=item.sampleid).aggregate(Sum('count'))
                        myNode2 = {
                            'title': 'Sample: ' + item.sample_name + '; Reads: ' + str(reads['count__sum']),
                            'id': item.sampleid,
                            'field': field,
                            'value': values[j],
                            'table': 'user',
                            'tooltip': 'Project: ' + item.projectid.project_name,
                            'hideCheckbox': True,
                            'isFolder': False
                        }
                        myNode1['children'].append(myNode2)
                    myNode.append(myNode1)

        res = simplejson.dumps(myNode, encoding="Latin-1")
        return HttpResponse(res, content_type='application/json')


def getTaxaTree(request):
    samples = Sample.objects.all()
    samples.query = pickle.loads(request.session['selected_samples'])
    selected = samples.values_list('sampleid')

    filterList = []
    for sample in selected:
        reads = Profile.objects.filter(sampleid=sample).aggregate(Sum('count'))
        if reads['count__sum'] is not None:
            filterList.append(sample[0])

    filtered = Sample.objects.filter(sampleid__in=filterList).values_list('sampleid')
    selected_taxa = Profile.objects.filter(sampleid_id__in=filtered)

    myTree = {'title': 'Taxa Name', 'tooltip': 'root', 'isFolder': False, 'hideCheckbox': True, 'expand': True, 'children': []}

    kingdoms = Kingdom.objects.all().filter(kingdomid__in=selected_taxa.values_list('kingdomid').distinct()).order_by('kingdomName')

    for kingdom in kingdoms:
        myNode = {
            'title': kingdom.kingdomName,
            'id': kingdom.kingdomid,
            'tooltip': "Kingdom",
            'isFolder': True,
            'expand': False,
            'children': []
        }
        phylas = kingdom.phyla_set.filter(phylaid__in=selected_taxa.values_list('phylaid').distinct()).order_by('phylaName')
        for phyla in phylas:
            myNode1 = {
                'title': phyla.phylaName,
                'id': phyla.phylaid,
                'tooltip': "Phyla",
                'isFolder': True,
                'isLazy': True
            }
            myNode['children'].append(myNode1)
        myTree['children'].append(myNode)

    # Convert result list to a JSON string
    res = simplejson.dumps(myTree, encoding="Latin-1")

    # Support for the JSONP protocol.
    response_dict = {}
    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')

    response_dict = {}
    response_dict.update({'children': myTree})
    return HttpResponse(response_dict, content_type='application/javascript')


def getTaxaTreeChildren(request):
    samples = Sample.objects.all()
    samples.query = pickle.loads(request.session['selected_samples'])
    selected = samples.values_list('sampleid')

    filterList = []
    for sample in selected:
        reads = Profile.objects.filter(sampleid=sample).aggregate(Sum('count'))
        if reads['count__sum'] is not None:
            filterList.append(sample[0])

    filtered = Sample.objects.filter(sampleid__in=filterList).values_list('sampleid')
    selected_taxa = Profile.objects.filter(sampleid_id__in=filtered)

    if request.is_ajax():
        taxa = request.GET["tooltip"]
        id = request.GET["id"]

        nodes = []
        if taxa == 'Phyla':
            qs = Class.objects.filter(classid__in=selected_taxa.values_list('classid').distinct()).filter(**{'phylaid': id}).order_by('className')
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
            qs = Order.objects.filter(orderid__in=selected_taxa.values_list('orderid').distinct()).filter(**{'classid': id}).order_by('orderName')
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
            qs = Family.objects.filter(familyid__in=selected_taxa.values_list('familyid').distinct()).filter(**{'orderid': id}).order_by('familyName')
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
            qs = Genus.objects.filter(genusid__in=selected_taxa.values_list('genusid').distinct()).filter(**{'familyid': id}).order_by('genusName')
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
            qs = Species.objects.filter(speciesid__in=selected_taxa.values_list('speciesid').distinct()).filter(**{'genusid': id}).order_by('speciesName')
            for item in qs:
                myNode = {
                    'title': item.speciesName,
                    'id': item.speciesid,
                    'tooltip': "Species",
                    'isLazy': True
                }
                nodes.append(myNode)

        res = simplejson.dumps(nodes, encoding="Latin-1")
        return HttpResponse(res, content_type='application/json')


def makeUpdateTree(request):
    myTree = {'title': 'All Uploads', 'isFolder': True, 'expand': True, 'hideCheckbox': True, 'children': []}

    projects = Project.objects.none()
    if request.user.is_superuser:
        projects = Project.objects.all()
    if request.user.is_authenticated():
        path_list = Reference.objects.filter(Q(author=request.user)).values_list('projectid_id')
        projects = Project.objects.all().filter( Q(projectid__in=path_list) )

    for project in projects:
        myNode = {
            'title': "Project: " + str(project.project_name),
            'tooltip': project.project_desc,
            'isFolder': True,
            'hideCheckbox': True,
            'children': []
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
    res = simplejson.dumps(myTree, encoding="Latin-1")

    # Support for the JSONP protocol.
    response_dict = {}
    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')

    response_dict = {}
    response_dict.update({'children': myTree})
    return HttpResponse(response_dict, content_type='application/javascript')


def makeReproTree(request):
    myTree = {'title': 'All Uploads', 'isFolder': True, 'expand': True, 'hideCheckbox': True, 'children': []}

    projects = Project.objects.none()
    if request.user.is_superuser:
        projects = Project.objects.all().filter(reference__raw=True)
    if request.user.is_authenticated():
        path_list = Reference.objects.filter(Q(author=request.user)).values_list('projectid_id')
        projects = Project.objects.all().filter( Q(projectid__in=path_list) ).filter(reference__raw=True)

    for project in projects:
        myNode = {
            'title': "Project: " + str(project.project_name),
            'tooltip': project.project_desc,
            'isFolder': True,
            'children': []
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

            myNode2={
                'title': ref.alignDB,
                'tooltip': 'Current align file...',
                'hideCheckbox': True,
                'isFolder': False
            }
            myNode1['children'].append(myNode2)

            myNode3={
                'title': ref.templateDB,
                'tooltip': 'Current template file...',
                'hideCheckbox': True,
                'isFolder': False
            }
            myNode1['children'].append(myNode3)

            myNode4={
                'title': ref.taxonomyDB,
                'tooltip': 'Current taxonomy file...',
                'hideCheckbox': True,
                'isFolder': False
            }
            myNode1['children'].append(myNode4)
            myNode['children'].append(myNode1)

        myTree['children'].append(myNode)

    # Convert result list to a JSON string
    res = simplejson.dumps(myTree, encoding="Latin-1")

    # Support for the JSONP protocol.
    response_dict = {}
    if 'callback' in request.GET:
        response_dict = request.GET['callback'] + "(" + res + ")"
        return HttpResponse(response_dict, content_type='application/json')

    response_dict = {}
    response_dict.update({'children': myTree})
    return HttpResponse(response_dict, content_type='application/javascript')


