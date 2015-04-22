import pickle
import operator
import simplejson
from django.http import HttpResponse
from django.db.models import Q, Sum
from models import Project, Sample
from models import Kingdom, Class, Order, Family, Genus, Species, Profile


def getProjectTree(request):
    myTree = {'title': 'All Projects', 'isFolder': True, 'expand': True, 'hideCheckbox': True, 'children': []}

    projects = Project.objects.all()

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
            if reads['count__sum'] is not None:
                myNode = {
                    'title': sample.sample_name + ' (' + str(reads['count__sum']) + ')',
                    'tooltip': sample.title,
                    'id': sample.sampleid,
                    'isFolder': False
                }
                nodes.append(myNode)

        res = simplejson.dumps(nodes, encoding="Latin-1")
        return HttpResponse(res, content_type='application/json')


def getSampleCatTree(request):
    myTree = {'title': 'Meta Data: Categorical', 'id': 'root', 'tooltip': 'root', 'isFolder': False,  'hideCheckbox': True, 'expand': True, 'children': []}
    mimark = {'title': 'MIMARKs', 'id': 'mimark', 'tooltip': 'Category', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    collect = {'title': 'Sample Collection', 'id': 'collect', 'tooltip': 'Category', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    soil_class = {'title': 'Soil Classification', 'id': 'soil_class', 'tooltip': 'Category', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    management = {'title': 'Management', 'id': 'management', 'tooltip': 'Category', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    user = {'title': 'User-defined', 'id': 'user', 'tooltip': 'Category', 'isFolder': True,  'hideCheckbox': True, 'children': []}

    list = ['sample_name', 'organism', 'seq_method', 'collection_date', 'biome', 'feature', 'geo_loc_country', 'geo_loc_state', 'geo_loc_city', 'geo_loc_farm', 'geo_loc_plot', 'material']
    for i in range(len(list)):
        myNode = {'title': list[i], 'isFolder': True, 'tooltip': 'Field', 'isLazy': True, 'children': []}
        mimark['children'].append(myNode)

    list = ['depth', 'pool_dna_extracts', 'samp_collection_device', 'sieving', 'storage_cond']
    for i in range(len(list)):
        myNode = {'title': list[i], 'isFolder': True, 'tooltip': 'Field', 'isLazy': True, 'children': []}
        collect['children'].append(myNode)

    list = ['drainage_class', 'fao_class', 'horizon', 'local_class', 'profile_position', 'slope_aspect', 'soil_type', 'texture_class']
    for i in range(len(list)):
        myNode = {'title': list[i], 'isFolder': True, 'tooltip': 'Field', 'isLazy': True, 'children': []}
        soil_class['children'].append(myNode)

    list = ['agrochem_amendments', 'agrochem_amendments_desc', 'biological_amendments', 'biological_amendments_desc', 'cover_crop', 'crop_rotation', 'cur_land_use', 'cur_vegetation', 'cur_crop', 'cur_cultivar', 'organic', 'previous_land_use', 'soil_amendments', 'soil_amendments_desc', 'tillage']
    for i in range(len(list)):
        myNode = {'title': list[i], 'isFolder': True, 'tooltip': 'Field', 'isLazy': True, 'children': []}
        management['children'].append(myNode)

    list = ['usr_cat1', 'usr_cat2', 'usr_cat3', 'usr_cat4', 'usr_cat5', 'usr_cat6']
    for i in range(len(list)):
        myNode = {'title': list[i], 'isFolder': True, 'tooltip': 'Field', 'isLazy': True, 'children': []}
        user['children'].append(myNode)

    myTree['children'].append(mimark)
    myTree['children'].append(collect)
    myTree['children'].append(soil_class)
    myTree['children'].append(management)
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

    if request.is_ajax():
        field = request.GET["field"]
        mimark = ['sample_name', 'organism', 'seq_method', 'collection_date', 'biome', 'feature', 'geo_loc_country', 'geo_loc_state', 'geo_loc_city',  'geo_loc_farm', 'geo_loc_plot', 'material']
        collect = ['depth', 'pool_dna_extracts', 'samp_collection_device', 'sieving', 'storage_cond']
        soil_class = ['drainage_class', 'fao_class', 'horizon', 'local_class', 'profile_position', 'slope_aspect', 'soil_type', 'texture_class']
        management = ['agrochem_amendments', 'agrochem_amendments_desc', 'biological_amendments', 'biological_amendments_desc', 'cover_crop', 'crop_rotation', 'cur_land_use', 'cur_vegetation', 'cur_crop', 'cur_cultivar', 'organic', 'previous_land_use', 'soil_amendments', 'soil_amendments_desc', 'tillage']
        user = ['usr_cat1', 'usr_cat2', 'usr_cat3', 'usr_cat4', 'usr_cat5', 'usr_cat6']

        myNode = []
        if field in mimark:
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            values = Sample.objects.values_list(field, flat='True').filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).distinct().order_by(field)
            for j in range(len(values)):
                myNode1 = {
                    'title': values[j],
                    'id': field,
                    'tooltip': 'Value',
                    'isFolder': True,
                    'children': []
                }
                args_list = []
                args_list.append(Q(**{field: values[j]}))
                items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).order_by('sample_name')
                for item in items:
                    myNode2 = {
                        'title': 'Sample: ' + item.sample_name,
                        'id': item.sampleid,
                        'tooltip': 'Project: ' + item.projectid.project_name,
                        'hideCheckbox': True,
                        'isFolder': False
                    }
                    myNode1['children'].append(myNode2)
                myNode.append(myNode1)

        elif field in collect:
            table_field = 'collect__' + field
            exclude_list = []
            exclude_list.append(Q(**{table_field: 'null'}))
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).distinct().order_by(table_field)
            for j in range(len(values)):
                myNode1 = {
                    'title': values[j],
                    'id': field,
                    'tooltip': 'Value',
                    'isFolder': True,
                    'children': []
                }
                args_list = []
                args_list.append(Q(**{table_field: values[j]}))
                items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).order_by('sample_name')
                for item in items:
                    myNode2 = {
                        'title': 'Sample: ' + item.sample_name,
                        'id': item.sampleid,
                        'tooltip': 'Project: ' + item.projectid.project_name,
                        'hideCheckbox': True,
                        'isFolder': False
                    }
                    myNode1['children'].append(myNode2)
                myNode.append(myNode1)

        elif field in soil_class:
            table_field = 'soil_class__' + field
            exclude_list = []
            exclude_list.append(Q(**{table_field: 'null'}))
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).distinct().order_by(table_field)
            for j in range(len(values)):
                myNode1 = {
                    'title': values[j],
                    'id': field,
                    'tooltip': 'Value',
                    'isFolder': True,
                    'children': []
                }
                args_list = []
                args_list.append(Q(**{table_field: values[j]}))
                items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).order_by('sample_name')
                for item in items:
                    myNode2 = {
                        'title': 'Sample: ' + item.sample_name,
                        'id': item.sampleid,
                        'tooltip': 'Project: ' + item.projectid.project_name,
                        'hideCheckbox': True,
                        'isFolder': False
                    }
                    myNode1['children'].append(myNode2)
                myNode.append(myNode1)

        elif field in management:
            table_field = 'management__' + field
            exclude_list = []
            exclude_list.append(Q(**{table_field: 'null'}))
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).distinct().order_by(table_field)
            for j in range(len(values)):
                myNode1 = {
                    'title': values[j],
                    'id': field,
                    'tooltip': 'Value',
                    'isFolder': True,
                    'children': []
                }
                args_list = []
                args_list.append(Q(**{table_field: values[j]}))
                items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).order_by('sample_name')
                for item in items:
                    myNode2 = {
                        'title': 'Sample: ' + item.sample_name,
                        'id': item.sampleid,
                        'tooltip': 'Project: ' + item.projectid.project_name,
                        'hideCheckbox': True,
                        'isFolder': False
                    }
                    myNode1['children'].append(myNode2)
                myNode.append(myNode1)


        elif field in user:
            table_field = 'user__' + field
            exclude_list = []
            exclude_list.append(Q(**{table_field: 'null'}))
            values = Sample.objects.values_list(table_field, flat='True').filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).distinct().order_by(table_field)
            for j in range(len(values)):
                myNode1 = {
                    'title': values[j],
                    'id': field,
                    'tooltip': 'Value',
                    'isFolder': True,
                    'children': []
                }
                args_list = []
                args_list.append(Q(**{table_field: values[j]}))
                items = Sample.objects.filter(reduce(operator.or_, args_list)).filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).order_by('sample_name')
                for item in items:
                    myNode2 = {
                        'title': 'Sample: ' + item.sample_name,
                        'id': item.sampleid,
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

    myTree = {'title': 'Meta Data: Quantitative', 'tooltip': 'root', 'isFolder': False,  'hideCheckbox': True, 'expand': True, 'children': []}
    mimark = {'title': 'MIMARKs', 'tooltip': 'MIMARKs', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    collect = {'title': 'Sample Collection', 'tooltip': 'Sample Collection', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    climate = {'title': 'Climate', 'tooltip': 'Climate', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    soil_class = {'title': 'Soil Classification', 'tooltip': 'Soil Classification', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    soil_nutrient = {'title': 'Soil pH / Nutrients', 'tooltip': 'Soil Nutrient', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    microbial = {'title': 'Microbial Biomass', 'tooltip': 'Microbial Biomass', 'isFolder': True,  'hideCheckbox': True, 'children': []}
    user = {'title': 'User-defined', 'tooltip': 'User_defined', 'isFolder': True,  'hideCheckbox': True, 'children': []}

    list = ['latitude', 'longitude', 'elevation']
    for i in range(len(list)):
        myNode = {'title': list[i], 'tooltip': 'mimark', 'isFolder': True, 'isLazy': True, 'children': []}
        mimark['children'].append(myNode)

    list = ['samp_size', 'samp_weight_dna_ext']
    for i in range(len(list)):
        myNode = {'title': list[i], 'tooltip': 'collect', 'isFolder': True, 'isLazy': True, 'children': []}
        collect['children'].append(myNode)

    list = ['annual_season_precpt', 'annual_season_temp']
    for i in range(len(list)):
        myNode = {'title': list[i], 'tooltip': 'climate', 'isFolder': True, 'isLazy': True, 'children': []}
        climate['children'].append(myNode)

    list = ['bulk_density', 'porosity', 'slope_gradient', 'water_content_soil']
    for i in range(len(list)):
        myNode = {'title': list[i], 'tooltip': 'soil_class', 'isFolder': True, 'isLazy': True, 'children': []}
        soil_class['children'].append(myNode)

    list = ['pH', 'EC', 'tot_C', 'tot_OM', 'tot_N', 'NO3_N', 'NH4_N', 'P', 'K', 'S', 'Zn', 'Fe', 'Cu', 'Mn', 'Ca', 'Mg', 'Na', 'B']
    for i in range(len(list)):
        myNode = {'title': list[i], 'tooltip': 'soil_nutrient', 'isFolder': True, 'isLazy': True, 'children': []}
        soil_nutrient['children'].append(myNode)

    list = ['rRNA_copies', 'microbial_biomass_C', 'microbial_biomass_N', 'microbial_respiration']
    for i in range(len(list)):
        myNode = {'title': list[i], 'tooltip': 'microbial', 'isFolder': True, 'isLazy': True, 'children': []}
        microbial['children'].append(myNode)

    list = ['usr_quant1', 'usr_quant2', 'usr_quant3', 'usr_quant4', 'usr_quant5', 'usr_quant6']
    for i in range(len(list)):
        myNode = {'title': list[i], 'tooltip': 'user', 'isFolder': True, 'isLazy': True, 'children': []}
        user['children'].append(myNode)

    myTree['children'].append(mimark)
    myTree['children'].append(collect)
    myTree['children'].append(climate)
    myTree['children'].append(soil_class)
    myTree['children'].append(soil_nutrient)
    myTree['children'].append(microbial)
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

    if request.is_ajax():
        field = request.GET["field"]
        mimark = ['latitude', 'longitude', 'elevation']
        collect = ['samp_size', 'samp_weight_dna_ext']
        climate = ['annual_season_precpt', 'annual_season_temp']
        soil_class = ['bulk_density', 'porosity', 'slope_gradient', 'water_content_soil']
        soil_nutrient = ['pH', 'EC', 'tot_C', 'tot_OM', 'tot_N', 'NO3_N', 'NH4_N', 'P', 'K', 'S', 'Zn', 'Fe', 'Cu', 'Mn', 'Ca', 'Mg', 'Na', 'B']
        microbial = ['rRNA_copies', 'microbial_biomass_C', 'microbial_biomass_N', 'microbial_respiration']
        user = ['usr_quant1', 'usr_quant2', 'usr_quant3', 'usr_quant4', 'usr_quant5', 'usr_quant6']

        myNode = []
        if field in mimark:
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            items = Sample.objects.filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).order_by('sample_name')
            for item in items:
                myNode1 = {
                    'title': 'Sample: ' + item.sample_name,
                    'id': item.sampleid,
                    'tooltip': 'Project: ' + item.projectid.project_name,
                    'hideCheckbox': True,
                    'isFolder': False
                }
                myNode.append(myNode1)

        elif field in collect:
            table_field = 'collect__' + field
            exclude_list = []
            exclude_list.append(Q(**{table_field: 'null'}))
            items = Sample.objects.filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).order_by('sample_name')
            for item in items:
                myNode1 = {
                    'title': 'Sample: ' + item.sample_name,
                    'id': item.sampleid,
                    'tooltip': 'Project: ' + item.projectid.project_name,
                    'hideCheckbox': True,
                    'isFolder': False
                }
                myNode.append(myNode1)

        elif field in climate:
            table_field = 'climate__' + field
            exclude_list = []
            exclude_list.append(Q(**{table_field: 'null'}))
            items = Sample.objects.filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).order_by('sample_name')
            for item in items:
                myNode1 = {
                    'title': 'Sample: ' + item.sample_name,
                    'id': item.sampleid,
                    'tooltip': 'Project: ' + item.projectid.project_name,
                    'hideCheckbox': True,
                    'isFolder': False
                }
                myNode.append(myNode1)

        elif field in soil_class:
            table_field = 'soil_class__' + field
            exclude_list = []
            exclude_list.append(Q(**{table_field: 'null'}))
            items = Sample.objects.filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).order_by('sample_name')
            for item in items:
                myNode1 = {
                    'title': 'Sample: ' + item.sample_name,
                    'id': item.sampleid,
                    'tooltip': 'Project: ' + item.projectid.project_name,
                    'hideCheckbox': True,
                    'isFolder': False
                }
                myNode.append(myNode1)

        elif field in soil_nutrient:
            table_field = 'soil_nutrient__' + field
            exclude_list = []
            exclude_list.append(Q(**{table_field: 'null'}))
            items = Sample.objects.filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).order_by('sample_name')
            for item in items:
                myNode1 = {
                    'title': 'Sample: ' + item.sample_name,
                    'id': item.sampleid,
                    'tooltip': 'Project: ' + item.projectid.project_name,
                    'hideCheckbox': True,
                    'isFolder': False
                }
                myNode.append(myNode1)

        elif field in microbial:
            table_field = 'microbial__' + field
            exclude_list = []
            exclude_list.append(Q(**{table_field: 'null'}))
            items = Sample.objects.filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).order_by('sample_name')
            for item in items:
                myNode1 = {
                    'title': 'Sample: ' + item.sample_name,
                    'id': item.sampleid,
                    'tooltip': 'Project: ' + item.projectid.project_name,
                    'hideCheckbox': True,
                    'isFolder': False
                }
                myNode.append(myNode1)

        elif field in user:
            table_field = 'user__' + field
            exclude_list = []
            exclude_list.append(Q(**{table_field: 'null'}))
            items = Sample.objects.filter(sampleid__in=selected).exclude(reduce(operator.or_, exclude_list)).order_by('sample_name')
            for item in items:
                myNode1 = {
                    'title': 'Sample: ' + item.sample_name,
                    'id': item.sampleid,
                    'tooltip': 'Project: ' + item.projectid.project_name,
                    'hideCheckbox': True,
                    'isFolder': False
                }
                myNode.append(myNode1)

        res = simplejson.dumps(myNode, encoding="Latin-1")
        return HttpResponse(res, content_type='application/json')


def getTaxaTree(request):
    samples = Sample.objects.all()
    samples.query = pickle.loads(request.session['selected_samples'])
    selected = samples.values_list('sampleid')
    selected_taxa = Profile.objects.filter(sampleid_id__in=selected)

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
    selected_taxa = Profile.objects.filter(sampleid_id__in=selected)

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


