import multiprocessing as mp
import numpy as np
import operator
import pandas as pd
from database.models import Kingdom, Phyla, Class, Order, Family, Genus, Species
from django.db.models import Q
from numpy import *
from numpy.random import mtrand
from scipy.spatial.distance import *
import math
from pyper import *


def catDiffAbundDF(qs1, metaDict):
    sampleTableList = ['sample_name', 'organism', 'seq_method', 'collection_date', 'biome', 'feature', 'geo_loc_country', 'geo_loc_state', 'geo_loc_city', 'geo_loc_farm', 'geo_loc_plot', 'material']
    soilTableList = ['depth', 'pool_dna_extracts', 'samp_collection_device', 'sieving', 'storage_cond', 'drainage_class', 'fao_class', 'horizon', 'local_class', 'profile_position', 'slope_aspect', 'soil_type', 'texture_class', 'agrochem_amendments', 'agrochem_amendments_desc', 'biological_amendments', 'biological_amendments_desc', 'cover_crop', 'crop_rotation', 'cur_land_use', 'cur_vegetation', 'cur_crop', 'cur_cultivar', 'organic', 'previous_land_use', 'soil_amendments', 'soil_amendments_desc', 'tillage']
    human_gutTableList = ['age', 'body_mass_index', 'body_product', 'chem_administration', 'diet', 'disease', 'ethnicity', 'family_relationship', 'grastointest_disord', 'genotype', 'height', 'host_body_temp', 'host_subject_id', 'ihmc_medication_code', 'last_meal', 'liver_disord', 'medic_hist_perform', 'nose_throat_disord', 'occupation', 'organism_count', 'oxy_stat_samp', 'perturbation', 'phenotype', 'pulse', 'rel_to_oxygen', 'samp_collect_device', 'samp_mat_process', 'sap_salinity', 'samp_size', 'samp_store_loc', 'samp_store_temp', 'sex', 'special_diet', 'temp', 'tissue', 'tot_mass', 'user_defined']
    microbialTableList = ['alkalinity', 'alkyl_diethers', 'altitude', 'aminopept_act', 'ammonium', 'bacteria_carb_prod', 'biomass', 'bishomohopanol', 'bromide', 'calcium', 'carb_nitro_ratio', 'chem_administration', 'chloride', 'chlorophyll', 'diether_lipids', 'diss_carb_dioxide', 'diss_hydrogen', 'diss_inorg_carb', 'diss_org_carb', 'diss_org_nitro', 'diss_oxygen', 'glucosidase_act', 'magnesium', 'mean_frict_vel', 'mean_peak_frict_vel', 'methane', 'n_alkanes', 'nitrate', 'nitrite', 'nitro', 'org_carb', 'org_matter', 'org_nitro', 'organism_count', 'oxy_stat_samp', 'part_org_carb', 'perturbation', 'petroleum_hydrocarb', 'ph', 'phaeopigments', 'phosphate', 'phosplipid_fatt_acid', 'potassium', 'pressure', 'redox_potential', 'rel_to_oxygen', 'salinity', 'samp_collect_device', 'samp_mat_process', 'samp_size', 'samp_store_dur', 'samp_store_loc', 'samp_store_temp', 'silicate', 'sodium', 'sulfate', 'sulfide', 'temp', 'tot_carb', 'tot_nitro', 'tot_org_carb', 'turbidity', 'water_content', 'user_defined']
    human_associatedTableList = ['age', 'amniotic_fluid_color', 'blood_blood_disord', 'body_mass_index', 'body_product', 'chem_administration', 'diet', 'diet_last_six_month', 'disease', 'drug_usage', 'ethnicity', 'family_relationship', 'fetal_health_stat', 'genotype', 'gestation_state', 'height', 'hiv_stat', 'host_body_temp', 'host_subject_id', 'ihmc_medication_code', 'kidney_disord', 'last_meal', 'maternal_health_stat', 'medic_hist_perform', 'nose_throat_disord', 'occupation', 'organism_count', 'oxy_stat_samp', 'perturbation', 'pet_farm_animal', 'phenotype', 'pulmonary_disord', 'pulse', 'rel_to_oxygen', 'samp_collect_device', 'samp_mat_process', 'samp_salinity', 'samp_size', 'samp_sotre_dur', 'samp_store_loc', 'samp_store_temp', 'sex', 'smoker', 'study_complt_stat', 'temp', 'tissue', 'tot_mass', 'travel_out_six_month', 'twin_sibling', 'urine_collect_meth', 'urogenit_tract_disor', 'weight_loss_3_month', 'user_defined']
    airTableList = ['barometric_press', 'carb_dioxide', 'carb_monoxide', 'chem_administration', 'elev', 'humidity', 'methane', 'organism_count', 'oxy_stat_samp', 'oxygen', 'perturbation', 'pollutants', 'rel_to_oxygen', 'resp_part_matter', 'samp_collect_device', 'samp_mat_process', 'samp_salinity', 'samp_size', 'samp_store_dur', 'samp_store_loc', 'samp_sotre_temp', 'solar_irradiance', 'temp', 'ventilation_rate', 'ventiliation_type', 'volatile_org_comp', 'wind_direction', 'wind_speed', 'user_defined']
    usrTableList = ['usr_cat1', 'usr_cat2', 'usr_cat3', 'usr_cat4', 'usr_cat5', 'usr_cat6']

    metaDF = pd.DataFrame()
    for key in metaDict:
        value = metaDict[key]
        args_list = []
        field_list = []

        if key in sampleTableList:
            field_list.append('sampleid')
            if key != 'sample_name':
                field_list.append('sample_name')
            field = str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key in soilTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'soil__' + str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key in human_gutTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'human_gut__' + str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key in human_associatedTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'human_associated__' + str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key in microbialTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'microbial__' + str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key in usrTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'user__' + str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

    return metaDF


