import multiprocessing as mp
import numpy as np
import operator
import pandas as pd
from database.models import Kingdom, Phyla, Class, Order, Family, Genus, Species, OTU_01, OTU_03
from django.db.models import Q
from numpy import *
from numpy.random import mtrand
from scipy.spatial.distance import *
import math
from pyper import *


def catUnivMetaDF(qs1, metaDict):
    sampleTableList = ['sample_name', 'organism', 'seq_method', 'collection_date', 'biome', 'feature', 'geo_loc_country', 'geo_loc_state', 'geo_loc_city', 'geo_loc_farm', 'geo_loc_plot', 'material']
    soilTableList = ['depth', 'pool_dna_extracts', 'samp_collection_device', 'sieving', 'storage_cond', 'drainage_class', 'fao_class', 'horizon', 'local_class', 'profile_position', 'slope_aspect', 'soil_type', 'texture_class', 'agrochem_amendments', 'agrochem_amendments_desc', 'biological_amendments', 'biological_amendments_desc', 'cover_crop', 'crop_rotation', 'cur_land_use', 'cur_vegetation', 'cur_crop', 'cur_cultivar', 'organic', 'previous_land_use', 'soil_amendments', 'soil_amendments_desc', 'tillage']
    human_gutTableList = ['age', 'body_mass_index', 'body_product', 'chem_administration', 'diet', 'disease', 'ethnicity', 'family_relationship', 'grastointest_disord', 'genotype', 'height', 'host_body_temp', 'host_subject_id', 'ihmc_medication_code', 'last_meal', 'liver_disord', 'medic_hist_perform', 'nose_throat_disord', 'occupation', 'organism_count', 'oxy_stat_samp', 'perturbation', 'phenotype', 'pulse', 'rel_to_oxygen', 'samp_collect_device', 'samp_mat_process', 'sap_salinity', 'samp_size', 'samp_store_loc', 'samp_store_temp', 'sex', 'special_diet', 'temp', 'tissue', 'tot_mass', 'user_defined']
    microbialTableList = ['alkalinity', 'alkyl_diethers', 'altitude', 'aminopept_act', 'ammonium', 'bacteria_carb_prod', 'biomass', 'bishomohopanol', 'bromide', 'calcium', 'carb_nitro_ratio', 'chem_administration', 'chloride', 'chlorophyll', 'diether_lipids', 'diss_carb_dioxide', 'diss_hydrogen', 'diss_inorg_carb', 'diss_org_carb', 'diss_org_nitro', 'diss_oxygen', 'glucosidase_act', 'magnesium', 'mean_frict_vel', 'mean_peak_frict_vel', 'methane', 'n_alkanes', 'nitrate', 'nitrite', 'nitro', 'org_carb', 'org_matter', 'org_nitro', 'organism_count', 'oxy_stat_samp', 'part_org_carb', 'perturbation', 'petroleum_hydrocarb', 'ph', 'phaeopigments', 'phosphate', 'phosplipid_fatt_acid', 'potassium', 'pressure', 'redox_potential', 'rel_to_oxygen', 'salinity', 'samp_collect_device', 'samp_mat_process', 'samp_size', 'samp_store_dur', 'samp_store_loc', 'samp_store_temp', 'silicate', 'sodium', 'sulfate', 'sulfide', 'temp', 'tot_carb', 'tot_nitro', 'tot_org_carb', 'turbidity', 'water_content', 'user_defined']
    waterTableList = ['alkalinity', 'alkyl_diethers', 'altitude', 'aminopept_act', 'ammonium', 'atmospheric_data', 'bac_prod', 'bac_resp', 'bacteria_carb_prod', 'biomass', 'bishomohopanol', 'bromide', 'calcium', 'carb_nitro_ratio', 'chem_administration', 'chloride', 'chlorophyll', 'conduc', 'density', 'diether_lipids', 'diss_carb_dioxide', 'diss_hydrogen', 'diss_inorg_carb', 'diss_inorg_nitro', 'diss_inorg_phosp', 'diss_org_carb', 'diss_org_nitro', 'diss_oxygen', 'down_par', 'elev', 'fluor', 'glucosidase_act', 'light_intensity', 'magnesium', 'mean_frict_vel', 'mean_peak_frict_vel', 'n_alkanes', 'nitrate', 'nitrite', 'nitro', 'org_carb', 'org_matter', 'org_nitro', 'organism_count', 'oxy_stat_samp', 'part_org_carb', 'part_org_nitro', 'perturbation', 'pretroleum_hydrocarb', 'ph', 'phaeopigments', 'phosphate', 'phosplipid_fatt_acid', 'photon_flux', 'potassium', 'pressure', 'primary_prod', 'redox_potential', 'rel_to_oxygen', 'samp_mat_process', 'samp_salinity', 'samp_size', 'samp_store_dur', 'samp_store_loc', 'samp_store_temp', 'samp_vol_we_dna_ext', 'silicate', 'sodium', 'soluble_react_phosp', 'source_material_id', 'sulfate', 'sulfide', 'suspen_part_matter', 'temp', 'tidal_stage', 'tot_depth_water_col', 'tot_diss_nitro', 'tot_inorg_nitro', 'tot_nitro', 'tot_part_carb', 'tot_phosp', 'water_current', 'user_defined']
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

        elif key in airTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'air__' + str(key)
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

        elif key in waterTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'water__' + str(key)
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


def quantUnivMetaDF(qs1, metaDict):
    metaDF = pd.DataFrame()
    final_fieldList = []
    for key in metaDict:
        value = metaDict[key]
        field_list = []

        if key == 'mimark':
            field_list.append('sampleid')
            field_list.append('sample_name')
            field_list.append(value)
            final_fieldList.append(value)
            exclude_list = []
            exclude_list.append(Q(**{value: 'null'}))
            qs2 = qs1.values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            if metaDF.empty:
                metaDF = tempDF
            else:
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key == 'soil':
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'soil__' + str(value)
            field_list.append(field)
            final_fieldList.append(value)
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: value}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key == 'air':
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'air__' + str(value)
            field_list.append(field)
            final_fieldList.append(value)
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: value}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key == 'water':
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'water__' + str(value)
            field_list.append(field)
            final_fieldList.append(value)
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: value}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key == 'human_associated':
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'human_associated__' + str(value)
            field_list.append(field)
            final_fieldList.append(value)
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: value}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key == 'microbial':
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'microbial__' + str(value)
            field_list.append(field)
            final_fieldList.append(value)
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: value}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key == 'human_gut':
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'human_gut__' + str(value)
            field_list.append(field)
            final_fieldList.append(value)
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: value}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key == 'user':
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'user__' + str(value)
            field_list.append(field)
            final_fieldList.append(value)
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: value}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

    return metaDF


def normalizeUniv(df, taxaDict, mySet, meth, reads, metaDF):
    df2 = df.reset_index()
    taxaID = ['kingdomid', 'phylaid', 'classid', 'orderid', 'familyid', 'genusid', 'speciesid', 'otuid1', 'otuid3']

    countDF = pd.DataFrame()
    DESeq_error = 'no'
    if meth == 1 or meth == 4:
        countDF = df2.reset_index(drop=True)

    elif meth == 2 or meth == 3:
        if reads >= 0:
            countDF = df2[taxaID].reset_index(drop=True)
            manager = mp.Manager()
            d = manager.dict()

            numcore = mp.cpu_count()-1 or 1
            processes = [mp.Process(target=weightedProb, args=(x, numcore, reads, mySet, df, meth, d)) for x in range(numcore)]

            for p in processes:
                p.start()
            for p in processes:
                p.join()

            for key, value in d.items():
                countDF[key] = value

        elif reads < 0:
            countDF = df2.reset_index(drop=True)

    elif meth == 5:
        countDF = df2[taxaID].reset_index(drop=True)
        r = R(RCMD="R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
        df3 = df2.drop(taxaID, axis=1)
        r.assign("count", df3)
        r.assign("metaDF", metaDF)
        r("trt <- factor(metaDF$merge)")

        r("library(DESeq2)")
        r("colData <- data.frame(row.names=colnames(count), trt=trt)")
        r("dds <- DESeqDataSetFromMatrix(countData=count, colData=colData, design= ~ trt)")
        r("dds <- estimateSizeFactors(dds)")
        pycds = r.get("sizeFactors(dds)")
        colList = df3.columns.tolist()

        found = False
        if pycds is list:
            for thing in pycds:
                if str(thing) == "None":
                    found = True
        else:
            if pycds is None:
                found = True

        if not found:
            DESeq_error = 'no'
            cdsDF = pd.DataFrame(r.get("counts(dds, normalize=TRUE)"), columns=[colList])
            countDF[colList] = cdsDF[colList]
        else:
            DESeq_error = 'yes'
            countDF = df2.reset_index(drop=True)

    relabundDF = pd.DataFrame(countDF[taxaID])
    binaryDF = pd.DataFrame(countDF[taxaID])
    diversityDF = pd.DataFrame(countDF[taxaID])
    for i in mySet:
        relabundDF[i] = countDF[i].div(countDF[i].sum(), axis=0)
        binaryDF[i] = countDF[i].apply(lambda x: 1 if x != 0 else 0)
        diversityDF[i] = relabundDF[i].apply(lambda x: -1 * x * math.log(x) if x > 0 else 0)

    rowsList = []
    namesDF = pd.DataFrame()
    normDF = pd.DataFrame()
    for key in taxaDict:
        taxaList = taxaDict[key]

        if isinstance(taxaList, unicode):
            if key == 'Kingdom':
                qs1 = Kingdom.objects.filter(kingdomid=taxaList).values('kingdomid', 'kingdomName')
                namesDF = pd.DataFrame.from_records(qs1, columns=['kingdomid', 'kingdomName'])
                namesDF.rename(columns={'kingdomid': 'taxa_id', 'kingdomName': 'taxa_name'}, inplace=True)
            elif key == 'Phyla':
                qs1 = Phyla.objects.filter(phylaid=taxaList).values('phylaid', 'phylaName')
                namesDF = pd.DataFrame.from_records(qs1, columns=['phylaid', 'phylaName'])
                namesDF.rename(columns={'phylaid': 'taxa_id', 'phylaName': 'taxa_name'}, inplace=True)
            elif key == 'Class':
                qs1 = Class.objects.filter(classid=taxaList).values('classid', 'className')
                namesDF = pd.DataFrame.from_records(qs1, columns=['classid', 'className'])
                namesDF.rename(columns={'classid': 'taxa_id', 'className': 'taxa_name'}, inplace=True)
            elif key == 'Order':
                qs1 = Order.objects.filter(orderid=taxaList).values('orderid', 'orderName')
                namesDF = pd.DataFrame.from_records(qs1, columns=['orderid', 'orderName'])
                namesDF.rename(columns={'orderid': 'taxa_id', 'orderName': 'taxa_name'}, inplace=True)
            elif key == 'Family':
                qs1 = Family.objects.filter(familyid=taxaList).values('familyid', 'familyName')
                namesDF = pd.DataFrame.from_records(qs1, columns=['familyid', 'familyName'])
                namesDF.rename(columns={'familyid': 'taxa_id', 'familyName': 'taxa_name'}, inplace=True)
            elif key == 'Genus':
                qs1 = Genus.objects.filter(genusid=taxaList).values('genusid', 'genusName')
                namesDF = pd.DataFrame.from_records(qs1, columns=['genusid', 'genusName'])
                namesDF.rename(columns={'genusid': 'taxa_id', 'genusName': 'taxa_name'}, inplace=True)
            elif key == 'Species':
                qs1 = Species.objects.filter(speciesid=taxaList).values('speciesid', 'speciesName')
                namesDF = pd.DataFrame.from_records(qs1, columns=['speciesid', 'speciesName'])
                namesDF.rename(columns={'speciesid': 'taxa_id', 'speciesName': 'taxa_name'}, inplace=True)

            elif key == 'OTU_0.01':
                qs1 = OTU_01.objects.filter(otuid1=taxaList).values('otuid1')
                namesDF = pd.DataFrame.from_records(qs1, columns=['otuid1'])
                namesDF.rename(columns={'otuid1': 'taxa_id'}, inplace=True)
            elif key == 'OTU_0.03':
                qs1 = OTU_03.objects.filter(otuid3=taxaList).values('otuid3')
                namesDF = pd.DataFrame.from_records(qs1, columns=['otuid3'])
                namesDF.rename(columns={'otuid3': 'taxa_id'}, inplace=True)

        else:
            if key == 'Kingdom':
                qs1 = Kingdom.objects.filter(kingdomid__in=taxaList).values('kingdomid', 'kingdomName')
                namesDF = pd.DataFrame.from_records(qs1, columns=['kingdomid', 'kingdomName'])
                namesDF.rename(columns={'kingdomid': 'taxa_id', 'kingdomName': 'taxa_name'}, inplace=True)
            elif key == 'Phyla':
                qs1 = Phyla.objects.filter(phylaid__in=taxaList).values('phylaid', 'phylaName')
                namesDF = pd.DataFrame.from_records(qs1, columns=['phylaid', 'phylaName'])
                namesDF.rename(columns={'phylaid': 'taxa_id', 'phylaName': 'taxa_name'}, inplace=True)
            elif key == 'Class':
                qs1 = Class.objects.filter(classid__in=taxaList).values('classid', 'className')
                namesDF = pd.DataFrame.from_records(qs1, columns=['classid', 'className'])
                namesDF.rename(columns={'classid': 'taxa_id', 'className': 'taxa_name'}, inplace=True)
            elif key == 'Order':
                qs1 = Order.objects.filter(orderid__in=taxaList).values('orderid', 'orderName')
                namesDF = pd.DataFrame.from_records(qs1, columns=['orderid', 'orderName'])
                namesDF.rename(columns={'orderid': 'taxa_id', 'orderName': 'taxa_name'}, inplace=True)
            elif key == 'Family':
                qs1 = Family.objects.filter(familyid__in=taxaList).values('familyid', 'familyName')
                namesDF = pd.DataFrame.from_records(qs1, columns=['familyid', 'familyName'])
                namesDF.rename(columns={'familyid': 'taxa_id', 'familyName': 'taxa_name'}, inplace=True)
            elif key == 'Genus':
                qs1 = Genus.objects.filter(genusid__in=taxaList).values('genusid', 'genusName')
                namesDF = pd.DataFrame.from_records(qs1, columns=['genusid', 'genusName'])
                namesDF.rename(columns={'genusid': 'taxa_id', 'genusName': 'taxa_name'}, inplace=True)
            elif key == 'Species':
                qs1 = Species.objects.filter(speciesid__in=taxaList).values('speciesid', 'speciesName')
                namesDF = pd.DataFrame.from_records(qs1, columns=['speciesid', 'speciesName'])
                namesDF.rename(columns={'speciesid': 'taxa_id', 'speciesName': 'taxa_name'}, inplace=True)

            elif key == 'OTU_0.01':
                qs1 = OTU_01.objects.filter(otuid1=taxaList).values('otuid1')
                namesDF = pd.DataFrame.from_records(qs1, columns=['otuid1'])
                namesDF.rename(columns={'otuid1': 'taxa_id', 'speciesName': 'taxa_name'}, inplace=True)
            elif key == 'OTU_0.03':
                qs1 = OTU_03.objects.filter(otuid3=taxaList).values('otuid3')
                namesDF = pd.DataFrame.from_records(qs1, columns=['otuid3'])
                namesDF.rename(columns={'otuid3': 'taxa_id'}, inplace=True)

        if key == 'Kingdom':
            rank = 'Kingdom'
            field = 'kingdomid'
        elif key == 'Phyla':
            rank = 'Phyla'
            field = 'phylaid'
        elif key == 'Class':
            rank = 'Class'
            field = 'classid'
        elif key == 'Order':
            rank = 'Order'
            field = 'orderid'
        elif key == 'Family':
            rank = 'Family'
            field = 'familyid'
        elif key == 'Genus':
            rank = 'Genus'
            field = 'genusid'
        elif key == 'Species':
            rank = 'Species'
            field = 'speciesid'

        elif key == 'OTU_0.01':
            rank = 'OTU_0.01'
            field = 'otuid1'
        elif key == 'OTU_0.03':
            rank = 'OTU_0.03'
            field = 'otuid3'

        for i in mySet:
            if meth == 4:
                groupAbund = relabundDF.groupby(field)[i].sum()
            else:
                groupAbund = countDF.groupby(field)[i].sum()
            groupRich = binaryDF.groupby(field)[i].sum()
            groupDiversity = diversityDF.groupby(field)[i].sum()
            if isinstance(taxaList, unicode):
                myDict = {}
                myDict['sampleid'] = i
                myDict['rank'] = rank
                myDict['taxa_id'] = taxaList
                myDict['abund'] = groupAbund[taxaList]
                myDict['rich'] = groupRich[taxaList]
                myDict['diversity'] = groupDiversity[taxaList]
                rowsList.append(myDict)
            else:
                for j in taxaList:
                    myDict = {}
                    myDict['sampleid'] = i
                    myDict['rank'] = rank
                    myDict['taxa_id'] = j
                    myDict['abund'] = groupAbund[j]
                    myDict['rich'] = groupRich[j]
                    myDict['diversity'] = groupDiversity[j]
                    rowsList.append(myDict)
        DF1 = pd.DataFrame(rowsList, columns=['sampleid', 'rank', 'taxa_id', 'abund', 'rich', 'diversity'])
        DF1 = DF1.merge(namesDF, on='taxa_id', how='outer')
        otupres = False
        for key in taxaDict:
            if key == 'OTU_0.01' or 'OTU_0.03':
                otupres = True
        if otupres:
            DF1 = DF1[['sampleid', 'rank', 'taxa_id', 'abund', 'rich', 'diversity']]
        else:
            DF1 = DF1[['sampleid', 'rank', 'taxa_id', 'taxa_name', 'abund', 'rich', 'diversity']]

        if normDF.empty:
            normDF = DF1
        else:
            normDF = normDF.append(DF1)
    return normDF, DESeq_error


def weightedProb(x, cores, reads, mySet, df, meth, d):
    high = mySet.__len__()
    set = mySet[x:high:cores]

    for i in set:
        arr = asarray(df[i])
        cols = shape(arr)
        sample = arr.astype(dtype=np.float64)

        if meth == 3:
            prob = (sample + 0.1) / (sample.sum() + cols[0] * 0.1)
        else:
            prob = sample / sample.sum()

        temp = np.zeros(cols)
        for n in range(reads):
            sub = np.random.mtrand.choice(range(sample.size), size=1, replace=False, p=prob)
            temp2 = np.zeros(cols)
            np.put(temp2, sub, 1)
            temp = np.core.umath.add(temp, temp2)
        d[i] = temp

