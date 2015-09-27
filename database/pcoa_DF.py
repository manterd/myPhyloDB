import operator
import pandas as pd
from django.db.models import Q
import numpy as np
from numpy import *
from numpy.random import mtrand
from scipy.spatial.distance import *
import multiprocessing as mp
from pyper import *


def catPCoAMetaDF(qs1, metaDict):
    sampleTableList = ['sample_name', 'organism', 'material', 'seq_platform', 'seq_gene', 'seq_gene_region', 'seq_for_primer', 'seq_rev_primer', 'collection_date', 'biome', 'feature', 'geo_loc_country', 'geo_loc_state', 'geo_loc_city', 'geo_loc_farm', 'geo_loc_plot']
    human_associatedTableList = ['samp_collect_device', 'samp_mat_process', 'samp_size', 'samp_store_temp', 'samp_store_dur', 'samp_type', 'samp_location', 'samp_temp', 'samp_ph', 'samp_oxy_stat', 'samp_salinity', 'host_subject_id', 'host_age', 'host_pulse', 'host_gender', 'host_ethnicity', 'host_height', 'host_weight', 'host_bmi', 'host_weight_loss_3_month', 'host_body_temp', 'host_occupation', 'pet_farm_animal', 'smoker', 'diet_type', 'diet_duration', 'diet_frequency', 'diet_last_six_month', 'last_meal', 'medic_hist_perform', 'disease_type', 'disease_location', 'disease_duration', 'organism_count', 'tumor_location', 'tumor_mass', 'tumor_stage', 'drug_usage', 'drug_type', 'drug_duration', 'drug_frequency', 'perturbation', 'pert_type', 'pert_duration', 'pert_frequency', 'fetal_health_stat', 'amniotic_fluid_color', 'gestation_stat', 'maternal_health_stat']
    soilTableList = ['samp_collection_device', 'samp_size', 'samp_depth', 'sieve_size', 'storage_cond', 'samp_weight_dna_ext', 'pool_dna_extracts', 'fao_class', 'local_class', 'texture_class', 'porosity', 'profile_position', 'slope_aspect', 'slope_gradient', 'bulk_density', 'drainage_class', 'water_content_soil', 'cur_land_use', 'cur_vegetation', 'cur_crop', 'cur_cultivar', 'crop_rotation', 'cover_crop', 'fert_amendment_class', 'fert_placement', 'fert_type', 'fert_tot_amount', 'fert_N_tot_amount', 'fert_P_tot_amount', 'fert_K_tot_amount', 'irrigation_type', 'irrigation_tot_amount', 'residue_removal', 'residue_growth_stage', 'residue_removal_percent', 'tillage_event', 'tillage_event_depth', 'amend1_class', 'amend1_active_ingredient', 'amend1_tot_amount', 'amend2_class', 'amend2_active_ingredient', 'amend2_tot_amount', 'amend3_class', 'amend3_active_ingredient', 'amend3_tot_amount', 'rRNA_copies', 'microbial_biomass_C', 'microbial_biomass_N', 'microbial_respiration', 'soil_pH', 'soil_EC', 'soil_C', 'soil_OM', 'soil_N', 'soil_NO3_N', 'soil_NH4_N', 'soil_P', 'soil_K', 'soil_S', 'soil_Zn', 'soil_Fe', 'soil_Cu', 'soil_Mn', 'soil_Ca', 'soil_Mg', 'soil_Na', 'soil_B', 'plant_C', 'plant_N', 'plant_P', 'plant_K', 'plant_Ca', 'plant_Mg', 'plant_S', 'plant_Na', 'plant_Cl', 'plant_Al', 'plant_B', 'plant_Cu', 'plant_Fe', 'plant_Mn', 'plant_Zn', 'crop_tot_biomass_fw', 'crop_tot_biomass_dw', 'crop_tot_above_biomass_fw', 'crop_tot_above_biomass_dw', 'crop_tot_below_biomass_fw', 'crop_tot_below_biomass_dw', 'harv_fraction', 'harv_fresh_weight', 'harv_dry_weight', 'ghg_chamber_placement', 'ghg_N2O', 'ghg_CO2', 'ghg_NH4']
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
            qs2 = qs1.values(*field_list).filter(reduce(operator.or_, args_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list).dropna()
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
            qs2 = qs1.values(*field_list).filter(reduce(operator.or_, args_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list).dropna()
            tempDF.rename(columns={field: key}, inplace=True)
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
            qs2 = qs1.values(*field_list).filter(reduce(operator.or_, args_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list).dropna()
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key in usrTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'userdefined__' + str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            qs2 = qs1.values(*field_list).filter(reduce(operator.or_, args_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list).dropna()
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

    return metaDF


def quantPCoAMetaDF(qs1, metaDict):
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
            qs2 = qs1.values(*field_list)
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list).dropna()
            if metaDF.empty:
                metaDF = tempDF
                metaDF[value] = metaDF[value].astype(float)
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')
                metaDF[value] = metaDF[value].astype(float)

        elif key == 'human_associated':
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'human_associated__' + str(value)
            field_list.append(field)
            final_fieldList.append(value)
            qs2 = qs1.values(*field_list)
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list).dropna()
            tempDF.rename(columns={field: value}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
                metaDF[value] = metaDF[value].astype(float)
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')
                metaDF[value] = metaDF[value].astype(float)

        elif key == 'soil':
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'soil__' + str(value)
            field_list.append(field)
            final_fieldList.append(value)
            qs2 = qs1.values(*field_list)
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list).dropna()
            tempDF.rename(columns={field: value}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
                metaDF[value] = metaDF[value].astype(float)
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')
                metaDF[value] = metaDF[value].astype(float)

        elif key == 'user':
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'userdefined__' + str(value)
            field_list.append(field)
            final_fieldList.append(value)
            qs2 = qs1.values(*field_list)
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list).dropna()
            tempDF.rename(columns={field: value}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
                metaDF[value] = metaDF[value].astype(float)
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')
                metaDF[value] = metaDF[value].astype(float)

    return metaDF


def normalizePCoA(df, taxaLevel, mySet, meth, reads, metaDF):
    df2 = df.reset_index()

    taxaID = ''
    if taxaLevel == 0:
        taxaID = 'kingdomid'
    elif taxaLevel == 1:
        taxaID = 'phylaid'
    elif taxaLevel == 2:
        taxaID = 'classid'
    elif taxaLevel == 3:
        taxaID = 'orderid'
    elif taxaLevel == 4:
        taxaID = 'familyid'
    elif taxaLevel == 5:
        taxaID = 'genusid'
    elif taxaLevel == 6:
        taxaID = 'speciesid'

    normDF = pd.DataFrame()
    DESeq_error = ''
    if meth == 1:
        normDF = df2.reset_index(drop=True)

    elif meth == 2 or meth == 3:
        if reads >= 0:
            normDF[taxaID] = df2[taxaID].reset_index(drop=True)

            manager = mp.Manager()
            d = manager.dict()

            numcore = mp.cpu_count()-1 or 1
            processes = [mp.Process(target=weightedProb, args=(x, numcore, reads, mySet, df, meth, d)) for x in range(numcore)]

            for p in processes:
                p.start()
            for p in processes:
                p.join()

            for key, value in d.items():
                normDF[key] = value

        elif reads < 0:
            normDF[taxaID] = df2.reset_index(drop=True)

    elif meth == 4:
        normDF[taxaID] = df2[taxaID].reset_index(drop=True)
        for i in mySet:
            normDF[i] = df2[i].div(df2[i].sum(), axis=0)

    elif meth == 5:
        if os.name == 'nt':
            r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
        else:
            r = R(RCMD="R/R-Linux/bin/R")
        df3 = df2.drop(taxaID, axis=1)
        r.assign("count", df3)
        r.assign("metaDF", metaDF)
        r("trt <- factor(metaDF$merge)")

        r("library(DESeq2)")
        r("colData <- data.frame(row.names=colnames(count), trt=trt)")
        r("dds <- DESeqDataSetFromMatrix(countData=count, colData=colData, design= ~ trt)")
        r("dds <- estimateSizeFactors(dds)")
        pycds = r.get("sizeFactors(dds)")

        if pycds is not None:
            DESeq_error = 'no'
            colList = df3.columns.tolist()
            indList = df2[taxaID].tolist()
            normDF = pd.DataFrame(r.get("counts(dds, normalize=TRUE)"), columns=[colList])
            normDF[taxaID] = indList
        elif pycds is None:
            DESeq_error = 'yes'
            sizeFactor = []
            for i in mySet:
                sizeFactor.append(1)
            r.assign("sizeFactor", sizeFactor)
            r("dds$sizeFactor <- sizeFactor")
            colList = df3.columns.tolist()
            indList = df2[taxaID].tolist()
            normDF = pd.DataFrame(r.get("counts(dds, normalize=TRUE)"), columns=[colList])
            normDF[taxaID] = indList

    elif meth == 6:
        if os.name == 'nt':
            r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
        else:
            r = R(RCMD="R/R-Linux/bin/R")
        df3 = df2.drop(taxaID, axis=1)
        r.assign("count", df3)
        r.assign("metaDF", metaDF)
        r("trt <- factor(metaDF$merge)")

        r("library(DESeq2)")
        r("colData <- data.frame(row.names=colnames(count), trt=trt)")
        r("dds <- DESeqDataSetFromMatrix(countData=count, colData=colData, design= ~ trt)")
        r("dds <- estimateSizeFactors(dds)")
        pycds = r.get("sizeFactors(dds)")

        if pycds is not None:
            DESeq_error = 'no'
            r("vsd <- varianceStabilizingTransformation(dds)")
            colList = df3.columns.tolist()
            indList = df2[taxaID].tolist()
            normDF = pd.DataFrame(r.get("assay(vsd)"), columns=[colList])
            normDF[taxaID] = indList
        elif pycds is None:
            DESeq_error = 'yes'
            sizeFactor = []
            for i in mySet:
                sizeFactor.append(1)
            r.assign("sizeFactor", sizeFactor)
            r("dds$sizeFactor <- sizeFactor")
            r("vsd <- varianceStabilizingTransformation(dds)")
            colList = df3.columns.tolist()
            indList = df2[taxaID].tolist()
            normDF = pd.DataFrame(r.get("assay(vsd)"), columns=[colList])
            normDF[taxaID] = indList

    normDF.set_index(taxaID, inplace=True)
    finalDF = normDF.transpose()
    return finalDF, DESeq_error


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
