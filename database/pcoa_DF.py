import operator
import pandas as pd
from django.db.models import Q
import numpy as np
from numpy import *
from numpy.random import mtrand
from scipy.spatial.distance import *
import multiprocessing as mp
import math
import os
from pyper import *

def catPCoAMetaDF(qs1, metaDict):
    sampleTableList = ['sample_name', 'organism', 'seq_method', 'collection_date', 'biome', 'feature', 'geo_loc_country', 'geo_loc_state', 'geo_loc_city', 'geo_loc_farm', 'geo_loc_plot', 'material']
    collectTableList = ['depth', 'pool_dna_extracts', 'samp_collection_device', 'sieving', 'storage_cond']
    soil_classTableList = ['drainage_class', 'fao_class', 'horizon', 'local_class', 'profile_position', 'slope_aspect', 'soil_type', 'texture_class']
    mgtTableList = ['agrochem_amendments', 'agrochem_amendments_desc', 'biological_amendments', 'biological_amendments_desc', 'cover_crop', 'crop_rotation', 'cur_land_use', 'cur_vegetation', 'cur_crop', 'cur_cultivar', 'organic', 'previous_land_use', 'soil_amendments', 'soil_amendments_desc', 'tillage']
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
            qs2 = qs1.values(*field_list).filter(reduce(operator.or_, args_list)).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key in collectTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'collect__' + str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.values(*field_list).filter(reduce(operator.or_, args_list)).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key in soil_classTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'soil_class__' + str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.values(*field_list).filter(reduce(operator.or_, args_list)).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key in mgtTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'management__' + str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.values(*field_list).filter(reduce(operator.or_, args_list)).exclude(reduce(operator.or_, exclude_list))
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
            qs2 = qs1.values(*field_list).filter(reduce(operator.or_, args_list)).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
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
            exclude_list = []
            exclude_list.append(Q(**{value: 'null'}))
            qs2 = qs1.values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            if metaDF.empty:
                metaDF = tempDF
                metaDF[value] = metaDF[value].astype(float)
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')
                metaDF[value] = metaDF[value].astype(float)


        elif key == 'collect':
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'collect__' + str(value)
            field_list.append(field)
            final_fieldList.append(value)
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: value}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
                metaDF[value] = metaDF[value].astype(float)
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')
                metaDF[value] = metaDF[value].astype(float)

        elif key == 'climate':
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'climate__' + str(value)
            field_list.append(field)
            final_fieldList.append(value)
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: value}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
                metaDF[value] = metaDF[value].astype(float)
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')
                metaDF[value] = metaDF[value].astype(float)

        elif key == 'soil_class':
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'soil_class__' + str(value)
            field_list.append(field)
            final_fieldList.append(value)
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: value}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
                metaDF[value] = metaDF[value].astype(float)
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')
                metaDF[value] = metaDF[value].astype(float)

        elif key == 'soil_nutrient':
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'soil_nutrient__' + str(value)
            field_list.append(field)
            final_fieldList.append(value)
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: value}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
                metaDF[value] = metaDF[value].astype(float)
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')
                metaDF[value] = metaDF[value].astype(float)

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
                metaDF[value] = metaDF[value].astype(float)
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')
                metaDF[value] = metaDF[value].astype(float)

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
        # Create CountDataSet for DESeq
        r = R(RCMD="R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
        df3 = df2.drop(taxaID, axis=1)
        r.assign("countTable", df3)
        r.assign("metaDF", metaDF)
        r("condition <- factor(metaDF$merge)")
        r("library(DESeq)")
        r("cds <- newCountDataSet(countTable, condition)")
        r("cds <- estimateSizeFactors(cds)")
        pycds = r.get("sizeFactors(cds)")

        found = 0
        for thing in pycds:
            if str(thing) == "None":
                found += 1

        if found == 0:
            DESeq_error = 'no'
            colList = df3.columns.tolist()
            indList = df2[taxaID].tolist()
            normDF = pd.DataFrame(r.get("counts(cds, normalize=TRUE)"), columns=[colList])
            normDF[taxaID] = indList
        else:
            DESeq_error = 'yes'
            sizeFactor = []
            for i in mySet:
                sizeFactor.append(1)
            r.assign("sizeFactor", sizeFactor)
            r("cds$sizeFactor <- sizeFactor")
            colList = df3.columns.tolist()
            indList = df2[taxaID].tolist()
            normDF = pd.DataFrame(r.get("counts(cds, normalize=TRUE)"), columns=[colList])
            normDF[taxaID] = indList

    elif meth == 6:
        # Create CountDataSet for DESeq
        r = R(RCMD="R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
        df3 = df2.drop(taxaID, axis=1)
        r.assign("countTable", df3)
        r.assign("metaDF", metaDF)
        r("condition <- factor(metaDF$merge)")
        r("library(DESeq)")
        r("cds <- newCountDataSet(countTable, condition)")
        r("cds <- estimateSizeFactors(cds)")
        pycds = r.get("sizeFactors(cds)")

        found = 0
        for thing in pycds:
            if str(thing) == "None":
                found += 1

        if found == 0:
            DESeq_error = 'no'
            r("cds <- estimateDispersions(cds, method='blind', fitType='local')")
            r("vsd <- getVarianceStabilizedData(cds)")
            colList = df3.columns.tolist()
            indList = df2[taxaID].tolist()
            normDF = pd.DataFrame(r.get("vsd"), columns=[colList])
            normDF[taxaID] = indList
        else:
            DESeq_error = 'yes'
            sizeFactor = []
            for i in mySet:
                sizeFactor.append(1)
            r.assign("sizeFactor", sizeFactor)
            r("cds$sizeFactor <- sizeFactor")
            r("cds <- estimateDispersions(cds, method='blind', fitType='local')")
            r("vsd <- getVarianceStabilizedData(cds)")
            colList = df3.columns.tolist()
            indList = df2[taxaID].tolist()
            normDF = pd.DataFrame(r.get("vsd"), columns=[colList])
            normDF[taxaID] = indList

    normDF.set_index(taxaID, inplace=True)
    finalDF = normDF.transpose()
    return finalDF, DESeq_error


def weightedProb(x, cores, reads, mySet, df, meth, d):
    high = len(mySet)
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
