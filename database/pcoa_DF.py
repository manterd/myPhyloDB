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


def normalizePCoA(df, taxaLevel, mySet, reads, depvar):
    df2 = df.reset_index()
    taxaID = ['kingdomid', 'phylaid', 'classid', 'orderid', 'familyid', 'genusid', 'speciesid']

    field = ''
    rank = ''
    if taxaLevel == 1:
        rank = 'Kingdom'
        field = 'kingdomid'
    elif taxaLevel == 2:
        rank = 'Phyla'
        field = 'phylaid'
    elif taxaLevel == 3:
        rank = 'Class'
        field = 'classid'
    elif taxaLevel == 4:
        rank = 'Order'
        field = 'orderid'
    elif taxaLevel == 5:
        rank = 'Family'
        field = 'familyid'
    elif taxaLevel == 6:
        rank = 'Genus'
        field = 'genusid'
    elif taxaLevel == 7:
        rank = 'Species'
        field = 'speciesid'

    countDF = pd.DataFrame()
    if reads >= 0:
        countDF = df2[taxaID].reset_index(drop=True)
        manager = mp.Manager()
        d = manager.dict()

        if os.name == 'nt':
            numcore = 1
        else:
            numcore = mp.cpu_count()-1 or 1

        if depvar == 2:
            processes = [mp.Process(target=weightedProb1, args=(x, numcore, reads, mySet, df, d)) for x in range(numcore)]
        else:
            processes = [mp.Process(target=weightedProb2, args=(x, numcore, reads, mySet, df, d)) for x in range(numcore)]

        for p in processes:
            p.start()
        for p in processes:
            p.join()

        for key, value in d.items():
            countDF[key] = value

    elif reads < 0:
        countDF = df2.reset_index(drop=True)

    relabundDF = pd.DataFrame(countDF[taxaID])
    binaryDF = pd.DataFrame(countDF[taxaID])
    diversityDF = pd.DataFrame(countDF[taxaID])
    for i in mySet:
        relabundDF[i] = countDF[i].div(countDF[i].sum(), axis=0)
        binaryDF[i] = countDF[i].apply(lambda x: 1 if x != 0 else 0)
        diversityDF[i] = relabundDF[i].apply(lambda x: -1 * x * math.log(x) if x > 0 else 0)

    normDF = pd.DataFrame(columns=['sampleid', 'rank', 'count', 'rel_abund', 'rich', 'diversity'])
    for i in mySet:
        tmpDF = pd.DataFrame()
        tmpDF['count'] = countDF.groupby(field)[i].sum()
        tmpDF['rel_abund'] = relabundDF.groupby(field)[i].sum()
        tmpDF['rich'] = binaryDF.groupby(field)[i].sum()
        tmpDF['diversity'] = diversityDF.groupby(field)[i].sum()
        tmpDF['sampleid'] = i
        tmpDF['rank'] = rank
        tmpDF.reset_index(inplace=True)
        normDF = normDF.append(tmpDF)
    normDF.rename(columns={field: 'taxaid'}, inplace=True)

    return normDF


def weightedProb1(x, cores, reads, mySet, df, d):
    high = len(mySet)
    set = mySet[x:high:cores]

    for i in set:
        arr = asarray(df[i])
        cols = shape(arr)
        sample = arr.astype(dtype=np.float64)
        if sample.sum() <= reads:
            prob = (sample + 0.1) / (sample.sum() + cols[0] * 0.1)
        else:
            prob = sample / sample.sum()
        d[i] = prob


def weightedProb2(x, cores, reads, mySet, df, d):
    high = len(mySet)
    set = mySet[x:high:cores]

    for i in set:
        arr = asarray(df[i])
        cols = shape(arr)
        sample = arr.astype(dtype=np.float64)

        if sample.sum() <= reads:
            prob = (sample + 0.1) / (sample.sum() + cols[0] * 0.1)
        else:
            prob = sample / sample.sum()

        temp = np.zeros(cols)
        for n in range(reads):
            sub = np.random.mtrand.choice(range(sample.size), size=1, replace=True, p=prob)
            temp2 = np.zeros(cols)
            np.put(temp2, sub, 1)
            temp = np.core.umath.add(temp, temp2)
        d[i] = temp
