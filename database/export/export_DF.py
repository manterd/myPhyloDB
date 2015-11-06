import math
import multiprocessing as mp
from numpy import *
import numpy as np
from numpy.random import mtrand
import pandas as pd
from pyper import *
from scipy.spatial.distance import *

from database.models import Kingdom, Phyla, Class, Order, Family, Genus, Species, Sample, Human_Associated, Soil, UserDefined


def UnivMetaDF(idDict):
    sampleTableList = Sample._meta.get_all_field_names()
    human_associatedTableList = Human_Associated._meta.get_all_field_names()
    soilTableList = Soil._meta.get_all_field_names()
    usrTableList = UserDefined._meta.get_all_field_names()

    #idList = []
    #fieldList = []
    #for key in idDict:
    #    fieldList.append(key)
    #    idList.extend(idDict[key])

    metaDF = pd.DataFrame()
    index = 0
    for key in idDict:
        tempDF = pd.DataFrame()
        mySet = idDict[key]

        if key in sampleTableList:
            if key == 'sample_name':
                fields = ['sampleid', key]
            else:
                fields = ['sampleid', 'sample_name', key]
            qs2 = Sample.objects.filter(sampleid__in=mySet).values(*fields)
            tempDF = pd.DataFrame.from_records(qs2, columns=fields).dropna()

        elif key in human_associatedTableList:
            fields = ['sampleid', 'sampleid__sample_name', key]
            qs2 = Human_Associated.objects.filter(sampleid__in=mySet).values(*fields)
            tempDF = pd.DataFrame.from_records(qs2, columns=fields).dropna()
            tempDF.rename(columns={'sampleid__sample_name': 'sample_name'}, inplace=True)

        elif key in soilTableList:
            fields = ['sampleid', 'sampleid__sample_name', key]
            qs2 = Soil.objects.filter(sampleid__in=mySet).values(*fields)
            tempDF = pd.DataFrame.from_records(qs2, columns=fields).dropna()
            tempDF.rename(columns={'sampleid__sample_name': 'sample_name'}, inplace=True)

        elif key in usrTableList:
            fields = ['sampleid', 'sampleid__sample_name', key]
            qs2 = UserDefined.objects.filter(sampleid__in=mySet).values(*fields)
            tempDF = pd.DataFrame.from_records(qs2, columns=fields).dropna()
            tempDF.rename(columns={'sampleid__sample_name': 'sample_name'}, inplace=True)

        tempDF.set_index('sampleid', inplace=True)

        if index == 0:
            metaDF = tempDF
        else:
            cols_to_use = tempDF.columns - metaDF.columns
            metaDF = pd.merge(metaDF, tempDF[cols_to_use], left_index=True, right_index=True)
        index += 1

    return metaDF


def normalizeUniv(df, taxaDict, mySet, meth, reads, metaDF, iters):
    df2 = df.reset_index()
    taxaID = ['kingdomid', 'phylaid', 'classid', 'orderid', 'familyid', 'genusid', 'speciesid']

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
            processes = [mp.Process(target=weightedProb, args=(x, numcore, reads, iters, mySet, df, meth, d)) for x in range(numcore)]

            for p in processes:
                p.start()
            for p in processes:
                p.join()

            for key, value in d.items():
                countDF[key] = value/iters

        elif reads < 0:
            countDF = df2.reset_index(drop=True)

    elif meth == 5:
        countDF = df2[taxaID].reset_index(drop=True)
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
        DF1 = DF1[['sampleid', 'rank', 'taxa_id', 'taxa_name', 'abund', 'rich', 'diversity']]

        if normDF.empty:
            normDF = DF1
        else:
            normDF = normDF.append(DF1)
    return normDF, DESeq_error


def weightedProb(x, cores, reads, iters, mySet, df, meth, d):
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
        for n in range(reads*iters):
            sub = np.random.mtrand.choice(range(sample.size), size=1, replace=False, p=prob)
            temp2 = np.zeros(cols)
            np.put(temp2, sub, 1)
            temp = np.core.umath.add(temp, temp2)
        d[i] = temp

