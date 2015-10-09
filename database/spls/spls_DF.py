import numpy as np
from numpy.random import mtrand
import pandas as pd
from pyper import *
from scipy.spatial.distance import *
import multiprocessing as mp

from database.models import Sample, Human_Associated, Soil, UserDefined


def SPLSMetaDF(idDict):
    sampleTableList = Sample._meta.get_all_field_names()
    human_associatedTableList = Human_Associated._meta.get_all_field_names()
    soilTableList = Soil._meta.get_all_field_names()
    usrTableList = UserDefined._meta.get_all_field_names()

    idList = []
    fieldList = []
    for key in idDict:
        fieldList.append(key)
        idList.extend(idDict[key])

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


def normalizeSPLS(df, taxaLevel, mySet, meth, reads, metaDF, iters):
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
            processes = [mp.Process(target=weightedProb, args=(x, numcore, reads, iters, mySet, df, meth, d)) for x in range(numcore)]

            for p in processes:
                p.start()
            for p in processes:
                p.join()

            for key, value in d.items():
                normDF[key] = value/iters

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


def weightedProb(x, cores, reads, iters, mySet, df, meth, d):
    high = mySet.__len__()
    set = mySet[x:high:cores]

    for i in set:
        arr = np.asarray(df[i])
        cols = np.shape(arr)
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
