import math
import multiprocessing as mp
from numpy import *
import numpy as np
from numpy.random import mtrand
import pandas as pd
from pyper import *
from scipy.spatial.distance import *
import threading

from database.models import Sample, Air, Human_Associated, Microbial, Soil, Water, UserDefined
from database.models import Species


def UnivMetaDF(sampleList):
    tableNames = ['project', 'reference', 'sample', 'air', 'human_associated', 'microbial', 'soil', 'water', 'userdefined', 'profile']
    idList = ['projectid', 'refid', u'id']

    your_fields = Sample._meta.local_fields
    sampleTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in sampleTableList:
            sampleTableList.remove(x)

    your_fields = Air._meta.local_fields
    airTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in airTableList:
            airTableList.remove(x)
    for x in idList:
        if x in airTableList:
            airTableList.remove(x)

    your_fields = Human_Associated._meta.local_fields
    humanAssocTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in humanAssocTableList:
            humanAssocTableList.remove(x)
    for x in idList:
        if x in humanAssocTableList:
            humanAssocTableList.remove(x)

    your_fields = Microbial._meta.local_fields
    microbialTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in microbialTableList:
            microbialTableList.remove(x)
    for x in idList:
        if x in microbialTableList:
            microbialTableList.remove(x)

    your_fields = Soil._meta.local_fields
    soilTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in soilTableList:
            soilTableList.remove(x)
    for x in idList:
        if x in soilTableList:
            soilTableList.remove(x)

    your_fields = Water._meta.local_fields
    waterTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in waterTableList:
            waterTableList.remove(x)
    for x in idList:
        if x in waterTableList:
            waterTableList.remove(x)

    your_fields = UserDefined._meta.local_fields
    usrTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in usrTableList:
            usrTableList.remove(x)
    for x in idList:
        if x in usrTableList:
            usrTableList.remove(x)

    metaDF = pd.DataFrame(list(Sample.objects.all().filter(sampleid__in=sampleList).values(*sampleTableList)))
    metaDF.set_index('sampleid', drop=True, inplace=True)

    tempDF = pd.DataFrame(list(Air.objects.all().filter(sampleid_id__in=sampleList).values(*airTableList)))
    if not tempDF.empty:
        tempDF.set_index('sampleid', drop=True, inplace=True)
        metaDF = pd.merge(metaDF, tempDF, left_index=True, right_index=True)

    tempDF = pd.DataFrame(list(Human_Associated.objects.all().filter(sampleid_id__in=sampleList).values(*humanAssocTableList)))
    if not tempDF.empty:
        tempDF.set_index('sampleid', drop=True, inplace=True)
        metaDF = pd.merge(metaDF, tempDF, left_index=True, right_index=True)

    tempDF = pd.DataFrame(list(Microbial.objects.all().filter(sampleid_id__in=sampleList).values(*microbialTableList)))
    if not tempDF.empty:
        tempDF.set_index('sampleid', drop=True, inplace=True)
        metaDF = pd.merge(metaDF, tempDF, left_index=True, right_index=True)

    tempDF = pd.DataFrame(list(Soil.objects.all().filter(sampleid_id__in=sampleList).values(*soilTableList)))
    if not tempDF.empty:
        tempDF.set_index('sampleid', drop=True, inplace=True)
        metaDF = pd.merge(metaDF, tempDF, left_index=True, right_index=True)

    tempDF = pd.DataFrame(list(Water.objects.all().filter(sampleid_id__in=sampleList).values(*waterTableList)))
    if not tempDF.empty:
        tempDF.set_index('sampleid', drop=True, inplace=True)
        metaDF = pd.merge(metaDF, tempDF, left_index=True, right_index=True)

    tempDF = pd.DataFrame(list(UserDefined.objects.all().filter(sampleid_id__in=sampleList).values(*usrTableList)))
    if not tempDF.empty:
        tempDF.set_index('sampleid', drop=True, inplace=True)
        metaDF = pd.merge(metaDF, tempDF, left_index=True, right_index=True)

    return metaDF


def normalizeUniv(df, taxaDict, mySet, meth, reads, metaDF, iters, Proc):
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

            if Proc < 1:
                Proc = 1

            if os.name == 'nt':
                numcore = 1
                processes = [threading.Thread(target=weightedProb, args=(x, numcore, reads, iters, mySet, df, meth, d,)) for x in range(numcore)]
            else:
                numcore = min(mp.cpu_count(), Proc)
                processes = [mp.Process(target=weightedProb, args=(x, numcore, reads, iters, mySet, df, meth, d,)) for x in range(numcore)]

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

        #R adds Xs to colnames that begin with a number, need to replace
        r("rows <- rownames(metaDF)")
        r("colnames(count) <- rows")

        r("trt <- factor(metaDF$merge)")

        r("library(DESeq2)")
        r("dds <- DESeqDataSetFromMatrix(countData = count, colData = metaDF, design = ~ sample_name)")

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

    field = 'speciesid'
    taxaList = taxaDict['Species']

    qs1 = Species.objects.filter(speciesid__in=taxaList).values('kingdomid__kingdomName', 'phylaid__phylaName', 'classid__className', 'orderid__orderName', 'familyid__familyName', 'genusid__genusName', 'speciesName', 'kingdomid', 'phylaid', 'classid', 'orderid', 'familyid', 'genusid', 'speciesid')
    namesDF = pd.DataFrame.from_records(qs1, columns=['kingdomid__kingdomName', 'phylaid__phylaName', 'classid__className', 'orderid__orderName', 'familyid__familyName', 'genusid__genusName', 'speciesName', 'kingdomid', 'phylaid', 'classid', 'orderid', 'familyid', 'genusid', 'speciesid'])
    namesDF.rename(columns={'kingdomid__kingdomName': 'kingdomName', 'phylaid__phylaName': 'phylaName', 'classid__className' : 'className', 'orderid__orderName' : 'orderName', 'familyid__familyName' : 'familyName', 'genusid__genusName' : 'genusName'}, inplace=True)

    normDF = pd.DataFrame()
    for item in mySet:
        rowsList = []
        groupRelAbund = relabundDF.groupby(field)[item].sum()
        groupAbund = countDF.groupby(field)[item].sum()
        groupRich = binaryDF.groupby(field)[item].sum()
        groupDiversity = diversityDF.groupby(field)[item].sum()

        if isinstance(taxaList, unicode):
            myDict = {}
            myDict['sampleid'] = item
            myDict['taxa_id'] = taxaList
            myDict['rel_abund'] = groupRelAbund[taxaList]
            myDict['abund'] = groupAbund[taxaList]
            myDict['rich'] = groupRich[taxaList]
            myDict['diversity'] = groupDiversity[taxaList]
            rowsList.append(myDict)
        else:
            for taxa in taxaList:
                myDict = {}
                myDict['sampleid'] = item
                myDict['taxa_id'] = taxa
                myDict['rel_abund'] = groupRelAbund[taxa]
                myDict['abund'] = groupAbund[taxa]
                myDict['rich'] = groupRich[taxa]
                myDict['diversity'] = groupDiversity[taxa]
                rowsList.append(myDict)

        DF1 = pd.DataFrame(rowsList, columns=['sampleid', 'taxa_id', 'rel_abund', 'abund', 'rich', 'diversity'])

        DF1.rename(columns={'taxa_id': 'speciesid'}, inplace=True)
        DF1 = DF1.merge(namesDF, on='speciesid', how='outer')

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
        cols = np.ma.shape(arr)
        otus = cols[0]
        sample = arr.astype(dtype=np.float64)
        if meth == 3:
            prob = (sample + 0.1) / (sample.sum() + otus * 0.1)
        else:
            prob = sample / sample.sum()

        temp = np.zeros((otus,), dtype=np.int)
        for n in range(reads*iters):
            if meth == 3:
                sub = np.random.mtrand.choice(range(sample.size), size=1, replace=False, p=prob)
            else:
                sub = np.random.mtrand.choice(range(sample.size), size=1, replace=True, p=prob)
            temp2 = np.zeros((otus,), dtype=np.int)
            np.put(temp2, sub, [1])
            temp = np.core.umath.add(temp, temp2)
        d[i] = temp

