import datetime
from django.http import HttpResponse
from django.db.models import Sum
from django_pandas.io import read_frame
import logging
import multiprocessing as mp
from numpy import *
import numpy as np
from numpy.random import mtrand
import pandas as pd
import pickle
from pyper import *
import simplejson
import threading

from database.models import Sample, Air, Human_Associated, Microbial, Soil, Water, UserDefined
from database.models import Species, Profile
from database.utils import taxaProfileDF, stoppableThread


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}
stop0 = True
stops = {}
thread0 = stoppableThread()
res = ''
LOG_FILENAME = 'error_log.txt'


def statusNorm(request):
    global base, stage, time1, time2, TimeDiff
    if request.is_ajax():
        RID = request.GET["all"]
        time2[RID] = time.time()
        try:
            TimeDiff[RID] = time2[RID] - time1[RID]
        except:
            TimeDiff[RID] = 0
        myDict = {}
        try:
            stage[RID] = str(base[RID]) + '&#10;Analysis has been running for %.1f seconds' % TimeDiff[RID]
        except:
            stage[RID] = '<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
        myDict['stage'] = stage[RID]
        json_data = simplejson.dumps(myDict, encoding="Latin-1")
        return HttpResponse(json_data, content_type='application/json')


def removeRIDNorm(request):
    global base, stage, time1, time2, TimeDiff
    try:
        if request.is_ajax():
            RID = request.GET["all"]
            base.pop(RID, None)
            stage.pop(RID, None)
            time1.pop(RID, None)
            time2.pop(RID, None)
            TimeDiff.pop(RID, None)
            stops.pop(RID, None)
            return True
        else:
            return False
    except:
        return False


def stopNorm(request):
    global thread0, stops, stop0, res
    if request.is_ajax():
        RID = request.GET["all"]
        stops[RID] = True
        stop0 = True
        thread0.terminate()
        thread0.join()
        removeRIDNorm(request)

        res = ''
        myDict = {}
        myDict['error'] = 'none'
        myDict['message'] = 'Your analysis has been stopped!'
        stop = simplejson.dumps(myDict)
        return HttpResponse(stop, content_type='application/json')


def getNormCatData(request):
    global res, thread0, stop0
    if request.is_ajax():
        stop0 = False
        thread0 = stoppableThread(target=loopCat, args=(request,))
        thread0.start()
        thread0.join()
        removeRIDNorm(request)
        return HttpResponse(res, content_type='application/json')


def loopCat(request):
    global base, stage, time1, TimeDiff, res, stops, stop0
    try:
        while True:
            if request.is_ajax():
                # Get variables from web page
                allJson = request.GET["all"]
                all = simplejson.loads(allJson)

                RID = str(all["RID"])
                stops[RID] = False

                time1[RID] = time.time()  # Moved these down here so RID is available
                base[RID] = 'Step 1 of 4: Querying database...'

                NormMeth = int(all["NormMeth"])
                Proc = int(all["Proc"])

                remove = int(all["Remove"])
                cutoff = int(all["Cutoff"])
                Iters = int(all["Iters"])
                NormVal = all["NormVal"]
                size_on = int(all["MinSize"])

                # Get selected samples from cookie and query database for sample info
                samples = Sample.objects.all()
                samples.query = pickle.loads(request.session['selected_samples'])
                selected = samples.values_list('sampleid')
                qs1 = Sample.objects.all().filter(sampleid__in=selected)

                # Generate a list of sequence reads per sample and filter samples if minimum samplesize
                countList = []
                subList = []
                size = 1
                if size_on == 1:
                    size = int(all["MinVal"])
                    for sample in qs1:
                        total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                        if total['count__sum'] >= size:
                            countList.append(total['count__sum'])
                            subList.append(sample.sampleid)
                    qs1 = qs1.filter(sampleid__in=subList)
                else:
                    for sample in qs1:
                        total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                        countList.append(total['count__sum'])

                if not countList:
                    myDict = {}
                    myDict['error'] = "Error with Normalization!\nYour minimum sample has caused all samples to be removed!"
                    res = simplejson.dumps(myDict)
                    return None

                tabular_on = int(all['Tabular'])
                biom_on = int(all['Biome'])

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                # Calculate min/median/max of sequence reads for rarefaction
                NormReads = 0
                if not NormMeth == 1:
                    minSize = int(min(countList))
                    medianSize = int(np.median(np.array(countList)))
                    maxSize = int(max(countList))

                    if NormVal == "min":
                        NormReads = minSize
                    elif NormVal == "median":
                        NormReads = medianSize
                    elif NormVal == "max":
                        NormReads = maxSize
                    elif NormVal == "none":
                        NormReads = -1
                    else:
                        NormReads = int(all["NormVal"])

                result = ''
                result += 'Data Normalization:\n'
                newList = []
                if NormMeth == 2:
                    for sample in qs1:
                        total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                        if total['count__sum'] is not None and int(total['count__sum']) >= NormReads:
                            id = sample.sampleid
                            newList.append(id)
                else:
                    for sample in qs1:
                        total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                        if total['count__sum'] is not None:
                            id = sample.sampleid
                            newList.append(id)

                if not newList:
                    myDict = {}
                    myDict['error'] = "Error with Normalization!\nYour sub-sample size has caused all samples to be removed!"
                    res = simplejson.dumps(myDict)
                    return None

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                metaDF = UnivMetaDF(newList)

                # remove emptycols
                metaDF.replace('nan', np.nan, inplace=True)
                metaDF.dropna(axis=1, how='all', inplace=True)
                metaDF.dropna(axis=0, how='all', inplace=True)
                lenB, col = metaDF.shape

                normRem = len(selected) - lenB

                result += str(lenB) + ' selected sample(s) were included in the final analysis...\n'
                if normRem > 0:
                    result += str(normRem) + ' sample(s) did not met the desired normalization criteria...\n'
                result += '\n'

                # Create unique list of samples in meta dataframe (may be different than selected samples)
                myList = metaDF.index.values.tolist()

                # Create dataframe with all taxa/count data by sample
                taxaDF = taxaProfileDF(myList)

                # Select only the taxa of interest if user used the selectAll button
                taxaDict = {}
                qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('speciesid', flat='True').distinct()
                taxaDict['Species'] = qs3

                base[RID] = 'Step 1 of 4: Querying database...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 2 of 4: Normalizing data...'

                normDF, DESeq_error = normalizeUniv(taxaDF, taxaDict, myList, NormMeth, NormReads, metaDF, Iters, Proc, RID)

                normDF['rel_abund'] = normDF['rel_abund'].astype(float)
                normDF['abund'] = normDF['abund'].round(0).astype(int)
                normDF['rich'] = normDF['rich'].round(0).astype(int)
                normDF['diversity'] = normDF['diversity'].astype(float)

                if remove == 1:
                    grouped = normDF.groupby('speciesid')
                    goodIDs = []
                    for name, group in grouped:
                        if group['abund'].sum() > cutoff:
                            goodIDs.append(name)
                    normDF = normDF.loc[normDF['speciesid'].isin(goodIDs)]

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDict = {}
                if NormMeth == 1:
                    result += 'No normalization was performed...\n'
                elif NormMeth == 2 or NormMeth == 3:
                    result += 'Data were rarefied to ' + str(NormReads) + ' sequence reads with ' + str(Iters) + ' iteration(s)...\n'
                elif NormMeth == 4:
                    result += 'Data were normalized by the total number of sequence reads...\n'
                elif NormMeth == 5 and DESeq_error == 'no':
                    result += 'Data were normalized by DESeq2...\n'
                elif NormMeth == 5 and DESeq_error == 'yes':
                    result += 'DESeq2 cannot run estimateSizeFactors...\n'
                    result += 'Analysis was run without normalization...\n'
                    result += 'To try again, please select fewer samples or another normalization method...\n'

                if size_on == 1:
                    result += "Samples with fewer than " + str(size) + " reads were removed from your analysis...\n"
                else:
                    result += "No minimum samples size was applied...\n"

                if remove == 1:
                    result += "Species with fewer than " + str(cutoff) + " read(s) were removed from your analysis\n"
                else:
                    result += "No minimum species size was applied...\n"

                result += '===============================================\n\n\n'
                finalDict['text'] = result

                normDF.set_index('sampleid', inplace=True)
                finalDF = pd.merge(metaDF, normDF, left_index=True, right_index=True)

                if 'rRNA_copies' in finalDF.columns:
                    finalDF['abund_16S'] = finalDF['rel_abund'] * finalDF['rRNA_copies']
                else:
                    finalDF['abund_16S'] = 0

                if NormMeth == 4:
                    finalDF.drop('abund', axis=1, inplace=True)
                    finalDF.rename(columns={'rel_abund': 'abund'}, inplace=True)
                    finalDF[['abund', 'rich', 'diversity']] = finalDF[['abund', 'rich', 'diversity']].astype(float)
                    finalDF['abund_16S'] = finalDF['abund_16S'].astype(int)

                else:
                    finalDF.drop('rel_abund', axis=1, inplace=True)
                    finalDF[['rich', 'diversity']] = finalDF[['rich', 'diversity']].astype(float)
                    finalDF[['abund', 'abund_16S']] = finalDF[['abund', 'abund_16S']].astype(int)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDF.reset_index(drop=False, inplace=True)
                finalDF.rename(columns={'index': 'sampleid'}, inplace=True)

                #re-order finalDF
                metaDFList = list(metaDF.columns.values)
                removeList = ['projectid', 'refid', 'sample_name']
                for x in removeList:
                    if x in metaDFList:
                        metaDFList.remove(x)
                metaDFList = ['projectid', 'refid', 'sampleid', 'sample_name'] + metaDFList
                metaDFList = metaDFList + ['kingdomid', 'kingdomName', 'phylaid', 'phylaName', 'classid', 'className', 'orderid', 'orderName', 'familyid', 'familyName', 'genusid', 'genusName', 'speciesid', 'speciesName', 'abund', 'abund_16S', 'rich', 'diversity']
                finalDF = finalDF[metaDFList]

                # save location info to session
                myDir = 'database/norm/saved/'
                name = request.user
                ip = request.META.get('REMOTE_ADDR')
                path = str(myDir) + 'savedDF_' + str(name) + '_' + str(ip) + '.pkl'
                request.session['savedDF'] = pickle.dumps(path)
                request.session['NormMeth'] = NormMeth

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 2 of 4: Normalizing data...done!'
                base[RID] = 'Step 3 of 4: Formatting biome data...'

                biome = {}

                nameList = []
                for i in myList:
                    nameList.append({"id": str(i), "metadata": metaDF.loc[i].to_dict()})

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[RID]:
                        res = ''
                        return None
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                # get list of lists with abundances
                taxaOnlyDF = finalDF.loc[:, ['sampleid', 'kingdomName', 'phylaName', 'className', 'orderName', 'familyName', 'genusName', 'speciesName', 'speciesid', 'abund']]
                taxaOnlyDF = taxaOnlyDF.pivot(index='speciesid', columns='sampleid', values='abund')
                dataList = taxaOnlyDF.values.tolist()

                # get list of taxa
                namesDF = finalDF.loc[:, ['sampleid', 'speciesid']]
                namesDF['taxa'] = finalDF.loc[:, ['kingdomName', 'phylaName', 'className', 'orderName', 'familyName', 'genusName', 'speciesName']].values.tolist()
                namesDF = namesDF.pivot(index='speciesid', columns='sampleid', values='taxa')

                taxaList = []
                for index, row in namesDF.iterrows():
                    metaDict = {}
                    metaDict['taxonomy'] = row[0]
                    taxaList.append({"id": index, "metadata": metaDict})

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[RID]:
                        res = ''
                        return None
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                biome['format'] = 'Biological Observation Matrix 0.9.1-dev'
                biome['format_url'] = 'http://biom-format.org/documentation/format_versions/biom-1.0.html'
                biome['type'] = 'OTU table'
                biome['generated_by'] = 'myPhyloDB'
                biome['date'] = str(datetime.datetime.now())
                biome['matrix_type'] = 'dense'
                biome['matrix_element_type'] = 'float'
                biome['rows'] = taxaList
                biome['columns'] = nameList
                biome['data'] = dataList

                if biom_on == 1:
                    biome_json = simplejson.dumps(biome, ensure_ascii=True, indent=4, sort_keys=True)
                    finalDict['biome'] = str(biome_json)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 3 of 4: Formatting biome data...done!'
                base[RID] = 'Step 4 of 4: Formatting result table...'

                if tabular_on == 1:
                    data = finalDF.values.tolist()
                    cols = finalDF.columns.values.tolist()
                    colList = []
                    for item in cols:
                        colDict = {}
                        colDict['title'] = item
                        colList.append(colDict)
                    finalDict['data'] = data
                    finalDict['columns'] = colList

                #base[RID] = 'Step 4 of 5: Formatting result table...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                return None

    except:
        if not stop0:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with Normalization!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
        return None


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


def normalizeUniv(df, taxaDict, mySet, meth, reads, metaDF, iters, Proc, RID):
    global res, stops
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
                processes = [threading.Thread(target=weightedProb, args=(x, numcore, reads, iters, mySet, df, meth, d, RID,)) for x in range(numcore)]
            else:
                numcore = min(mp.cpu_count(), Proc)
                processes = [mp.Process(target=weightedProb, args=(x, numcore, reads, iters, mySet, df, meth, d, RID,)) for x in range(numcore)]

            for p in processes:
                p.start()

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            for p in processes:
                p.join()

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            for key, value in d.items():
                countDF[key] = value/iters

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

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

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[RID]:
            res = ''
            return None
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    field = 'speciesid'
    taxaList = taxaDict['Species']

    qs = Species.objects.filter(speciesid__in=taxaList)
    namesDF = read_frame(qs, fieldnames=['kingdomid__kingdomName', 'phylaid__phylaName', 'classid__className', 'orderid__orderName', 'familyid__familyName', 'genusid__genusName', 'speciesName', 'kingdomid__kingdomid', 'phylaid__phylaid', 'classid__classid', 'orderid__orderid', 'familyid__familyid', 'genusid__genusid', 'speciesid'])
    namesDF.rename(columns={'kingdomid__kingdomid': 'kingdomid', 'phylaid__phylaid': 'phylaid', 'classid__classid' : 'classid', 'orderid__orderid' : 'orderid', 'familyid__familyid' : 'familyid', 'genusid__genusid' : 'genusid'}, inplace=True)
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

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[RID]:
            res = ''
            return None
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        DF1 = pd.DataFrame(rowsList, columns=['sampleid', 'taxa_id', 'rel_abund', 'abund', 'rich', 'diversity'])

        DF1.rename(columns={'taxa_id': 'speciesid'}, inplace=True)
        DF1 = DF1.merge(namesDF, on='speciesid', how='outer')

        if normDF.empty:
            normDF = DF1
        else:
            normDF = normDF.append(DF1)

    return normDF, DESeq_error


def weightedProb(x, cores, reads, iters, mySet, df, meth, d, RID):
    global res, stops
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

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[RID]:
            res = ''
            return None
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        temp = np.zeros((otus,), dtype=np.int)
        for n in range(reads*iters):
            if meth == 3:
                sub = np.random.mtrand.choice(range(sample.size), size=1, replace=False, p=prob)
            else:
                sub = np.random.mtrand.choice(range(sample.size), size=1, replace=True, p=prob)
            temp2 = np.zeros((otus,), dtype=np.int)
            np.put(temp2, sub, [1])
            temp = np.core.umath.add(temp, temp2)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[RID]:
                res = ''
                return None
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        d[i] = temp



