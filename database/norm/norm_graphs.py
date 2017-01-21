import datetime
from django.http import HttpResponse
from django.db.models import Sum
from django_pandas.io import read_frame
import logging
from numpy import *
import numpy as np
from numpy.random.mtrand import RandomState
import pandas as pd
import pickle
from pyper import *
import simplejson

from database.models import Sample, Air, Human_Associated, Microbial, Soil, Water, UserDefined
from database.models import OTU_97, Profile
from database.utils import taxaProfileDF
import database.queue


curSamples = {}
totSamples = {}
LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getNorm(request, RID, stopList, PID):
    try:
        if request.is_ajax():
            # Get variables from web page
            allJson = request.body.split('&')[0]
            all = simplejson.loads(allJson)
            database.queue.setBase(RID, 'Step 1 of 6: Querying database...')

            NormMeth = int(all["NormMeth"])

            remove = int(all["Remove"])
            cutoff = int(all["Cutoff"])
            Iters = int(all["Iters"])
            Lambda = float(all["Lambda"])
            NormVal = all["NormVal"]
            size_on = int(all["MinSize"])

            # Get selected samples from user's folder and query database for sample info
            myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
            path = str(myDir) + 'usr_sel_samples.pkl'
            with open(path, 'rb') as f:
                selected = pickle.load(f)
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
                    if total['count__sum'] is not None:
                        countList.append(total['count__sum'])
                        subList.append(sample.sampleid)
                qs1 = qs1.filter(sampleid__in=subList)

            if not countList:
                myDict = {}
                myDict['error'] = "Error with Normalization!\nYour minimum sample has caused all samples to be removed!"
                res = simplejson.dumps(myDict)
                return HttpResponse(res, content_type='application/json')

            tabular_on = int(all['Tabular'])
            biom_on = int(all['Biome'])

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
                    if int(total['count__sum']) >= NormReads:
                        id = sample.sampleid
                        newList.append(id)
            else:
                for sample in qs1:
                    total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                    if int(total['count__sum']) > 0:
                        id = sample.sampleid
                        newList.append(id)

            if not newList:
                myDict = {}
                myDict['error'] = "Error with Normalization!\nYour sub-sample size has caused all samples to be removed!"
                res = simplejson.dumps(myDict)
                return HttpResponse(res, content_type='application/json')

            metaDF = UnivMetaDF(newList, RID, stopList, PID)

            # remove emptycols
            metaDF.replace('nan', np.nan, inplace=True)
            metaDF.dropna(axis=1, how='all', inplace=True)
            metaDF.dropna(axis=0, how='all', inplace=True)
            lenB, col = metaDF.shape
            normRem = len(selected) - lenB

            result += str(lenB) + ' selected sample(s) were included in the final analysis...\n'
            result += str(normRem) + ' sample(s) did not met the desired normalization criteria...\n'
            result += '\n'

            # Create unique list of samples in meta dataframe (may be different than selected samples)
            myList = metaDF.index.values.tolist()

            # Create dataframe with all taxa/count data by sample
            taxaDF = taxaProfileDF(myList)

            # Select only the taxa of interest if user used the selectAll button
            taxaDict = {}
            qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('otuid', flat='True').distinct()
            taxaDict['OTU_97'] = qs3

            database.queue.setBase(RID, 'Step 1 of 6: Querying database...done!')

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stopList[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            database.queue.setBase(RID, 'Step 2 of 6: Sub-sampling data...')

            normDF, DESeq_error = normalizeUniv(taxaDF, taxaDict, myList, NormMeth, NormReads, metaDF, Iters, Lambda, RID, stopList, PID)

            database.queue.setBase(RID, 'Step 4 of 6: Calculating indices...')

            normDF['rel_abund'] = normDF['rel_abund'].astype(float)
            normDF['abund'] = normDF['abund'].round(0).astype(int)
            normDF['rich'] = normDF['rich'].round(0).astype(int)
            normDF['diversity'] = normDF['diversity'].astype(float)

            if remove == 1:
                grouped = normDF.groupby('otuid')
                goodIDs = []
                for name, group in grouped:
                    if group['abund'].sum() > cutoff:
                        goodIDs.append(name)
                normDF = normDF.loc[normDF['otuid'].isin(goodIDs)]

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
                result += "Phylotypes with fewer than " + str(cutoff) + " read(s) were removed from your analysis\n"
            else:
                result += "No minimum otu size was applied...\n"

            result += '===============================================\n\n\n'
            finalDict['text'] = result

            normDF.set_index('sampleid', inplace=True)
            finalDF = pd.merge(metaDF, normDF, left_index=True, right_index=True)

            if 'rRNA_copies' in finalDF.columns:
                finalDF['abund_16S'] = finalDF['rel_abund'] * finalDF['rRNA_copies'] / 1000
                finalDF.fillna(0, inplace=True)
            else:
                finalDF['abund_16S'] = 0

            finalDF[['rel_abund', 'rich', 'diversity']] = finalDF[['rel_abund', 'rich', 'diversity']].astype(float)
            finalDF[['abund', 'abund_16S']] = finalDF[['abund', 'abund_16S']].astype(int)

            finalDF.reset_index(drop=False, inplace=True)
            finalDF.rename(columns={'index': 'sampleid'}, inplace=True)

            #re-order finalDF
            metaDFList = list(metaDF.columns.values)
            removeList = ['projectid', 'refid', 'sample_name']
            for x in removeList:
                if x in metaDFList:
                    metaDFList.remove(x)
            metaDFList = ['projectid', 'refid', 'sampleid', 'sample_name'] + metaDFList
            metaDFList = metaDFList + ['kingdomid', 'kingdomName', 'phylaid', 'phylaName', 'classid', 'className', 'orderid', 'orderName', 'familyid', 'familyName', 'genusid', 'genusName', 'speciesid', 'speciesName', 'otuid', 'otuName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']
            finalDF = finalDF[metaDFList]

            # save location info to session and save in temp/norm
            myDir = 'myPhyloDB/media/temp/norm/'
            path = str(myDir) + str(RID) + '.pkl'
            request.session['savedDF'] = pickle.dumps(path)
            request.session['NormMeth'] = NormMeth

            if not os.path.exists(myDir):
                os.makedirs(myDir)
            finalDF.to_pickle(path)

            # save file to users temp/ folder
            myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
            path = str(myDir) + 'usr_norm_data.csv'

            if not os.path.exists(myDir):
                os.makedirs(myDir)
            finalDF.to_csv(path, sep=',')

            database.queue.setBase(RID, 'Step 4 of 6: Calculating indices...done')

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stopList[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            database.queue.setBase(RID, 'Step 5 of 6: Formatting biome data...')

            biome = {}

            nameList = []
            for i in myList:
                nameList.append({"id": str(i), "metadata": metaDF.loc[i].to_dict()})

            # get list of lists with abundances
            taxaOnlyDF = finalDF.loc[:, ['sampleid', 'kingdomName', 'phylaName', 'className', 'orderName', 'familyName', 'genusName', 'speciesName', 'otuName', 'otuid', 'abund']]
            taxaOnlyDF = taxaOnlyDF.pivot(index='otuid', columns='sampleid', values='abund')
            dataList = taxaOnlyDF.values.tolist()

            # get list of taxa
            namesDF = finalDF.loc[:, ['sampleid', 'otuid']]
            namesDF['taxa'] = finalDF.loc[:, ['kingdomName', 'phylaName', 'className', 'orderName', 'familyName', 'genusName', 'speciesName', 'otuName']].values.tolist()
            namesDF = namesDF.pivot(index='otuid', columns='sampleid', values='taxa')

            taxaList = []
            for index, row in namesDF.iterrows():
                metaDict = {}
                metaDict['taxonomy'] = row[0]
                taxaList.append({"id": index, "metadata": metaDict})

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

            myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
            path = str(myDir) + 'usr_norm_data.biom'
            with open(path, 'w') as outfile:
                simplejson.dump(biome, outfile, ensure_ascii=True, indent=4, sort_keys=True)

            database.queue.setBase(RID, 'Step 5 of 6: Formatting biome data...done!')

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stopList[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            database.queue.setBase(RID, 'Step 6 of 6: Formatting result table...')

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

            database.queue.setBase(RID, 'Step 6 of 6: Formatting result table...done!')

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stopList[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            finalDict['error'] = 'none'
            res = simplejson.dumps(finalDict)
            return HttpResponse(res, content_type='application/json')

    except Exception as e:
        if not stopList[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "There was an error during your normalization:\nError: " + str(e.message) + "\nTimestamp: " + str(datetime.datetime.now())
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def UnivMetaDF(sampleList, RID, stopList, PID):
    tableNames = ['project', 'reference', 'sample', 'air', 'human_associated', 'microbial', 'soil', 'water', 'userdefined', 'profile']
    idList = ['projectid', 'refid', u'id']

    your_fields = Sample._meta.local_fields
    sampleTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in sampleTableList:
            sampleTableList.remove(x)

    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
    if stopList[PID] == RID:
        res = ''
        return HttpResponse(res, content_type='application/json')
    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    your_fields = Air._meta.local_fields
    airTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in airTableList:
            airTableList.remove(x)
    for x in idList:
        if x in airTableList:
            airTableList.remove(x)

    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
    if stopList[PID] == RID:
        res = ''
        return HttpResponse(res, content_type='application/json')
    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    your_fields = Human_Associated._meta.local_fields
    humanAssocTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in humanAssocTableList:
            humanAssocTableList.remove(x)
    for x in idList:
        if x in humanAssocTableList:
            humanAssocTableList.remove(x)

    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
    if stopList[PID] == RID:
        res = ''
        return HttpResponse(res, content_type='application/json')
    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    your_fields = Microbial._meta.local_fields
    microbialTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in microbialTableList:
            microbialTableList.remove(x)
    for x in idList:
        if x in microbialTableList:
            microbialTableList.remove(x)

    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
    if stopList[PID] == RID:
        res = ''
        return HttpResponse(res, content_type='application/json')
    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    your_fields = Soil._meta.local_fields
    soilTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in soilTableList:
            soilTableList.remove(x)
    for x in idList:
        if x in soilTableList:
            soilTableList.remove(x)

    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
    if stopList[PID] == RID:
        res = ''
        return HttpResponse(res, content_type='application/json')
    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    your_fields = Water._meta.local_fields
    waterTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in waterTableList:
            waterTableList.remove(x)
    for x in idList:
        if x in waterTableList:
            waterTableList.remove(x)

    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
    if stopList[PID] == RID:
        res = ''
        return HttpResponse(res, content_type='application/json')
    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

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


def normalizeUniv(df, taxaDict, mySet, meth, reads, metaDF, iters, Lambda, RID, stopList, PID):
    global curSamples, totSamples
    df2 = df.reset_index()
    taxaID = ['kingdomid', 'phylaid', 'classid', 'orderid', 'familyid', 'genusid', 'speciesid', 'otuid']

    countDF = pd.DataFrame()
    DESeq_error = 'no'
    if meth == 1 or meth == 4:
        countDF = df2.reset_index(drop=True)

    elif meth == 2:
        if reads >= 0:
            countDF = df2[taxaID].reset_index(drop=True)

            curSamples[RID] = 0
            totSamples[RID] = len(mySet)

            myArr = df2[mySet].as_matrix()
            probArr = myArr.T
            finalArr = rarefaction_remove(probArr, RID, reads=reads, iters=iters)
            for i in xrange(len(mySet)):
                countDF[mySet[i]] = finalArr[i]

            curSamples[RID] = 0

        elif reads < 0:
            countDF = df2.reset_index(drop=True)

    elif meth == 3:
        if reads >= 0:
            countDF = df2[taxaID].reset_index(drop=True)

            curSamples[RID] = 0
            totSamples[RID] = len(mySet)

            myArr = df2[mySet].as_matrix()
            probArr = myArr.T
            finalArr = rarefaction_keep(probArr, RID, reads=reads, iters=iters, myLambda=Lambda)
            for i in xrange(len(mySet)):
                countDF[mySet[i]] = finalArr[i]

            curSamples[RID] = 0

        elif reads < 0:
            countDF = df2.reset_index(drop=True)

    elif meth == 5:
        countDF = df2[taxaID].reset_index(drop=True)
        if os.name == 'nt':
            r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
        else:
            r = R(RCMD="R/R-Linux/bin/R")

        database.queue.setBase(RID, 'Verifying R packages...missing packages are being installed')

        r("list.of.packages <- c('DESeq2')")
        r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
        r("if (length(new.packages)) source('http://bioconductor.org/biocLite.R')")
        print r("if (length(new.packages)) biocLite(new.packages)")

        database.queue.setBase(RID, 'Step 2 of 6: Sub-sampling data...')

        print r("library(DESeq2)")

        df3 = df2.drop(taxaID, axis=1)

        r.assign("count", df3)
        r.assign("metaDF", metaDF)

        r("rows <- rownames(metaDF)")
        r("colnames(count) <- rows")

        r("trt <- factor(metaDF$merge)")

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

    database.queue.setBase(RID, 'Step 2 of 6: Sub-sampling data...done!')
    database.queue.setBase(RID, 'Step 3 of 6: Tabulating data...')

    relabundDF = pd.DataFrame(countDF[taxaID])
    binaryDF = pd.DataFrame(countDF[taxaID])
    diversityDF = pd.DataFrame(countDF[taxaID])
    for i in mySet:
        relabundDF[i] = countDF[i].div(countDF[i].sum(), axis=0)
        binaryDF[i] = countDF[i].apply(lambda x: 1 if x != 0 else 0)
        diversityDF[i] = relabundDF[i].apply(lambda x: -1 * x * math.log(x) if x > 0 else 0)

    field = 'otuid'
    taxaList = taxaDict['OTU_97']

    qs = OTU_97.objects.filter(otuid__in=taxaList)
    namesDF = read_frame(qs, fieldnames=['kingdomid__kingdomName', 'phylaid__phylaName', 'classid__className', 'orderid__orderName', 'familyid__familyName', 'genusid__genusName', 'speciesid__speciesName', 'otuName', 'kingdomid__kingdomid', 'phylaid__phylaid', 'classid__classid', 'orderid__orderid', 'familyid__familyid', 'genusid__genusid', 'speciesid__speciesid', 'otuid'])
    namesDF.rename(columns={'kingdomid__kingdomid': 'kingdomid', 'phylaid__phylaid': 'phylaid', 'classid__classid' : 'classid', 'orderid__orderid' : 'orderid', 'familyid__familyid' : 'familyid', 'genusid__genusid' : 'genusid', 'speciesid__speciesid' : 'speciesid'}, inplace=True)
    namesDF.rename(columns={'kingdomid__kingdomName': 'kingdomName', 'phylaid__phylaName': 'phylaName', 'classid__className' : 'className', 'orderid__orderName' : 'orderName', 'familyid__familyName' : 'familyName', 'genusid__genusName' : 'genusName', 'speciesid__speciesName' : 'speciesName'}, inplace=True)

    curSamples[RID] = 0
    totSamples[RID] = len(mySet)

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
                if stopList[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stopList[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        DF1 = pd.DataFrame(rowsList, columns=['sampleid', 'taxa_id', 'rel_abund', 'abund', 'rich', 'diversity'])

        DF1.rename(columns={'taxa_id': 'otuid'}, inplace=True)
        DF1 = DF1.merge(namesDF, on='otuid', how='outer')

        if normDF.empty:
            normDF = DF1
        else:
            normDF = normDF.append(DF1)

        curSamples[RID] += 1
        database.queue.setBase(RID, 'Step 3 of 6: Tabulating data...\nTabulating is complete for ' + str(curSamples[RID]) + ' out of ' + str(totSamples[RID]) + ' samples')

    curSamples[RID] = 0

    return normDF, DESeq_error


def rarefaction_remove(M, RID, reads=0, iters=0):
    global curSamples, totSamples
    nsamp = M.shape[0]

    Mrarefied = np.empty_like(M)
    for i in range(nsamp):
        counts = M[i]
        nz = counts.nonzero()[0]
        unpacked = np.concatenate([np.repeat(np.array(j,), counts[j]) for j in nz])
        myArr = np.zeros(len(counts), dtype=int)
        for n in xrange(iters):
            permuted = np.random.permutation(unpacked)[:reads]
            binArr = np.zeros(len(counts), dtype=int)
            for p in permuted:
                binArr[p] += 1
            if n == 0:
                myArr = binArr
            else:
                myArr = np.vstack((myArr, binArr))
        Mrarefied[i] = np.mean(myArr, axis=0)
        curSamples[RID] += 1
        database.queue.setBase(RID, 'Step 2 of 6: Sub-sampling data...\nSub-sampling is complete for ' + str(curSamples[RID]) + ' out of ' + str(totSamples[RID]) + ' samples')
    return Mrarefied


def rarefaction_keep(M, RID, reads=0, iters=0, myLambda=0.1):
    global curSamples, totSamples
    noccur = np.sum(M, axis=1)  # number of occurrences for each sample
    nvar = M.shape[1]  # number of variables
    nsamp = M.shape[0]  # number of samples

    Mrarefied = np.empty_like(M)
    for i in range(nsamp):
        p = (M[i] + myLambda) / (float(noccur[i]) + nvar * myLambda)
        myArr = np.zeros(nvar)
        for n in xrange(iters):
            prng = RandomState()
            choice = prng.choice(nvar, size=reads, replace=True, p=p)
            binArr = np.bincount(choice, minlength=nvar)
            if n == 0:
                myArr = binArr
            else:
                myArr = np.vstack((myArr, binArr))

        Mrarefied[i] = np.mean(myArr, axis=0)
        curSamples[RID] += 1
        database.queue.setBase(RID, 'Step 2 of 6: Sub-sampling data...\nSub-sampling is complete for ' + str(curSamples[RID]) + ' out of ' + str(totSamples[RID]) + ' samples')
    return Mrarefied


def getTab(request):
    if request.is_ajax():
        myDict = {}
        myDir = '/myPhyloDB/media/usr_temp/' + str(request.user) + '/'
        fileName = str(myDir) + 'usr_norm_data.csv'  # this file for curNorm
        myDict['name'] = str(fileName)
        res = simplejson.dumps(myDict)
        return HttpResponse(res, content_type='application/json')


def getBiom(request):
    if request.is_ajax():
        myDict = {}
        myDir = '/myPhyloDB/media/usr_temp/' + str(request.user) + '/'
        fileName = str(myDir) + 'usr_norm_data.biom'
        myDict['name'] = str(fileName)
        res = simplejson.dumps(myDict)
        return HttpResponse(res, content_type='application/json')

