import datetime
from django.http import HttpResponse
from django.db.models import Sum
import logging
import numpy as np
import pandas as pd
import pickle
from pyper import *
import simplejson

from database.norm.norm_DF import UnivMetaDF, normalizeUniv
from database.models import Sample, Profile
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
            # Get selected samples from cookie and query database for sample info
            samples = Sample.objects.all()
            samples.query = pickle.loads(request.session['selected_samples'])
            selected = samples.values_list('sampleid')
            qs1 = Sample.objects.all().filter(sampleid__in=selected)

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
                mtColsDF = metaDF.replace(to_replace='None', value=np.nan)
                mtColsDF.replace(to_replace='nan', value=np.nan, inplace=True)
                mtColsDF.dropna(axis=1, how='all', inplace=True)
                mtColsList = mtColsDF.columns.values.tolist()
                metaDF = metaDF[mtColsList]

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

                normDF, DESeq_error = normalizeUniv(taxaDF, taxaDict, myList, NormMeth, NormReads, metaDF, Iters, Proc)

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

                finalDF['abund_16S'] = finalDF['rel_abund'] * finalDF['rRNA_copies']

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
                    res_table = finalDF.to_html(classes="table display")
                    res_table = res_table.replace('border="1"', 'border="0"')
                    finalDict['res_table'] = str(res_table)

                base[RID] = 'Step 4 of 4: Formatting result table...done!'

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




