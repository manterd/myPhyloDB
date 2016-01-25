import datetime
from django.http import HttpResponse
import logging
import pandas as pd
import pickle
import simplejson
import numpy as np
from pyper import *

from database.utils import multidict, stoppableThread
from database.models import Kingdom, Phyla, Class, Order, Family, Genus, Species


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}
stop2 = False
stops = {}
thread2 = stoppableThread()
res = ''
LOG_FILENAME = 'error_log.txt'


def updateDiffAbund(request):
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
            stage[RID] = str(base[RID]) + '<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
        except:
            stage[RID] = '<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
        myDict['stage'] = stage[RID]
        json_data = simplejson.dumps(myDict, encoding="Latin-1")
        return HttpResponse(json_data, content_type='application/json')


def removeRIDDIFF(request):
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


def stopDiffAbund(request):
    global thread2, stops, stop2, res
    if request.is_ajax():
        RID = request.GET["all"]
        stops[RID] = True
        stop2 = True
        thread2.terminate()
        thread2.join()
        removeRIDDIFF(request)

        res = ''
        myDict = {}
        myDict['error'] = 'none'
        myDict['message'] = 'Your analysis has been stopped!'
        stop = simplejson.dumps(myDict)
        return HttpResponse(stop, content_type='application/json')


def getDiffAbund(request):
    global res, thread2, stop2
    if request.is_ajax():
        stop2 = False
        thread2 = stoppableThread(target=loopCat, args=(request,))
        thread2.start()
        thread2.join()
        removeRIDDIFF(request)
        return HttpResponse(res, content_type='application/json')


def loopCat(request):
    global res, base, stage, time1, TimeDiff, stops, stop2
    try:
        while True:
            if request.is_ajax():
                # Get variables from web page
                allJson = request.GET["all"]
                all = simplejson.loads(allJson)

                RID = str(all["RID"])
                stops[RID] = False

                time1[RID] = time.time()  # Moved these down here so RID is available
                base[RID] = 'Step 1 of 5: Selecting your chosen meta-variables...'

                path = pickle.loads(request.session['savedDF'])
                savedDF = pd.read_pickle(path)

                selectAll = int(all["selectAll"])
                DepVar = int(all["DepVar"])

                result = ''
                if selectAll == 1:
                    result += 'Taxa level: Kingdom' + '\n'
                elif selectAll == 2:
                    result += 'Taxa level: Phyla' + '\n'
                elif selectAll == 3:
                    result += 'Taxa level: Class' + '\n'
                elif selectAll == 4:
                    result += 'Taxa level: Order' + '\n'
                elif selectAll == 5:
                    result += 'Taxa level: Family' + '\n'
                elif selectAll == 6:
                    result += 'Taxa level: Genus' + '\n'
                elif selectAll == 7:
                    result += 'Taxa level: Species' + '\n'

                # Select samples and meta-variables from savedDF
                metaValsCat = all['metaValsCat']
                metaIDsCat = all['metaIDsCat']

                metaDictCat = {}
                catFields = []
                catValues = []
                if metaValsCat:
                    metaDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaValsCat)
                    for key in sorted(metaDictCat):
                        catFields.append(key)
                        catValues.extend(metaDictCat[key])

                catFields_edit = []
                removed = []
                for i in metaDictCat:
                    levels = len(set(metaDictCat[i]))
                    if levels > 1:
                        catFields_edit.append(i)
                    else:
                        removed.append(i)

                if not catFields_edit:
                    myDict = {}
                    myDict['error'] = "Selected meta data only has one level.\nPlease select a different variable(s)."
                    res = simplejson.dumps(myDict)
                    return None

                catSampleIDs = []
                if metaIDsCat:
                    idDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaIDsCat)
                    for key in sorted(idDictCat):
                        catSampleIDs.extend(idDictCat[key])

                # Removes samples (rows) that are not in our samplelist
                metaDF = savedDF.loc[savedDF['sampleid'].isin(catSampleIDs)]

                if metaDictCat:
                    for key in metaDictCat:
                        metaDF = metaDF.loc[metaDF[key].isin(metaDictCat[key])]

                # Create combined metadata column - DiffAbund only
                if len(catFields_edit) > 1:
                    for index, row in metaDF.iterrows():
                        metaDF.ix[index, 'merge'] = "; ".join(row[catFields_edit])
                else:
                    metaDF['merge'] = metaDF[catFields_edit[0]]

                # command is unique to DiffAbund
                metaDF = metaDF.loc[:, ['merge']]

                result += 'Categorical variables selected by user: ' + ", ".join(catFields) + '\n'
                result += 'Categorical variables removed (contains only 1 level): ' + ", ".join(removed) + '\n'
                result += '===============================================\n'

                base[RID] = 'Step 1 of 5: Selecting your chosen meta-variables...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 2 of 5: Selecting your chosen taxa...'

                # get selected taxa fro each rank selected in the tree
                taxaDF = pd.DataFrame(columns=['sampleid', 'rank', 'taxa_id', 'taxa_name', 'abund', 'abund_16S', 'rich', 'diversity'])

                if selectAll == 1:
                    taxaDF = savedDF.loc[:, ['sampleid', 'kingdomid', 'kingdomName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'kingdomid': 'taxa_id', 'kingdomName': 'taxa_name'}, inplace=True)
                    taxaDF.loc[:, 'rank'] = 'Kingdom'
                elif selectAll == 2:
                    taxaDF = savedDF.loc[:, ['sampleid', 'phylaid', 'phylaName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'phylaid': 'taxa_id', 'phylaName': 'taxa_name'}, inplace=True)
                    taxaDF.loc[:, 'rank'] = 'Phyla'
                elif selectAll == 3:
                    taxaDF = savedDF.loc[:, ['sampleid', 'classid', 'className', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'classid': 'taxa_id', 'className': 'taxa_name'}, inplace=True)
                    taxaDF.loc[:, 'rank'] = 'Class'
                elif selectAll == 4:
                    taxaDF = savedDF.loc[:, ['sampleid', 'orderid', 'orderName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'orderid': 'taxa_id', 'orderName': 'taxa_name'}, inplace=True)
                    taxaDF.loc[:, 'rank'] = 'Order'
                elif selectAll == 5:
                    taxaDF = savedDF.loc[:, ['sampleid', 'familyid', 'familyName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'familyid': 'taxa_id', 'familyName': 'taxa_name'}, inplace=True)
                    taxaDF.loc[:, 'rank'] = 'Family'
                elif selectAll == 6:
                    taxaDF = savedDF.loc[:, ['sampleid', 'genusid', 'genusName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'genusid': 'taxa_id', 'genusName': 'taxa_name'}, inplace=True)
                    taxaDF.loc[:, 'rank'] = 'Genus'
                elif selectAll == 7:
                    taxaDF = savedDF.loc[:, ['sampleid', 'speciesid', 'speciesName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'speciesid': 'taxa_id', 'speciesName': 'taxa_name'}, inplace=True)
                    taxaDF.loc[:, 'rank'] = 'Species'

                finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')
                
                wantedList = ['sampleid', 'merge', 'rank', 'taxa_name', 'taxa_id']
                finalDF = finalDF.groupby(wantedList)[['abund', 'abund_16S']].sum()

                finalDF.reset_index(drop=False, inplace=True)

                taxaSums = finalDF.groupby('taxa_id').sum()
                goodList = taxaSums[taxaSums['abund'] > 0].index.values.tolist()
                finalDF = finalDF.loc[finalDF['taxa_id'].isin(goodList)]

                count_rDF = pd.DataFrame()
                if DepVar == 1:
                    count_rDF = finalDF.pivot(index='taxa_id', columns='sampleid', values='abund')
                elif DepVar == 4:
                    count_rDF = finalDF.pivot(index='taxa_id', columns='sampleid', values='abund_16S')

                meta_rDF = finalDF.drop_duplicates(subset='sampleid', keep='last')
                wantedList = ['sampleid', 'merge']
                meta_rDF = meta_rDF[wantedList]
                meta_rDF.set_index('sampleid', drop=True, inplace=True)

                base[RID] = 'Step 2 of 5: Selecting your chosen taxa...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 3 of 5: Performing statistical test...'

                finalDict = {}
                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                r.assign("metaDF", meta_rDF)
                r("trt <- factor(metaDF$merge)")

                r.assign("count", count_rDF)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                DESeq_error = ''
                r("library(DESeq2)")
                r("colData <- data.frame(row.names=colnames(count), trt=trt)")
                r("dds <- DESeqDataSetFromMatrix(countData=count, colData=colData, design= ~ trt)")

                r("dds <- estimateSizeFactors(dds)")
                pycds = r.get("sizeFactors(dds)")

                if pycds is not None:
                    DESeq_error = 'no'
                    r("dds <- estimateDispersions(dds)")
                    r("dds <- nbinomWaldTest(dds)")

                elif pycds is None:
                    DESeq_error = 'yes'
                    r("sizeFactor <- rep(1, length(trt))")
                    r("dds$sizeFactor <- sizeFactor")
                    r("dds <- estimateDispersions(dds)")
                    r("dds <- nbinomWaldTest(dds)")

                if DESeq_error == 'no':
                    result += 'Data were normalized by DESeq2...\n'
                elif DESeq_error == 'yes':
                    result += 'DESeq2 cannot run estimateSizeFactors...\n'
                    result += 'Analysis was run without size normalization...\n'
                    result += 'To try again, select a different sample combination or increase the minimum sample size...\n'
                result += '===============================================\n\n\n'

                mergeList = meta_rDF['merge'].tolist()
                mergeSet = list(set(mergeList))

                finalDF = pd.DataFrame()
                for i, val in enumerate(mergeSet):
                    start = i + 1
                    stop = int(len(mergeSet))
                    for j in range(start, stop):
                        if i != j:
                            r.assign("trt1", mergeSet[i])
                            r.assign("trt2", mergeSet[j])
                            r("res <- results(dds, contrast=c('trt', trt1, trt2))")
                            r("baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,dds$trt==trt1, drop=FALSE])")
                            r("baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,dds$trt==trt2, drop=FALSE])")
                            r("df <- data.frame(id=rownames(res), baseMean=res$baseMean, baseMeanA=baseMeanA, baseMeanB=baseMeanB, log2FoldChange=res$log2FoldChange, stderr=res$lfcSE, stat=res$stat, pval=res$pvalue, padj=res$padj)")
                            nbinom_res = r.get("df")

                            names = []
                            try:
                                for item in nbinom_res["id"]:
                                    names.append(findTaxa(item))
                            except Exception as e:
                                print("Failed at loop, "+str(e.args))

                            nbinom_res['Taxa Name'] = names
                            nbinom_res.rename(columns={'id': 'Taxa ID'}, inplace=True)
                            stuff = ['Taxa ID', 'Taxa Name', ' baseMean ', ' baseMeanA ', ' baseMeanB ', ' log2FoldChange ', ' stderr ', ' stat ', ' pval ', ' padj ']
                            nbinom_res = nbinom_res.reindex(columns=stuff)

                            iterationName = str(mergeSet[i]) + ' vs ' + str(mergeSet[j])

                            nbinom_res['Comparison'] = iterationName

                            nbinom_res.rename(columns={' baseMean ': 'baseMean'}, inplace=True)
                            nbinom_res.rename(columns={' baseMeanA ': 'baseMeanA'}, inplace=True)
                            nbinom_res.rename(columns={' baseMeanB ': 'baseMeanB'}, inplace=True)
                            nbinom_res.rename(columns={' log2FoldChange ': 'log2FoldChange'}, inplace=True)
                            nbinom_res.rename(columns={' stderr ': 'StdErr'}, inplace=True)
                            nbinom_res.rename(columns={' stat ': 'Stat'}, inplace=True)
                            nbinom_res.rename(columns={' pval ': 'p-value'}, inplace=True)
                            nbinom_res.rename(columns={' padj ': 'p-adjusted'}, inplace=True)
                            nbinom_res[['p-value', 'p-adjusted']].astype(float)

                            finalDF = pd.concat([finalDF, nbinom_res])

                            base[RID] = 'Step 3 of 5: Performing statistical test...' + str(iterationName) + ' is done!'

                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                            if stops[RID]:
                                return None
                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[RID]:
                        res = ''
                        return None
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 3 of 5: Performing statistical test...done!'
                base[RID] = 'Step 4 of 5: Formatting graph data for display...'

                seriesList = []
                xAxisDict = {}
                yAxisDict = {}

                '''DESeq2 flags genes with abnormal distribtuions and sets the p-value to NaN
                We will replace NaN with 1, so they can be converted to floats, but still not
                be interpreted as significant'''
                #finalDF.fillna(1, inplace=True)

                grouped = finalDF.groupby('Comparison')

                listOfShapes = ['circle', 'square', 'triangle', 'triangle-down', 'diamond',]
                shapeIterator = 0

                FdrVal = float(all['FdrVal'])
                for name, group in grouped:
                    nosigDF = group[group["p-adjusted"] > FdrVal]
                    nosigData = nosigDF[["baseMean", "log2FoldChange"]].values.astype(np.float).tolist()
                    sigDF = group[group["p-adjusted"] <= FdrVal]
                    sigData = sigDF[["baseMean", "log2FoldChange"]].values.astype(np.float).tolist()

                    seriesDict = {}
                    seriesDict['name'] = "NotSig: " + str(name)
                    seriesDict['data'] = nosigData
                    markerDict = {}
                    markerDict['symbol'] = listOfShapes[shapeIterator]
                    seriesDict['marker'] = markerDict
                    seriesList.append(seriesDict)

                    seriesDict = {}
                    seriesDict['name'] = "Sig: " + str(name)
                    seriesDict['data'] = sigData
                    markerDict = {}
                    markerDict['symbol'] = listOfShapes[shapeIterator]
                    seriesDict['marker'] = markerDict
                    seriesList.append(seriesDict)

                    shapeIterator += 1
                    if shapeIterator >= len(listOfShapes):
                        shapeIterator = 0

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[RID]:
                        res = ''
                        return None
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                xTitle = {}
                xTitle['text'] = "baseMean"
                xAxisDict['title'] = xTitle
                xAxisDict['type'] = 'logarithmic'

                yTitle = {}
                yTitle['text'] = "log2FoldChange"
                yAxisDict['title'] = yTitle
                yAxisDict['type'] = 'linear'

                finalDict['series'] = seriesList
                finalDict['xAxis'] = xAxisDict
                finalDict['yAxis'] = yAxisDict

                base[RID] = 'Step 4 of 5: Formatting graph data for display...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 5 of 5:  Formatting nbinomTest results for display...'

                finalDF = finalDF[['Comparison', 'Taxa ID', 'Taxa Name', 'baseMean', 'baseMeanA', 'baseMeanB', 'log2FoldChange', 'StdErr', 'Stat', 'p-value', 'p-adjusted']]
                res_table = finalDF.to_html(classes="table display")
                res_table = res_table.replace('border="1"', 'border="0"')
                finalDict['res_table'] = str(res_table)
                finalDict['text'] = result

                base[RID] = 'Step 5 of 5: Formatting nbinomTest results for display...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                return None

    except:
        if not stop2:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with Differential Abundance!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
        return None


def findTaxa(id):
    taxa = ""
    try:
        temp = Kingdom.objects.filter(kingdomid=id)
        taxa += temp[0].kingdomName
    except:
        # not kingdom, try next one
        try:
            temp = Phyla.objects.filter(phylaid=id)
            taxa += temp[0].phylaName
        except:
            try:
                temp = Class.objects.filter(classid=id)
                taxa += temp[0].className
            except:
                try:
                    temp = Order.objects.filter(orderid=id)
                    taxa += temp[0].orderName
                except:
                    try:
                        temp = Family.objects.filter(familyid=id)
                        taxa += temp[0].familyName
                    except:
                        try:
                            temp = Genus.objects.filter(genusid=id)
                            taxa += temp[0].genusName
                        except:
                            try:
                                temp = Species.objects.filter(speciesid=id)
                                taxa += temp[0].speciesName
                            except:
                                # not found, error!
                                print("Could not find taxa for "+str(id))
    return taxa