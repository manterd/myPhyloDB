#import ast
import datetime
#from django import db
from django.http import HttpResponse
import logging
#import multiprocessing as mp
import numpy as np
import pandas as pd
import pickle
from pyper import *
import simplejson
#import threading
import time

#from database.models import PICRUSt
#from database.models import ko_lvl1, ko_lvl2, ko_lvl3, ko_entry
#from database.models import nz_lvl1, nz_lvl2, nz_lvl3, nz_lvl4, nz_entry
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

                button3 = int(all['button3'])
                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])

                result = ''
                if button3 == 1:
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
                elif button3 == 2:
                    if keggAll == 1:
                        result += 'KEGG Pathway level: 1' + '\n'
                    elif keggAll == 2:
                        result += 'KEGG Pathway level: 2' + '\n'
                    elif keggAll == 3:
                        result += 'KEGG Pathway level: 3' + '\n'
                elif button3 == 3:
                    if nzAll == 1:
                        result += 'KEGG Enzyme level: 1' + '\n'
                    elif nzAll == 2:
                        result += 'KEGG Enzyme level: 2' + '\n'
                    elif nzAll == 3:
                        result += 'KEGG Enzyme level: 3' + '\n'
                    elif nzAll == 4:
                        result += 'KEGG Enzyme level: 4' + '\n'
                    elif keggAll == 5:
                        result += 'KEGG Enzyme level: GIBBs' + '\n'
                    elif keggAll == 6:
                        result += 'KEGG Enzyme level: Nitrogen cycle' + '\n'

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
                tempDF = savedDF.loc[savedDF['sampleid'].isin(catSampleIDs)]

                if metaDictCat:
                    for key in metaDictCat:
                        tempDF = tempDF.loc[tempDF[key].isin(metaDictCat[key])]

                result += 'Categorical variables selected by user: ' + ", ".join(catFields) + '\n'
                result += 'Categorical variables removed from analysis (contains only 1 level): ' + ", ".join(removed) + '\n'
                result += '\n===============================================\n'

                base[RID] = 'Step 1 of 5: Selecting your chosen meta-variables...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 2 of 5: Selecting your chosen taxa or KEGG level...'

                # get selected taxa fro each rank selected in the tree
                DepVar = 1
                finalDF = pd.DataFrame()
                if button3 == 1:
                    DepVar = int(all["DepVar_taxa"])
                    finalDF = getTaxaDF(selectAll, tempDF, catFields_edit, DepVar, RID)

                if button3 == 2:
                    DepVar = int(all["DepVar_kegg"])
                    finalDF = getKeggDF(keggAll, tempDF, catFields_edit, DepVar, RID)

                if button3 == 3:
                    DepVar = int(all["DepVar_nz"])
                    finalDF = getNZDF(nzAll, tempDF, catFields_edit, DepVar, RID)

                count_rDF = pd.DataFrame()
                if DepVar == 1:
                    finalDF['abund'] = finalDF['abund'].round(0).astype(int)
                    count_rDF = finalDF.pivot(index='rank_id', columns='sampleid', values='abund')
                elif DepVar == 4:
                    finalDF['abund_16S'] = finalDF['abund_16S'].round(0).astype(int)
                    count_rDF = finalDF.pivot(index='rank_id', columns='sampleid', values='abund_16S')

                temp_rDF = savedDF.drop_duplicates(subset='sampleid', take_last=True)

                # Create combined metadata column - DiffAbund only
                meta_rDF = temp_rDF.copy()
                if len(catFields) > 1:
                    for index, row in temp_rDF.iterrows():
                       meta_rDF.loc[index, 'merge'] = "; ".join(row[catFields])
                else:
                    meta_rDF.loc[:, 'merge'] = temp_rDF.loc[:, catFields[0]]

                wantedList = ['sampleid', 'merge', 'sample_name']
                meta_rDF = meta_rDF.loc[:, wantedList]
                meta_rDF.set_index('sampleid', drop=True, inplace=True)

                base[RID] = 'Step 2 of 5: Selecting your chosen taxa or KEGG level...done'

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
                result += '\n===============================================\n'

                if DepVar == 1:
                    result += 'Dependent Variable: Abundance' + '\n'
                elif DepVar == 2:
                    result += 'Dependent Variable: Species Richness' + '\n'
                elif DepVar == 3:
                    result += 'Dependent Variable: Species Diversity' + '\n'
                elif DepVar == 4:
                    result += 'Dependent Variable: Abundance (rRNA gene copies)' + '\n'
                result += '\n===============================================\n\n\n'

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
                            if button3 == 1:
                                for item in nbinom_res["id"]:
                                    match = findTaxa(item)
                                    if match == 'Not found':
                                        names.append(item)
                                    else:
                                        names.append(match)

                            if button3 == 2:
                                for item in nbinom_res["id"]:
                                    match = findKEGG(item)
                                    if match == 'Not found':
                                        names.append(item)
                                    else:
                                        names.append(match)

                            if button3 == 3:
                                for item in nbinom_res["id"]:
                                    match = findNZ(item)
                                    if match == 'Not found':
                                        names.append(item)
                                    else:
                                        names.append(match)

                            nbinom_res.loc[:, 'Rank Name'] = names
                            nbinom_res.rename(columns={'id': 'Rank ID'}, inplace=True)
                            stuff = ['Rank ID', 'Rank Name', ' baseMean ', ' baseMeanA ', ' baseMeanB ', ' log2FoldChange ', ' stderr ', ' stat ', ' pval ', ' padj ']
                            nbinom_res = nbinom_res.reindex(columns=stuff)

                            iterationName = str(mergeSet[i]) + ' vs ' + str(mergeSet[j])
                            nbinom_res.loc[:, 'Comparison'] = iterationName

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

                grouped = finalDF.groupby('Comparison')

                listOfShapes = ['circle', 'square', 'triangle', 'triangle-down', 'diamond',]
                shapeIterator = 0

                FdrVal = float(all['FdrVal'])
                for name, group in grouped:
                    nosigDF = group[group["p-adjusted"] > FdrVal]
                    nosigData = []
                    for index, row in nosigDF.iterrows():
                        dataDict = {}
                        dataDict['name'] = row['Rank Name']
                        dataDict['x'] = float(row['baseMean'])
                        dataDict['y'] = float(row['log2FoldChange'])
                        nosigData.append(dataDict)

                    seriesDict = {}
                    seriesDict['name'] = "NotSig: " + str(name)
                    seriesDict['data'] = nosigData
                    markerDict = {}
                    markerDict['symbol'] = listOfShapes[shapeIterator]
                    seriesDict['marker'] = markerDict
                    seriesList.append(seriesDict)

                    sigDF = group[group["p-adjusted"] <= FdrVal]
                    sigData = []
                    for index, row in sigDF.iterrows():
                        dataDict = {}
                        dataDict['name'] = row['Rank Name']
                        dataDict['x'] = float(row['baseMean'])
                        dataDict['y'] = float(row['log2FoldChange'])
                        sigData.append(dataDict)

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
                xTitle['style'] = {'color': 'black', 'fontSize': '18px', 'fontWeight': 'bold'}
                xAxisDict['title'] = xTitle
                xAxisDict['type'] = 'logarithmic'

                yTitle = {}
                yTitle['text'] = "log2FoldChange"
                yTitle['style'] = {'color': 'black', 'fontSize': '18px', 'fontWeight': 'bold'}
                yAxisDict['title'] = yTitle
                yAxisDict['type'] = 'linear'

                styleDict = {'style': {'color': 'black', 'fontSize': '14px'}}
                xAxisDict['labels'] = styleDict
                yAxisDict['labels'] = styleDict

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

                finalDF = finalDF[['Comparison', 'Rank ID', 'Rank Name', 'baseMean', 'baseMeanA', 'baseMeanB', 'log2FoldChange', 'StdErr', 'Stat', 'p-value', 'p-adjusted']]
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
                                return 'Not found'
    return taxa


def findKEGG(id):
    taxa = ""
    try:
        temp = ko_lvl1.objects.using('picrust').filter(ko_lvl1_id=id)
        taxa += temp[0].ko_lvl1_name
    except:
        try:
            temp = ko_lvl2.objects.using('picrust').filter(ko_lvl2_id=id)
            taxa += temp[0].ko_lvl2_name
        except:
            try:
                temp = ko_lvl3.objects.using('picrust').filter(ko_lvl3_id=id)
                taxa += temp[0].ko_lvl3_name
            except:
                try:
                    temp = ko_entry.objects.using('picrust').filter(ko_lvl4_id=id)
                    taxa += temp[0].ko_name
                except:
                    return 'Not found'
    return taxa


def findNZ(id):
    taxa = ""
    try:
        temp = nz_lvl1.objects.using('picrust').filter(nz_lvl1_id=id)
        taxa += temp[0].nz_lvl1_name
    except:
        try:
            temp = nz_lvl2.objects.using('picrust').filter(nz_lvl2_id=id)
            taxa += temp[0].nz_lvl2_name
        except:
            try:
                temp = nz_lvl3.objects.using('picrust').filter(nz_lvl3_id=id)
                taxa += temp[0].nz_lvl3_name
            except:
                try:
                    temp = nz_lvl4.objects.using('picrust').filter(nz_lvl4_id=id)
                    taxa += temp[0].nz_lvl4_name
                except:
                    try:
                        temp = nz_entry.objects.using('picrust').filter(nz_lvl5_id=id)
                        taxa += temp[0].nz_name
                    except:
                        return 'Not found'
    return taxa


def getTaxaDF(selectAll, tempDF, allFields, DepVar, RID):
    global base, stops, stop2, res
    try:
        base[RID] = 'Step 2 of 4: Selecting your chosen taxa or KEGG level...'
        taxaDF = pd.DataFrame(columns=['sampleid', 'rank', 'rank_id', 'rank_name', 'abund', 'abund_16S'])

        if selectAll == 2:
            taxaDF = tempDF.loc[:, ['sampleid', 'phylaid', 'phylaName', 'abund', 'abund_16S']]
            taxaDF.rename(columns={'phylaid': 'rank_id', 'phylaName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Phyla'
        elif selectAll == 3:
            taxaDF = tempDF.loc[:, ['sampleid', 'classid', 'className', 'abund', 'abund_16S']]
            taxaDF.rename(columns={'classid': 'rank_id', 'className': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Class'
        elif selectAll == 4:
            taxaDF = tempDF.loc[:, ['sampleid', 'orderid', 'orderName', 'abund', 'abund_16S']]
            taxaDF.rename(columns={'orderid': 'rank_id', 'orderName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Order'
        elif selectAll == 5:
            taxaDF = tempDF.loc[:, ['sampleid', 'familyid', 'familyName', 'abund', 'abund_16S']]
            taxaDF.rename(columns={'familyid': 'rank_id', 'familyName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Family'
        elif selectAll == 6:
            taxaDF = tempDF.loc[:, ['sampleid', 'genusid', 'genusName', 'abund', 'abund_16S']]
            taxaDF.rename(columns={'genusid': 'rank_id', 'genusName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Genus'
        elif selectAll == 7:
            taxaDF = tempDF.loc[:, ['sampleid', 'speciesid', 'speciesName', 'abund', 'abund_16S']]
            taxaDF.rename(columns={'speciesid': 'rank_id', 'speciesName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Species'

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[RID]:
            res = ''
            return None
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        # Create combined metadata column - DiffAbund only
        if len(allFields) > 1:
            for index, row in tempDF.iterrows():
                tempDF.loc[index, 'merge'] = "; ".join(row.loc[index, allFields])
        else:
            tempDF.loc[:, 'merge'] = tempDF.loc[:, allFields[0]]

        metaDF = tempDF.loc[:, ['merge']]
        finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')

        wantedList = ['sampleid', 'merge', 'rank', 'rank_name', 'rank_id']
        if DepVar == 1:
            finalDF = finalDF.groupby(wantedList)[['abund']].sum()
        elif DepVar == 4:
            finalDF = finalDF.groupby(wantedList)[['abund_16S']].sum()

        finalDF.reset_index(drop=False, inplace=True)
        return finalDF

    except:
        if not stop2:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with Differential Abundance!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)

