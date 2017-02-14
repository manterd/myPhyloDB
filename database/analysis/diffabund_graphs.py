import datetime
from django.http import HttpResponse
import logging
import numpy as np
import pandas as pd
from pyper import *
import json

import database.utils
import database.queue


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getDiffAbund(request, stops, RID, PID):
    try:
        while True:
            if request.is_ajax():
                # Get variables from web page
                allJson = request.body.split('&')[0]
                all = json.loads(allJson)
                database.queue.setBase(RID, 'Step 1 of 6: Reading normalized data file...')

                database.queue.setBase(RID, 'Step 2 of 6: Selecting your chosen meta-variables...')

                treeType = int(all['treeType'])
                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])

                result = ''
                if treeType == 1:
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
                    elif selectAll == 9:
                        result += 'Taxa level: OTU_99' + '\n'
                elif treeType == 2:
                    if keggAll == 1:
                        result += 'KEGG Pathway level: 1' + '\n'
                    elif keggAll == 2:
                        result += 'KEGG Pathway level: 2' + '\n'
                    elif keggAll == 3:
                        result += 'KEGG Pathway level: 3' + '\n'
                elif treeType == 3:
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
                metaValsQuant = []
                metaIDsQuant = []

                treeType = int(all['treeType'])
                DepVar = int(all["DepVar"])

                # Create meta-variable DataFrame, final sample list, final category and quantitative field lists based on tree selections
                savedDF, metaDF, finalSampleIDs, catFields, remCatFields, quantFields, catValues, quantValues = database.utils.getMetaDF(request.user, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, DepVar)
                allFields = catFields + quantFields

                # round data to fix normalization type issues
                savedDF['abund'] = savedDF['abund'].round(0).astype(int)
                savedDF['rel_abund'] = savedDF['rel_abund'].round(0).astype(int)
                savedDF['abund_16S'] = savedDF['abund_16S'].round(0).astype(int)
                savedDF['rich'] = savedDF['rich'].round(0).astype(int)
                savedDF['diversity'] = savedDF['diversity'].round(0).astype(int)

                if not catFields:
                    error = "Selected categorical variable(s) contain only one level.\nPlease select different variable(s)."
                    myDict = {'error': error}
                    res = json.dumps(myDict)
                    return HttpResponse(res, content_type='application/json')

                if not finalSampleIDs:
                    error = "No valid samples were contained in your final dataset.\nPlease select different variable(s)."
                    myDict = {'error': error}
                    res = json.dumps(myDict)
                    return HttpResponse(res, content_type='application/json')

                result = ''
                result += 'Categorical variables selected by user: ' + ", ".join(catFields + remCatFields) + '\n'
                result += 'Categorical variables not included in the statistical analysis (contains only 1 level): ' + ", ".join(remCatFields) + '\n'
                result += 'Quantitative variables selected by user: ' + ", ".join(quantFields) + '\n'
                result += '===============================================\n\n'

                database.queue.setBase(RID, 'Step 2 of 6: Selecting your chosen meta-variables...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 3 of 6: Selecting your chosen taxa or KEGG level...')
                # status for stage 2 could have a loop counter or equivalent, currently the longest stage by far with
                # no proper indicator of how long it will run

                # filter otus based on user settings
                remUnclass = all['remUnclass']
                remZeroes = all['remZeroes']
                perZeroes = int(all['perZeroes'])
                filterData = all['filterData']
                filterPer = int(all['filterPer'])
                filterMeth = int(all['filterMeth'])
                mapTaxa = 'no'

                finalDF = pd.DataFrame()
                if treeType == 1:
                    if selectAll != 8:
                        filteredDF = database.utils.filterDF(savedDF, DepVar, selectAll, remUnclass, remZeroes, perZeroes, filterData, filterPer, filterMeth)
                    else:
                        filteredDF = savedDF.copy()

                    finalDF, missingList = database.utils.getTaxaDF(selectAll, '', filteredDF, metaDF, allFields, DepVar, RID, stops, PID)

                    if selectAll == 8:
                        result += '\nThe following PGPRs were not detected: ' + ", ".join(missingList) + '\n'
                        result += '===============================================\n'

                if treeType == 2:
                    finalDF, allDF = database.utils.getKeggDF(keggAll, '', savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)

                if treeType == 3:
                    finalDF, allDF = database.utils.getNZDF(nzAll, '', savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)

                if finalDF.empty:
                    error = "Selected taxa were not found in your selected samples."
                    myDict = {'error': error}
                    res = json.dumps(myDict)
                    return HttpResponse(res, content_type='application/json')

                # make sure column types are correct
                finalDF[catFields] = finalDF[catFields].astype(str)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/diffabund/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                count_rDF = pd.DataFrame()
                if DepVar == 0:
                    finalDF['abund'] = finalDF['abund'].round(0).astype(int)
                    count_rDF = finalDF.pivot(index='rank_id', columns='sampleid', values='abund')
                elif DepVar == 4:
                    finalDF['abund_16S'] = finalDF['abund_16S'].round(0).astype(int)
                    count_rDF = finalDF.pivot(index='rank_id', columns='sampleid', values='abund_16S')

                count_rDF.fillna(0, inplace=True)

                # Create combined metadata column - DiffAbund only
                if len(catFields) > 1:
                    for index, row in metaDF.iterrows():
                       metaDF.loc[index, 'merge'] = "; ".join(row[catFields])
                else:
                    metaDF.loc[:, 'merge'] = metaDF.loc[:, catFields[0]]

                database.queue.setBase(RID, 'Step 3 of 6: Selecting your chosen taxa or KEGG level...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 4 of 6: Performing statistical test...')

                finalDict = {}
                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                database.queue.setBase(RID, 'Verifying R packages...missing packages are being installed')

                r("list.of.packages <- c('DESeq2')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                r("if (length(new.packages)) source('http://bioconductor.org/biocLite.R')")
                print r("if (length(new.packages)) biocLite(new.packages)")

                database.queue.setBase(RID, 'Step 4 of 6: Performing statistical test...')

                print r("library(DESeq2)")

                r.assign("metaDF", metaDF)
                r("trt <- factor(metaDF$merge)")

                r.assign("count", count_rDF)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                if DepVar == 0:
                    result += 'Dependent Variable: Abundance' + '\n'
                elif DepVar == 4:
                    result += 'Dependent Variable: Total Abundance' + '\n'
                result += '\n===============================================\n\n\n'

                r("colData <- data.frame(row.names=colnames(count), trt=trt)")
                r("dds <- DESeqDataSetFromMatrix(countData=count, colData=colData, design = ~ trt)")

                r("sizeFactor <- rep(1, length(trt))")
                r("dds$sizeFactor <- sizeFactor")
                r("dds <- estimateDispersions(dds)")
                r("dds <- nbinomWaldTest(dds)")

                mergeList = metaDF['merge'].tolist()
                mergeSet = list(set(mergeList))
                nbinom_res = pd.DataFrame()
                for i, val in enumerate(mergeSet):
                    start = i + 1
                    stop = int(len(mergeSet))
                    for j in range(start, stop):
                        if i != j:
                            r.assign("trt1", mergeSet[i])
                            r.assign("trt2", mergeSet[j])

                            # round values, convert to int
                            r("res <- results(dds, contrast=c('trt', trt1, trt2))")
                            r("baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,dds$trt==trt1, drop=FALSE])")
                            r("baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,dds$trt==trt2, drop=FALSE])")
                            r("df <- data.frame(rank_id=rownames(res), baseMean=res$baseMean, baseMeanA=baseMeanA, baseMeanB=baseMeanB, log2FoldChange=-res$log2FoldChange, stderr=res$lfcSE, stat=res$stat, pval=res$pvalue, padj=res$padj)")
                            nbinom_res = r.get("df")

                            if nbinom_res is None:
                                myDict = {'error': "DESeq failed!\nPlease try a different data combination."}
                                res = json.dumps(myDict)
                                return HttpResponse(res, content_type='application/json')

                            # remove taxa that failed (i.e., both trts are zero or log2FoldChange is NaN)
                            nbinom_res = nbinom_res.loc[pd.notnull(nbinom_res[' log2FoldChange '])]

                            if treeType == 1:
                                idList = database.utils.getFullTaxonomy(list(nbinom_res.rank_id.unique()))
                                nbinom_res['Taxonomy'] = nbinom_res['rank_id'].map(idList)
                            elif treeType == 2:
                                idList = database.utils.getFullKO(list(nbinom_res.rank_id.unique()))
                                nbinom_res['Taxonomy'] = nbinom_res['rank_id'].map(idList)
                            elif treeType == 3:
                                idList = database.utils.getFullNZ(list(nbinom_res.rank_id.unique()))
                                nbinom_res['Taxonomy'] = nbinom_res['rank_id'].map(idList)

                            nbinom_res.rename(columns={'rank_id': 'Rank ID'}, inplace=True)

                            iterationName = str(mergeSet[i]) + ' vs ' + str(mergeSet[j])
                            nbinom_res.insert(1, 'Comparison', iterationName)
                            nbinom_res.rename(columns={' baseMean ': 'baseMean'}, inplace=True)
                            nbinom_res.rename(columns={' baseMeanA ': 'baseMeanA'}, inplace=True)
                            nbinom_res.rename(columns={' baseMeanB ': 'baseMeanB'}, inplace=True)
                            nbinom_res.rename(columns={' log2FoldChange ': 'log2FoldChange'}, inplace=True)
                            nbinom_res.rename(columns={' stderr ': 'StdErr'}, inplace=True)
                            nbinom_res.rename(columns={' stat ': 'Stat'}, inplace=True)
                            nbinom_res.rename(columns={' pval ': 'p-value'}, inplace=True)
                            nbinom_res.rename(columns={' padj ': 'p-adjusted'}, inplace=True)

                            database.queue.setBase(RID, 'Step 4 of 6: Performing statistical test...' + str(iterationName) + ' is done!')

                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                            if stops[PID] == RID:
                                res = ''
                                return HttpResponse(res, content_type='application/json')
                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 4 of 6: Performing statistical test...done!')
                database.queue.setBase(RID, 'Step 5 of 6: Formatting graph data for display...')

                seriesList = []
                xAxisDict = {}
                yAxisDict = {}

                grouped = nbinom_res.groupby('Comparison')

                listOfShapes = ['circle', 'square', 'triangle', 'triangle-down', 'diamond',]
                shapeIterator = 0

                FdrVal = float(all['FdrVal'])
                for name, group in grouped:
                    nosigDF = group[group["p-adjusted"] > FdrVal]
                    nosigData = []
                    for index, row in nosigDF.iterrows():
                        dataDict = {}
                        dataDict['name'] = 'ID: ' + row['Rank ID']
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
                        dataDict['name'] = row['Rank ID']
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
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                xTitle = {}
                xTitle['text'] = "baseMean"
                xTitle['style'] = {'fontSize': '18px', 'fontWeight': 'bold'}
                xAxisDict['title'] = xTitle
                xAxisDict['type'] = 'logarithmic'

                yTitle = {}
                yTitle['text'] = "log2FoldChange"
                yTitle['style'] = {'fontSize': '18px', 'fontWeight': 'bold'}
                yAxisDict['title'] = yTitle
                yAxisDict['type'] = 'linear'

                styleDict = {'style': {'fontSize': '14px'}}
                xAxisDict['labels'] = styleDict
                yAxisDict['labels'] = styleDict

                finalDict['series'] = seriesList
                finalDict['xAxis'] = xAxisDict
                finalDict['yAxis'] = yAxisDict

                database.queue.setBase(RID, 'Step 5 of 6: Formatting graph data for display...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 6 of 6:  Formatting nbinomTest results for display...')

                nbinom_res.replace(to_replace='N/A', value=np.nan, inplace=True)
                nbinom_res.dropna(axis=1, how='all', inplace=True)
                res_table = nbinom_res.to_html(classes="table display")
                res_table = res_table.replace('border="1"', 'border="0"')
                finalDict['res_table'] = str(res_table)

                finalDict['text'] = result

                database.queue.setBase(RID, 'Step 6 of 6: Formatting nbinomTest results for display...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDict['error'] = 'none'
                res = json.dumps(finalDict)

                return HttpResponse(res, content_type='application/json')

    except Exception as e:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "There was an error during your analysis:\nError: " + str(e.message) + "\nTimestamp: " + str(datetime.datetime.now())
            res = json.dumps(myDict)
            return HttpResponse(res, content_type='application/json')



