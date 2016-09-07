import datetime
from django.http import HttpResponse
import logging
import pandas as pd
from pyper import *
import simplejson

from database.utils import multidict
from database.utils_kegg import getTaxaDF, getKeggDF, getNZDF
from database.utils_kegg import getFullTaxonomy, getFullKO, getFullNZ, insertTaxaInfo
import database.queue


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getDiffAbund(request, stops, RID, PID):
    try:
        while True:
            if request.is_ajax():
                # Get variables from web page
                allJson = request.body.split('&')[0]
                all = simplejson.loads(allJson)
                database.queue.setBase(RID, 'Step 1 of 5: Selecting your chosen meta-variables...')
                myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
                path = str(myDir) + 'usr_norm_data.csv'

                with open(path, 'rb') as f:
                    savedDF = pd.read_csv(f, index_col=0, sep=',')

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
                    return HttpResponse(res, content_type='application/json')

                catSampleIDs = []
                if metaIDsCat:
                    idDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaIDsCat)
                    for key in sorted(idDictCat):
                        catSampleIDs.extend(idDictCat[key])

                # Removes samples (rows) that are not in our samplelist
                tempDF = savedDF.loc[savedDF['sampleid'].isin(catSampleIDs)]

                # make sure column types are correct
                tempDF[catFields_edit] = tempDF[catFields_edit].astype(str)

                if metaDictCat:
                    for key in metaDictCat:
                        tempDF = tempDF.loc[tempDF[key].isin(metaDictCat[key])]

                metaDF = tempDF[catFields_edit]

                result += 'Categorical variables selected by user: ' + ", ".join(catFields) + '\n'
                result += 'Categorical variables removed from analysis (contains only 1 level): ' + ", ".join(removed) + '\n'
                result += '===============================================\n\n'

                button3 = int(all['button3'])
                DepVar = 1
                if button3 == 1:
                    DepVar = int(all["DepVar_taxa"])
                elif button3 == 2:
                    DepVar = int(all["DepVar_kegg"])
                elif button3 == 3:
                    DepVar = int(all["DepVar_nz"])

                if DepVar == 4:
                    savedDF = savedDF.loc[savedDF['abund_16S'] != 0]
                    rows, cols = savedDF.shape
                    if rows < 1:
                        myDict = {'error': "Error: no qPCR or 'rRNA gene copies' data were found for this dataset"}
                        res = simplejson.dumps(myDict)
                        return HttpResponse(res, content_type='application/json')

                    finalSampleList = pd.unique(savedDF.sampleid.ravel().tolist())
                    remSampleList = list(set(catSampleIDs) - set(finalSampleList))

                    result += str(len(remSampleList)) + " samples were removed from analysis (missing 'rRNA gene copies' data)\n"
                    result += '===============================================\n\n'

                database.queue.setBase(RID, 'Step 1 of 5: Selecting your chosen meta-variables...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 2 of 5: Selecting your chosen taxa or KEGG level...')

                # get selected taxa fro each rank selected in the tree
                finalDF = pd.DataFrame()
                if button3 == 1:
                    finalDF, missingList = getTaxaDF('abund', selectAll, '', savedDF, metaDF, catFields_edit, DepVar, RID, stops, PID)
                    if selectAll == 8:
                        result += '\nThe following PGPRs were not detected: ' + ", ".join(missingList) + '\n'
                        result += '===============================================\n'

                if button3 == 2:
                    finalDF = getKeggDF('abund', keggAll, '', savedDF, tempDF, catFields_edit, DepVar, RID, stops, PID)

                if button3 == 3:
                    finalDF = getNZDF('abund', nzAll, '', savedDF, tempDF, catFields_edit, DepVar, RID, stops, PID)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/diffabund/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                count_rDF = pd.DataFrame()
                if DepVar == 1:
                    finalDF['abund'] = finalDF['abund'].round(0).astype(int)
                    count_rDF = finalDF.pivot(index='rank_id', columns='sampleid', values='abund')
                elif DepVar == 4:
                    finalDF['abund_16S'] = finalDF['abund_16S'].round(0).astype(int)
                    count_rDF = finalDF.pivot(index='rank_id', columns='sampleid', values='abund_16S')

                temp_rDF = savedDF.drop_duplicates(subset='sampleid', take_last=True)

                # Removes samples (rows) that are not in our samplelist
                temp_rDF = temp_rDF.loc[temp_rDF['sampleid'].isin(catSampleIDs)]

                # make sure column types are correct
                temp_rDF[catFields_edit] = temp_rDF[catFields_edit].astype(str)

                if metaDictCat:
                    for key in metaDictCat:
                        temp_rDF = temp_rDF.loc[temp_rDF[key].isin(metaDictCat[key])]

                # Create combined metadata column - DiffAbund only
                meta_rDF = temp_rDF.copy()
                if len(catFields) > 1:
                    for index, row in temp_rDF.iterrows():
                       meta_rDF.loc[index, 'merge'] = "; ".join(row[catFields_edit])
                else:
                    meta_rDF.loc[:, 'merge'] = temp_rDF.loc[:, catFields_edit[0]]

                wantedList = ['sampleid', 'merge', 'sample_name']
                meta_rDF = meta_rDF.loc[:, wantedList]
                meta_rDF.set_index('sampleid', drop=True, inplace=True)
                database.queue.setBase(RID, 'Step 2 of 5: Selecting your chosen taxa or KEGG level...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 3 of 5: Performing statistical test...')

                finalDict = {}
                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                r.assign("metaDF", meta_rDF)
                r("trt <- factor(metaDF$merge)")

                r.assign("count", count_rDF)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                r("library(DESeq2)")
                r("colData <- data.frame(row.names=colnames(count), trt=trt)")
                r("dds <- DESeqDataSetFromMatrix(countData=count, colData=colData, design = ~ trt)")

                r("sizeFactor <- rep(1, length(trt))")
                r("dds$sizeFactor <- sizeFactor")
                r("dds <- estimateDispersions(dds)")
                r("dds <- nbinomWaldTest(dds)")

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
                            # run DESeq here?
                            r("res <- results(dds, contrast=c('trt', trt1, trt2))")
                            r("baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,dds$trt==trt1, drop=FALSE])")
                            r("baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,dds$trt==trt2, drop=FALSE])")
                            r("df <- data.frame(rank_id=rownames(res), baseMean=res$baseMean, baseMeanA=baseMeanA, baseMeanB=baseMeanB, log2FoldChange=-res$log2FoldChange, stderr=res$lfcSE, stat=res$stat, pval=res$pvalue, padj=res$padj)")
                            nbinom_res = r.get("df")

                            # remove taxa that failed (i.e., both trts are zero or log2FoldChange is NaN)
                            nbinom_res = nbinom_res.loc[pd.notnull(nbinom_res[' log2FoldChange '])]
                            # /\ AttributeError: 'NoneType' object has no attribute 'loc' /\

                            if button3 == 1:
                                zipped = getFullTaxonomy(nbinom_res['rank_id'])
                                insertTaxaInfo(button3, zipped, nbinom_res, pos=1)
                            elif button3 == 2:
                                zipped = getFullKO(nbinom_res['rank_id'])
                                insertTaxaInfo(button3, zipped, nbinom_res, pos=1)
                            elif button3 == 3:
                                zipped = getFullNZ(nbinom_res['rank_id'])
                                insertTaxaInfo(button3, zipped, nbinom_res, pos=1)

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
                            finalDF = pd.concat([finalDF, nbinom_res])

                            database.queue.setBase(RID, 'Step 3 of 5: Performing statistical test...' + str(iterationName) + ' is done!')

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

                database.queue.setBase(RID, 'Step 3 of 5: Performing statistical test...done!')
                database.queue.setBase(RID, 'Step 4 of 5: Formatting graph data for display...')

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
                        dataDict['name'] = row['Rank ID']

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

                database.queue.setBase(RID, 'Step 4 of 5: Formatting graph data for display...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 5 of 5:  Formatting nbinomTest results for display...')

                res_table = finalDF.to_html(classes="table display")
                res_table = res_table.replace('border="1"', 'border="0"')
                finalDict['res_table'] = str(res_table)
                finalDict['text'] = result

                database.queue.setBase(RID, 'Step 5 of 5: Formatting nbinomTest results for display...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)

                return HttpResponse(res, content_type='application/json')

    except Exception as e:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "There was an error during your analysis:\nError: " + str(e.message) + "\nTimestamp: " + str(datetime.datetime.now())
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')



