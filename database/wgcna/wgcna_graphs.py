import datetime
from django.http import HttpResponse
import logging
from natsort import natsorted
import pandas as pd
from pyper import *
from PyPDF2 import PdfFileReader, PdfFileMerger
import simplejson
import sys

from database.utils import multidict
from database.utils_kegg import getTaxaDF, getKeggDF, getNZDF
from database.utils_kegg import getFullTaxonomy, getFullKO, getFullNZ, insertTaxaInfo
import database.queue


reload(sys)
sys.setdefaultencoding('utf8')

LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getWGCNA(request, stops, RID, PID):
    try:
        while True:
            if request.is_ajax():
                allJson = request.body.split('&')[0]
                all = simplejson.loads(allJson)
                database.queue.setBase(RID, 'Step 1 of 6: Selecting your chosen meta-variables...')
                myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
                path = str(myDir) + 'usr_norm_data.csv'

                with open(path, 'rb') as f:
                    savedDF = pd.read_csv(f, index_col=0, sep=',')

                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])

                result = ''
                button3 = int(all['button3'])
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
                metaValsQuant = all['metaValsQuant']
                metaIDsQuant = all['metaIDsQuant']

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

                catSampleIDs = []
                if metaIDsCat:
                    idDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaIDsCat)
                    for key in sorted(idDictCat):
                        if idDictCat[key] not in catSampleIDs:
                            catSampleIDs.extend(idDictCat[key])

                metaDictQuant = {}
                quantFields = []
                quantValues = []
                if metaValsQuant:
                    metaDictQuant = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaValsQuant)
                    for key in sorted(metaDictQuant):
                        quantFields.append(key)
                        quantValues.extend(metaDictQuant[key])

                quantSampleIDs = []
                if metaIDsQuant:
                    idDictQuant = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaIDsQuant)
                    for key in sorted(idDictQuant):
                        if idDictQuant[key] not in quantSampleIDs:
                            quantSampleIDs.extend(idDictQuant[key])

                allSampleIDs = list(set(catSampleIDs) | set(quantSampleIDs))
                allFields = catFields_edit + quantFields

                # Removes samples (rows) that are not in our samplelist
                metaDF = savedDF.drop_duplicates(subset='sampleid', take_last=True)
                if allSampleIDs:
                    metaDF = metaDF.loc[metaDF['sampleid'].isin(allSampleIDs)]

                # make sure column types are correct
                metaDF[catFields_edit] = metaDF[catFields_edit].astype(str)
                metaDF[quantFields] = metaDF[quantFields].astype(float)

                if metaDictCat:
                    for key in metaDictCat:
                        metaDF = metaDF.loc[metaDF[key].isin(metaDictCat[key])]

                if metaDictQuant:
                    for key in metaDictQuant:
                        valueList = [float(x) for x in metaDictQuant[key]]
                        metaDF = metaDF.loc[metaDF[key].isin(valueList)]

                finalSampleList = metaDF.sampleid.tolist()
                wantedList = allFields + ['sampleid', 'sample_name']
                metaDF = metaDF[wantedList]
                metaDF.set_index('sampleid', drop=True, inplace=True)

                if allFields:
                    result += 'Categorical variables selected by user: ' + ", ".join(catFields) + '\n'
                    result += 'Categorical variables removed from analysis (contains only 1 level): ' + ", ".join(removed) + '\n'
                    result += 'Quantitative variables selected by user: ' + ", ".join(quantFields) + '\n'
                else:
                    result += 'No meta variables were selected\n'
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

                    remSampleList = list(set(catSampleIDs) - set(finalSampleList))

                    result += str(len(remSampleList)) + " samples were removed from analysis (missing 'rRNA gene copies' data)\n"
                    result += '===============================================\n\n'

                # settings block
                networkType = all['networkType']
                corType = all['corType']
                maxPOutliers = int(all['maxPOutliers'])
                deepSplit = float(all['deepSplit'])
                detectCutHeight = float(all['detectCutHeight'])
                minModuleSize = int(all['minModuleSize'])
                reassignThreshold = float(all['reassignThreshold'])
                minCoreKME = float(all['minCoreKME'])
                minCoreKMESize = int(all['minCoreKMESize'])
                minKMEtoStay = float(all['minKMEtoStay'])
                mergeCutHeight = float(all['mergeCutHeight'])
                minEdge = float(all['minEdge'])
                minKME = float(all['minKME'])
                maxNGenes = int(all['maxNGenes'])
                graphLayout = all['graphLayout']

                result += 'Settings used in WGCNA calculation:\n'
                result += '\nTopological Overlap:\n'
                result += 'networkType: ' + str(networkType) + '\n'
                result += 'corType: ' + str(corType) + '\n'
                result += 'maxPOutliers: ' + str(maxPOutliers) + '\n'

                result += '\nBasic Tree Cut:\n'
                result += 'deepSplit: ' + str(deepSplit) + '\n'
                result += 'detectCutHeight: ' + str(detectCutHeight) + '\n'
                result += 'minModuleSize: ' + str(minModuleSize) + '\n'

                result += '\nGene Reassignment:\n'
                result += 'reassignThreshold: ' + str(reassignThreshold) + '\n'
                result += 'minCoreKME: ' + str(minCoreKME) + '\n'
                result += 'minCoreKMESize: ' + str(minCoreKMESize) + '\n'
                result += 'minKMEtoStay: ' + str(minKMEtoStay) + '\n'

                result += '\nModule Merging:\n'
                result += 'mergeCutHeight: ' + str(mergeCutHeight) + '\n'

                result += '\nNetwork Graph:\n'
                result += 'minEdge: ' + str(minEdge) + '\n'
                result += 'minKME: ' + str(minKME) + '\n'
                result += 'maxNGenes: ' + str(maxNGenes) + '\n'
                result += 'graphLayout: ' + str(graphLayout) + '\n'

                database.queue.setBase(RID, 'Step 1 of 6: Selecting your chosen meta-variables...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 2 of 6: Selecting your chosen taxa or KEGG level...')

                finalDF = pd.DataFrame()
                if button3 == 1:
                    finalDF, missingList = getTaxaDF('rel_abund', selectAll, '', savedDF, metaDF, catFields_edit,DepVar, RID, stops, PID)
                    if selectAll == 8:
                        result += '\nThe following PGPRs were not detected: ' + ", ".join(missingList) + '\n'
                        result += '===============================================\n'

                if button3 == 2:
                    finalDF = getKeggDF('rel_abund', keggAll, '', savedDF, metaDF, catFields_edit, DepVar, RID, stops, PID)

                if button3 == 3:
                    finalDF = getNZDF('rel_abund', nzAll, '', savedDF, metaDF, catFields_edit, DepVar, RID, stops, PID)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/wgcna/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                count_rDF = pd.DataFrame()
                if DepVar == 1:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='rel_abund')
                elif DepVar == 2:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='rich')
                elif DepVar == 3:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='diversity')
                elif DepVar == 4:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='abund_16S')

                database.queue.setBase(RID, 'Step 2 of 6: Selecting your chosen taxa...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 3 of 6: WGCNA analysis...')

                finalDict = {}
                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                r("library(ggplot2)")
                r("library(reshape2)")
                r("library(WGCNA)")
                r("library(cluster)")
                r("library(ggplus)")
                r("allowWGCNAThreads()")
                r("options(stringAsFactors=FALSE)")
                r("pdf_counter <- 1")

                path = 'myPhyloDB/media/temp/wgcna/Rplots/%s' % RID
                if not os.path.exists(path):
                    os.makedirs(path)

                r.assign("path", path)
                r.assign("RID", RID)

                r.assign("datExpr", count_rDF)
                r.assign("geneID", count_rDF.columns.values.tolist())
                r.assign("sampleID", metaDF.index.values.tolist())
                r("datExpr[] <- lapply(datExpr, as.numeric)")
                r("names(datExpr) <- geneID")

                r("options(strngsAsFactors=FALSE)")
                r("nGenes <- ncol(datExpr)")
                r("nSamples <- nrow(datExpr)")

                # data cleaning
                outString = str(r("gsg <- goodSamplesGenes(datExpr, minNSamples=3, minNGenes=3, verbose=3)"))

                error = "\nOutput from WGCNA...goodSampleGenes\n\n"
                lines = outString.split('\n')
                for line in lines[1:]:
                        error += str(line) + '\n'
                error += '===============================================\n'
                r("collectGarbage()")

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #


                gsg = r.get("gsg")
                if not gsg:
                    finalDict['text'] = error
                    finalDict['error'] = "Network construction failed!\nPlease check the 'Test Results' section for more information."
                    res = simplejson.dumps(finalDict)
                    return HttpResponse(res, content_type='application/json')

                r(" if (!gsg$allOK) { \
                        if ( sum(!gsg$goodGenes) > 0) { \
                            datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes] \
                        } \
                    } \
                ")

                # R settings assigns
                r.assign("networkType", networkType)
                r.assign("corType", corType)
                r.assign("maxPOutliers", maxPOutliers)
                r.assign("deepSplit", deepSplit)
                r.assign("detectCutHeight", detectCutHeight)
                r.assign("minModuleSize", minModuleSize)
                r.assign("reassignThreshold", reassignThreshold)
                r.assign("minCoreKME", minCoreKME)
                r.assign("minCoreKMESize", minCoreKMESize)
                r.assign("minKMEtoStay", minKMEtoStay)
                r.assign("mergeCutHeight", mergeCutHeight)
                r.assign("minEdge", minEdge)
                r.assign("minKME", minKME)
                r.assign("maxNGenes", maxNGenes)
                r.assign("graphLayout", graphLayout)


                # Set soft-thresholding power based on number of samples
                nSamples, col = metaDF.shape
                if networkType == 'unsigned':
                    if nSamples < 20:
                        myPower = 10
                    elif nSamples >= 20 and nSamples < 30:
                        myPower = 9
                    elif nSamples >= 30 and nSamples < 40:
                        myPower = 8
                    elif nSamples >= 40 and nSamples < 60:
                        myPower = 7
                    else:
                        myPower = 6
                else:
                    if nSamples < 20:
                        myPower = 20
                    elif nSamples >= 20 and nSamples < 30:
                        myPower = 18
                    elif nSamples >= 30 and nSamples < 40:
                        myPower = 16
                    elif nSamples >= 40 and nSamples < 60:
                        myPower = 14
                    else:
                        myPower = 12
                r.assign("myPower", myPower)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                #WGCNA calculation
                outString = str(r("net <- blockwiseModules(datExpr, power=myPower, \
                    networkType=networkType, minModuleSize=minModuleSize, \
                    reassignThreshold=reassignThreshold, mergeCutHeight=mergeCutHeight, numericLabels=TRUE, \
                    pamRespectsDendro=FALSE, saveTOMs=TRUE, saveTOMFileBase=paste(path, 'TOM', sep='/'), \
                    minNSamples=3, minNGenes=3, verbose=3, maxBlockSize=nGenes, corType=corType, \
                    maxPOutliers=maxPOutliers, deepSplit=deepSplit, detectCutHeight=detectCutHeight, \
                    minCoreKME=minCoreKME, minCoreKMESize=minCoreKMESize, minKMEtoStay=minKMEtoStay) \
                "))

                error = "\nOutput from WGCNA...blockwiseModules\n\n"
                lines = outString.split('\n')
                for line in lines[1:]:
                    error += line
                    error += '\n'
                error += '===============================================\n'
                r("collectGarbage()")

                net = r.get("net")
                if not net:
                    finalDict['text'] = error
                    finalDict['error'] = "Network construction failed!\nPlease check the 'Test Results' section for more information."
                    res = simplejson.dumps(finalDict)
                    return HttpResponse(res, content_type='application/json')

                ### create our own custom color vector
                r('BaseColors = c("turquoise", "blue", "brown", "yellow", "green", "red", \
                    "pink", "magenta", "purple", "greenyellow", "tan", "salmon", "cyan", \
                    "midnightblue", "lightcyan", "grey60", "lightgreen", "lightyellow", \
                    "royalblue", "darkred", "darkgreen", "darkturquoise", "darkgrey", \
                    "orange", "darkorange", "skyblue", "saddlebrown", "steelblue", \
                    "paleturquoise", "violet", "darkolivegreen", "darkmagenta")')

                r('Rcolors <- colors()[-grep("grey", colors())]')
                r('Rcolors <- Rcolors[-grep("gray", Rcolors)]')
                r('Rcolors <- Rcolors[-grep("black", Rcolors)]')
                r('Rcolors <- Rcolors[-grep("white", Rcolors)]')
                r('InBase <- match(BaseColors, Rcolors)')
                r('ExtraColors <- Rcolors[-c(InBase[!is.na(InBase)])]')
                r('ExtraColors <- sample(ExtraColors)')
                r('colorSeq <- c(BaseColors, ExtraColors)')

                #Reformat colors
                r("mergedColors <- labels2colors(net$colors, colorSeq=colorSeq)")
                r("unmergedColors <- labels2colors(net$unmergedColors, colorSeq=colorSeq)")
                r("colorlevels <- unique(mergedColors)")
                r("consTree <- net$dendrograms[[1]]")
                r("datME <- net$MEs")

                # rename MEs with color not index
                r("MEColors = labels2colors(as.numeric(substring(names(datME), 3)), colorSeq=colorSeq)")
                r("rownames(datME) <- rownames(datExpr)")
                r("names(datME) <- paste('ME.', MEColors, sep='')")

                # Get datExprframe with list of genes by module
                r("rank_id <- names(datExpr)")
                r("levels <- unique(net$colors)")

                r("df <- as.data.frame(rank_id)")
                r("for (i in 1:length(levels) ) { \
                    modBool <- (net$colors==levels[i]); \
                    color <- labels2colors(levels[i], colorSeq=colorSeq); \
                    df[,color] <- modBool; \
                }")

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                # Fundamental network concepts
                r("ADJ1 <- abs(cor(datExpr, use='p'))^6")
                r("ADJ1")
                r("netStats <- fundamentalNetworkConcepts(ADJ1)")
                r("connectDF <- as.data.frame(t(rbind(netStats$Connectivity, \
                    netStats$ScaledConnectivity, netStats$ClusterCoef, netStats$MAR)))")
                r("names(connectDF) <- c('Whole Network Connectivity', 'Scaled Connectivity', 'Clustering Coefficient', 'Maximum Adjacency Ratio')")
                r("connectDF$Module <- mergedColors")
                r("connectDF <- connectDF[,c(5,1,2,3,4)]")

                moduleDF = r.get("connectDF")
                index = r.get("rownames(connectDF)")
                moduleDF = pd.DataFrame(moduleDF)
                moduleDF.rename(columns={' Whole Network Connectivity ': 'Whole Network Connectivity'}, inplace=True)
                moduleDF.insert(0, 'rank_id', index)

                zipped = []  # JUMP  Potentially run an if to check for all select level
                if button3 == 1:
                    zipped = getFullTaxonomy(moduleDF['rank_id'])  # had selectAll AND moduleDF
                elif button3 == 2:
                    zipped = getFullKO(moduleDF['rank_id'])
                elif button3 == 3:
                    zipped = getFullNZ(moduleDF['rank_id'])


                if button3 == 1:
                    # removed split based on select level as getFullTaxonomy returns a full set
                    k, p, c, o, f, g, s = map(None, *zipped)
                    moduleDF.insert(1, 'Kingdom', k)
                    moduleDF.insert(2, 'Phyla', p)
                    moduleDF.insert(3, 'Class', c)
                    moduleDF.insert(4, 'Order', o)
                    moduleDF.insert(5, 'Family', f)
                    moduleDF.insert(6, 'Genus', g)
                    moduleDF.insert(7, 'Species', s)
                if button3 == 2:
                    # same as taxa
                    L1, L2, L3 = map(None, *zipped)
                    moduleDF.insert(1, 'Level_1', L1)
                    moduleDF.insert(2, 'Level_2', L2)
                    moduleDF.insert(3, 'Level_3', L3)
                if button3 == 3:
                    # same as taxa
                    L1, L2, L3, L4 = map(None, *zipped)
                    moduleDF.insert(1, 'Level_1', L1)
                    moduleDF.insert(2, 'Level_2', L2)
                    moduleDF.insert(3, 'Level_3', L3)
                    moduleDF.insert(4, 'Level_4', L4)

                    moduleDF.fillna(value='N/A', inplace=True)

                dist_table = moduleDF.to_html(classes="table display")
                dist_table = dist_table.replace('border="1"', 'border="0"')
                finalDict['dist_table'] = str(dist_table)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                # kME
                r("datKME <- signedKME(datExpr, datME, outputColumnName='MM')")
                r("rownames(datKME) <- names(datExpr)")
                kmeDF = r.get("datKME")
                kmeDF = pd.DataFrame(kmeDF)
                index = r.get("names(datExpr)")
                kmeDF.insert(0, 'rank_id', index)

                zipped = []
                if button3 == 1:
                    zipped = getFullTaxonomy(list(kmeDF['rank_id']))
                    insertTaxaInfo(button3, zipped, kmeDF, pos=1)
                elif button3 == 2:
                    zipped = getFullKO(list(kmeDF['rank_id']))
                    insertTaxaInfo(button3, zipped, kmeDF, pos=1)
                elif button3 == 3:
                    zipped = getFullNZ(list(kmeDF['rank_id']))
                    insertTaxaInfo(button3, zipped, kmeDF, pos=1)

                kme_table = kmeDF.to_html(classes="table display")
                kme_table = kme_table.replace('border="1"', 'border="0"')
                finalDict['kme_table'] = str(kme_table)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                # Create Eigengene Table
                r.assign("meta", metaDF)
                r("aovDF <- data.frame(datME)")
                r("names(aovDF) <- names(datME)")
                r("rownames(aovDF) <- rownames(datME)")

                r("aovDF <- cbind(meta, aovDF)")
                aovDF = r.get("aovDF")
                indices = r.get("rownames(aovDF)")
                eigDF = pd.DataFrame(aovDF)
                eigDF['sampleid'] = indices
                cols = list(eigDF)
                cols.insert(0, cols.pop(cols.index('sampleid')))
                eigDF = eigDF.ix[:, cols]
                res_table = eigDF.to_html(classes="table display")
                res_table = res_table.replace('border="1"', 'border="0"')
                finalDict['res_table'] = str(res_table)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                # plot dendrogram
                r("pdf(paste(path, '/wgcna_temp', pdf_counter, '.pdf', sep=''), height=6, width=5+0.01*nGenes)")
                r("pdf_counter <- pdf_counter + 1")
                r("plotDendroAndColors(net$dendrograms[[1]], \
                    cbind(mergedColors[net$blockGenes[[1]]], unmergedColors[net$blockGenes[[1]]]), \
                    c('Merged', 'Unmerged'), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, \
                    guideHang=0.05, main='Taxa Dendrogram')")
                r("dev.off()")

                #relating the modules
                r("pdf(paste(path, '/wgcna_temp', pdf_counter, '.pdf', sep=''), height=6, width=5+0.01*nGenes)")
                r("pdf_counter <- pdf_counter + 1")
                nModules = r.get("colorlevels")
                if len(nModules) > 2:
                    r("dissimME <- (1-t(cor(datME, method='p')))/2")
                    r("hclustdatME <- hclust(as.dist(dissimME), method='average')")
                    r("plot(hclustdatME, main='Eigengene Dendrogram')")
                else:
                    r("plot(0:10, type='n', axes=FALSE, bty='n', xlab='', ylab='', \
                        main='Eigengene Dendrogram\n\nInsufficient data to generate (n < 3)!') \
                    ")
                r("dev.off()")

                # Visualize the TOM Similarity Matrix
                r("load(paste(path, 'TOM-block.1.RData', sep='/'))")
                r("TOM.mat <- as.matrix(TOM)")
                r("diss = (1-TOM.mat)^6")
                r("diag(diss)=NA")
                r("pdf(paste(path, '/wgcna_temp', pdf_counter, '.pdf', sep=''), height=3+0.01*nGenes, width=5+0.01*nGenes)")
                r("pdf_counter <- pdf_counter + 1")
                r("TOMplot(diss, net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], main='TOM Dissimilarity Matrix')")
                r("dev.off()")

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                r("probes <- names(datExpr)")

                # Filter genes based on kME [only consider positive KMEs -- otherwise use abs(datKME)]
                r("modKME <- (datKME > minKME)")
                r("rownames(modKME) <- names(datExpr)")
                r("modKME <- modKME * 1")
                r("filter1 <- as.vector((rowSums(modKME)) > 0)")

                # Filter genes based on kME Rank
                r("rankKME <- apply(-abs(datKME), 2, rank, ties.method='min')")
                r("rownames(rankKME) <- names(datExpr)")
                r("rankKME <- (rankKME <= maxNGenes)")
                r("rankKME <- rankKME * 1")
                r("filter2 <- as.vector((rowSums(rankKME) > 0))")

                r("filtDF <- as.data.frame(cbind(filter1, filter2))")
                r("filtDF <- filtDF * 1")

                r("filter <- as.vector((rowSums(filtDF)==2))")

                # Create neural network
                r("load(paste(path, 'TOM-block.1.RData', sep='/'))")
                r("TOM.mat <- as.matrix(TOM)")

                r.assign("annotDF", moduleDF)
                r("annotDF$key <- paste(net$colors, '-', rownames(annotDF), sep='')")

                if button3 == 1:
                    if selectAll == 2:
                        moduleDF["name"] = moduleDF[['Kingdom', 'Phyla']].apply(lambda x: ';'.join(x), axis=1)
                        r.assign("annotDF$name", moduleDF['name'])
                    elif selectAll == 3:
                        moduleDF['name'] = moduleDF[['Phyla', 'Class']].apply(lambda x: ';'.join(x), axis=1)
                        r.assign("annotDF$name", moduleDF['name'])
                    elif selectAll == 4:
                        moduleDF['name'] = moduleDF[['Phyla', 'Class', 'Order']].apply(lambda x: ';'.join(x), axis=1)
                        r.assign("annotDF$name", moduleDF['name'])
                    elif selectAll == 5:
                        moduleDF['name'] = moduleDF[['Phyla', 'Class', 'Order', 'Family']].apply(lambda x: ';'.join(x), axis=1)
                        r.assign("annotDF$name", moduleDF['name'])
                    elif selectAll == 6:
                        moduleDF['name'] = moduleDF[['Phyla', 'Class', 'Order', 'Family', 'Genus']].apply(lambda x: ';'.join(x), axis=1)
                        r.assign("annotDF$name", moduleDF['name'])
                    elif selectAll == 7:
                        moduleDF['name'] = moduleDF[['Phyla', 'Class', 'Order', 'Family', 'Genus', 'Species']].apply(lambda x: ';'.join(x), axis=1)
                        r.assign("annotDF$name", moduleDF['name'])
                elif button3 == 2:
                    if keggAll == 1:
                        moduleDF["name"] = moduleDF['Level_1']
                        r.assign("annotDF$name", moduleDF['name'])
                    elif keggAll == 2:
                        moduleDF['name'] = moduleDF[['Level_1', 'Level_2']].apply(lambda x: ';'.join(x), axis=1)
                        r.assign("annotDF$name", moduleDF['name'])
                    elif keggAll == 3:
                        moduleDF['name'] = moduleDF[['Level_1', 'Level_2', 'Level_3']].apply(lambda x: ';'.join(x), axis=1)
                        r.assign("annotDF$name", moduleDF['name'])
                elif button3 == 3:
                    if nzAll == 1:
                        moduleDF["name"] = moduleDF['Level_1']
                        r.assign("annotDF$name", moduleDF['name'])
                    elif nzAll == 2:
                        moduleDF['name'] = moduleDF[['Level_1', 'Level_2']].apply(lambda x: ';'.join(x), axis=1)
                        r.assign("annotDF$name", moduleDF['name'])
                    elif nzAll == 3:
                        moduleDF['name'] = moduleDF[['Level_1', 'Level_2', 'Level_3']].apply(lambda x: ';'.join(x), axis=1)
                        r.assign("annotDF$name", moduleDF['name'])
                    elif nzAll >= 4:
                        moduleDF['name'] = moduleDF[['Level_1', 'Level_2', 'Level_3', 'Level_4']].apply(lambda x: ';'.join(x), axis=1)
                        r.assign("annotDF$name", moduleDF['name'])

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 3 of 6: WGCNA analysis...done!')
                database.queue.setBase(RID, 'Step 4 of 6: Creating network graph...')

                r("library('igraph')")

                ### now apply the kME filter created above
                r("modProbes <- probes[filter]")
                r("TOM.mat <- as.matrix(TOM)")
                r("modTOM <- TOM.mat[filter,filter]")
                r("modAnnot <- annotDF[filter,]")
                r("modColors <- mergedColors[filter]")

                ### create network in Cytoscape format
                r("cyt = exportNetworkToCytoscape(modTOM, \
                    weighted=TRUE, \
                    threshold=minEdge, \
                    nodeNames = modProbes, \
                    altNodeNames = modAnnot[,'key'], \
                    nodeAttr = modColors, \
                )")

                ### check to make sure there is data
                nGenes = r.get("nrow(cyt$nodeData)")
                if not nGenes > 0:
                    finalDict['text'] = error
                    finalDict['error'] = "All genes were removed from your network graph!\nPlease check your Network Graph settings."
                    res = simplejson.dumps(finalDict)
                    return HttpResponse(res, content_type='application/json')

                ### get node data
                r("edgeDF <- as.data.frame(cyt$edgeData)")
                r("edgeDF <- edgeDF[,(1:3)]")
                r("names(edgeDF) <- c('from', 'to', 'weight')")
                r("edgeDF$type <- 'mention'")

                ### get node data
                r("nodeDF <- as.data.frame(cyt$nodeData)")

                # Create graph
                r("nnet <- graph.data.frame(edgeDF, nodeDF, directed=T)")

                r("layout <- paste('layout', graphLayout, sep='.')")
                r("l <- do.call(layout, list(nnet))")
                r('path')
                r("pdf(paste(path, '/wgcna_temp', pdf_counter, '.pdf', sep=''), height=8, width=8)")
                r("pdf_counter <- pdf_counter + 1")
                r("p <- plot(nnet, \
                    main='Eigengene Network', \
                    edge.color='gray', \
                    edge.arrow.size=0, \
                    edge.width=edgeDF$weight / max(edgeDF$weight) * 5, \
                    edge.label=NA, \
                    edge.label.cex=0.7, \
                    edge.label.color='blue', \
                    vertex.size=14, \
                    vertex.label=paste(nodeDF$altName), \
                    vertex.label.cex=0.7, \
                    vertex.label.color='black', \
                    vertex.frame.color=adjustcolor('black', alpha=0.75), \
                    vertex.color=adjustcolor(nodeDF$nodeAttr, alpha=0.75), \
                    layout=l, \
                )")
                r("dev.off()")

                ### create Legend
                r("modAnnot <- modAnnot[modAnnot$rank_id %in% nodeDF$nodeName,]")
                r("legendDF <- data.frame(key=modAnnot$key, module=modAnnot$Module, rank_id=modAnnot$rank_id, modAnnot$name)")
                legendDF = r.get("legendDF")
                legend = legendDF.to_html(classes="table display")
                legend = legend.replace('border="1"', 'border="0"')
                finalDict['legend'] = str(legend)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 4 of 6: Creating network graph...done!')
                database.queue.setBase(RID, 'Step 5 of 6: Analyzing eigengenes based on chosen meta-variables...')

                #Finding modules that relate to a trait (quantitative variable)
                if quantFields:
                    if len(catFields_edit) == 1:
                        r.assign("meta", metaDF)
                        r.assign("quantFields", quantFields)
                        r.assign("catFields_edit", catFields_edit[0])
                        r("levels = levels(meta[,paste(catFields_edit)])")
                        levels = r.get("levels")
                        for level in levels:
                            r.assign("level", level)
                            r("filter <- meta[,paste(catFields_edit)]==level")
                            r("datME.mod <- as.matrix(datME)[filter,]")
                            r("meta.mod <- as.matrix(meta[,paste(quantFields)])[filter,]")
                            r("nSamples.mod <- nrow(meta.mod)")
                            r("moduleTraitCor <- cor(datME.mod, meta.mod, use='p')")
                            r("moduleTraitPvalue <- corPvalueFisher(moduleTraitCor, nSamples.mod)")
                            r("textMatrix <- paste(signif(moduleTraitCor, 2), '\n(', signif(moduleTraitPvalue, 1), ')', sep='')")
                            r("dim(textMatrix) <- dim(moduleTraitCor)")
                            r("nrow <- as.integer(dim(textMatrix)[1])")
                            r("ncol <- as.integer(dim(textMatrix)[2])")
                            r("pdf(paste(path, '/wgcna_temp', pdf_counter, '.pdf', sep=''), height=2+0.2*nrow, width=4+0.2*ncol)")
                            r("pdf_counter <- pdf_counter + 1")
                            r("par(mar=c(6, 8.8, 3, 2.2))")
                            r("labeledHeatmap( \
                                Matrix = moduleTraitCor, \
                                xLabels = names(meta[,paste(quantFields)]), \
                                yLabels = names(datME), \
                                ySymbols = names(datME), \
                                colorLabels = FALSE, \
                                colors = blueWhiteRed(50), \
                                textMatrix = textMatrix, \
                                setStdMargins = FALSE, \
                                cex.text = 0.5, \
                                zlim = c(-1,1), \
                                main = paste('Eigengene--Trait Correlations\n', level, sep=''), \
                                )")
                            r("dev.off()")

                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                            if stops[PID] == RID:
                                res = ''
                                return HttpResponse(res, content_type='application/json')
                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                        # Now All data
                        r("moduleTraitCor <- cor(datME, meta[,paste(quantFields)], use='p')")
                        r("moduleTraitPvalue <- corPvalueFisher(moduleTraitCor, nSamples)")
                        r("textMatrix <- paste(signif(moduleTraitCor, 2), '\n(', signif(moduleTraitPvalue, 1), ')', sep='')")
                        r("dim(textMatrix) <- dim(moduleTraitCor)")
                        r("nrow <- as.integer(dim(textMatrix)[1])")
                        r("ncol <- as.integer(dim(textMatrix)[2])")
                        r("pdf(paste(path, '/wgcna_temp', pdf_counter, '.pdf', sep=''), height=2+0.2*nrow, width=4+0.2*ncol)")
                        r("pdf_counter <- pdf_counter + 1")
                        r("par(mar=c(6, 8.8, 3, 2.2))")
                        r("labeledHeatmap( \
                            Matrix = moduleTraitCor, \
                            xLabels = names(meta[,paste(quantFields)]), \
                            yLabels = names(datME), \
                            ySymbols = names(datME), \
                            colorLabels = FALSE, \
                            colors = blueWhiteRed(50), \
                            textMatrix = textMatrix, \
                            setStdMargins = FALSE, \
                            cex.text = 0.5, \
                            zlim = c(-1,1), \
                            main = paste('Eigengene--Trait Correlations\nAll Data'), \
                        )")
                        r("dev.off()")

                    elif len(catFields_edit) > 1:
                        for index, row in metaDF.iterrows():
                           metaDF.loc[index, 'merge'] = "; ".join(row[catFields_edit])
                        r.assign("meta", metaDF)
                        r.assign("quantFields", quantFields)
                        r("levels = levels(meta$merge)")
                        levels = r.get("levels")
                        for level in levels:
                            r.assign("level", level)
                            r("filter <- meta[,'merge']==level")
                            r("datME.mod <- as.matrix(datME)[filter,]")
                            r("meta.mod <- as.matrix(meta[,paste(quantFields)])[filter,]")
                            r("nSamples.mod <- nrow(meta.mod)")
                            r("moduleTraitCor <- cor(datME.mod, meta.mod, use='p')")
                            r("moduleTraitPvalue <- corPvalueFisher(moduleTraitCor, nSamples.mod)")
                            r("textMatrix <- paste(signif(moduleTraitCor, 2), '\n(', signif(moduleTraitPvalue, 1), ')', sep='')")
                            r("dim(textMatrix) <- dim(moduleTraitCor)")
                            r("nrow <- as.integer(dim(textMatrix)[1])")
                            r("ncol <- as.integer(dim(textMatrix)[2])")
                            r("pdf(paste(path, '/wgcna_temp', pdf_counter, '.pdf', sep=''), height=2+0.2*nrow, width=4+0.2*ncol)")
                            r("pdf_counter <- pdf_counter + 1")
                            r("par(mar=c(6, 8.8, 3, 2.2))")
                            r("labeledHeatmap( \
                                Matrix = moduleTraitCor, \
                                xLabels = names(meta[,paste(quantFields)]), \
                                yLabels = names(datME), \
                                ySymbols = names(datME), \
                                colorLabels = FALSE, \
                                colors = blueWhiteRed(50), \
                                textMatrix = textMatrix, \
                                setStdMargins = FALSE, \
                                cex.text = 0.5, \
                                zlim = c(-1,1), \
                                main = paste('Eigengene--Trait Correlations\n', level, sep=''), \
                                )")
                            r("dev.off()")

                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                            if stops[PID] == RID:
                                res = ''
                                return HttpResponse(res, content_type='application/json')
                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                        # Now All data
                        r("moduleTraitCor <- cor(datME, meta[,paste(quantFields)], use='p')")
                        r("moduleTraitPvalue <- corPvalueFisher(moduleTraitCor, nSamples)")
                        r("moduleTraitCor")
                        r("moduleTraitPvalue")
                        r("textMatrix <- paste(signif(moduleTraitCor, 2), '\n(', signif(moduleTraitPvalue, 1), ')', sep='')")
                        r("dim(textMatrix) <- dim(moduleTraitCor)")
                        r("nrow <- as.integer(dim(textMatrix)[1])")
                        r("ncol <- as.integer(dim(textMatrix)[2])")
                        r("pdf(paste(path, '/wgcna_temp', pdf_counter, '.pdf', sep=''), height=2+0.2*nrow, width=4+0.2*ncol)")
                        r("pdf_counter <- pdf_counter + 1")
                        r("par(mar=c(6, 8.8, 3, 2.2))")
                        r("labeledHeatmap( \
                            Matrix = moduleTraitCor, \
                            xLabels = names(meta[,paste(quantFields)]), \
                            yLabels = names(datME), \
                            ySymbols = names(datME), \
                            colorLabels = FALSE, \
                            colors = blueWhiteRed(50), \
                            textMatrix = textMatrix, \
                            setStdMargins = FALSE, \
                            cex.text = 0.5, \
                            zlim = c(-1,1), \
                            main = paste('Eigengene--Trait Correlations\nAll Data'), \
                        )")
                        r("dev.off()")

                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                        if stops[PID] == RID:
                            res = ''
                            return HttpResponse(res, content_type='application/json')
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    else:
                        r.assign("meta", metaDF)
                        r.assign("quantFields", quantFields)
                        r("moduleTraitCor <- cor(datME, meta[,paste(quantFields)], use='p')")
                        r("moduleTraitPvalue <- corPvalueFisher(moduleTraitCor, nSamples)")
                        r("textMatrix <- paste(signif(moduleTraitCor, 2), '\n(', signif(moduleTraitPvalue, 1), ')', sep='')")
                        r("dim(textMatrix) <- dim(moduleTraitCor)")
                        r("nrow <- as.integer(dim(textMatrix)[1])")
                        r("ncol <- as.integer(dim(textMatrix)[2])")
                        r("pdf(paste(path, '/wgcna_temp', pdf_counter, '.pdf', sep=''), height=2+0.2*nrow, width=4+0.2*ncol)")
                        r("pdf_counter <- pdf_counter + 1")
                        r("par(mar=c(6, 8.8, 3, 2.2))")
                        r("labeledHeatmap( \
                            Matrix = moduleTraitCor, \
                            xLabels = names(meta[,paste(quantFields)]), \
                            yLabels = names(datME), \
                            ySymbols = names(datME), \
                            colorLabels = FALSE, \
                            colors = blueWhiteRed(50), \
                            textMatrix = textMatrix, \
                            setStdMargins = FALSE, \
                            cex.text = 0.5, \
                            zlim = c(-1,1), \
                            main = paste('Eigengene--Trait Correlations\nAll Data'), \
                        )")
                        r("dev.off()")

                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                        if stops[PID] == RID:
                            res = ''
                            return HttpResponse(res, content_type='application/json')
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                else:
                    if len(catFields_edit) == 1:
                        r.assign("meta", metaDF)
                        trtStr = ' * '.join(catFields_edit)
                        r.assign("trtStr", trtStr)

                        r("modules <- unique(mergedColors)")
                        r("modules <- paste('ME', modules, sep='.')")
                        r.assign("cols", catFields_edit)
                        r("aovDF$x <- aovDF[,paste(cols)]")
                        r("dat.m <- melt(aovDF, id.vars='x', measure.vars=modules)")

                        r("h <- as.integer(length(modules)/2)")
                        r("if (h <= 1) { h = 4 } else { h = 8 }")
                        r("w <- as.integer(length(modules)/2)")
                        r("if (w == 0) { w = 4 } else { w = 8 }")

                        r("pdf(paste(path, '/wgcna_temp', pdf_counter, '.pdf', sep=''), height=h, width=w)")
                        r("pdf_counter <- pdf_counter + 1")
                        r("p <- ggplot(dat.m, aes(x=x,  y=value, fill=x))")
                        r("p <- p + geom_boxplot(notch=FALSE)")
                        r("p <- p + theme(axis.text.x=element_text(angle=45, hjust=1, size=12), strip.text.x=element_text(face='bold'), strip.text.y=element_text(face='bold'))")
                        r("p <- p + labs(x='', y='Module Eigengenes', title='Eigengene Boxplots')")
                        r("p <- p + scale_y_continuous(limits=c(-1, 1))")
                        r("p <- p + scale_fill_brewer(palette='Set1', name='Legend')")
                        r("facet_multiple(plot=p, facets='variable', ncol=2, nrow=2, scale='free')")
                        r("dev.off()")

                        test = str(r(" for (i in 1:length(modules)) { \
                            module <- modules[i]; \
                            cmd <- paste('fit <- aov(', module, '~ ', trtStr, ', data=aovDF)', sep=''); \
                            eval(parse(text=cmd)); \
                            print(paste('Module:', module)); \
                            print(summary(fit)); \
                            print(TukeyHSD(fit)); \
                         }"))

                        result += '===============================================\n'
                        result += "ANOVA Results:\n\n"
                        lines = test.split('\n')
                        for line in lines[1:]:
                            result += line
                            result += '\n'
                        result += '===============================================\n'

                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                        if stops[PID] == RID:
                            res = ''
                            return HttpResponse(res, content_type='application/json')
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    elif len(catFields_edit) > 1:
                        r.assign("meta", metaDF)
                        trtStr = ' * '.join(catFields_edit)
                        r.assign("trtStr", trtStr)

                        r("modules <- unique(mergedColors)")
                        r("modules <- paste('ME', modules, sep='.')")
                        r.assign("cols", catFields_edit)
                        r("aovDF$x <- apply(aovDF[,paste(cols)], 1, paste, collapse='; ')")
                        r("dat.m <- melt(aovDF, id.vars='x', measure.vars=modules)")

                        r("h <- as.integer((length(modules)/2))")
                        r("if (h <= 1) {h = 4 } else { h = 8 }")
                        r("w <- as.integer((length(modules)/2))")
                        r("if (w == 0) {w = 4 } else { w = 8 }")

                        r("pdf(paste(path, '/wgcna_temp', pdf_counter, '.pdf', sep=''), height=h, width=w)")
                        r("pdf_counter <- pdf_counter + 1")
                        r("p <- ggplot(dat.m, aes(x=x,  y=value, fill=x))")
                        r("p <- p + geom_boxplot(notch=FALSE)")
                        r("p <- p + theme(axis.text.x = element_text(angle=45, hjust=1, size=12), strip.text.x=element_text(face='bold'), strip.text.y=element_text(face='bold'))")
                        r("p <- p + labs(x='', y='Module Eigengenes', title='Eigengene Boxplots')")
                        r("p <- p + scale_y_continuous(limits=c(-1, 1))")
                        r("p <- p + scale_fill_brewer(palette='Set1', name='Legend')")
                        r("facet_multiple(plot=p, facets='variable', ncol=2, nrow=2, scale='free')")

                        r("dev.off()")

                        test = str(r(" for (i in 1:length(modules)) { \
                            module <- modules[i]; \
                            cmd <- paste('fit <- aov(', module, '~ ', trtStr, ', data=aovDF)', sep=''); \
                            eval(parse(text=cmd)); \
                            print(paste('Module:', module)); \
                            print(summary(fit)); \
                            print(TukeyHSD(fit)); \
                         }"))

                        result += '===============================================\n'
                        result += "ANOVA Results:\n\n"
                        lines = test.split('\n')
                        for line in lines[1:]:
                            result += line
                            result += '\n'
                        result += '===============================================\n'

                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                        if stops[PID] == RID:
                            res = ''
                            return HttpResponse(res, content_type='application/json')
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 5 of 6: Analyzing eigengenes based on chosen meta-variables...done!')
                database.queue.setBase(RID, 'Step 6 of 6: Pooling pdf files for display...!')

                # Combining Pdf files
                finalFile = 'myPhyloDB/media/temp/wgcna/Rplots/' + str(RID) + '/wgcna_final.pdf'

                pdf_files = [f for f in os.listdir(path) if f.endswith("pdf")]
                pdf_files = natsorted(pdf_files, key=lambda y: y.lower())

                merger = PdfFileMerger()
                for filename in pdf_files:
                    merger.append(PdfFileReader(os.path.join(path, filename), 'rb'))

                merger.write(finalFile)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDict['text'] = result
                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                return HttpResponse(res, content_type='application/json')

    except Exception as e:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict['error'] = "There was an error during your analysis:\nError: " + str(e.message) + "\nTimestamp: " + str(datetime.datetime.now())
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')
