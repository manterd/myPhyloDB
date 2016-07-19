import ast
import datetime
import math
from django import db
from django.http import HttpResponse
from django_pandas.io import read_frame
import logging
from natsort import natsorted
import multiprocessing as mp
import numpy as np
import pandas as pd
from pyper import *
from PyPDF2 import PdfFileReader, PdfFileMerger
import simplejson
import shutil
import sys
import threading

from database.models import PICRUSt
from database.models import Phyla, Class, Order, Family, Genus, Species
from database.models import ko_lvl1, ko_lvl2, ko_lvl3, ko_entry
from database.models import nz_lvl1, nz_lvl2, nz_lvl3, nz_lvl4, nz_entry
from database.utils import multidict
import database.queue
from config import local_cfg


reload(sys)
sys.setdefaultencoding('utf8')

base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}

LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def statusWGCNA(request):
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
            if TimeDiff[RID] == 0:
                stage[RID] = 'Analysis has been placed in queue, there are '+str(database.queue.stat(RID))+' others in front of you.'
            else:
                stage[RID] = str(base[RID]) + '<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
        except:
            if TimeDiff[RID] == 0:
                stage[RID] = 'In queue'
            else:
                stage[RID] = '<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
        myDict['stage'] = stage[RID]
        json_data = simplejson.dumps(myDict, encoding="Latin-1")
        return HttpResponse(json_data, content_type='application/json')


def removeRIDWGCNA(RID):
    global base, stage, time1, time2, TimeDiff
    try:
        base.pop(RID, None)
        stage.pop(RID, None)
        time1.pop(RID, None)
        time2.pop(RID, None)
        TimeDiff.pop(RID, None)
        return True
    except:
        return False


def getWGCNA(request, stops, RID, PID):
    global base, stage, time1, TimeDiff
    try:
        while True:
            if request.is_ajax():
                allJson = request.body
                all = simplejson.loads(allJson)

                time1[RID] = time.time()
                base[RID] = 'Step 1 of 6: Selecting your chosen meta-variables...'

                myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
                path = str(myDir) + 'usr_norm_data.csv'

                with open(path, 'rb') as f:
                    savedDF = pd.read_csv(f, index_col=0, sep='\t')

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
                        catSampleIDs.extend(idDictCat[key])

                metaValsQuant = all['metaValsQuant']
                metaIDsQuant = all['metaIDsQuant']

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
                        quantSampleIDs.extend(idDictQuant[key])

                allSampleIDs = catSampleIDs + quantSampleIDs
                allFields = catFields_edit + quantFields

                # Removes samples (rows) that are not in our samplelist
                if allSampleIDs:
                    tempDF = savedDF.loc[savedDF['sampleid'].isin(allSampleIDs)]
                else:
                    tempDF = savedDF

                if metaDictCat:
                    for key in metaDictCat:
                        tempDF = tempDF.loc[tempDF[key].isin(metaDictCat[key])]

                if metaDictQuant:
                    for key in metaDictQuant:
                        valueList = [float(x) for x in metaDictQuant[key]]
                        tempDF = tempDF.loc[tempDF[key].isin(valueList)]

                if allFields:
                    wantedList = allFields + ['sampleid']
                else:
                    wantedList = ['sampleid']
                metaDF = tempDF[wantedList]

                if allFields:
                    result += 'Categorical variables selected by user: ' + ", ".join(catFields) + '\n'
                    result += 'Categorical variables removed from analysis (contains only 1 level): ' + ", ".join(removed) + '\n'
                    result += 'Quantitative variables selected by user: ' + ", ".join(quantFields) + '\n'
                else:
                    result += 'No meta variables were selected\n'
                result += '===============================================\n'

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

                base[RID] = 'Step 1 of 6: Selecting your chosen meta-variables...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 2 of 6: Selecting your chosen taxa or KEGG level...'

                DepVar = 1
                finalDF = pd.DataFrame()
                if button3 == 1:
                    DepVar = int(all["DepVar_taxa"])
                    finalDF = getTaxaDF(selectAll, savedDF, metaDF, DepVar, RID, stops, PID)

                if button3 == 2:
                    DepVar = int(all["DepVar_kegg"])
                    finalDF = getKeggDF(keggAll, savedDF, tempDF, DepVar, RID, stops, PID)

                if button3 == 3:
                    DepVar = int(all["DepVar_nz"])
                    finalDF = getNZDF(nzAll, savedDF, tempDF, DepVar, RID, stops, PID)

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

                meta_rDF = savedDF.drop_duplicates(subset='sampleid', take_last=True)

                # Removes samples (rows) that are not in our samplelist
                meta_rDF = meta_rDF.loc[meta_rDF['sampleid'].isin(allSampleIDs)]

                if metaDictCat:
                    for key in metaDictCat:
                        meta_rDF = meta_rDF.loc[meta_rDF[key].isin(metaDictCat[key])]

                if metaDictQuant:
                    for key in metaDictQuant:
                        valueList = [float(x) for x in metaDictQuant[key]]
                        meta_rDF = meta_rDF.loc[meta_rDF[key].isin(valueList)]

                wantedList = ['sampleid', 'sample_name'] + allFields
                meta_rDF = meta_rDF[wantedList]
                meta_rDF.set_index('sampleid', drop=True, inplace=True)

                base[RID] = 'Step 2 of 6: Selecting your chosen taxa...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 3 of 6: WGCNA analysis...'

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
                r.assign("sampleID", meta_rDF.index.values.tolist())
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
                nSamples, col = meta_rDF.shape
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

                zipped = []
                if button3 == 1:
                    zipped = getFullTaxonomy(selectAll, moduleDF['rank_id'])
                elif button3 == 2:
                    zipped = getFullKO(keggAll, moduleDF['rank_id'])
                elif button3 == 3:
                    zipped = getFullNZ(nzAll, moduleDF['rank_id'])

                if button3 == 1:
                    if selectAll == 2:
                        k, p = map(None, *zipped)
                        moduleDF.insert(1, 'Kingdom', k)
                        moduleDF.insert(2, 'Phyla', p)
                    elif selectAll == 3:
                        k, p, c = map(None, *zipped)
                        moduleDF.insert(1, 'Kingdom', k)
                        moduleDF.insert(2, 'Phyla', p)
                        moduleDF.insert(3, 'Class', c)
                    elif selectAll == 4:
                        k, p, c, o = map(None, *zipped)
                        moduleDF.insert(1, 'Kingdom', k)
                        moduleDF.insert(2, 'Phyla', p)
                        moduleDF.insert(3, 'Class', c)
                        moduleDF.insert(4, 'Order', o)
                    elif selectAll == 5:
                        k, p, c, o, f = map(None, *zipped)
                        moduleDF.insert(1, 'Kingdom', k)
                        moduleDF.insert(2, 'Phyla', p)
                        moduleDF.insert(3, 'Class', c)
                        moduleDF.insert(4, 'Order', o)
                        moduleDF.insert(5, 'Family', f)
                    elif selectAll == 6:
                        k, p, c, o, f, g = map(None, *zipped)
                        moduleDF.insert(1, 'Kingdom', k)
                        moduleDF.insert(2, 'Phyla', p)
                        moduleDF.insert(3, 'Class', c)
                        moduleDF.insert(4, 'Order', o)
                        moduleDF.insert(5, 'Family', f)
                        moduleDF.insert(6, 'Genus', g)
                    elif selectAll == 7:
                        k, p, c, o, f, g, s = map(None, *zipped)
                        moduleDF.insert(1, 'Kingdom', k)
                        moduleDF.insert(2, 'Phyla', p)
                        moduleDF.insert(3, 'Class', c)
                        moduleDF.insert(4, 'Order', o)
                        moduleDF.insert(5, 'Family', f)
                        moduleDF.insert(6, 'Genus', g)
                        moduleDF.insert(7, 'Species', s)
                if button3 == 2:
                    if keggAll == 1:
                        L1 = [x[0] for x in zipped]
                        moduleDF.insert(1, 'Level_1', L1)
                    if keggAll == 2:
                        L1, L2 = map(None, *zipped)
                        moduleDF.insert(1, 'Level_1', L1)
                        moduleDF.insert(2, 'Level_2', L2)
                    if keggAll == 3:
                        L1, L2, L3 = map(None, *zipped)
                        moduleDF.insert(1, 'Level_1', L1)
                        moduleDF.insert(2, 'Level_2', L2)
                        moduleDF.insert(3, 'Level_3', L3)
                if button3 == 3:
                    if nzAll == 1:
                        L1 = [x[0] for x in zipped]
                        moduleDF.insert(1, 'Level_1', L1)
                    if nzAll == 2:
                        L1, L2 = map(None, *zipped)
                        moduleDF.insert(1, 'Level_1', L1)
                        moduleDF.insert(2, 'Level_2', L2)
                    if nzAll == 3:
                        L1, L2, L3 = map(None, *zipped)
                        moduleDF.insert(1, 'Level_1', L1)
                        moduleDF.insert(2, 'Level_2', L2)
                        moduleDF.insert(3, 'Level_3', L3)
                    if nzAll == 4:
                        L1, L2, L3, L4 = map(None, *zipped)
                        moduleDF.insert(1, 'Level_1', L1)
                        moduleDF.insert(2, 'Level_2', L2)
                        moduleDF.insert(3, 'Level_3', L3)
                        moduleDF.insert(4, 'Level_4', L4)
                    if nzAll == 5:
                        L1, L2, L3, L4 = map(None, *zipped)
                        moduleDF.insert(1, 'Level_1', L1)
                        moduleDF.insert(2, 'Level_2', L2)
                        moduleDF.insert(3, 'Level_3', L3)
                        moduleDF.insert(4, 'Level_4', L4)
                    if nzAll == 6:
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
                    zipped = getFullTaxonomy(selectAll, kmeDF['rank_id'])
                elif button3 == 2:
                    zipped = getFullKO(keggAll, kmeDF['rank_id'])
                elif button3 == 3:
                    zipped = getFullNZ(nzAll, kmeDF['rank_id'])

                if button3 == 1:
                    if selectAll == 2:
                        k, p = map(None, *zipped)
                        kmeDF.insert(1, 'Kingdom', k)
                        kmeDF.insert(2, 'Phyla', p)
                    elif selectAll == 3:
                        k, p, c = map(None, *zipped)
                        kmeDF.insert(1, 'Kingdom', k)
                        kmeDF.insert(2, 'Phyla', p)
                        kmeDF.insert(3, 'Class', c)
                    elif selectAll == 4:
                        k, p, c, o = map(None, *zipped)
                        kmeDF.insert(1, 'Kingdom', k)
                        kmeDF.insert(2, 'Phyla', p)
                        kmeDF.insert(3, 'Class', c)
                        kmeDF.insert(4, 'Order', o)
                    elif selectAll == 5:
                        k, p, c, o, f = map(None, *zipped)
                        kmeDF.insert(1, 'Kingdom', k)
                        kmeDF.insert(2, 'Phyla', p)
                        kmeDF.insert(3, 'Class', c)
                        kmeDF.insert(4, 'Order', o)
                        kmeDF.insert(5, 'Family', f)
                    elif selectAll == 6:
                        k, p, c, o, f, g = map(None, *zipped)
                        kmeDF.insert(1, 'Kingdom', k)
                        kmeDF.insert(2, 'Phyla', p)
                        kmeDF.insert(3, 'Class', c)
                        kmeDF.insert(4, 'Order', o)
                        kmeDF.insert(5, 'Family', f)
                        kmeDF.insert(6, 'Genus', g)
                    elif selectAll == 7:
                        k, p, c, o, f, g, s = map(None, *zipped)
                        kmeDF.insert(1, 'Kingdom', k)
                        kmeDF.insert(2, 'Phyla', p)
                        kmeDF.insert(3, 'Class', c)
                        kmeDF.insert(4, 'Order', o)
                        kmeDF.insert(5, 'Family', f)
                        kmeDF.insert(6, 'Genus', g)
                        kmeDF.insert(7, 'Species', s)
                if button3 == 2:
                    if keggAll == 1:
                        L1 = [i[0] for i in zipped]
                        kmeDF.insert(1, 'Level_1', L1)
                    if keggAll == 2:
                        L1, L2 = map(None, *zipped)
                        kmeDF.insert(1, 'Level_1', L1)
                        kmeDF.insert(2, 'Level_2', L2)
                    if keggAll == 3:
                        L1, L2, L3 = map(None, *zipped)
                        kmeDF.insert(1, 'Level_1', L1)
                        kmeDF.insert(2, 'Level_2', L2)
                        kmeDF.insert(3, 'Level_3', L3)
                if button3 == 3:
                    if nzAll == 1:
                        L1 = [i[0] for i in zipped]
                        kmeDF.insert(1, 'Level_1', L1)
                    if nzAll == 2:
                        L1, L2 = map(None, *zipped)
                        kmeDF.insert(1, 'Level_1', L1)
                        kmeDF.insert(2, 'Level_2', L2)
                    if nzAll == 3:
                        L1, L2, L3 = map(None, *zipped)
                        kmeDF.insert(1, 'Level_1', L1)
                        kmeDF.insert(2, 'Level_2', L2)
                        kmeDF.insert(3, 'Level_3', L3)
                    if nzAll == 4:
                        L1, L2, L3, L4 = map(None, *zipped)
                        kmeDF.insert(1, 'Level_1', L1)
                        kmeDF.insert(2, 'Level_2', L2)
                        kmeDF.insert(3, 'Level_3', L3)
                        kmeDF.insert(4, 'Level_4', L4)
                    if nzAll == 5:
                        L1, L2, L3, L4 = map(None, *zipped)
                        kmeDF.insert(1, 'Level_1', L1)
                        kmeDF.insert(2, 'Level_2', L2)
                        kmeDF.insert(3, 'Level_3', L3)
                        kmeDF.insert(4, 'Level_4', L4)
                    if nzAll == 6:
                        L1, L2, L3, L4 = map(None, *zipped)
                        kmeDF.insert(1, 'Level_1', L1)
                        kmeDF.insert(2, 'Level_2', L2)
                        kmeDF.insert(3, 'Level_3', L3)
                        kmeDF.insert(4, 'Level_4', L4)

                    kmeDF.fillna(value='N/A', inplace=True)

                kme_table = kmeDF.to_html(classes="table display")
                kme_table = kme_table.replace('border="1"', 'border="0"')
                finalDict['kme_table'] = str(kme_table)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                # Create Eigengene Table
                r.assign("meta", meta_rDF)
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
                r("pdf(paste(path, '/wgcna_temp', pdf_counter, '.pdf', sep=''), height=6, width=5+0.025*nGenes)")
                r("pdf_counter <- pdf_counter + 1")
                r("plotDendroAndColors(net$dendrograms[[1]], \
                    cbind(mergedColors[net$blockGenes[[1]]], unmergedColors[net$blockGenes[[1]]]), \
                    c('Merged', 'Unmerged'), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, \
                    guideHang=0.05, main='Taxa Dendrogram')")
                r("dev.off()")

                #relating the modules
                r("pdf(paste(path, '/wgcna_temp', pdf_counter, '.pdf', sep=''), height=6, width=5+0.025*nGenes)")
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
                r("pdf(paste(path, '/wgcna_temp', pdf_counter, '.pdf', sep=''), height=3+0.025*nGenes, width=5+0.025*nGenes)")
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

                base[RID] = 'Step 3 of 6: WGCNA analysis...done!'
                base[RID] = 'Step 4 of 6: Creating network graph...'

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

                base[RID] = 'Step 4 of 6: Creating network graph...done!'
                base[RID] = 'Step 5 of 6: Analyzing eigengenes based on chosen meta-variables...'

                #Finding modules that relate to a trait (quantitative variable)
                if quantFields:
                    if len(catFields_edit) == 1:
                        r.assign("meta", meta_rDF)
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
                        for index, row in meta_rDF.iterrows():
                           meta_rDF.loc[index, 'merge'] = "; ".join(row[catFields_edit])
                        r.assign("meta", meta_rDF)
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
                        r.assign("meta", meta_rDF)
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
                        r.assign("meta", meta_rDF)
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
                        r.assign("meta", meta_rDF)
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

                base[RID] = 'Step 5 of 6: Analyzing eigengenes based on chosen meta-variables...done!'
                base[RID] = 'Step 6 of 6: Pooling pdf files for display...!'

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
                removeRIDWGCNA(RID)
                return HttpResponse(res, content_type='application/json')

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "Error with WGCNA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = simplejson.dumps(myDict)
            removeRIDWGCNA(RID)
            return HttpResponse(res, content_type='application/json')


def getTaxaDF(selectAll, savedDF, metaDF, DepVar, RID, stops, PID):
    global base
    try:
        base[RID] = 'Step 2 of 4: Selecting your chosen taxa or KEGG level...'
        taxaDF = pd.DataFrame(columns=['sampleid', 'rank', 'rank_id', 'rank_name', 'rel_abund', 'abund_16S'])

        if selectAll == 2:
            taxaDF = savedDF.loc[:, ['sampleid', 'phylaid', 'phylaName', 'rel_abund', 'abund_16S']]
            taxaDF.rename(columns={'phylaid': 'rank_id', 'phylaName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Phyla'
        elif selectAll == 3:
            taxaDF = savedDF.loc[:, ['sampleid', 'classid', 'className', 'rel_abund', 'abund_16S']]
            taxaDF.rename(columns={'classid': 'rank_id', 'className': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Class'
        elif selectAll == 4:
            taxaDF = savedDF.loc[:, ['sampleid', 'orderid', 'orderName', 'rel_abund', 'abund_16S']]
            taxaDF.rename(columns={'orderid': 'rank_id', 'orderName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Order'
        elif selectAll == 5:
            taxaDF = savedDF.loc[:, ['sampleid', 'familyid', 'familyName', 'rel_abund', 'abund_16S']]
            taxaDF.rename(columns={'familyid': 'rank_id', 'familyName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Family'
        elif selectAll == 6:
            taxaDF = savedDF.loc[:, ['sampleid', 'genusid', 'genusName', 'rel_abund', 'abund_16S']]
            taxaDF.rename(columns={'genusid': 'rank_id', 'genusName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Genus'
        elif selectAll == 7:
            taxaDF = savedDF.loc[:, ['sampleid', 'speciesid', 'speciesName', 'rel_abund', 'abund_16S']]
            taxaDF.rename(columns={'speciesid': 'rank_id', 'speciesName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Species'

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        taxaDF.drop('sampleid', axis=1, inplace=True)
        finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')

        wantedList = ['sampleid', 'rank', 'rank_name', 'rank_id']
        if DepVar == 1:
            finalDF = finalDF.groupby(wantedList)[['rel_abund']].sum()
        elif DepVar == 4:
            finalDF = finalDF.groupby(wantedList)[['abund_16S']].sum()

        finalDF.reset_index(drop=False, inplace=True)
        return finalDF

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with WGCNA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getKeggDF(keggAll, savedDF, tempDF, DepVar, RID, stops, PID):
    global base
    try:
        base[RID] = 'Step 2 of 8: Selecting your chosen taxa or KEGG level...'
        koDict = {}
        if keggAll == 1:
            keys = ko_lvl1.objects.using('picrust').values_list('ko_lvl1_id', flat=True)
            for key in keys:
                koList = ko_entry.objects.using('picrust').filter(ko_lvl1_id_id=key).values_list('ko_orthology', flat=True)
                if koList:
                    koDict[key] = koList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif keggAll == 2:
            keys = ko_lvl2.objects.using('picrust').values_list('ko_lvl2_id', flat=True)
            for key in keys:
                koList = ko_entry.objects.using('picrust').filter(ko_lvl2_id_id=key).values_list('ko_orthology', flat=True)
                if koList:
                    koDict[key] = koList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif keggAll == 3:
            keys = ko_lvl3.objects.using('picrust').values_list('ko_lvl3_id', flat=True)
            for key in keys:
                koList = ko_entry.objects.using('picrust').filter(ko_lvl3_id_id=key).values_list('ko_orthology', flat=True)
                if koList:
                    koDict[key] = koList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        if stops[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        # create sample and species lists based on meta data selection
        wanted = ['sampleid', 'speciesid', 'rel_abund', 'abund_16S']
        profileDF = tempDF.loc[:, wanted]
        profileDF.set_index('speciesid', inplace=True)

        # get PICRUSt data for species
        speciesList = pd.unique(profileDF.index.ravel().tolist())
        qs = PICRUSt.objects.using('picrust').filter(speciesid__in=speciesList)
        picrustDF = read_frame(qs, fieldnames=['speciesid__speciesid', 'geneCount'])
        picrustDF.set_index('speciesid__speciesid', inplace=True)

        path = 'myPhyloDB/media/temp/wgcna/' + str(RID)
        if not os.path.exists(path):
            os.makedirs(path)

        maxCPU = mp.cpu_count()
        maxCPU /= 3
        maxCPU = math.trunc(maxCPU)
        if maxCPU < 1:
            maxCPU = 1

        if os.name == 'nt':
            numcore = 1
            listDF = np.array_split(picrustDF, numcore)
            processes = [threading.Thread(target=sumStuff, args=(listDF[x], koDict, RID, x, stops, PID)) for x in xrange(numcore)]
        else:
            numcore = min(local_cfg.usr_numcore, maxCPU)
            listDF = np.array_split(picrustDF, numcore)
            processes = [threading.Thread(target=sumStuff, args=(listDF[x], koDict, RID, x, stops, PID)) for x in xrange(numcore)]

        for p in processes:
            p.start()

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        for p in processes:
            p.join()

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        levelList = []
        for key in koDict:
            levelList.append(key)

        picrustDF = pd.DataFrame()
        for i in xrange(numcore):
            path = 'myPhyloDB/media/temp/wgcna/'+str(RID)+'/file%d.temp' % i
            frame = pd.read_csv(path)
            picrustDF = picrustDF.append(frame, ignore_index=True)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        shutil.rmtree('myPhyloDB/media/temp/wgcna/'+str(RID))
        picrustDF.set_index('speciesid', inplace=True)
        picrustDF[picrustDF > 0] = 1

        # merge to get final gene counts for all selected samples
        taxaDF = pd.merge(profileDF, picrustDF, left_index=True, right_index=True, how='inner')

        for level in levelList:
            if DepVar == 1:
                taxaDF[level] = taxaDF['rel_abund'] * taxaDF[level]
            elif DepVar == 4:
                taxaDF[level] = taxaDF['abund_16S'] * taxaDF[level]

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        taxaDF = taxaDF.groupby('sampleid')[levelList].agg('sum')
        taxaDF.reset_index(drop=False, inplace=True)

        if DepVar == 1:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rel_abund')
        elif DepVar == 4:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund_16S')

        wanted = ['sampleid']
        metaDF = savedDF.loc[:, wanted]
        metaDF.set_index('sampleid', drop=True, inplace=True)
        grouped = metaDF.groupby(level=0)
        metaDF = grouped.last()

        taxaDF.set_index('sampleid', drop=True, inplace=True)
        finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')
        finalDF.reset_index(drop=False, inplace=True)

        finalDF['rank'] = ''
        finalDF['rank_name'] = ''
        for index, row in finalDF.iterrows():
            if ko_lvl1.objects.using('picrust').filter(ko_lvl1_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl1'
                finalDF.loc[index, 'rank_name'] = ko_lvl1.objects.using('picrust').get(ko_lvl1_id=row['rank_id']).ko_lvl1_name
            elif ko_lvl2.objects.using('picrust').filter(ko_lvl2_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl2'
                finalDF.loc[index, 'rank_name'] = ko_lvl2.objects.using('picrust').get(ko_lvl2_id=row['rank_id']).ko_lvl2_name
            elif ko_lvl3.objects.using('picrust').filter(ko_lvl3_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl3'
                finalDF.loc[index, 'rank_name'] = ko_lvl3.objects.using('picrust').get(ko_lvl3_id=row['rank_id']).ko_lvl3_name
            elif ko_entry.objects.using('picrust').filter(ko_lvl4_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl4'
                finalDF.loc[index, 'rank_name'] = ko_entry.objects.using('picrust').get(ko_lvl4_id=row['rank_id']).ko_name

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        return finalDF

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with WGCNA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getNZDF(nzAll, savedDF, tempDF, DepVar, RID, stops, PID):
    global base
    try:
        nzDict = {}
        if nzAll == 1:
            keys = nz_lvl1.objects.using('picrust').values_list('nz_lvl1_id', flat=True)
            for key in keys:
                nzList = nz_entry.objects.using('picrust').filter(nz_lvl1_id_id=key).values_list('nz_orthology', flat=True)
                if nzList:
                    nzDict[key] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 2:
            keys = nz_lvl2.objects.using('picrust').values_list('nz_lvl2_id', flat=True)
            for key in keys:
                nzList = nz_entry.objects.using('picrust').filter(nz_lvl2_id_id=key).values_list('nz_orthology', flat=True)
                if nzList:
                    nzDict[key] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 3:
            keys = nz_lvl3.objects.using('picrust').values_list('nz_lvl3_id', flat=True)
            for key in keys:
                nzList = nz_entry.objects.using('picrust').filter(nz_lvl3_id_id=key).values_list('nz_orthology', flat=True)
                if nzList:
                    nzDict[key] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 4:
            keys = nz_lvl4.objects.using('picrust').values_list('nz_lvl4_id', flat=True)
            for key in keys:
                nzList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=key).values_list('nz_orthology', flat=True)
                if nzList:
                    nzDict[key] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 5:
            # 1.18.6.1  nitrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.18.6.1  nitrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.3.3.11  pyrroloquinoline-quinone synthase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.3.3.11  pyrroloquinoline-quinone synthase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.4.99.5  glycine dehydrogenase (cyanide-forming)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.4.99.5  glycine dehydrogenase (cyanide-forming)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.1.1.76  (S,S)-butanediol dehydrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.1.1.76  (S,S)-butanediol dehydrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.14  chitinase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='3.2.1.14  chitinase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 4.1.1.74  indolepyruvate decarboxylase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='4.1.1.74  indolepyruvate decarboxylase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 6.3.2.39  aerobactin synthase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='6.3.2.39  aerobactin synthase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.4  cellulase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.4  cellulase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.91  cellulose 1,4-beta-cellobiosidase (non-reducing end)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.91  cellulose 1,4-beta-cellobiosidase (non-reducing end)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.21  beta-glucosidase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.21  beta-glucosidase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.8  endo-1,4-beta-xylanase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.8  endo-1,4-beta-xylanase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.37  xylan 1,4-beta-xylosidase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.37  xylan 1,4-beta-xylosidase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.5.1.4  amidase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.5.1.4  amidase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.5.1.5  urease
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.5.1.5  urease').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.1.3.1  alkaline phosphatase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.1.3.1  alkaline phosphatase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.1.3.2  acid phosphatase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.1.3.2  acid phosphatase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.1.6.1  arylsulfatase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.1.6.1  arylsulfatase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

        elif nzAll == 6:
            # 1.18.6.1  nitrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.18.6.1  nitrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.5.1.5  urease
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='3.5.1.5  urease').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.14.99.39  ammonia monooxygenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.14.99.39  ammonia monooxygenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.6  hydroxylamine dehydrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.6  hydroxylamine dehydrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.99.4  nitrate reductase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.99.4  nitrate reductase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.1  nitrite reductase (NO-forming)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.1  nitrite reductase (NO-forming)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.5  nitric oxide reductase (cytochrome c)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.5  nitric oxide reductase (cytochrome c)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.4  nitrous-oxide reductase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.4  nitrous-oxide reductase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

        # create sample and species lists based on meta data selection
        wanted = ['sampleid', 'speciesid', 'rel_abund', 'abund_16S']
        profileDF = tempDF.loc[:, wanted]
        profileDF.set_index('speciesid', inplace=True)

        # get PICRUSt data for species
        speciesList = pd.unique(profileDF.index.ravel().tolist())
        qs = PICRUSt.objects.using('picrust').filter(speciesid__in=speciesList)
        picrustDF = read_frame(qs, fieldnames=['speciesid__speciesid', 'geneCount'])
        picrustDF.set_index('speciesid__speciesid', inplace=True)

        path = 'myPhyloDB/media/temp/wgcna/' + str(RID)
        if not os.path.exists(path):
            os.makedirs(path)

        maxCPU = mp.cpu_count()
        maxCPU /= 3
        maxCPU = math.trunc(maxCPU)
        if maxCPU < 1:
            maxCPU = 1

        if os.name == 'nt':
            numcore = 1
            listDF = np.array_split(picrustDF, numcore)
            processes = [threading.Thread(target=sumStuff, args=(listDF[x], nzDict, RID, x, stops, PID)) for x in xrange(numcore)]
        else:
            numcore = min(local_cfg.usr_numcore, maxCPU)
            listDF = np.array_split(picrustDF, numcore)
            processes = [threading.Thread(target=sumStuff, args=(listDF[x], nzDict, RID, x, stops, PID)) for x in xrange(numcore)]

        for p in processes:
            p.start()

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        for p in processes:
            p.join()

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        levelList = []
        for key in nzDict:
            levelList.append(key)

        picrustDF = pd.DataFrame()
        for i in xrange(numcore):
            path = 'myPhyloDB/media/temp/wgcna/'+str(RID)+'/file%d.temp' % i
            frame = pd.read_csv(path)
            picrustDF = picrustDF.append(frame, ignore_index=True)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        shutil.rmtree('myPhyloDB/media/temp/wgcna/'+str(RID))
        picrustDF.set_index('speciesid', inplace=True)
        picrustDF[picrustDF > 0] = 1

        # merge to get final gene counts for all selected samples
        taxaDF = pd.merge(profileDF, picrustDF, left_index=True, right_index=True, how='inner')
        for level in levelList:
            if DepVar == 1:
                taxaDF[level] = taxaDF['rel_abund'] * taxaDF[level]
            elif DepVar == 4:
                taxaDF[level] = taxaDF['abund_16S'] * taxaDF[level]

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        taxaDF = taxaDF.groupby('sampleid')[levelList].agg('sum')
        taxaDF.reset_index(drop=False, inplace=True)

        if DepVar == 1:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rel_abund')
        elif DepVar == 4:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund_16S')

        wanted = ['sampleid']
        metaDF = savedDF.loc[:, wanted]
        metaDF.set_index('sampleid', drop=True, inplace=True)
        grouped = metaDF.groupby(level=0)
        metaDF = grouped.last()

        taxaDF.set_index('sampleid', drop=True, inplace=True)
        finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')
        finalDF.reset_index(drop=False, inplace=True)

        finalDF['rank'] = ''
        finalDF['rank_name'] = ''
        for index, row in finalDF.iterrows():
            if nz_lvl1.objects.using('picrust').filter(nz_lvl1_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl1'
                finalDF.loc[index, 'rank_name'] = nz_lvl1.objects.using('picrust').get(nz_lvl1_id=row['rank_id']).nz_lvl1_name
            elif nz_lvl2.objects.using('picrust').filter(nz_lvl2_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl2'
                finalDF.loc[index, 'rank_name'] = nz_lvl2.objects.using('picrust').get(nz_lvl2_id=row['rank_id']).nz_lvl2_name
            elif nz_lvl3.objects.using('picrust').filter(nz_lvl3_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl3'
                finalDF.loc[index, 'rank_name'] = nz_lvl3.objects.using('picrust').get(nz_lvl3_id=row['rank_id']).nz_lvl3_name
            elif nz_lvl4.objects.using('picrust').filter(nz_lvl4_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl4'
                finalDF.loc[index, 'rank_name'] = nz_lvl4.objects.using('picrust').get(nz_lvl4_id=row['rank_id']).nz_lvl4_name
            elif nz_entry.objects.using('picrust').filter(nz_lvl5_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl5'
                finalDF.loc[index, 'rank_name'] = nz_entry.objects.using('picrust').get(nz_lvl5_id=row['rank_id']).nz_name

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        return finalDF

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with WGCNA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def sumStuff(slice, koDict, RID, num, stops, PID):
    global base
    db.close_old_connections()

    f = open('myPhyloDB/media/temp/wgcna/'+str(RID)+'/file'+str(num)+".temp", 'w')

    keyList = []
    for key in koDict:
        keyList.append(key)

    f.write('speciesid,'+",".join(keyList)+'\n')

    for index, row in slice.iterrows():
        d = ast.literal_eval(row['geneCount'])

        f.write(str(index)+',')
        sumList = []
        for key in koDict:
            sum = 0.0
            myList = koDict[key]
            for k in myList:
                if k in d:
                    sum += d[k]
            sumList.append(sum)

        f.write(','.join(map(str, sumList)))
        f.write('\n')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    f.close()


def removeWGCNAFiles(request):
    if request.is_ajax():
        RID = request.GET["all"]

        rmPath = 'myPhyloDB/media/temp/wgcna/Rplots/%s' % RID
        shutil.rmtree(rmPath, ignore_errors=True)

        file = "myPhyloDB/media/temp/wgcna/" + str(RID) + ".pkl"
        if os.path.exists(file):
            os.remove(file)

        file = "myPhyloDB/media/temp/wgcna/" + str(RID) + ".csv"
        if os.path.exists(file):
            os.remove(file)

        return HttpResponse()


def getTabWGCNA(request):
    if request.is_ajax():
        RID = request.GET["all"]
        myDir = 'myPhyloDB/media/temp/wgcna/'
        fileName = str(myDir) + str(RID) + '.pkl'
        savedDF = pd.read_pickle(fileName)

        myDir = 'myPhyloDB/media/temp/wgcna/'
        fileName = str(myDir) + str(RID) + '.csv'
        savedDF.to_csv(fileName)

        myDict = {}
        myDir = '/myPhyloDB/media/temp/wgcna/'
        fileName = str(myDir) + str(RID) + '.csv'
        myDict['name'] = str(fileName)
        res = simplejson.dumps(myDict)

        return HttpResponse(res, content_type='application/json')


def getFullTaxonomy(level, id):
    record = []

    if level == 2:
        record = Phyla.objects.all().filter(phylaid__in=id).values_list('kingdomid_id__kingdomName', 'phylaName')
    elif level == 3:
        record = Class.objects.all().filter(classid__in=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'className')
    elif level == 4:
        record = Order.objects.all().filter(orderid__in=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderName')
    elif level == 5:
        record = Family.objects.all().filter(familyid__in=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyName')
    elif level == 6:
        record = Genus.objects.all().filter(genusid__in=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyid_id__familyName', 'genusName')
    elif level == 7:
        record = Species.objects.all().filter(speciesid__in=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyid_id__familyName', 'genusid_id__genusName', 'speciesName')

    return record


def getFullKO(level, id):
    record = []

    if level == 1:
        record = ko_lvl1.objects.using('picrust').all().filter(ko_lvl1_id__in=id).values_list('ko_lvl1_name')
    elif level == 2:
        record = ko_lvl2.objects.using('picrust').all().filter(ko_lvl2_id__in=id).values_list('ko_lvl1_id_id__ko_lvl1_name', 'ko_lvl2_name')
    elif level == 3:
        record = ko_lvl3.objects.using('picrust').all().filter(ko_lvl3_id__in=id).values_list('ko_lvl1_id_id__ko_lvl1_name', 'ko_lvl2_id_id__ko_lvl2_name', 'ko_lvl3_name')

    return record


def getFullNZ(level, id):
    record = []

    if level == 1:
        record = nz_lvl1.objects.using('picrust').all().filter(nz_lvl1_id__in=id).values_list('nz_lvl1_name')
    elif level == 2:
        record = nz_lvl2.objects.using('picrust').all().filter(nz_lvl2_id__in=id).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_name')
    elif level == 3:
        record = nz_lvl3.objects.using('picrust').all().filter(nz_lvl3_id__in=id).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_name')
    elif level == 4:
        record = nz_lvl4.objects.using('picrust').all().filter(nz_lvl4_id__in=id).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_id_id__nz_lvl3_name', 'nz_lvl4_name')
    elif level == 5:
        for item in id:
            if nz_lvl3.objects.using('picrust').all().filter(nz_lvl3_id=item).exists():
                qs = nz_lvl3.objects.using('picrust').all().filter(nz_lvl3_id=item).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_name')
                record.extend(qs)
            elif nz_lvl4.objects.using('picrust').all().filter(nz_lvl4_id=item).exists():
                qs = nz_lvl4.objects.using('picrust').all().filter(nz_lvl4_id=item).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_id_id__nz_lvl3_name', 'nz_lvl4_name')
                record.extend(qs)
            elif nz_entry.objects.using('picrust').all().filter(nz_lvl5_id=item).exists():
                qs = nz_entry.objects.using('picrust').all().filter(nz_lvl5_id=item).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_id_id__nz_lvl3_name', 'nz_lvl4_id_id__nz_lvl4_name')
                record.extend(qs)
    elif level == 6:
        record = nz_lvl4.objects.using('picrust').all().filter(nz_lvl4_id__in=id).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_id_id__nz_lvl3_name', 'nz_lvl4_name')

    return record
