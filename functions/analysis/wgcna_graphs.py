import datetime
from django.http import HttpResponse
import logging
import numpy as np
from natsort import natsorted
import pandas as pd
from pyper import *
from PyPDF2 import PdfFileReader, PdfFileMerger
import json
import sys

import functions


reload(sys)
sys.setdefaultencoding('utf8')

LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getWGCNA(request, stops, RID, PID):
    try:
        while True:
            if request.is_ajax():
                allJson = request.body.split('&')[0]
                all = json.loads(allJson)
                functions.setBase(RID, 'Step 1 of 7: Reading normalized data file...')

                functions.setBase(RID, 'Step 2 of 7: Selecting your chosen meta-variables...')
                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])

                # Select samples and meta-variables from savedDF
                metaValsCat = all['metaValsCat']
                metaIDsCat = all['metaIDsCat']
                metaValsQuant = all['metaValsQuant']
                metaIDsQuant = all['metaIDsQuant']

                treeType = int(all['treeType'])
                DepVar = int(all["DepVar"])

                # Create meta-variable DataFrame, final sample list, final category and quantitative field lists based on tree selections
                savedDF, metaDF, finalSampleIDs, catFields, remCatFields, quantFields, catValues, quantValues = functions.getMetaDF(request.user, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, DepVar)
                allFields = catFields + quantFields

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

                result += 'Categorical variables selected by user: ' + ", ".join(catFields + remCatFields) + '\n'
                result += 'Categorical variables not included in the statistical analysis (contains only 1 level): ' + ", ".join(remCatFields) + '\n'
                result += 'Quantitative variables selected by user: ' + ", ".join(quantFields) + '\n'
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

                functions.setBase(RID, 'Step 2 of 7: Selecting your chosen meta-variables...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 3 of 7: Selecting your chosen taxa or KEGG level...')

                # filter phylotypes based on user settings
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
                        filteredDF = functions.filterDF(savedDF, DepVar, selectAll, remUnclass, remZeroes, perZeroes, filterData, filterPer, filterMeth)
                    else:
                        filteredDF = savedDF.copy()

                    finalDF, missingList = functions.getTaxaDF(selectAll, '', filteredDF, metaDF, allFields, DepVar, RID, stops, PID)

                    if selectAll == 8:
                        result += '\nThe following PGPRs were not detected: ' + ", ".join(missingList) + '\n'
                        result += '===============================================\n'

                if treeType == 2:
                    finalDF, allDF = functions.getKeggDF(keggAll, '', savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)

                if treeType == 3:
                    finalDF, allDF = functions.getNZDF(nzAll, '', savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)

                if finalDF.empty:
                    error = "Selected taxa were not found in your selected samples."
                    myDict = {'error': error}
                    res = json.dumps(myDict)
                    return HttpResponse(res, content_type='application/json')

                # make sure column types are correct
                finalDF[catFields] = finalDF[catFields].astype(str)
                finalDF[quantFields] = finalDF[quantFields].astype(float)

                # transform Y, if requested
                transform = int(all["transform"])
                finalDF = functions.transformDF(transform, DepVar, finalDF)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/wgcna/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                count_rDF = pd.DataFrame()
                if DepVar == 0:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='abund')
                elif DepVar == 1:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='rel_abund')
                elif DepVar == 2:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='rich')
                elif DepVar == 3:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='diversity')
                elif DepVar == 4:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='abund_16S')

                count_rDF.fillna(0, inplace=True)

                functions.setBase(RID, 'Step 3 of 7: Selecting your chosen taxa...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 4 of 7: WGCNA analysis...')

                finalDict = {}
                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                functions.setBase(RID, 'Verifying R packages...missing packages are being installed')

                # R packages from biocLite
                r("list.of.packages <- c('XML')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                r("if (length(new.packages)) source('http://bioconductor.org/biocLite.R')")
                r("if (length(new.packages)) biocLite(new.packages)")

                # R packages from biocLite
                r("list.of.packages <- c('impute', 'preprocessCore', 'GO.db', 'AnnotationDbi')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                r("if (length(new.packages)) source('http://bioconductor.org/biocLite.R')")
                r("if (length(new.packages)) biocLite(new.packages)")

                # R packages from cran
                r("list.of.packages <- c('devtools', 'ggplot2', 'reshape2', 'WGCNA', 'cluster', 'igraph')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                print r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

                # R packages from github
                r("list.of.packages <- c('ggplus')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                print r("if (length(new.packages)) library('devtools')")
                print r("if (length(new.packages)) install_github('guiastrennec/ggplus')")

                functions.setBase(RID, 'Step 4 of 7: WGCNA analysis...')

                # Load R libraries
                print r("library(ggplot2)")
                print r("library(reshape2)")
                print r("library(WGCNA)")
                print r("library(cluster)")
                print r("library(ggplus)")
                print r("library(igraph)")

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
                metaDF.set_index('sampleid', inplace=True)
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
                    res = json.dumps(finalDict)
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
                    res = json.dumps(finalDict)
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

                if treeType == 1:
                    idList = functions.getFullTaxonomy(list(moduleDF.rank_id.unique()))
                    moduleDF['Taxonomy'] = moduleDF['rank_id'].map(idList)
                elif treeType == 2:
                    idList = functions.getFullKO(list(moduleDF.rank_id.unique()))
                    moduleDF['Pathway'] = moduleDF['rank_id'].map(idList)
                elif treeType == 3:
                    idList = functions.getFullNZ(list(moduleDF.rank_id.unique()))
                    moduleDF['Enzyme'] = moduleDF['rank_id'].map(idList)

                moduleDF.replace(to_replace='N/A', value=np.nan, inplace=True)
                moduleDF.dropna(axis=1, how='all', inplace=True)
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

                if treeType == 1:
                    idList = functions.getFullTaxonomy(list(kmeDF.rank_id.unique()))
                    kmeDF['Taxonomy'] = kmeDF['rank_id'].map(idList)
                elif treeType == 2:
                    idList = functions.getFullKO(list(kmeDF.rank_id.unique()))
                    kmeDF['Taxonomy'] = kmeDF['rank_id'].map(idList)
                elif treeType == 3:
                    idList = functions.getFullNZ(list(kmeDF.rank_id.unique()))
                    kmeDF['Taxonomy'] = kmeDF['rank_id'].map(idList)

                kmeDF.replace(to_replace='N/A', value=np.nan, inplace=True)
                kmeDF.dropna(axis=1, how='all', inplace=True)
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
                r("rankKME <- apply(-abs(datKME), 2, rank, ties.method='random')")
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

                if treeType == 1:
                    r.assign("annotDF$name", moduleDF['Taxonomy'])
                elif treeType == 2:
                    r.assign("annotDF$name", moduleDF['Pathway'])
                elif treeType == 3:
                    r.assign("annotDF$name", moduleDF['Enzyme'])

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 4 of 7: WGCNA analysis...done!')
                functions.setBase(RID, 'Step 5 of 7: Creating network graph...')

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
                    res = json.dumps(finalDict)
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
                r("pdf(paste(path, '/wgcna_temp', pdf_counter, '.pdf', sep=''), height=8, width=8)")
                r("pdf_counter <- pdf_counter + 1")
                r("p <- plot(nnet, \
                    main='Eigengene Network', \
                    edge.color='gray', \
                    edge.arrow.size=0, \
                    edge.width=edgeDF$weight / max(edgeDF$weight) * 10, \
                    edge.label=NA, \
                    edge.label.cex=0.7, \
                    edge.label.color='blue', \
                    vertex.size=14, \
                    vertex.label=paste(nodeDF$altName), \
                    vertex.label.cex=0.7, \
                    vertex.label.color='black', \
                    vertex.frame.color=adjustcolor('gray', alpha=0.95), \
                    vertex.color=adjustcolor(nodeDF$nodeAttr, alpha=0.95), \
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

                functions.setBase(RID, 'Step 5 of 7: Creating network graph...done!')
                functions.setBase(RID, 'Step 6 of 7: Analyzing eigengenes based on chosen meta-variables...')

                #Finding modules that relate to a trait (quantitative variable)
                if quantFields:
                    if len(catFields) == 1:
                        r.assign("meta", metaDF)
                        r.assign("quantFields", quantFields)
                        r.assign("catFields", catFields[0])
                        r("levels = levels(meta[,paste(catFields)])")
                        levels = r.get("levels")
                        for level in levels:
                            r.assign("level", level)
                            r("filter <- meta[,paste(catFields)]==level")
                            r("datME.mod <- as.matrix(datME)[filter,]")
                            r("meta.mod <- as.matrix(meta[,paste(quantFields)])[filter,]")
                            r("nSamples.mod <- nrow(meta.mod)")
                            r("moduleTraitCor <- cor(datME.mod, meta.mod, use='p')")
                            r("moduleTraitPvalue <- corPvalueFisher(moduleTraitCor, nSamples.mod)")
                            r("textMatrix <- paste(signif(moduleTraitCor, 3), '\n(', signif(moduleTraitPvalue, 3), ')', sep='')")
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
                        r("textMatrix <- paste(signif(moduleTraitCor, 3), '\n(', signif(moduleTraitPvalue, 3), ')', sep='')")
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

                    elif len(catFields) > 1:
                        for index, row in metaDF.iterrows():
                           metaDF.loc[index, 'merge'] = "; ".join(row[catFields])
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
                            r("textMatrix <- paste(signif(moduleTraitCor, 3), '\n(', signif(moduleTraitPvalue, 3), ')', sep='')")
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
                        r("textMatrix <- paste(signif(moduleTraitCor, 3), '\n(', signif(moduleTraitPvalue, 3), ')', sep='')")
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
                        r("textMatrix <- paste(signif(moduleTraitCor, 3), '\n(', signif(moduleTraitPvalue, 3), ')', sep='')")
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
                    if len(catFields) == 1:
                        r.assign("meta", metaDF)
                        trtStr = ' * '.join(catFields)
                        r.assign("trtStr", trtStr)

                        r("modules <- unique(mergedColors)")
                        r("modules <- paste('ME', modules, sep='.')")
                        r.assign("cols", catFields)
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

                    elif len(catFields) > 1:
                        r.assign("meta", metaDF)
                        trtStr = ' * '.join(catFields)
                        r.assign("trtStr", trtStr)

                        r("modules <- unique(mergedColors)")
                        r("modules <- paste('ME', modules, sep='.')")
                        r.assign("cols", catFields)
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

                functions.setBase(RID, 'Step 6 of 7: Analyzing eigengenes based on chosen meta-variables...done!')
                functions.setBase(RID, 'Step 7 of 7: Pooling pdf files for display...!')

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
