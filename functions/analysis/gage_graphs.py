import datetime
from django.http import HttpResponse
import logging
from natsort import natsorted
import pandas as pd
from PyPDF2 import PdfFileReader, PdfFileMerger
from pyper import *
import json

from database.models import ko_entry

import functions


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getGAGE(request, stops, RID, PID):
    try:
        while True:
            if request.is_ajax():
                # Get variables from web page
                allJson = request.body.split('&')[0]
                all = json.loads(allJson)
                functions.setBase(RID, 'Step 1 of 6: Reading normalized data file...')

                # Select samples and meta-variables from savedDF
                functions.setBase(RID, 'Step 2 of 6: Selecting your chosen meta-variables...')
                metaValsCat = all['metaValsCat']
                metaIDsCat = all['metaIDsCat']
                metaValsQuant = []
                metaIDsQuant = []
                
                DepVar = int(all["DepVar"])

                # Create meta-variable DataFrame, final sample list, final category and quantitative field lists based on tree selections
                savedDF, metaDF, finalSampleIDs, catFields, remCatFields, quantFields, catValues, quantValues = functions.getMetaDF(request.user, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, DepVar, levelDep=True)

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
                
                functions.setBase(RID, 'Step 2 of 6: Selecting your chosen meta-variables...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 3 of 6: Mapping phylotypes to KEGG pathways...')

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                functions.setBase(RID, 'Verifying R packages...missing packages are being installed')

                # R packages from biocLite
                r("list.of.packages <- c('gage', 'edgeR', 'pathview')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                r("if (length(new.packages)) source('http://bioconductor.org/biocLite.R')")
                print r("if (length(new.packages)) biocLite(new.packages, type='source', suppressUpdate=T, dependencies=T)")

                # R packages from cran
                r("list.of.packages <- c('png', 'grid')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                print r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

                functions.setBase(RID, 'Step 3 of 6: Mapping phylotypes to KEGG pathways...')

                print r("library(gage)")
                print r("library(edgeR)")
                print r("library(pathview)")
                print r("library(png)")
                print r("library(grid)")

                keggString = all["kegg"]
                keggDict = json.JSONDecoder(object_pairs_hook=functions.multidict).decode(keggString)
                nameList = []
                for value in keggDict.itervalues():
                    if isinstance(value, list):
                        nameList.extend(value)
                    else:
                        nameList.append(value)

                # Enable this only if you want to update gage data and pathways
                '''
                r("list.of.packages <- c('gageData')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                r("if (length(new.packages)) source('http://bioconductor.org/biocLite.R')")
                print r("if (length(new.packages)) biocLite(new.packages)")

                r('library(gageData)')
                r("data(kegg.sets.ko)")
                r("save(kegg.sets.ko, file='myPhyloDB/media/kegg/kegg.sets.ko.RData')")
                r("for (name in names(kegg.sets.ko)) { \
                    id = substr(name, 3, 7); \
                    download.kegg(pathway.id=id, species='ko', kegg.dir='myPhyloDB/media/kegg/pathways', file.type=c('xml', 'png')) \
                    } \
                ")
                '''

                r("load('myPhyloDB/media/kegg/kegg.sets.ko.RData')")

                keggDict = {}
                r("selPaths <- vector()")
                for i in nameList:
                    pathStr = i.split('[PATH:')[1].split(']')[0]
                    r.assign("pathStr", pathStr)
                    r("selPath <- kegg.sets.ko[grepl(paste(pathStr), names(kegg.sets.ko))]")
                    key = r.get("names(selPath)")
                    value = r.get("selPath$ko")
                    keggDict[key] = value.tolist()
                    r("selPaths <- append(selPaths, names(selPath))")

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                keggAll = 4
                mapTaxa = 'no'
                finalDF, junk = functions.getKeggDF(keggAll, '', savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)

                # make sure column types are correct
                finalDF[catFields] = finalDF[catFields].astype(str)

                functions.setBase(RID, 'Step 4 of 6: Performing GAGE analysis...')

                # save location info to session
                myDir = 'myPhyloDB/media/temp/gage/'
                if not os.path.exists(myDir):
                    os.makedirs(myDir)

                path = str(myDir) + str(RID) + '.biom'
                functions.imploding_panda(path, 2, DepVar, finalSampleIDs, metaDF, finalDF)

                count_rDF = pd.DataFrame()
                if DepVar == 0:
                    count_rDF = finalDF.pivot(index='rank_id', columns='sampleid', values='abund')
                elif DepVar == 4:
                    count_rDF = finalDF.pivot(index='rank_id', columns='sampleid', values='abund_16S')

                # need to change rank_id to kegg orthologies for gage analysis
                count_rDF.reset_index(drop=False, inplace=True)
                count_rDF.rename(columns={'index': 'rank_id'}, inplace=True)
                idList = count_rDF.rank_id.tolist()
                idDict = {}
                for id in idList:
                    entry = ko_entry.objects.using('picrust').get(ko_lvl4_id=id).ko_orthology
                    idDict[id] = entry
                count_rDF['ko'] = count_rDF['rank_id'].map(idDict)
                count_rDF.drop('rank_id', axis=1, inplace=True)
                count_rDF.drop_duplicates(take_last=True, inplace=True)  # remove dups - KOs mapped to multiple pathways
                count_rDF.set_index('ko', drop=True, inplace=True)

                # Create combined metadata column
                if len(catFields) > 1:
                    for index, row in metaDF.iterrows():
                        metaDF.loc[index, 'merge'] = ".".join(row[catFields])
                else:
                    metaDF.loc[:, 'merge'] = metaDF.loc[:, catFields[0]]

                wantedList = ['merge', 'sampleid', 'sample_name']
                metaDF = metaDF.loc[:, wantedList]

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDict = {}
                metaDF.sort(columns='sampleid', inplace=True)
                r.assign("metaDF", metaDF)
                r("trt <- factor(metaDF$merge)")

                r.assign("count", count_rDF)
                r.assign("sampleIDs", count_rDF.columns.values.tolist())
                r("names(count) <- sampleIDs")

                r('e <- DGEList(counts=count)')
                r('e <- calcNormFactors(e, method="none")')

                r('design <- model.matrix(~ 0 + trt)')
                r('trtLevels <- levels(trt)')
                r('colnames(design) <- trtLevels')

                r('e <- estimateGLMCommonDisp(e, design)')
                r('e <- estimateGLMTrendedDisp(e, design)')
                r('e <- estimateGLMTagwiseDisp(e, design)')
                r('fit <- glmFit(e, design)')
                fit = r.get('fit')

                if not fit:
                    error = "edgeR failed!\nUsually this is caused by one or more taxa having a negative disperion.\nTry filtering your data to remove problematic taxa (e.g. remove phylotypes with 50% or more zeros)."
                    myDict = {'error': error}
                    res = json.dumps(myDict)
                    return HttpResponse(res, content_type='application/json')

                if DepVar == 0:
                    result += 'Dependent Variable: Abundance' + '\n'
                elif DepVar == 4:
                    result += 'Dependent Variable: Total Abundance' + '\n'
                result += '\n===============================================\n\n\n'

                levels = list(set(metaDF['merge'].tolist()))
                levels = natsorted(levels, key=lambda y: y.lower())

                r("pdf_counter <- 1")

                path = os.path.join('myPhyloDB', 'media', 'temp', 'gage', 'Rplots', RID)
                if not os.path.exists(path):
                    os.makedirs(path)

                r.assign("path", path)
                r("setwd(path)")
                r("options(width=5000)")
                r.assign("RID", RID)

                gageDF = pd.DataFrame(columns=['comparison', 'pathway', ' p.geomean ', ' stat.mean ', ' p.val ', ' q.val ', ' set.size '])
                diffDF = pd.DataFrame(columns=['comparison', 'kegg', ' baseMean ', ' baseMeanA ', ' baseMeanB ', ' logFC ', ' logCPM ', ' LR ', ' pval ', ' FDR '])

                mergeList = metaDF['merge'].tolist()
                mergeSet = list(set(mergeList))
                for i in xrange(len(levels)-1):
                    for j in xrange(i+1, len(levels)):
                        trt1 = levels[i]
                        trt2 = levels[j]
                        r.assign("trt1", trt1)
                        r.assign("trt2", trt2)

                        r('contVec <- sprintf("%s-%s", trt1, trt2)')
                        r('cont.matrix= makeContrasts(contVec, levels=design)')
                        r('lrt <- glmLRT(fit, contrast=cont.matrix)')
                        r("res <- as.data.frame(topTags(lrt, n=nrow(lrt$table)))")
                        r('res <- res[ order(row.names(res)), ]')
                        taxaIDs = r.get("row.names(res)")

                        r("change <- -res$logFC")
                        r("names(change) <- row.names(res)")

                        baseMean = count_rDF.mean(axis=1)
                        baseMean = baseMean.loc[baseMean.index.isin(taxaIDs)]

                        listA = metaDF[metaDF['merge'] == mergeSet[i]].sampleid.tolist()
                        baseMeanA = count_rDF[listA].mean(axis=1)
                        baseMeanA = baseMeanA.loc[baseMeanA.index.isin(taxaIDs)]

                        listB = metaDF[metaDF['merge'] == mergeSet[j]].sampleid.tolist()
                        baseMeanB = count_rDF[listB].mean(axis=1)
                        baseMeanB = baseMeanB.loc[baseMeanB.index.isin(taxaIDs)]

                        r.assign("baseMean", baseMean)
                        r.assign("baseMeanA", baseMeanA)
                        r.assign("baseMeanB", baseMeanB)

                        r('baseMean <- baseMean[ order(as.numeric(row.names(baseMean))), ]')
                        r('baseMeanA <- baseMeanA[ order(as.numeric(row.names(baseMeanA))), ]')
                        r('baseMeanB <- baseMeanB[ order(as.numeric(row.names(baseMeanB))), ]')

                        # output DiffAbund to DataTable
                        r("df <- data.frame(kegg=row.names(res), baseMean=baseMean, baseMeanA=baseMeanA, \
                            baseMeanB=baseMeanB, logFC=-res$logFC, logCPM=res$logCPM, \
                            LR=res$LR, pval=res$PValue, FDR=res$FDR) \
                        ")

                        nbinom_res = r.get("df")
                        nbinom_res.fillna(value=1.0, inplace=True)

                        if nbinom_res is None:
                            myDict = {'error': "edgeR failed!\nPlease try a different data combination."}
                            res = json.dumps(myDict)
                            return HttpResponse(res, content_type='application/json')

                        comparison = str(trt1) + ' vs. ' + str(trt2)
                        nbinom_res.insert(0, 'comparison', comparison)
                        diffDF = diffDF.append(nbinom_res, ignore_index=True)

                        ### GAGE analysis on all pathways...
                        r("gage.res <- gage(change, gsets=kegg.sets.ko, species='ko', same.dir=FALSE)")
                        r("df2 <- data.frame(pathway=row.names(gage.res$greater), p.geomean=gage.res$greater[, 1], stat.mean=gage.res$greater[, 2], \
                            p.val=gage.res$greater[, 3], q.val=gage.res$greater[, 4], \
                            set.size=gage.res$greater[, 5])")

                        compDF = r.get("df2")
                        compDF.insert(0, 'comparison', comparison)
                        gageDF = gageDF.append(compDF, ignore_index=True)

                        ### Get data way for pathview
                        # merge sign and sig to get vector (1=sig. positive, 0=not sig., -1=sig. negative)
                        r("binary <- change / abs(change)")
                        r("sig <- as.vector((res$PValue <= 0.05))")
                        r("sig <- sig * 1")

                        r("sig <- sig * binary")
                        r("names(sig) <- row.names(res)")

                        for key in keggDict.iterkeys():
                            r.assign("pathway", key)
                            r("pid <- substr(pathway, start=1, stop=7)")
                            r("pv <- pathview(gene.data=sig, pathway.id=pid, species='ko', kegg.dir='../../../../kegg/pathways', \
                                kegg.native=T,  multi.state=F, same.layer=T, low='red', mid='gray', high='green')")

                            # convert to pdf
                            r("pdf(paste('gage_temp', pdf_counter, '.pdf', sep=''))")
                            r("plot.new()")
                            r("pngRaster <- readPNG(paste(pid, 'pathview.png', sep='.'))")
                            r("grid.raster(pngRaster)")
                            r("mtext(paste(trt1, ' vs ', trt2, sep=''), side=3, line=3, col='blue')")
                            r("dev.off()")
                            r("pdf_counter <- pdf_counter + 1")

                        functions.setBase(RID, 'Step 4 of 6: Performing GAGE Analysis...\nComparison: ' + str(trt1) + ' vs ' + str(trt2) + ' is done!')

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

                functions.setBase(RID, 'Step 5 of 6: Pooling pdf files for display...')

                # Combining Pdf files
                finalFile = 'myPhyloDB/media/temp/gage/Rplots/' + str(RID) + '/gage_final.pdf'
                pdf_files = [f for f in os.listdir(path) if f.endswith("pdf")]
                if pdf_files:
                    pdf_files = natsorted(pdf_files, key=lambda y: y.lower())

                    merger = PdfFileMerger()
                    for filename in pdf_files:
                        merger.append(PdfFileReader(os.path.join(path, filename), 'rb'))

                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                        if stops[PID] == RID:
                            res = ''
                            return HttpResponse(res, content_type='application/json')
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    merger.write(finalFile)

                functions.setBase(RID, 'Step 6 of 6: Formatting result tables...this may take several minutes')
                #Export tables to html
                gage_table = gageDF.to_html(classes="table display")
                gage_table = gage_table.replace('border="1"', 'border="0"')
                finalDict['gage_table'] = str(gage_table)

                diff_table = diffDF.to_html(classes="table display")
                diff_table = diff_table.replace('border="1"', 'border="0"')
                finalDict['diff_table'] = str(diff_table)

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
