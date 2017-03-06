import collections
import datetime
from django.http import HttpResponse
import logging
import pandas as pd
from pyper import *
import json

import functions


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getCorr(request, stops, RID, PID):
    try:
        while True:
                allJson = request.body.split('&')[0]
                all = json.loads(allJson)
                functions.setBase(RID, 'Step 1 of 6: Reading normalized data file...')

                functions.setBase(RID, 'Step 2 of 6: Selecting your chosen meta-variables...')
                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])

                # Select samples and meta-variables from savedDF
                metaValsCat = []
                metaIDsCat = []
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

                functions.setBase(RID, 'Step 2 of 6: Selecting your chosen meta-variables...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 3 of 6: Selecting your chosen taxa or KEGG level...')

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
                finalDF[quantFields] = finalDF[quantFields].astype(float)

                # transform Y, if requested
                transform = int(all["transform"])
                finalDF = functions.transformDF(transform, DepVar, finalDF)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/corr/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                functions.setBase(RID, 'Step 3 of 6: Selecting your chosen taxa or KEGG level...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 4 of 6: Calculating Correlations Matrix...')

                if DepVar == 0:
                    result += 'Dependent Variable: Abundance' + '\n'
                elif DepVar == 1:
                    result += 'Dependent Variable: Relative Abundance' + '\n'
                elif DepVar == 2:
                    result += 'Dependent Variable: OTU Richness' + '\n'
                elif DepVar == 3:
                    result += 'Dependent Variable: OTU Diversity' + '\n'
                elif DepVar == 4:
                    result += 'Dependent Variable: Total Abundance' + '\n'
                result += '\n===============================================\n'

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

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                functions.setBase(RID, 'Verifying R packages...missing packages are being installed')

                # R packages from cran
                r("list.of.packages <- c('corrplot', 'RColorBrewer', 'WGCNA')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                print r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

                functions.setBase(RID, 'Step 4 of 6: Calculating correlation matirx...')

                print r("library(corrplot)")
                print r("library(WGCNA)")
                print r("library(RColorBrewer)")

                idList = count_rDF.columns.values.tolist()
                namesList = []
                if treeType == 1:
                    namesDict = functions.getFullTaxonomy(idList)
                    od = collections.OrderedDict(sorted(namesDict.items()))
                    namesList = od.values()
                    namesList = [i.split('|')[-1] for i in namesList]
                elif treeType == 2:
                    namesDict = functions.getFullKO(idList)
                    od = collections.OrderedDict(sorted(namesDict.items()))
                    namesList = od.values()
                    namesList = [i.split('|')[-1] for i in namesList]
                elif treeType == 3:
                    if nzAll < 5:
                        namesDict = functions.getFullNZ(idList)
                        od = collections.OrderedDict(sorted(namesDict.items()))
                        namesList = od.values()
                        namesList = [i.split('|')[-1] for i in namesList]
                    else:
                        namesList = [i.split(':', 1)[0] for i in idList]

                count_rDF.sort_index(axis=0, inplace=True)
                metaDF.sort('sampleid', inplace=True)

                r.assign("X", count_rDF)
                r('X <- X * 1.0')
                r.assign("names", namesList)
                r("colnames(X) <- names")

                if quantFields:
                    r.assign("Y", metaDF[quantFields])
                    r('M <- cor(as.matrix(Y), as.matrix(X))')
                    r('M[is.na(M)] <- 0')
                    r("M.p <- corPvalueFisher(M, nrow(X))")
                else:
                    r('M <- cor(X)')
                    r('M[is.na(M)] <- 0')
                    r("M.p <- corPvalueFisher(M, nrow(X))")

                path = "myPhyloDB/media/temp/corr/Rplots/" + str(RID) + ".corr.pdf"
                if os.path.exists(path):
                    os.remove(path)

                if not os.path.exists('myPhyloDB/media/temp/corr/Rplots'):
                    os.makedirs('myPhyloDB/media/temp/corr/Rplots')

                row, col = count_rDF.shape
                width = 2 + col*0.15
                if quantFields:
                    height = 4 + len(quantFields)*0.1
                else:
                    height = width
                file = "pdf('myPhyloDB/media/temp/corr/Rplots/" + str(RID) + ".corr.pdf', height=" + str(height) + ", width=" + str(width) + ")"
                r.assign("cmd", file)
                r("eval(parse(text=cmd))")

                if quantFields:
                    r('corrplot(M, p.mat=M.p, method="pie", sig.level=0.05, insig="blank", cl.length=5, tl.cex=0.7, cl.cex=0.5)')
                else:
                    r('corrplot(M, p.mat=M.p, method="pie", sig.level=0.05, insig="blank", cl.length=5, tl.cex=0.7, cl.cex=0.5, order="hclust")')

                r("dev.off()")

                functions.setBase(RID, 'Step 4 of 6: Calculating correlation matrix...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 5 of 6: Formatting correlation coefficient table...')

                finalDict = {}
                dat = r.get("M")
                if quantFields:
                    coeffsDF = pd.DataFrame(dat, columns=namesList, index=quantFields)
                else:
                    coeffsDF = pd.DataFrame(dat, columns=namesList, index=namesList)
                res_table = coeffsDF.to_html(classes="table display")
                res_table = res_table.replace('border="1"', 'border="0"')
                finalDict['coeff_table'] = str(res_table)

                dat = r.get("M.p")
                if quantFields:
                    coeffsDF = pd.DataFrame(dat, columns=namesList, index=quantFields)
                else:
                    coeffsDF = pd.DataFrame(dat, columns=namesList, index=namesList)
                res_table = coeffsDF.to_html(classes="table display")
                res_table = res_table.replace('border="1"', 'border="0"')
                finalDict['p_table'] = str(res_table)

                finalDict['text'] = result

                functions.setBase(RID, 'Step 6 of 6: Formatting graph data for display...done!')

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
