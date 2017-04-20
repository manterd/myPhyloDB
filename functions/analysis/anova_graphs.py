import datetime
from django.http import HttpResponse
import logging
import json
import numpy as np
import pandas as pd
from pyper import *
from scipy import stats

from database.models import Sample

import functions


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getCatUnivData(request, RID, stops, PID):
    try:
        while True:
            if request.is_ajax():
                # Get variables from web page
                allJson = request.body.split('&')[0]
                all = json.loads(allJson)
                functions.setBase(RID, 'Step 1 of 4: Selecting your chosen meta-variables...')

                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])
                sig_only = int(all["sig_only"])

                metaValsCat = all['metaValsCat']
                metaIDsCat = all['metaIDsCat']
                metaValsQuant = all['metaValsQuant']
                metaIDsQuant = all['metaIDsQuant']

                treeType = int(all['treeType'])
                DepVar = int(all["DepVar"])

                # Create meta-variable DataFrame, final sample list, final category and quantitative field lists based on tree selections
                savedDF, metaDF, finalSampleIDs, catFields, remCatFields, quantFields, catValues, quantValues = functions.getMetaDF(request.user, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, DepVar)
                allFields = catFields + quantFields

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

                functions.setBase(RID, 'Step 1 of 4: Selecting your chosen meta-variables...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 2 of 4: Selecting your chosen taxa or KEGG level...')

                # filter otus based on user settings
                remUnclass = all['remUnclass']
                remZeroes = all['remZeroes']
                perZeroes = int(all['perZeroes'])
                filterData = all['filterData']
                filterPer = int(all['filterPer'])
                filterMeth = int(all['filterMeth'])
                mapTaxa = all['map_taxa']

                finalDF = pd.DataFrame()
                allDF = pd.DataFrame()
                if treeType == 1:
                    if selectAll == 0 or selectAll == 8:
                        taxaString = all["taxa"]
                        taxaDict = json.JSONDecoder(object_pairs_hook=functions.multidict).decode(taxaString)
                        filteredDF = savedDF.copy()
                    else:
                        taxaDict = ''
                        filteredDF = functions.filterDF(savedDF, DepVar, selectAll, remUnclass, remZeroes, perZeroes, filterData, filterPer, filterMeth)
                    finalDF, missingList = functions.getTaxaDF(selectAll, taxaDict, filteredDF, metaDF, allFields, DepVar, RID, stops, PID)

                    if selectAll == 8:
                        result += '\nThe following PGPRs were not detected: ' + ", ".join(missingList) + '\n'
                        result += '===============================================\n'

                if treeType == 2:
                    keggDict = ''
                    if keggAll == 0:
                        keggString = all["kegg"]
                        keggDict = json.JSONDecoder(object_pairs_hook=functions.multidict).decode(keggString)
                    finalDF, allDF = functions.getKeggDF(keggAll, keggDict, savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)

                if treeType == 3:
                    keggDict = ''
                    if nzAll == 0:
                        keggString = all["nz"]
                        keggDict = json.JSONDecoder(object_pairs_hook=functions.multidict).decode(keggString)
                    finalDF, allDF = functions.getNZDF(nzAll, keggDict, savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)

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
                myDir = 'myPhyloDB/media/temp/anova/'
                if not os.path.exists(myDir):
                    os.makedirs(myDir)

                path = str(myDir) + str(RID) + '.biom'
                functions.imploding_panda(path, treeType, DepVar, finalSampleIDs, metaDF, finalDF)

                functions.setBase(RID, 'Step 2 of 4: Selecting your chosen taxa or KEGG level...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 3 of 4: Performing statistical test...')
                finalDict = {}
                seriesList = []
                xAxisDict = {}
                yAxisDict = {}

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                functions.setBase(RID, 'Verifying R packages...missing packages are being installed')

                # R packages from cran
                r("list.of.packages <- c('lsmeans', 'ggplot2', 'RColorBrewer', 'ggthemes')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                print r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

                functions.setBase(RID, 'Step 3 of 4: Performing statistical test...')

                print r("library(lsmeans)")
                print r("library(ggplot2)")
                print r("library(ggthemes)")
                print r("library(RColorBrewer)")
                print r('source("R/myFunctions/myFunctions.R")')

                # R graph
                r.assign('finalDF', finalDF)

                colorVal = all['colorVal']
                if colorVal != 'None':
                    r.assign('colorVal', colorVal)
                else:
                    r.assign('colorVal', 'rank_name')

                xVal = all['xVal']
                if xVal != 'None':
                    r.assign('xVal', xVal)
                else:
                    r.assign('xVal', catFields[0])

                gridVal_X = all['gridVal_X']
                r.assign('gridVal_X', gridVal_X)

                gridVal_Y = all['gridVal_Y']
                r.assign('gridVal_Y', gridVal_Y)

                if DepVar == 0:
                    r('DepVar <- "abund"')
                elif DepVar == 1:
                    r('DepVar <- "rel_abund"')
                elif DepVar == 2:
                    r('DepVar <- "rich"')
                elif DepVar == 3:
                    r('DepVar <- "diversity"')
                elif DepVar == 4:
                    r('DepVar <- "abund_16S"')

                if gridVal_X == 'None' and gridVal_Y == 'None':
                    r("gDF <- data.frame(x=finalDF[,paste(xVal)], y=finalDF[,paste(DepVar)], \
                        myFill=finalDF[,paste(colorVal)])")
                    r("gDF <- aggregate(gDF[, 'y'], list(gDF$x, gDF$myFill), mean)")
                    r("names(gDF) <- c('x', 'myFill', 'y')")
                elif gridVal_X != 'None' and gridVal_Y == 'None':
                    r("gDF <- data.frame(x=finalDF[,paste(xVal)], y=finalDF[,paste(DepVar)], \
                        gridVal_X=finalDF[,paste(gridVal_X)], \
                        myFill=finalDF[,paste(colorVal)])")
                    r("gDF <- aggregate(gDF[, 'y'], list(gDF$x, gDF$myFill, gDF$gridVal_X), mean)")
                    r("names(gDF) <- c('x', 'myFill', 'gridVal_X', 'y')")
                elif gridVal_X == 'None' and gridVal_Y != 'None':
                    r("gDF <- data.frame(x=finalDF[,paste(xVal)], y=finalDF[,paste(DepVar)], \
                        gridVal_Y=finalDF[,paste(gridVal_Y)], \
                        myFill=finalDF[,paste(colorVal)])")
                    r("gDF <- aggregate(gDF[, 'y'], list(gDF$x, gDF$myFill, gDF$gridVal_Y), mean)")
                    r("names(gDF) <- c('x', 'myFill', 'gridVal_Y', 'y')")
                elif gridVal_X != 'None' and gridVal_Y != 'None':
                    r("gDF <- data.frame(x=finalDF[,paste(xVal)], y=finalDF[,paste(DepVar)], \
                        gridVal_X=finalDF[,paste(gridVal_X)], gridVal_Y=finalDF[,paste(gridVal_Y)], \
                        myFill=finalDF[,paste(colorVal)])")
                    r("gDF <- aggregate(gDF[, 'y'], list(gDF$x, gDF$myFill, gDF$gridVal_X, gDF$gridVal_Y), mean)")
                    r("names(gDF) <- c('x', 'myFill', 'gridVal_X', 'gridVal_Y', 'y')")

                r("p <- ggplot(gDF, aes(x=x, y=y, fill=myFill))")
                r("p <- p + geom_bar(stat='identity', position='stack')")

                if gridVal_X != 'None' and gridVal_Y == 'None':
                    r("p <- p + facet_grid(. ~ gridVal_X)")
                    r("p <- p + theme(strip.text.x=element_text(size=10, colour='blue', angle=0))")
                elif gridVal_X == 'None' and gridVal_Y != 'None':
                    r("p <- p + facet_grid(gridVal_Y ~ .)")
                    r("p <- p + theme(strip.text.y=element_text(size=10, colour='blue', angle=90))")
                elif gridVal_X != 'None' and gridVal_Y != 'None':
                    r("p <- p + facet_grid(gridVal_Y ~ gridVal_X)")
                    r("p <- p + theme(strip.text.x=element_text(size=10, colour='blue', angle=0))")
                    r("p <- p + theme(strip.text.y=element_text(size=10, colour='blue', angle=90))")

                r("p <- p + theme(strip.text.x=element_text(size=10, colour='blue', angle=0))")
                r("p <- p + theme(strip.text.y=element_text(size=10, colour='blue', angle=90))")

                palette = all['palette']
                r.assign('palette', palette)
                if palette == 'gdocs':
                    r('pal <- gdocs_pal()(20)')
                elif palette == 'hc':
                    r('pal <- hc_pal()(10)')
                elif palette == 'Set1':
                    r('pal <- brewer.pal(8, "Set1")')
                elif palette == 'Set2':
                    r('pal <- brewer.pal(8, "Set2")')
                elif palette == 'Set3':
                    r('pal <- brewer.pal(12, "Set3")')
                elif palette == 'Paired':
                    r('pal <- brewer.pal(12, "Paired")')
                elif palette == 'Dark2':
                    r('pal <- brewer.pal(12, "Dark2")')
                elif palette == 'Accent':
                    r('pal <- brewer.pal(12, "Accent")')
                r('nColors <- length(pal)')

                r("p <- p + scale_fill_manual(values=rep(pal, ceiling(nlevels(gDF$myFill)/nColors)))")

                r("p <- p + theme(legend.text=element_text(size=7))")
                r("p <- p + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))")

                r("p <- p + theme(legend.title=element_blank())")
                r("p <- p + theme(legend.position='bottom')")
                r("p <- p + guides(fill=guide_legend(ncol=8))")

                if DepVar == 0:
                    r("p <- p + ylab('Abundance') + xlab('')")
                elif DepVar == 1:
                    r("p <- p + ylab('Relative Abundance') + xlab('')")
                elif DepVar == 2:
                    r("p <- p + ylab('OTU Richness') + xlab('')")
                elif DepVar == 3:
                    r("p <- p + ylab('OTU Diversity') + xlab('')")
                elif DepVar == 4:
                    r("p <- p + ylab('Total Abundance') + xlab('')")

                path = "myPhyloDB/media/temp/anova/Rplots"
                if not os.path.exists(path):
                    os.makedirs(path)

                r.assign("path", path)
                r.assign("RID", RID)
                r("file <- paste(path, '/', RID, '.anova.pdf', sep='')")
                r("p <- set_panel_size(p, height=unit(2.9, 'in'), width=unit(2.9, 'in'))")

                r("nlev <- nlevels(as.factor(gDF$gridVal_X))")

                r('if (nlev == 0) { \
                        myWidth <- 11 \
                    } else { \
                        myWidth <- 3*nlev+4 \
                }')

                r("nlev <- nlevels(as.factor(gDF$gridVal_Y))")
                r("nRow <- ceiling(nlevels(gDF$myFill)/8)")
                r('if (nlev == 0) { \
                        myHeight <- 8 \
                    } else { \
                        myHeight <- 2+(3*nlev)+(1*nRow) \
                }')

                r("ggsave(filename=file, plot=p, units='in', height=myHeight, width=myWidth, limitsize=F)")

                # group DataFrame by each taxa level selected
                grouped1 = finalDF.groupby(['rank_name', 'rank_id'])
                pValDict = {}
                counter = 1
                for name1, group1 in grouped1:
                    D = ''
                    r.assign("df", group1)
                    trtString = " * ".join(allFields)

                    if DepVar == 0:
                        anova_string = "fit <- aov(df$abund ~ " + str(trtString) + ", data=df)"
                        r.assign("cmd", anova_string)
                        r("eval(parse(text=cmd))")

                    elif DepVar == 1:
                        anova_string = "fit <- aov(df$rel_abund ~ " + str(trtString) + ", data=df)"
                        r.assign("cmd", anova_string)
                        r("eval(parse(text=cmd))")

                    elif DepVar == 2:
                        anova_string = "fit <- aov(df$rich ~ " + str(trtString) + ", data=df)"
                        r.assign("cmd", anova_string)
                        r("eval(parse(text=cmd))")

                    elif DepVar == 3:
                        anova_string = "fit <- aov(df$diversity ~ " + str(trtString) + ", data=df)"
                        r.assign("cmd", anova_string)
                        r("eval(parse(text=cmd))")

                    elif DepVar == 4:
                        anova_string = "fit <- aov(df$abund_16S ~ " + str(trtString) + ", data=df)"
                        r.assign("cmd", anova_string)
                        r("eval(parse(text=cmd))")

                    aov = r("summary(fit)")
                    pString = r("summary(fit)[[1]][['Pr(>F)']]")

                    tempStuff = pString.split(' ')
                    pList = []
                    for part in tempStuff:
                        try:
                            pList.append(float(part))
                        except Exception:
                            pass

                    if pList:
                        p_val = min(pList)
                        tempStuff = aov.split('\n')
                        for part in tempStuff:
                            if part != tempStuff[0]:
                                D += part + '\n'

                        fList = []
                        for part in tempStuff:
                            if part != tempStuff[0] and part != tempStuff[1]:
                                part = part.replace('\r', '')
                                part1 = part.split(' ')
                                if part1[0] == 'Residuals':
                                    break
                                if part1[0] not in fList:
                                    fList.append(part1[0])

                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                        if stops[PID] == RID:
                            res = ''
                            return HttpResponse(res, content_type='application/json')
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                        D += "\nLSmeans & Tukey's HSD post-hoc test:\n\n"

                        if len(quantFields) == 0:
                            for i in fList:
                                hsd_string = "lsm <- lsmeans(fit, list(pairwise ~ " + str(i) + "))"
                                r.assign("cmd", hsd_string)
                                r("eval(parse(text=cmd))")
                                r("options(width=5000)")
                                table = r("lsm")
                                tempStuff = table.split('\n')
                                for i in xrange(len(tempStuff)):
                                    if i > 0:
                                        D += tempStuff[i] + '\n'
                        else:
                            for i in fList:
                                if i not in quantFields:
                                    hsd_string = "lsm <- lsmeans(fit, list(pairwise ~ " + str(i) + "))"
                                    r.assign("cmd", hsd_string)
                                    r("eval(parse(text=cmd))")
                                    r("options(width=5000)")
                                    table = r("lsm")
                                    tempStuff = table.split('\n')
                                    for i in xrange(len(tempStuff)):
                                        if i > 0:
                                            D += tempStuff[i] + '\n'

                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                        if stops[PID] == RID:
                            res = ''
                            return HttpResponse(res, content_type='application/json')
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    else:
                        p_val = 1.0
                        D = 'ANOVA cannot be performed, please check that you have more than one treatment level and appropriate replication.\n'

                    pValDict[name1] = p_val

                    result += 'Name: ' + str(name1[0]) + '\n'
                    result += 'ID: ' + str(name1[1]) + '\n'
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

                    result += '\nANCOVA table:\n'
                    D = D.decode('utf-8')
                    result += D + '\n'
                    result += '===============================================\n'
                    result += '\n\n\n\n'

                    taxa_no = len(grouped1)
                    functions.setBase(RID, 'Step 3 of 4: Performing statistical test...taxa ' + str(counter) + ' of ' + str(taxa_no) + ' is complete!')
                    counter += 1

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 3 of 4: Performing statistical test...done!')
                functions.setBase(RID, 'Step 4 of 4: Formatting graph data for display...')

                grouped1 = finalDF.groupby(['rank_name', 'rank_id'])
                for name1, group1 in grouped1:
                    dataList = []
                    errorList = []
                    pValue = pValDict[name1]

                    if sig_only == 0:
                        if DepVar == 0:
                            mean = group1.groupby(catFields)['abund'].mean()
                            se = group1.groupby(catFields)['abund'].std()
                            se.fillna(0, inplace=True)
                            high = [x + y for x, y in zip(mean, se)]
                            low = [x - y for x, y in zip(mean, se)]
                            dataList = list(mean)
                            errorTuple = zip(low, high)
                            errorList = [list(elem) for elem in errorTuple]
                        elif DepVar == 1:
                            mean = group1.groupby(catFields)['rel_abund'].mean()
                            se = group1.groupby(catFields)['rel_abund'].std()
                            se.fillna(0, inplace=True)
                            high = [x + y for x, y in zip(mean, se)]
                            low = [x - y for x, y in zip(mean, se)]
                            dataList = list(mean)
                            errorTuple = zip(low, high)
                            errorList = [list(elem) for elem in errorTuple]
                        elif DepVar == 2:
                            mean = group1.groupby(catFields)['rich'].mean()
                            se = group1.groupby(catFields)['rich'].std()
                            se.fillna(0, inplace=True)
                            high = [x + y for x, y in zip(mean, se)]
                            low = [x - y for x, y in zip(mean, se)]
                            dataList = list(mean)
                            errorTuple = zip(low, high)
                            errorList = [list(elem) for elem in errorTuple]
                        elif DepVar == 3:
                            mean = group1.groupby(catFields)['diversity'].mean()
                            se = group1.groupby(catFields)['diversity'].std()
                            se.fillna(0, inplace=True)
                            high = [x + y for x, y in zip(mean, se)]
                            low = [x - y for x, y in zip(mean, se)]
                            dataList = list(mean)
                            errorTuple = zip(low, high)
                            errorList = [list(elem) for elem in errorTuple]
                        elif DepVar == 4:
                            mean = group1.groupby(catFields)['abund_16S'].mean()
                            se = group1.groupby(catFields)['abund_16S'].std()
                            se.fillna(0, inplace=True)
                            high = [x + y for x, y in zip(mean, se)]
                            low = [x - y for x, y in zip(mean, se)]
                            dataList = list(mean)
                            errorTuple = zip(low, high)
                            errorList = [list(elem) for elem in errorTuple]

                        seriesDict = {}
                        seriesDict['name'] = name1
                        seriesDict['type'] = 'column'
                        seriesDict['data'] = dataList
                        seriesList.append(seriesDict)

                        seriesDict = {}
                        seriesDict['name'] = name1
                        seriesDict['type'] = 'errorbar'
                        seriesDict['visible'] = False
                        seriesDict['data'] = errorList
                        seriesList.append(seriesDict)

                    elif sig_only == 1:
                        if pValue < 0.05:
                            if DepVar == 0:
                                mean = group1.groupby(catFields)['abund'].mean()
                                se = group1.groupby(catFields)['abund'].std()
                                se.fillna(0, inplace=True)
                                high = [x + y for x, y in zip(mean, se)]
                                low = [x - y for x, y in zip(mean, se)]
                                dataList = list(mean)
                                errorTuple = zip(low, high)
                                errorList = [list(elem) for elem in errorTuple]
                            elif DepVar == 1:
                                mean = group1.groupby(catFields)['rel_abund'].mean()
                                se = group1.groupby(catFields)['rel_abund'].std()
                                se.fillna(0, inplace=True)
                                high = [x + y for x, y in zip(mean, se)]
                                low = [x - y for x, y in zip(mean, se)]
                                dataList = list(mean)
                                errorTuple = zip(low, high)
                                errorList = [list(elem) for elem in errorTuple]
                            elif DepVar == 2:
                                mean = group1.groupby(catFields)['rich'].mean()
                                se = group1.groupby(catFields)['rich'].std()
                                se.fillna(0, inplace=True)
                                high = [x + y for x, y in zip(mean, se)]
                                low = [x - y for x, y in zip(mean, se)]
                                dataList = list(mean)
                                errorTuple = zip(low, high)
                                errorList = [list(elem) for elem in errorTuple]
                            elif DepVar == 3:
                                mean = group1.groupby(catFields)['diversity'].mean()
                                se = group1.groupby(catFields)['diversity'].std()
                                se.fillna(0, inplace=True)
                                high = [x + y for x, y in zip(mean, se)]
                                low = [x - y for x, y in zip(mean, se)]
                                dataList = list(mean)
                                errorTuple = zip(low, high)
                                errorList = [list(elem) for elem in errorTuple]
                            elif DepVar == 4:
                                mean = group1.groupby(catFields)['abund_16S'].mean()
                                se = group1.groupby(catFields)['abund_16S'].std()
                                se.fillna(0, inplace=True)
                                high = [x + y for x, y in zip(mean, se)]
                                low = [x - y for x, y in zip(mean, se)]
                                dataList = list(mean)
                                errorTuple = zip(low, high)
                                errorList = [list(elem) for elem in errorTuple]

                            seriesDict = {}
                            seriesDict['name'] = name1
                            seriesDict['type'] = 'column'
                            seriesDict['data'] = dataList
                            seriesList.append(seriesDict)

                            seriesDict = {}
                            seriesDict['name'] = name1
                            seriesDict['type'] = 'errorbar'
                            seriesDict['visible'] = False
                            seriesDict['data'] = errorList
                            seriesList.append(seriesDict)

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    '''
                    catFieldsList = []
                    for i in catFields:
                        catFieldsList.append(len(group1.groupby(i)))

                    catFields = [x for (y, x) in sorted(zip(catFieldsList, catFields))]
                    '''

                    if DepVar == 0:
                        grouped2 = group1.groupby(catFields)['abund'].mean()
                    elif DepVar == 1:
                        grouped2 = group1.groupby(catFields)['rel_abund'].mean()
                    elif DepVar == 4:
                        grouped2 = group1.groupby(catFields)['abund_16S'].mean()
                    elif DepVar == 2:
                        grouped2 = group1.groupby(catFields)['rich'].mean()
                    elif DepVar == 3:
                        grouped2 = group1.groupby(catFields)['diversity'].mean()
                    else:
                        raise Exception("Something went horribly wrong")

                    if catFields.__len__() == 1:
                        xAxisDict['categories'] = grouped2.index.values.tolist()
                    else:
                        g2indexvals = grouped2.index.values
                        level = g2indexvals[0].__len__()
                        labelTree = recLabels(g2indexvals, level)
                        xAxisDict['categories'] = labelTree['categories']

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                yTitle = {}
                if DepVar == 0:
                    yTitle['text'] = 'Abundance'
                elif DepVar == 1:
                    yTitle['text'] = 'Relative Abundance'
                elif DepVar == 2:
                    yTitle['text'] = 'OTU Richness'
                elif DepVar == 3:
                    yTitle['text'] = 'OTU Diversity'
                elif DepVar == 4:
                    yTitle['text'] = 'Total Abundance'
                yTitle['style'] = {'fontSize': '18px', 'fontWeight': 'bold'}

                if transform != 0:
                    tname = {
                        '1': "Ln", '2': "Log10", '3': "Sqrt", '4': "Logit", '5': "Arcsin"
                    }
                    yTitle['text'] = tname[str(transform)] + "(" + str(yTitle['text']) + ")"

                yAxisDict['title'] = yTitle

                xStyleDict = {'style': {'fontSize': '14px'}, 'rotation': 0}
                xAxisDict['labels'] = xStyleDict
                yStyleDict = {'style': {'fontSize': '14px'}}
                yAxisDict['labels'] = yStyleDict

                finalDict['series'] = seriesList
                finalDict['xAxis'] = xAxisDict
                finalDict['yAxis'] = yAxisDict
                finalDict['text'] = result

                if not seriesList:
                    finalDict['empty'] = 0
                else:
                    finalDict['empty'] = 1

                functions.setBase(RID, 'Step 4 of 4: Formatting graph data for display...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                # datatable of taxa mapped to selected kegg orthologies
                if not treeType == 1 and mapTaxa == 'yes':
                    myDir = 'myPhyloDB/media/temp/anova/'
                    fileName = str(myDir) + 'Mapped_Taxa.csv'
                    allDF.to_csv(fileName)

                finalDict['resType'] = 'res'
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


def recLabels(lists, level):
    if lists.__len__() == 0:
        return {}

    first = lists
    splitset = []
    for i in range(0, level):
        children = []
        parents = []
        for set in first:
            children.append(set[set.__len__()-1])
            parents.append(set[0:set.__len__()-1])
        first = parents
        splitset.append(children)
    return makeLabels(" ", splitset)


def makeLabels(name, list):
    retDict = {}
    if list.__len__() == 1:
        # final layer
        retDict['name'] = name
        retDict['categories'] = list[0]
        return retDict

    # change here
    children = []
    first = list[list.__len__()-1][0]
    iter = 0
    start = 0
    for stuff in list[list.__len__()-1]:
        if stuff != first:
            sublist = []
            for otherstuff in list[0:list.__len__()-1]:
                sublist.append(otherstuff[start:iter])
            children.append(makeLabels(first, sublist))
            first = stuff
            start = iter
        iter += 1

    # Repeat else condition at the end of the list
    sublist = []
    for otherstuff in list[0:list.__len__()-1]:
        sublist.append(otherstuff[start:iter])
    children.append(makeLabels(first, sublist))

    retDict['name'] = name
    retDict['categories'] = children
    return retDict


def getQuantUnivData(request, RID, stops, PID):
    try:
        while True:
            if request.is_ajax():
                # Get variables from web page
                allJson = request.body.split('&')[0]
                all = json.loads(allJson)

                functions.setBase(RID, 'Step 1 of 4: Selecting your chosen meta-variables...')
                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])
                sig_only = int(all["sig_only"])

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

                functions.setBase(RID, 'Step 1 of 4: Selecting your chosen meta-variables...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                functions.setBase(RID, 'Step 2 of 4: Selecting your chosen taxa or kegg level...')

                # filter otus based on user settings
                remUnclass = all['remUnclass']
                remZeroes = all['remZeroes']
                perZeroes = int(all['perZeroes'])
                filterData = all['filterData']
                filterPer = int(all['filterPer'])
                filterMeth = int(all['filterMeth'])
                mapTaxa = all['map_taxa']

                finalDF = pd.DataFrame()
                allDF = pd.DataFrame()
                if treeType == 1:
                    if selectAll == 0 or selectAll == 8:
                        taxaString = all["taxa"]
                        taxaDict = json.JSONDecoder(object_pairs_hook=functions.multidict).decode(taxaString)
                        filteredDF = savedDF.copy()
                    else:
                        taxaDict = ''
                        filteredDF = functions.filterDF(savedDF, DepVar, selectAll, remUnclass, remZeroes, perZeroes, filterData, filterPer, filterMeth)

                    finalDF, missingList = functions.getTaxaDF(selectAll, taxaDict, filteredDF, metaDF, allFields, DepVar, RID, stops, PID)

                    if selectAll == 8:
                        result += '\nThe following PGPRs were not detected: ' + ", ".join(missingList) + '\n'
                        result += '===============================================\n'

                if treeType == 2:
                    keggDict = ''
                    if keggAll == 0:
                        keggString = all["kegg"]
                        keggDict = json.JSONDecoder(object_pairs_hook=functions.multidict).decode(keggString)
                    finalDF, allDF = functions.getKeggDF(keggAll, keggDict, savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)

                if treeType == 3:
                    keggDict = ''
                    if nzAll == 0:
                        keggString = all["nz"]
                        keggDict = json.JSONDecoder(object_pairs_hook=functions.multidict).decode(keggString)
                    finalDF, allDF = functions.getNZDF(nzAll, keggDict, savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)

                # make sure column types are correct
                finalDF[catFields] = finalDF[catFields].astype(str)
                finalDF[quantFields] = finalDF[quantFields].astype(float)

                # transform Y, if requested
                transform = int(all["transform"])
                finalDF = functions.transformDF(transform, DepVar, finalDF)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/anova/'
                if not os.path.exists(myDir):
                    os.makedirs(myDir)

                path = str(myDir) + str(RID) + '.biom'
                functions.imploding_panda(path, treeType, DepVar, finalSampleIDs, metaDF, finalDF)

                functions.setBase(RID, 'Step 2 of 4: Selecting your chosen taxa or KEGG level...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                functions.setBase(RID, 'Step 3 of 4: Performing statistical test...!')

                finalDict = {}
                # group DataFrame by each taxa level selected
                shapes = ['circle', 'square', 'triangle', 'triangle-down', 'diamond']

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                functions.setBase(RID, 'Verifying R packages...missing packages are being installed')

                # R packages from cran
                r("list.of.packages <- c('ggplot2', 'RColorBrewer', 'ggthemes')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                print r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

                functions.setBase(RID, 'Step 3 of 4: Performing statistical test...')

                print r("library(ggplot2)")
                print r("library(ggthemes)")
                print r("library(RColorBrewer)")
                print r('source("R/myFunctions/myFunctions.R")')

                # R graph
                r.assign('finalDF', finalDF)

                colorVal = all['colorVal']
                if colorVal == 'None':
                    r("colorTrt <- c('All')")
                else:
                    r.assign("colorVal", colorVal)
                    r("colorTrt <- as.factor(finalDF[,paste(colorVal)])")

                r.assign('xVal', quantFields[0])

                gridVal_X = all['gridVal_X']
                if gridVal_X == 'None':
                    r("gridTrt_X <- c('All')")
                else:
                    r.assign("gridVal_X", gridVal_X)
                    r("gridTrt_X <- as.factor(finalDF[,paste(gridVal_X)])")

                gridVal_Y = all['gridVal_Y']
                if gridVal_Y == 'None':
                    r("gridTrt_Y <- c('All')")
                else:
                    r.assign("gridVal_Y", gridVal_Y)
                    r("gridTrt_Y <- as.factor(finalDF[,paste(gridVal_Y)])")

                shapeVal = all['shapeVal']
                if shapeVal == 'None':
                    r("shapeTrt <- c('All')")
                else:
                    r.assign("shapeVal", shapeVal)
                    r("shapeTrt <- as.factor(finalDF[,paste(shapeVal)])")

                if DepVar == 0:
                    r('DepVar <- "abund"')
                elif DepVar == 1:
                    r('DepVar <- "rel_abund"')
                elif DepVar == 2:
                    r('DepVar <- "rich"')
                elif DepVar == 3:
                    r('DepVar <- "diversity"')
                elif DepVar == 4:
                    r('DepVar <- "abund_16S"')

                r("gDF <- data.frame(x=finalDF[,paste(xVal)], y=finalDF[,paste(DepVar)], \
                    gridVal_X=gridTrt_X, gridVal_Y=gridTrt_Y, \
                    myColor=colorTrt, myShape=shapeTrt)")

                r("p <- ggplot(gDF, aes(x=x, y=y, fill=factor(myColor), shape=factor(myShape)) )")
                r("p <- p + geom_point(size=4)")

                if gridVal_X != 'None' and gridVal_Y == 'None':
                    r("p <- p + facet_grid(. ~ gridVal_X)")
                elif gridVal_X == 'None' and gridVal_Y != 'None':
                    r("p <- p + facet_grid(gridVal_Y ~ .)")
                elif gridVal_X != 'None' and gridVal_Y != 'None':
                    r("p <- p + facet_grid(gridVal_Y ~ gridVal_X)")

                r("p <- p + theme(strip.text.x=element_text(size=10, colour='blue', angle=0))")
                r("p <- p + theme(strip.text.y=element_text(size=10, colour='blue', angle=90))")

                palette = all['palette']
                r.assign('palette', palette)
                if palette == 'gdocs':
                    r('pal <- gdocs_pal()(20)')
                elif palette == 'hc':
                    r('pal <- hc_pal()(10)')
                elif palette == 'Set1':
                    r('pal <- brewer.pal(8, "Set1")')
                elif palette == 'Set2':
                    r('pal <- brewer.pal(8, "Set2")')
                elif palette == 'Set3':
                    r('pal <- brewer.pal(12, "Set3")')
                elif palette == 'Paired':
                    r('pal <- brewer.pal(12, "Paired")')
                elif palette == 'Dark2':
                    r('pal <- brewer.pal(12, "Dark2")')
                elif palette == 'Accent':
                    r('pal <- brewer.pal(12, "Accent")')
                r('nColors <- length(pal)')

                r('number <- nlevels(gDF$myColor)')
                r('colors <- rep(pal, length.out=number) ')
                r("p <- p + scale_fill_manual(name='', values=colors, guide=guide_legend(override.aes=list(shape=21)))")

                r('number <- nlevels(gDF$myShape)')
                r('shapes <- rep(c(21, 22, 23, 24, 25), length.out=number) ')
                r("p <- p + scale_shape_manual(name='', values=shapes)")

                r("p <- p + theme(legend.text=element_text(size=7))")
                r("p <- p + theme(legend.position='bottom')")

                r("my.formula <- y ~ x")
                r("p <- p + geom_smooth(method='lm', se=T, color='black', formula=my.formula)")

                if DepVar == 0:
                    r("p <- p + ylab('Abundance') + xlab(paste(xVal))")
                elif DepVar == 1:
                    r("p <- p + ylab('Relative Abundance') + xlab(paste(xVal))")
                elif DepVar == 2:
                    r("p <- p + ylab('OTU Richness') + xlab(paste(xVal))")
                elif DepVar == 3:
                    r("p <- p + ylab('OTU Diversity') + xlab(paste(xVal))")
                elif DepVar == 4:
                    r("p <- p + ylab('Total Abundance') + xlab(paste(xVal))")

                path = "myPhyloDB/media/temp/anova/Rplots"
                if not os.path.exists(path):
                    os.makedirs(path)

                r.assign("path", path)
                r.assign("RID", RID)
                r("file <- paste(path, '/', RID, '.anova.pdf', sep='')")
                r("p <- set_panel_size(p, height=unit(2.9, 'in'), width=unit(2.9, 'in'))")

                r("nlev <- nlevels(as.factor(gDF$gridVal_X))")
                r('if (nlev == 0) { \
                        myWidth <- 8 \
                    } else { \
                        myWidth <- 3*nlev+4, 50 \
                }')

                r("nlev <- nlevels(as.factor(gDF$gridVal_Y))")
                r('if (nlev == 0) { \
                        myHeight <- 8 \
                    } else { \
                        myHeight <- 3*nlev+4, 50 \
                }')

                r("ggsave(filename=file, plot=p, units='in', height=myHeight, width=myWidth, limitsize=F)")

                pValDict = {}
                counter = 1
                catLevels = len(set(catValues))
                grouped1 = finalDF.groupby(['rank_name', 'rank_id'])
                for name1, group1 in grouped1:
                    D = ''
                    r.assign("df", group1)

                    trtString = " * ".join(allFields)
                    if DepVar == 0:
                        anova_string = "fit <- lm(abund ~ " + str(trtString) + ", data=df)"
                        r.assign("cmd", anova_string)
                        r("eval(parse(text=cmd))")
                    elif DepVar == 1:
                        anova_string = "fit <- lm(rel_abund ~ " + str(trtString) + ", data=df)"
                        r.assign("cmd", anova_string)
                        r("eval(parse(text=cmd))")
                    elif DepVar == 2:
                        anova_string = "fit <- lm(rich ~ " + str(trtString) + ", data=df)"
                        r.assign("cmd", anova_string)
                        r("eval(parse(text=cmd))")
                    elif DepVar == 3:
                        anova_string = "fit <- lm(diversity ~ " + str(trtString) + ", data=df)"
                        r.assign("cmd", anova_string)
                        r("eval(parse(text=cmd))")
                    elif DepVar == 4:
                        anova_string = "fit <- lm(abund_16S ~ " + str(trtString) + ", data=df)"
                        r.assign("cmd", anova_string)
                        r("eval(parse(text=cmd))")

                    # calculate predicted scores (full model)
                    pred_string = "df$pred <- predict(fit, df)"
                    r.assign("cmd", pred_string)
                    r("eval(parse(text=cmd))")

                    r("options(width=5000)")
                    r("aov <- anova(fit)")
                    aov = r("aov")

                    tempStuff = aov.split('\n')
                    for i in xrange(len(tempStuff)):
                        if i >= 4:
                            D += tempStuff[i] + '\n'

                    fit = r("summary(fit)")
                    tempStuff = fit.split('\n')
                    for i in xrange(len(tempStuff)):
                        if i >= 8:
                            D += tempStuff[i] + '\n'
                    r("p_vals <- summary(fit)$coefficients[,4]")
                    p_vals = r.get("p_vals")

                    if not np.isnan(p_vals).any():
                        p_value = min(p_vals)
                        pValDict[name1] = p_value
                    else:
                        pValDict[name1] = np.nan

                    result += 'Name: ' + str(name1[0]) + '\n'
                    result += 'ID: ' + str(name1[1]) + '\n'
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

                    result += '\nANCOVA table:\n'
                    D = D.decode('utf-8')
                    result += D + '\n'
                    result += '===============================================\n'
                    result += '\n\n\n\n'

                    taxa_no = len(grouped1)
                    functions.setBase(RID, 'Step 3 of 4: Performing statistical test...taxa ' + str(counter) + ' of ' + str(taxa_no) + ' is complete!')
                    counter += 1

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                functions.setBase(RID, 'Step 3 of 4: Performing statistical test...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                functions.setBase(RID, 'Step 4 of 4: Formatting graph data for display...')

                finalDF['sample_name'] = ''
                for index, row in finalDF.iterrows():
                    val = Sample.objects.get(sampleid=row['sampleid']).sample_name
                    finalDF.loc[index, 'sample_name'] = val

                shapes_idx = 0
                seriesList = []
                grouped1 = finalDF.groupby(['rank_name', 'rank_id'])
                for name1, group1 in grouped1:
                    pValue = pValDict[name1]

                    if sig_only == 0:
                        if catLevels > 1:
                            grouped2 = group1.groupby(catFields)
                            for name2, group2 in grouped2:
                                dataList = []
                                x = []
                                y = []
                                if DepVar == 0:
                                    x = group2[quantFields[0]].values.astype(float).tolist()
                                    y = group2['abund'].values.astype(float).tolist()
                                elif DepVar == 1:
                                    x = group2[quantFields[0]].values.astype(float).tolist()
                                    y = group2['rel_abund'].values.astype(float).tolist()
                                elif DepVar == 2:
                                    x = group2[quantFields[0]].values.astype(float).tolist()
                                    y = group2['rich'].values.astype(float).tolist()
                                elif DepVar == 3:
                                    x = group2[quantFields[0]].values.astype(float).tolist()
                                    y = group2['diversity'].values.astype(float).tolist()
                                elif DepVar == 4:
                                    x = group2[quantFields[0]].values.astype(float).tolist()
                                    y = group2['abund_16S'].values.astype(float).tolist()

                                if DepVar == 0:
                                    for index, row in group2.iterrows():
                                        dataDict = {}
                                        dataDict['name'] = row['sample_name']
                                        dataDict['x'] = float(row[quantFields[0]])
                                        dataDict['y'] = float(row['abund'])
                                        dataList.append(dataDict)
                                elif DepVar == 1:
                                    for index, row in group2.iterrows():
                                        dataDict = {}
                                        dataDict['name'] = row['sample_name']
                                        dataDict['x'] = float(row[quantFields[0]])
                                        dataDict['y'] = float(row['rel_abund'])
                                        dataList.append(dataDict)
                                elif DepVar == 2:
                                    for index, row in group2.iterrows():
                                        dataDict = {}
                                        dataDict['name'] = row['sample_name']
                                        dataDict['x'] = float(row[quantFields[0]])
                                        dataDict['y'] = float(row['rich'])
                                        dataList.append(dataDict)
                                elif DepVar == 3:
                                    for index, row in group2.iterrows():
                                        dataDict = {}
                                        dataDict['name'] = row['sample_name']
                                        dataDict['x'] = float(row[quantFields[0]])
                                        dataDict['y'] = float(row['diversity'])
                                        dataList.append(dataDict)
                                elif DepVar == 4:
                                    for index, row in group2.iterrows():
                                        dataDict = {}
                                        dataDict['name'] = row['sample_name']
                                        dataDict['x'] = float(row[quantFields[0]])
                                        dataDict['y'] = float(row['abund_16S'])
                                        dataList.append(dataDict)

                                seriesDict = {}
                                seriesDict['turboThreshold'] = 0
                                seriesDict['type'] = 'scatter'
                                seriesDict['name'] = str(name1[0]) + ": " + str(name2)
                                seriesDict['data'] = dataList

                                markerDict = {}
                                markerDict['symbol'] = shapes[shapes_idx]
                                seriesDict['marker'] = markerDict
                                seriesDict['data'] = dataList

                                seriesList.append(seriesDict)

                                slp, inter, r_value, p, std_err = stats.linregress(x, y)
                                min_y = float(slp*min(x) + inter)
                                max_y = float(slp*max(x) + inter)
                                slope = "%0.3f" % slp
                                intercept = "%0.3f" % inter
                                r_sq = r_value * r_value
                                r_square = "%0.3f" % r_sq

                                regrList = []
                                regrList.append([float(min(x)), min_y])
                                regrList.append([float(max(x)), max_y])

                                regrDict = {}
                                regrDict['type'] = 'line'
                                sup2 = u"\u00B2"
                                regrDict['name'] = 'y = ' + str(slope) + 'x' + ' + ' + str(intercept) + '; R' + sup2 + ' = ' + str(r_square)
                                regrDict['data'] = regrList
                                regrDict['color'] = 'black'

                                markerDict = {}
                                markerDict['enabled'] = False
                                regrDict['marker'] = markerDict
                                seriesList.append(regrDict)

                                shapes_idx += 1
                                if shapes_idx >= len(shapes):
                                    shapes_idx = 0

                        else:   # if catLevel <=1
                            dataList = []
                            x = []
                            y = []
                            if DepVar == 0:
                                x = group1[quantFields[0]].values.astype(float).tolist()
                                y = group1['abund'].values.astype(float).tolist()
                            elif DepVar == 1:
                                x = group1[quantFields[0]].values.astype(float).tolist()
                                y = group1['rel_abund'].values.astype(float).tolist()
                            elif DepVar == 2:
                                x = group1[quantFields[0]].values.astype(float).tolist()
                                y = group1['rich'].values.astype(float).tolist()
                            elif DepVar == 3:
                                x = group1[quantFields[0]].values.astype(float).tolist()
                                y = group1['diversity'].values.astype(float).tolist()
                            elif DepVar == 4:
                                x = group1[quantFields[0]].values.astype(float).tolist()
                                y = group1['abund_16S'].values.astype(float).tolist()

                            if DepVar == 0:
                                for index, row in group1.iterrows():
                                    dataDict = {}
                                    dataDict['name'] = row['sample_name']
                                    dataDict['x'] = float(row[quantFields[0]])
                                    dataDict['y'] = float(row['abund'])
                                    dataList.append(dataDict)
                            elif DepVar == 1:
                                for index, row in group1.iterrows():
                                    dataDict = {}
                                    dataDict['name'] = row['sample_name']
                                    dataDict['x'] = float(row[quantFields[0]])
                                    dataDict['y'] = float(row['rel_abund'])
                                    dataList.append(dataDict)
                            elif DepVar == 2:
                                for index, row in group1.iterrows():
                                    dataDict = {}
                                    dataDict['name'] = row['sample_name']
                                    dataDict['x'] = float(row[quantFields[0]])
                                    dataDict['y'] = float(row['rich'])
                                    dataList.append(dataDict)
                            elif DepVar == 3:
                                for index, row in group1.iterrows():
                                    dataDict = {}
                                    dataDict['name'] = row['sample_name']
                                    dataDict['x'] = float(row[quantFields[0]])
                                    dataDict['y'] = float(row['diversity'])
                                    dataList.append(dataDict)
                            elif DepVar == 4:
                                for index, row in group1.iterrows():
                                    dataDict = {}
                                    dataDict['name'] = row['sample_name']
                                    dataDict['x'] = float(row[quantFields[0]])
                                    dataDict['y'] = float(row['abund_16S'])
                                    dataList.append(dataDict)

                            seriesDict = {}
                            seriesDict['turboThreshold'] = 0
                            seriesDict['type'] = 'scatter'
                            seriesDict['name'] = str(name1[0])
                            seriesDict['data'] = dataList

                            markerDict = {}
                            markerDict['symbol'] = shapes[shapes_idx]
                            seriesDict['marker'] = markerDict
                            seriesDict['data'] = dataList

                            seriesList.append(seriesDict)

                            slp, inter, r_value, p, std_err = stats.linregress(x, y)
                            min_y = float(slp*min(x) + inter)
                            max_y = float(slp*max(x) + inter)
                            slope = "%0.3f" % slp
                            intercept = "%0.3f" % inter
                            r_sq = r_value * r_value
                            r_square = "%0.3f" % r_sq

                            regrList = []
                            regrList.append([float(min(x)), min_y])
                            regrList.append([float(max(x)), max_y])

                            regrDict = {}
                            regrDict['type'] = 'line'
                            sup2 = u"\u00B2"
                            regrDict['name'] = 'y = ' + str(slope) + 'x' + ' + ' + str(intercept) + '; R' + sup2 + ' = ' + str(r_square)
                            regrDict['data'] = regrList
                            regrDict['color'] = 'black'

                            markerDict = {}
                            markerDict['enabled'] = False
                            regrDict['marker'] = markerDict
                            seriesList.append(regrDict)

                    elif sig_only == 1:
                        if pValue < 0.05:
                            if catLevels > 1:
                                grouped2 = group1.groupby(catFields)
                                for name2, group2 in grouped2:
                                    dataList = []
                                    x = []
                                    y = []
                                    if DepVar == 0:
                                        x = group2[quantFields[0]].values.astype(float).tolist()
                                        y = group2['abund'].values.astype(float).tolist()
                                    elif DepVar == 1:
                                        x = group2[quantFields[0]].values.astype(float).tolist()
                                        y = group2['rel_abund'].values.astype(float).tolist()
                                    elif DepVar == 2:
                                        x = group2[quantFields[0]].values.astype(float).tolist()
                                        y = group2['rich'].values.astype(float).tolist()
                                    elif DepVar == 3:
                                        x = group2[quantFields[0]].values.astype(float).tolist()
                                        y = group2['diversity'].values.astype(float).tolist()
                                    elif DepVar == 4:
                                        x = group2[quantFields[0]].values.astype(float).tolist()
                                        y = group2['abund_16S'].values.astype(float).tolist()

                                    if DepVar == 0:
                                        for index, row in group2.iterrows():
                                            dataDict = {}
                                            dataDict['name'] = row['sample_name']
                                            dataDict['x'] = float(row[quantFields[0]])
                                            dataDict['y'] = float(row['abund'])
                                            dataList.append(dataDict)
                                    elif DepVar == 1:
                                        for index, row in group2.iterrows():
                                            dataDict = {}
                                            dataDict['name'] = row['sample_name']
                                            dataDict['x'] = float(row[quantFields[0]])
                                            dataDict['y'] = float(row['rel_abund'])
                                            dataList.append(dataDict)
                                    elif DepVar == 2:
                                        for index, row in group2.iterrows():
                                            dataDict = {}
                                            dataDict['name'] = row['sample_name']
                                            dataDict['x'] = float(row[quantFields[0]])
                                            dataDict['y'] = float(row['rich'])
                                            dataList.append(dataDict)
                                    elif DepVar == 3:
                                        for index, row in group2.iterrows():
                                            dataDict = {}
                                            dataDict['name'] = row['sample_name']
                                            dataDict['x'] = float(row[quantFields[0]])
                                            dataDict['y'] = float(row['diversity'])
                                            dataList.append(dataDict)
                                    elif DepVar == 4:
                                        for index, row in group2.iterrows():
                                            dataDict = {}
                                            dataDict['name'] = row['sample_name']
                                            dataDict['x'] = float(row[quantFields[0]])
                                            dataDict['y'] = float(row['abund_16S'])
                                            dataList.append(dataDict)

                                    seriesDict = {}
                                    seriesDict['turboThreshold'] = 0
                                    seriesDict['type'] = 'scatter'
                                    seriesDict['name'] = str(name1[0]) + ": " + str(name2)
                                    seriesDict['data'] = dataList

                                    markerDict = {}
                                    markerDict['symbol'] = shapes[shapes_idx]
                                    seriesDict['marker'] = markerDict
                                    seriesDict['data'] = dataList

                                    seriesList.append(seriesDict)

                                    slp, inter, r_value, p, std_err = stats.linregress(x, y)
                                    min_y = float(slp*min(x) + inter)
                                    max_y = float(slp*max(x) + inter)
                                    slope = "%0.3f" % slp
                                    intercept = "%0.3f" % inter
                                    r_sq = r_value * r_value
                                    r_square = "%0.3f" % r_sq

                                    regrList = []
                                    regrList.append([float(min(x)), min_y])
                                    regrList.append([float(max(x)), max_y])

                                    regrDict = {}
                                    regrDict['type'] = 'line'
                                    regrDict['name'] = 'y = ' + str(slope) + 'x' + ' + ' + str(intercept) + '; R2 = ' + str(r_square)
                                    regrDict['data'] = regrList
                                    regrDict['color'] = 'black'

                                    markerDict = {}
                                    markerDict['enabled'] = False
                                    regrDict['marker'] = markerDict
                                    seriesList.append(regrDict)

                                    shapes_idx += 1
                                    if shapes_idx >= len(shapes):
                                        shapes_idx = 0

                            else:   # if catLevel <=1
                                dataList = []
                                x = []
                                y = []
                                if DepVar == 0:
                                    x = group1[quantFields[0]].values.astype(float).tolist()
                                    y = group1['abund'].values.astype(float).tolist()
                                elif DepVar == 1:
                                    x = group1[quantFields[0]].values.astype(float).tolist()
                                    y = group1['rel_abund'].values.astype(float).tolist()
                                elif DepVar == 2:
                                    x = group1[quantFields[0]].values.astype(float).tolist()
                                    y = group1['rich'].values.astype(float).tolist()
                                elif DepVar == 3:
                                    x = group1[quantFields[0]].values.astype(float).tolist()
                                    y = group1['diversity'].values.astype(float).tolist()
                                elif DepVar == 4:
                                    x = group1[quantFields[0]].values.astype(float).tolist()
                                    y = group1['abund_16S'].values.astype(float).tolist()

                                if DepVar == 0:
                                    for index, row in group1.iterrows():
                                        dataDict = {}
                                        dataDict['name'] = row['sample_name']
                                        dataDict['x'] = float(row[quantFields[0]])
                                        dataDict['y'] = float(row['abund'])
                                        dataList.append(dataDict)
                                elif DepVar == 1:
                                    for index, row in group1.iterrows():
                                        dataDict = {}
                                        dataDict['name'] = row['sample_name']
                                        dataDict['x'] = float(row[quantFields[0]])
                                        dataDict['y'] = float(row['rel_abund'])
                                        dataList.append(dataDict)
                                elif DepVar == 2:
                                    for index, row in group1.iterrows():
                                        dataDict = {}
                                        dataDict['name'] = row['sample_name']
                                        dataDict['x'] = float(row[quantFields[0]])
                                        dataDict['y'] = float(row['rich'])
                                        dataList.append(dataDict)
                                elif DepVar == 3:
                                    for index, row in group1.iterrows():
                                        dataDict = {}
                                        dataDict['name'] = row['sample_name']
                                        dataDict['x'] = float(row[quantFields[0]])
                                        dataDict['y'] = float(row['diversity'])
                                        dataList.append(dataDict)
                                elif DepVar == 4:
                                    for index, row in group1.iterrows():
                                        dataDict = {}
                                        dataDict['name'] = row['sample_name']
                                        dataDict['x'] = float(row[quantFields[0]])
                                        dataDict['y'] = float(row['abund_16S'])
                                        dataList.append(dataDict)

                                seriesDict = {}
                                seriesDict['turboThreshold'] = 0
                                seriesDict['type'] = 'scatter'
                                seriesDict['name'] = str(name1[0])
                                seriesDict['data'] = dataList

                                markerDict = {}
                                markerDict['symbol'] = shapes[shapes_idx]
                                seriesDict['marker'] = markerDict
                                seriesDict['data'] = dataList

                                seriesList.append(seriesDict)

                                slp, inter, r_value, p, std_err = stats.linregress(x, y)
                                min_y = float(slp*min(x) + inter)
                                max_y = float(slp*max(x) + inter)
                                slope = "%0.3f" % slp
                                intercept = "%0.3f" % inter
                                r_sq = r_value * r_value
                                r_square = "%0.3f" % r_sq

                                regrList = []
                                regrList.append([float(min(x)), min_y])
                                regrList.append([float(max(x)), max_y])

                                regrDict = {}
                                regrDict['type'] = 'line'
                                regrDict['name'] = 'y = ' + str(slope) + 'x' + ' + ' + str(intercept) + '; R2 = ' + str(r_square)
                                regrDict['data'] = regrList
                                regrDict['color'] = 'black'

                                markerDict = {}
                                markerDict['enabled'] = False
                                regrDict['marker'] = markerDict
                                seriesList.append(regrDict)

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                xAxisDict = {}
                xTitle = {}
                xTitle['text'] = quantFields[0]
                xTitle['style'] = {'fontSize': '18px', 'fontWeight': 'bold'}
                xAxisDict['title'] = xTitle

                yAxisDict = {}
                yTitle = {}
                if DepVar == 0:
                    yTitle['text'] = 'Abundance'
                elif DepVar == 1:
                    yTitle['text'] = 'Relative Abundance'
                elif DepVar == 2:
                    yTitle['text'] = 'OTU Richness'
                elif DepVar == 3:
                    yTitle['text'] = 'OTU Diversity'
                elif DepVar == 4:
                    yTitle['text'] = 'Total Abundance'
                yAxisDict['title'] = yTitle

                if transform != 0:
                    tname = {
                        '1': "Ln", '2': "Log10", '3': "Sqrt", '4': "Logit", '5': "Arcsin"
                    }
                    yTitle['text'] = tname[str(transform)] + "(" + str(yTitle['text']) + ")"

                yTitle['style'] = {'fontSize': '18px', 'fontWeight': 'bold'}
                yAxisDict['title'] = yTitle

                styleDict = {'style': {'fontSize': '14px'}}
                xAxisDict['labels'] = styleDict
                yAxisDict['labels'] = styleDict

                finalDict['series'] = seriesList
                finalDict['xAxis'] = xAxisDict
                finalDict['yAxis'] = yAxisDict

                if not seriesList:
                    finalDict['empty'] = 0
                else:
                    finalDict['empty'] = 1

                functions.setBase(RID, 'Step 4 of 4: Formatting graph data for display...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                # datatable of taxa mapped to selected kegg orthologies
                if not treeType == 1 and mapTaxa == 'yes':
                    records = allDF.values.tolist()
                    finalDict['taxData'] = json.dumps(records)
                    columns = allDF.columns.values.tolist()
                    finalDict['taxColumns'] = json.dumps(columns)

                finalDict['resType'] = 'res'
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


