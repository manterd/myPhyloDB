import datetime
from django.http import HttpResponse
import logging
import json
import numpy as np
import pandas as pd
from pyper import *
from scipy import stats

from database.models import Sample

import database.utils
import database.queue


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getCatUnivData(request, RID, stops, PID):
    try:
        while True:
            if request.is_ajax():
                # Get variables from web page
                allJson = request.body.split('&')[0]
                all = json.loads(allJson)
                database.queue.setBase(RID, 'Step 1 of 4: Selecting your chosen meta-variables...')

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
                savedDF, metaDF, finalSampleIDs, catFields, remCatFields, quantFields, catValues, quantValues = database.utils.getMetaDF(request.user, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, DepVar)
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

                database.queue.setBase(RID, 'Step 1 of 4: Selecting your chosen meta-variables...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 2 of 4: Selecting your chosen taxa or KEGG level...')

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
                        taxaDict = json.JSONDecoder(object_pairs_hook=database.utils.multidict).decode(taxaString)
                        filteredDF = savedDF.copy()
                    else:
                        taxaDict = ''
                        filteredDF = database.utils.filterDF(savedDF, DepVar, selectAll, remUnclass, remZeroes, perZeroes, filterData, filterPer, filterMeth)
                    finalDF, missingList = database.utils.getTaxaDF(selectAll, taxaDict, filteredDF, metaDF, allFields, DepVar, RID, stops, PID)

                    if selectAll == 8:
                        result += '\nThe following PGPRs were not detected: ' + ", ".join(missingList) + '\n'
                        result += '===============================================\n'

                if treeType == 2:
                    keggDict = ''
                    if keggAll == 0:
                        keggString = all["kegg"]
                        keggDict = json.JSONDecoder(object_pairs_hook=database.utils.multidict).decode(keggString)
                    finalDF, allDF = database.utils.getKeggDF(keggAll, keggDict, savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)

                if treeType == 3:
                    keggDict = ''
                    if nzAll == 0:
                        keggString = all["nz"]
                        keggDict = json.JSONDecoder(object_pairs_hook=database.utils.multidict).decode(keggString)
                    finalDF, allDF = database.utils.getNZDF(nzAll, keggDict, savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)

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
                finalDF = database.utils.transformDF(transform, DepVar, finalDF)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/anova/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                database.queue.setBase(RID, 'Step 2 of 4: Selecting your chosen taxa or KEGG level...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 3 of 4: Performing statistical test...')
                finalDict = {}
                seriesList = []
                xAxisDict = {}
                yAxisDict = {}

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                database.queue.setBase(RID, 'Verifying R packages...missing packages are being installed')

                # R packages from cran
                r("list.of.packages <- c('lsmeans')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                print r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

                database.queue.setBase(RID, 'Step 3 of 4: Performing statistical test...')

                print r("library(lsmeans)")

                # group DataFrame by each taxa level selected
                grouped1 = finalDF.groupby(['rank', 'rank_name', 'rank_id'])
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

                    result += 'Level: ' + str(name1[0]) + '\n'
                    result += 'Name: ' + str(name1[1]) + '\n'
                    result += 'ID: ' + str(name1[2]) + '\n'
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
                    database.queue.setBase(RID, 'Step 3 of 4: Performing statistical test...taxa ' + str(counter) + ' of ' + str(taxa_no) + ' is complete!')
                    counter += 1

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 3 of 4: Performing statistical test...done!')
                database.queue.setBase(RID, 'Step 4 of 4: Formatting graph data for display...')

                grouped1 = finalDF.groupby(['rank', 'rank_name', 'rank_id'])
                for name1, group1 in grouped1:
                    dataList = []
                    errorList = []
                    pValue = pValDict[name1]

                    if sig_only == 0:
                        if DepVar == 0:
                            mean = group1.groupby(catFields)['abund'].mean()
                            se = group1.groupby(catFields)['abund'].std() / np.sqrt(group1.groupby(catFields)['abund'].count())
                            se.fillna(0, inplace=True)
                            high = [x + y for x, y in zip(mean, se)]
                            low = [x - y for x, y in zip(mean, se)]
                            dataList = list(mean)
                            errorTuple = zip(low, high)
                            errorList = [list(elem) for elem in errorTuple]
                        elif DepVar == 1:
                            mean = group1.groupby(catFields)['rel_abund'].mean()
                            se = group1.groupby(catFields)['rel_abund'].std() / np.sqrt(group1.groupby(catFields)['rel_abund'].count())
                            se.fillna(0, inplace=True)
                            high = [x + y for x, y in zip(mean, se)]
                            low = [x - y for x, y in zip(mean, se)]
                            dataList = list(mean)
                            errorTuple = zip(low, high)
                            errorList = [list(elem) for elem in errorTuple]
                        elif DepVar == 2:
                            mean = group1.groupby(catFields)['rich'].mean()
                            se = group1.groupby(catFields)['rich'].std() / np.sqrt(group1.groupby(catFields)['rich'].count())
                            se.fillna(0, inplace=True)
                            high = [x + y for x, y in zip(mean, se)]
                            low = [x - y for x, y in zip(mean, se)]
                            dataList = list(mean)
                            errorTuple = zip(low, high)
                            errorList = [list(elem) for elem in errorTuple]
                        elif DepVar == 3:
                            mean = group1.groupby(catFields)['diversity'].mean()
                            se = group1.groupby(catFields)['diversity'].std() / np.sqrt(group1.groupby(catFields)['diversity'].count())
                            se.fillna(0, inplace=True)
                            high = [x + y for x, y in zip(mean, se)]
                            low = [x - y for x, y in zip(mean, se)]
                            dataList = list(mean)
                            errorTuple = zip(low, high)
                            errorList = [list(elem) for elem in errorTuple]
                        elif DepVar == 4:
                            mean = group1.groupby(catFields)['abund_16S'].mean()
                            se = group1.groupby(catFields)['abund_16S'].std() / np.sqrt(group1.groupby(catFields)['abund_16S'].count())
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
                                se = group1.groupby(catFields)['abund'].std() / np.sqrt(group1.groupby(catFields)['abund'].count())
                                se.fillna(0, inplace=True)
                                high = [x + y for x, y in zip(mean, se)]
                                low = [x - y for x, y in zip(mean, se)]
                                dataList = list(mean)
                                errorTuple = zip(low, high)
                                errorList = [list(elem) for elem in errorTuple]
                            elif DepVar == 1:
                                mean = group1.groupby(catFields)['rel_abund'].mean()
                                se = group1.groupby(catFields)['rel_abund'].std() / np.sqrt(group1.groupby(catFields)['rel_abund'].count())
                                se.fillna(0, inplace=True)
                                high = [x + y for x, y in zip(mean, se)]
                                low = [x - y for x, y in zip(mean, se)]
                                dataList = list(mean)
                                errorTuple = zip(low, high)
                                errorList = [list(elem) for elem in errorTuple]
                            elif DepVar == 2:
                                mean = group1.groupby(catFields)['rich'].mean()
                                se = group1.groupby(catFields)['rich'].std() / np.sqrt(group1.groupby(catFields)['rich'].count())
                                se.fillna(0, inplace=True)
                                high = [x + y for x, y in zip(mean, se)]
                                low = [x - y for x, y in zip(mean, se)]
                                dataList = list(mean)
                                errorTuple = zip(low, high)
                                errorList = [list(elem) for elem in errorTuple]
                            elif DepVar == 3:
                                mean = group1.groupby(catFields)['diversity'].mean()
                                se = group1.groupby(catFields)['diversity'].std() / np.sqrt(group1.groupby(catFields)['diversity'].count())
                                se.fillna(0, inplace=True)
                                high = [x + y for x, y in zip(mean, se)]
                                low = [x - y for x, y in zip(mean, se)]
                                dataList = list(mean)
                                errorTuple = zip(low, high)
                                errorList = [list(elem) for elem in errorTuple]
                            elif DepVar == 4:
                                mean = group1.groupby(catFields)['abund_16S'].mean()
                                se = group1.groupby(catFields)['abund_16S'].std() / np.sqrt(group1.groupby(catFields)['abund_16S'].count())
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

                    if catFields.__len__() == 1:
                        xAxisDict['categories'] = grouped2.index.values.tolist()
                    else:
                        g2indexvals = grouped2.index.values
                        level = g2indexvals[0].__len__()
                        labelTree = recLabels(g2indexvals, level)
                        xAxisDict['categories'] = [labelTree]

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

                database.queue.setBase(RID, 'Step 4 of 4: Formatting graph data for display...done!')

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

                database.queue.setBase(RID, 'Step 1 of 4: Selecting your chosen meta-variables...')
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
                savedDF, metaDF, finalSampleIDs, catFields, remCatFields, quantFields, catValues, quantValues = database.utils.getMetaDF(request.user, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, DepVar)
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

                database.queue.setBase(RID, 'Step 1 of 4: Selecting your chosen meta-variables...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                database.queue.setBase(RID, 'Step 2 of 4: Selecting your chosen taxa or kegg level...')

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
                        taxaDict = json.JSONDecoder(object_pairs_hook=database.utils.multidict).decode(taxaString)
                        filteredDF = savedDF.copy()
                    else:
                        taxaDict = ''
                        filteredDF = database.utils.filterDF(savedDF, DepVar, selectAll, remUnclass, remZeroes, perZeroes, filterData, filterPer, filterMeth)

                    finalDF, missingList = database.utils.getTaxaDF(selectAll, taxaDict, filteredDF, metaDF, allFields, DepVar, RID, stops, PID)

                    if selectAll == 8:
                        result += '\nThe following PGPRs were not detected: ' + ", ".join(missingList) + '\n'
                        result += '===============================================\n'

                if treeType == 2:
                    keggDict = ''
                    if keggAll == 0:
                        keggString = all["kegg"]
                        keggDict = json.JSONDecoder(object_pairs_hook=database.utils.multidict).decode(keggString)
                    finalDF, allDF = database.utils.getKeggDF(keggAll, keggDict, savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)

                if treeType == 3:
                    keggDict = ''
                    if nzAll == 0:
                        keggString = all["nz"]
                        keggDict = json.JSONDecoder(object_pairs_hook=database.utils.multidict).decode(keggString)
                    finalDF, allDF = database.utils.getNZDF(nzAll, keggDict, savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)

                # make sure column types are correct
                finalDF[catFields] = finalDF[catFields].astype(str)
                finalDF[quantFields] = finalDF[quantFields].astype(float)

                # transform Y, if requested
                transform = int(all["transform"])
                finalDF = database.utils.transformDF(transform, DepVar, finalDF)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/anova/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                database.queue.setBase(RID, 'Step 2 of 4: Selecting your chosen taxa or KEGG level...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                database.queue.setBase(RID, 'Step 3 of 4: Performing statistical test...!')
                finalDict = {}
                # group DataFrame by each taxa level selected
                shapes = ['circle', 'square', 'triangle', 'triangle-down', 'diamond']

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)
                pValDict = {}
                counter = 1
                catLevels = len(set(catValues))
                grouped1 = finalDF.groupby(['rank', 'rank_name', 'rank_id'])
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

                    result += 'Level: ' + str(name1[0]) + '\n'
                    result += 'Name: ' + str(name1[1]) + '\n'
                    result += 'ID: ' + str(name1[2]) + '\n'
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
                    database.queue.setBase(RID, 'Step 3 of 4: Performing statistical test...taxa ' + str(counter) + ' of ' + str(taxa_no) + ' is complete!')
                    counter += 1

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                database.queue.setBase(RID, 'Step 3 of 4: Performing statistical test...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                database.queue.setBase(RID, 'Step 4 of 4: Formatting graph data for display...')

                finalDF['sample_name'] = ''
                for index, row in finalDF.iterrows():
                    val = Sample.objects.get(sampleid=row['sampleid']).sample_name
                    finalDF.loc[index, 'sample_name'] = val

                shapes_idx = 0
                seriesList = []
                grouped1 = finalDF.groupby(['rank', 'rank_name', 'rank_id'])
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
                                seriesDict['name'] = str(name1[1]) + ": " + str(name2)
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
                            seriesDict['name'] = str(name1[1])
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
                                    seriesDict['name'] = str(name1[1]) + ": " + str(name2)
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
                                seriesDict['name'] = str(name1[1])
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

                database.queue.setBase(RID, 'Step 4 of 4: Formatting graph data for display...done!')

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


