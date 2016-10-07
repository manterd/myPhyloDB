import datetime
from django.http import HttpResponse
import logging
import numpy as np
import pandas as pd
from pyper import *
from scipy import stats
import simplejson

from database.models import Sample
from database.utils import multidict, getMetaDF, transformDF
from database.utils_kegg import getTaxaDF, getKeggDF, getNZDF, filterDF
import database.queue


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getCatUnivData(request, RID, stops, PID):
    try:
        while True:
            if request.is_ajax():
                # Get variables from web page
                allJson = request.body.split('&')[0]
                all = simplejson.loads(allJson)
                database.queue.setBase(RID, 'Step 1 of 4: Selecting your chosen meta-variables...')
                myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
                path = str(myDir) + 'usr_norm_data.csv'

                with open(path, 'rb') as f:
                    savedDF = pd.read_csv(f, index_col=0, sep=',')

                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])
                sig_only = int(all["sig_only"])

                metaValsCat = all['metaValsCat']
                metaIDsCat = all['metaIDsCat']
                metaValsQuant = all['metaValsQuant']
                metaIDsQuant = all['metaIDsQuant']

                button3 = int(all['button3'])
                DepVar = int(all["DepVar"])

                # Create meta-variable DataFrame, final sample list, final category and quantitative field lists based on tree selections
                savedDF, metaDF, finalSampleIDs, catFields, remCatFields, quantFields, catValues, quantValues = getMetaDF(savedDF, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, DepVar)
                allFields = catFields + quantFields

                if not catFields:
                    error = "Selected categorical variable(s) contain only one level.\nPlease select different variable(s)."
                    myDict = {'error': error}
                    res = simplejson.dumps(myDict)
                    return HttpResponse(res, content_type='application/json')

                if not finalSampleIDs:
                    error = "No valid samples were contained in your final dataset.\nPlease select different variable(s)."
                    myDict = {'error': error}
                    res = simplejson.dumps(myDict)
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

                finalDF = pd.DataFrame()
                if button3 == 1:
                    taxaString = all["taxa"]
                    taxaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(taxaString)
                    finalDF, missingList = getTaxaDF(selectAll, taxaDict, savedDF, metaDF, allFields, DepVar, RID, stops, PID)
                    if selectAll == 8:
                        result += '\nThe following PGPRs were not detected: ' + ", ".join(missingList) + '\n'
                        result += '===============================================\n'

                elif button3 == 2:
                    keggString = all["kegg"]
                    keggDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(keggString)
                    finalDF = getKeggDF(keggAll, keggDict, savedDF, metaDF, allFields, DepVar, RID, stops, PID)

                elif button3 == 3:
                    nzString = all["nz"]
                    nzDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(nzString)
                    finalDF = getNZDF(nzAll, nzDict, savedDF, metaDF, allFields, DepVar, RID, stops, PID)

                # filter phylotypes based on user settings
                remUnclass = all['remUnclass']
                remZeroes = all['remZeroes']
                perZeroes = int(all['perZeroes'])
                filterData = all['filterData']
                filterPer = int(all['filterPer'])
                filterMeth = int(all['filterMeth'])
                finalDF = filterDF(finalDF, remUnclass, remZeroes, perZeroes, filterData, filterPer, filterMeth)

                # make sure column types are correct
                finalDF[catFields] = finalDF[catFields].astype(str)
                finalDF[quantFields] = finalDF[quantFields].astype(float)

                # transform Y, if requested
                transform = int(all["transform"])
                finalDF = transformDF(transform, DepVar, finalDF)

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
                colors_idx = 0
                colors = [
                    "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                    "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                    "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                    "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                    "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                    "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                    "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                    "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                    "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                    "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
                    "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
                    "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
                    "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58",
                    "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
                    "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
                    "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5",
                    "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
                    "#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01",
                    "#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966",
                    "#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
                    "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#C6DC99", "#203B3C",
                    "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
                    "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
                    "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
                    "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F",
                    "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
                    "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
                    "#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00",
                    "#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
                    "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25",
                    "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"
                ]

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

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
                        except:
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
                        r("library(lsmeans)")
                        # r('TukeyHSD(fit)')

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
                        result += 'Dependent Variable: Species Richness' + '\n'
                    elif DepVar == 3:
                        result += 'Dependent Variable: Species Diversity' + '\n'
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
                    pValue = pValDict[name1]

                    if sig_only == 0:
                        if DepVar == 0:
                            grouped2 = group1.groupby(catFields)['abund'].mean()
                            dataList = list(grouped2)
                        elif DepVar == 1:
                            grouped2 = group1.groupby(catFields)['rel_abund'].mean()
                            dataList = list(grouped2)
                        elif DepVar == 2:
                            grouped2 = group1.groupby(catFields)['rich'].mean()
                            dataList = list(grouped2)
                        elif DepVar == 3:
                            grouped2 = group1.groupby(catFields)['diversity'].mean()
                            dataList = list(grouped2)
                        elif DepVar == 4:
                            grouped2 = group1.groupby(catFields)['abund_16S'].mean()
                            dataList = list(grouped2)

                        seriesDict = {}
                        seriesDict['name'] = name1
                        seriesDict['color'] = colors[colors_idx]
                        seriesDict['data'] = dataList
                        seriesList.append(seriesDict)

                    elif sig_only == 1:
                        if pValue < 0.05:
                            if DepVar == 0:
                                grouped2 = group1.groupby(catFields)['abund'].mean()
                                dataList = list(grouped2)
                            elif DepVar == 1:
                                grouped2 = group1.groupby(catFields)['rel_abund'].mean()
                                dataList = list(grouped2)
                            elif DepVar == 2:
                                grouped2 = group1.groupby(catFields)['rich'].mean()
                                dataList = list(grouped2)
                            elif DepVar == 3:
                                grouped2 = group1.groupby(catFields)['diversity'].mean()
                                dataList = list(grouped2)
                            elif DepVar == 4:
                                grouped2 = group1.groupby(catFields)['abund_16S'].mean()
                                dataList = list(grouped2)

                            seriesDict = {}
                            seriesDict['name'] = name1
                            seriesDict['color'] = colors[colors_idx]
                            seriesDict['data'] = dataList
                            seriesList.append(seriesDict)

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    colors_idx += 1
                    if colors_idx >= len(colors):
                        colors_idx = 0

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
                    yTitle['text'] = 'Species Richness'
                elif DepVar == 3:
                    yTitle['text'] = 'Species Diversity'
                elif DepVar == 4:
                    yTitle['text'] = 'Total Abundance'
                yTitle['style'] = {'color': 'black', 'fontSize': '18px', 'fontWeight': 'bold'}

                if transform != 0:
                    tname = {
                        '1': "Ln", '2': "Log10", '3': "Sqrt", '4': "Logit", '5': "Arcsin"
                    }
                    yTitle['text'] = tname[str(transform)] + "(" + yTitle['text'] + ")"

                yAxisDict['title'] = yTitle

                xStyleDict = {'style': {'color': 'black', 'fontSize': '14px'}, 'rotation': 0}
                xAxisDict['labels'] = xStyleDict
                yStyleDict = {'style': {'color': 'black', 'fontSize': '14px'}}
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


                finalDict['resType'] = 'res'
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
                all = simplejson.loads(allJson)

                database.queue.setBase(RID, 'Step 1 of 4: Selecting your chosen meta-variables...')
                myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
                path = str(myDir) + 'usr_norm_data.csv'

                with open(path, 'rb') as f:
                    savedDF = pd.read_csv(f, index_col=0, sep=',')

                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])
                sig_only = int(all["sig_only"])

                # Select samples and meta-variables from savedDF
                metaValsCat = all['metaValsCat']
                metaIDsCat = all['metaIDsCat']
                metaValsQuant = all['metaValsQuant']
                metaIDsQuant = all['metaIDsQuant']

                button3 = int(all['button3'])
                DepVar = int(all["DepVar"])

                # Create meta-variable DataFrame, final sample list, final category and quantitative field lists based on tree selections
                savedDF, metaDF, finalSampleIDs, catFields, remCatFields, quantFields, catValues, quantValues = getMetaDF(savedDF, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, DepVar)
                allFields = catFields + quantFields

                if not finalSampleIDs:
                    error = "No valid samples were contained in your final dataset.\nPlease select different variable(s)."
                    myDict = {'error': error}
                    res = simplejson.dumps(myDict)
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

                finalDF = pd.DataFrame()
                if button3 == 1:
                    taxaString = all["taxa"]
                    taxaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(taxaString)
                    finalDF, missingList = getTaxaDF(selectAll, taxaDict, savedDF, metaDF, allFields, DepVar, RID, stops, PID)
                    if selectAll == 8:
                        result += '\nThe following PGPRs were not detected: ' + ", ".join(missingList) + '\n'
                        result += '===============================================\n'

                if button3 == 2:
                    keggString = all["kegg"]
                    keggDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(keggString)
                    finalDF, missingList = getKeggDF(keggAll, keggDict, savedDF, metaDF, allFields, DepVar, RID, stops, PID)

                if button3 == 3:
                    nzString = all["nz"]
                    nzDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(nzString)
                    finalDF, missingList = getNZDF(nzAll, nzDict, savedDF, metaDF, allFields, DepVar, RID, stops, PID)

                # filter phylotypes based on user settings
                remUnclass = all['remUnclass']
                remZeroes = all['remZeroes']
                perZeroes = int(all['perZeroes'])
                filterData = all['filterData']
                filterPer = int(all['filterPer'])
                filterMeth = int(all['filterMeth'])
                finalDF = filterDF(finalDF, remUnclass, remZeroes, perZeroes, filterData, filterPer, filterMeth)

                # make sure column types are correct
                finalDF[catFields] = finalDF[catFields].astype(str)
                finalDF[quantFields] = finalDF[quantFields].astype(float)

                # transform Y, if requested
                transform = int(all["transform"])
                finalDF = transformDF(transform, DepVar, finalDF)

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
                colors = [
                    "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                    "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                    "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                    "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                    "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                    "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                    "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                    "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                    "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                    "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
                    "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
                    "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
                    "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58",
                    "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
                    "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
                    "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5",
                    "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
                    "#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01",
                    "#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966",
                    "#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
                    "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#C6DC99", "#203B3C",
                    "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
                    "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
                    "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
                    "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F",
                    "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
                    "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
                    "#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00",
                    "#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
                    "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25",
                    "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"
                ]

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
                        result += 'Dependent Variable: Species Richness' + '\n'
                    elif DepVar == 3:
                        result += 'Dependent Variable: Species Diversity' + '\n'
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

                colors_idx = 0
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
                                seriesDict['color'] = colors[colors_idx]

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
                                regrDict['color'] = colors[colors_idx]

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
                            seriesDict['color'] = colors[colors_idx]

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
                            regrDict['color'] = colors[colors_idx]

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
                                    seriesDict['color'] = colors[colors_idx]

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
                                    regrDict['color'] = colors[colors_idx]

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
                                    y = group1['rel_abund'].values.astype(float).tolist()
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
                                        dataDict['y'] = float(row['rel_abund'])
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
                                seriesDict['color'] = colors[colors_idx]

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
                                regrDict['color'] = colors[colors_idx]

                                markerDict = {}
                                markerDict['enabled'] = False
                                regrDict['marker'] = markerDict
                                seriesList.append(regrDict)

                    colors_idx += 1
                    if colors_idx >= len(colors):
                        colors_idx = 0

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                xAxisDict = {}
                xTitle = {}
                xTitle['text'] = quantFields[0]
                xTitle['style'] = {'color': 'black', 'fontSize': '18px', 'fontWeight': 'bold'}
                xAxisDict['title'] = xTitle

                yAxisDict = {}
                yTitle = {}
                if DepVar == 0:
                    yTitle['text'] = 'Abundance'
                elif DepVar == 1:
                    yTitle['text'] = 'Relative Abundance'
                elif DepVar == 2:
                    yTitle['text'] = 'Species Richness'
                elif DepVar == 3:
                    yTitle['text'] = 'Shannon Diversity'
                elif DepVar == 4:
                    yTitle['text'] = 'Total Abundance'
                yAxisDict['title'] = yTitle

                if transform != 0:
                    tname = {
                        '1': "Ln", '2': "Log10", '3': "Sqrt", '4': "Logit", '5': "Arcsin"
                    }
                    yTitle['text'] = tname[str(transform)] + "(" + yTitle['text'] + ")"

                yTitle['style'] = {'color': 'black', 'fontSize': '18px', 'fontWeight': 'bold'}
                yAxisDict['title'] = yTitle

                styleDict = {'style': {'color': 'black', 'fontSize': '14px'}}
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

                finalDict['resType'] = 'res'
                finalDict['text'] = result
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


