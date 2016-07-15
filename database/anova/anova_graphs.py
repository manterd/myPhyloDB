import ast
import datetime
from django import db
from django.http import HttpResponse
from django_pandas.io import read_frame
import multiprocessing as mp
import logging
import math
import numpy as np
import pandas as pd
import pickle
from pyper import *
from scipy import stats
import shutil
import simplejson
import threading

from database.models import Sample
from database.models import PICRUSt
from database.models import ko_lvl1, ko_lvl2, ko_lvl3, ko_entry
from database.models import nz_lvl1, nz_lvl2, nz_lvl3, nz_lvl4, nz_entry
from database.utils import multidict
import database.queue
from config import local_cfg


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}

LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def removeRIDANOVA(RID):
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


def statusANOVA(request):
    global base, stage, time1, time2, TimeDiff
    # if request.is_ajax():
    RID = request.GET["all"]
    time2[RID] = time.time()
    try:
        TimeDiff[RID] = time2[RID] - time1[RID]
    except:
        TimeDiff[RID] = 0
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
    myDict = {'stage': stage[RID]}
    json_data = simplejson.dumps(myDict, encoding="Latin-1")
    return HttpResponse(json_data, content_type='application/json')


def getCatUnivData(request, RID, stops, PID):
    global base, stage, time1, TimeDiff
    try:
        while True:
            if request.is_ajax():
                # Get variables from web page
                allJson = request.GET["all"]
                all = simplejson.loads(allJson)

                time1[RID] = time.time()  # Moved these down here so RID is available
                base[RID] = 'Step 1 of 4: Selecting your chosen meta-variables...'

                myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
                path = str(myDir) + 'usr_norm_data.csv'

                with open(path, 'rb') as f:
                    savedDF = pd.read_csv(f, index_col=0, sep='\t')

                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])
                sig_only = int(all["sig_only"])

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

                if not catFields_edit:
                    myDict = {}
                    myDict['error'] = "Selected meta data only has one level.\nPlease select different variable(s)."
                    res = simplejson.dumps(myDict)
                    return HttpResponse(res, content_type='application/json')

                catSampleIDs = []
                if metaIDsCat:
                    idDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaIDsCat)
                    for key in sorted(idDictCat):
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
                        quantSampleIDs.extend(idDictQuant[key])

                allSampleIDs = catSampleIDs + quantSampleIDs
                allFields = catFields_edit + quantFields

                # Removes samples (rows) that are not in our samplelist
                tempDF = savedDF.loc[savedDF['sampleid'].isin(allSampleIDs)]

                if metaDictCat:
                    for key in metaDictCat:
                        tempDF = tempDF.loc[tempDF[key].isin(metaDictCat[key])]

                if metaDictQuant:
                    for key in metaDictQuant:
                        valueList = [float(x) for x in metaDictQuant[key]]
                        tempDF = tempDF.loc[tempDF[key].isin(valueList)]

                metaDF = tempDF[allFields]

                result = ''
                result += 'Categorical variables selected by user: ' + ", ".join(catFields) + '\n'
                result += 'Categorical variables removed from analysis (contains only 1 level): ' + ", ".join(removed) + '\n'
                result += 'Quantitative variables selected by user: ' + ", ".join(quantFields) + '\n'
                result += '===============================================\n'

                base[RID] = 'Step 1 of 4: Selecting your chosen meta-variables...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 2 of 4: Selecting your chosen taxa or KEGG level...'
                button3 = int(all['button3'])

                DepVar = 1
                finalDF = pd.DataFrame()
                if button3 == 1:
                    DepVar = int(all["DepVar_taxa"])
                    taxaString = all["taxa"]
                    taxaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(taxaString)
                    finalDF = getTaxaDF(selectAll, taxaDict, savedDF, metaDF, allFields, DepVar, RID, stops, PID)

                if button3 == 2:
                    DepVar = int(all["DepVar_kegg"])
                    keggString = all["kegg"]
                    keggDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(keggString)
                    finalDF = getKeggDF(keggAll, keggDict, savedDF, tempDF, allFields, DepVar, RID, stops, PID)

                if button3 == 3:
                    DepVar = int(all["DepVar_nz"])
                    nzString = all["nz"]
                    nzDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(nzString)
                    finalDF = getNZDF(nzAll, nzDict, savedDF, tempDF, allFields, DepVar, RID, stops, PID)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/anova/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                base[RID] = 'Step 2 of 4: Selecting your chosen taxa or KEGG level...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 3 of 4: Performing statistical test...'
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

                # group DataFrame by each taxa level selected
                grouped1 = finalDF.groupby(['rank', 'rank_name', 'rank_id'])
                pValDict = {}
                counter = 1
                for name1, group1 in grouped1:
                    D = ''

                    if os.name == 'nt':
                        r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                    else:
                        r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                    r.assign("df", group1)
                    trtString = " * ".join(allFields)
                    if DepVar == 1:
                        anova_string = "fit <- aov(rel_abund ~ " + str(trtString) + ", data=df)"
                        r.assign("cmd", anova_string)
                        r("eval(parse(text=cmd))")
                    elif DepVar == 2:
                        anova_string = "fit <- aov(rich ~ " + str(trtString) + ", data=df)"
                        r.assign("cmd", anova_string)
                        r("eval(parse(text=cmd))")
                    elif DepVar == 3:
                        anova_string = "fit <- aov(diversity ~ " + str(trtString) + ", data=df)"
                        r.assign("cmd", anova_string)
                        r("eval(parse(text=cmd))")
                    elif DepVar == 4:
                        anova_string = "fit <- aov(abund_16S ~ " + str(trtString) + ", data=df)"
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
                            placeholder = ''

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
                    if DepVar == 1:
                        result += 'Dependent Variable: Relative Abundance' + '\n'
                    elif DepVar == 2:
                        result += 'Dependent Variable: Species Richness' + '\n'
                    elif DepVar == 3:
                        result += 'Dependent Variable: Species Diversity' + '\n'
                    elif DepVar == 4:
                        result += 'Dependent Variable: Total Abundance (rRNA gene copies)' + '\n'

                    result += '\nANCOVA table:\n'
                    D = D.decode('utf-8')
                    result += D + '\n'
                    result += '===============================================\n'
                    result += '\n\n\n\n'

                    taxa_no = len(grouped1)
                    base[RID] = 'Step 3 of 4: Performing statistical test...taxa ' + str(counter) + ' of ' + str(taxa_no) + ' is complete!'
                    counter += 1

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 3 of 4: Performing statistical test...done!'
                base[RID] = 'Step 4 of 4: Formatting graph data for display...'

                grouped1 = finalDF.groupby(['rank', 'rank_name', 'rank_id'])
                for name1, group1 in grouped1:
                    dataList = []
                    pValue = pValDict[name1]

                    if sig_only == 0:
                        if DepVar == 1:
                            grouped2 = group1.groupby(catFields_edit)['rel_abund'].mean()
                            dataList = list(grouped2)
                        elif DepVar == 2:
                            grouped2 = group1.groupby(catFields_edit)['rich'].mean()
                            dataList = list(grouped2)
                        elif DepVar == 3:
                            grouped2 = group1.groupby(catFields_edit)['diversity'].mean()
                            dataList = list(grouped2)
                        elif DepVar == 4:
                            grouped2 = group1.groupby(catFields_edit)['abund_16S'].mean()
                            dataList = list(grouped2)

                    elif sig_only == 1:
                        if pValue < 0.05:
                            if DepVar == 1:
                                grouped2 = group1.groupby(catFields_edit)['rel_abund'].mean()
                                dataList = list(grouped2)
                            elif DepVar == 2:
                                grouped2 = group1.groupby(catFields_edit)['rich'].mean()
                                dataList = list(grouped2)
                            elif DepVar == 3:
                                grouped2 = group1.groupby(catFields_edit)['diversity'].mean()
                                dataList = list(grouped2)
                            elif DepVar == 4:
                                grouped2 = group1.groupby(catFields_edit)['abund_16S'].mean()
                                dataList = list(grouped2)

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    seriesDict = {}
                    seriesDict['name'] = name1
                    seriesDict['color'] = colors[colors_idx]
                    seriesDict['data'] = dataList
                    seriesList.append(seriesDict)

                    colors_idx += 1
                    if colors_idx >= len(colors):
                        colors_idx = 0

                    if DepVar == 1:
                        grouped2 = group1.groupby(catFields_edit)['rel_abund'].mean()
                    elif DepVar == 4:
                        grouped2 = group1.groupby(catFields_edit)['abund_16S'].mean()
                    elif DepVar == 2:
                        grouped2 = group1.groupby(catFields_edit)['rich'].mean()
                    elif DepVar == 3:
                        grouped2 = group1.groupby(catFields_edit)['diversity'].mean()

                    if catFields_edit.__len__() == 1:
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
                if DepVar == 1:
                    yTitle['text'] = 'Relative Abundance'
                elif DepVar == 2:
                    yTitle['text'] = 'Species Richness'
                elif DepVar == 3:
                    yTitle['text'] = 'Species Diversity'
                elif DepVar == 4:
                    yTitle['text'] = 'Total Abundance (rRNA gene copies)'
                yTitle['style'] = {'color': 'black', 'fontSize': '18px', 'fontWeight': 'bold'}
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

                base[RID] = 'Step 4 of 4: Formatting graph data for display...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')

                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                removeRIDANOVA(RID)
                return HttpResponse(res, content_type='application/json')

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with ANcOVA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            removeRIDANOVA(RID)
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
    global base, stage, time1, TimeDiff
    try:
        while True:
            if request.is_ajax():
                # Get variables from web page
                allJson = request.GET["all"]
                all = simplejson.loads(allJson)

                time1[RID] = time.time()  # Moved these down here so RID is available
                base[RID] = 'Step 1 of 4: Selecting your chosen meta-variables...'

                path = pickle.loads(request.session['savedDF'])
                savedDF = pd.read_pickle(path)

                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])
                sig_only = int(all["sig_only"])

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
                        quantSampleIDs.extend(idDictQuant[key])

                allSampleIDs = catSampleIDs + quantSampleIDs
                allFields = catFields_edit + quantFields

                # Removes samples (rows) that are not in our samplelist
                tempDF = savedDF.loc[savedDF['sampleid'].isin(allSampleIDs)]

                if metaDictCat:
                    for key in metaDictCat:
                        tempDF = tempDF.loc[tempDF[key].isin(metaDictCat[key])]

                if metaDictQuant:
                    for key in metaDictCat:
                        tempDF = tempDF.loc[tempDF[key].isin(metaDictCat[key])]

                metaDF = tempDF[allFields]

                result = ''
                result += 'Categorical variables selected by user: ' + ", ".join(catFields) + '\n'
                result += 'Categorical variables removed from analysis (contains only 1 level): ' + ", ".join(removed) + '\n'
                result += 'Quantitative variables selected by user: ' + ", ".join(quantFields) + '\n'
                result += '===============================================\n'

                base[RID] = 'Step 1 of 4: Selecting your chosen meta-variables...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                base[RID] = 'Step 2 of 4: Selecting your chosen taxa or kegg level...'
                button3 = int(all['button3'])

                finalDF = pd.DataFrame()
                DepVar = 1
                if button3 == 1:
                    DepVar = int(all["DepVar_taxa"])
                    taxaString = all["taxa"]
                    taxaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(taxaString)
                    finalDF = getTaxaDF(selectAll, taxaDict, savedDF, metaDF, allFields, DepVar, RID, stops, PID)

                if button3 == 2:
                    DepVar = int(all["DepVar_kegg"])
                    keggString = all["kegg"]
                    keggDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(keggString)
                    finalDF = getKeggDF(keggAll, keggDict, savedDF, tempDF, allFields, DepVar, RID, stops, PID)

                if button3 == 3:
                    DepVar = int(all["DepVar_nz"])
                    nzString = all["nz"]
                    nzDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(nzString)
                    finalDF = getNZDF(nzAll, nzDict, savedDF, tempDF, allFields, DepVar, RID, stops, PID)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/anova/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                base[RID] = 'Step 2 of 4: Selecting your chosen taxa or KEGG level...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                base[RID] = 'Step 3 of 4: Performing statistical test...!'
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
                    if DepVar == 1:
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

                    result += 'Level: ' + str(name1[0]) + '\n'
                    result += 'Name: ' + str(name1[1]) + '\n'
                    result += 'ID: ' + str(name1[2]) + '\n'
                    if DepVar == 1:
                        result += 'Dependent Variable: Relative Abundance' + '\n'
                    elif DepVar == 2:
                        result += 'Dependent Variable: Species Richness' + '\n'
                    elif DepVar == 3:
                        result += 'Dependent Variable: Species Diversity' + '\n'
                    elif DepVar == 4:
                        result += 'Dependent Variable: Total Abundance (rRNA gene copies)' + '\n'

                    result += '\nANCOVA table:\n'
                    D = D.decode('utf-8')
                    result += D + '\n'
                    result += '===============================================\n'
                    result += '\n\n\n\n'

                    taxa_no = len(grouped1)
                    base[RID] = 'Step 3 of 4: Performing statistical test...taxa ' + str(counter) + ' of ' + str(taxa_no) + ' is complete!'
                    counter += 1

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                base[RID] = 'Step 3 of 4: Performing statistical test...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                base[RID] = 'Step 4 of 4: Formatting graph data for display...'


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
                            grouped2 = group1.groupby(catFields_edit)
                            for name2, group2 in grouped2:
                                dataList = []
                                x = []
                                y = []
                                if DepVar == 1:
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

                                if DepVar == 1:
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
                            if DepVar == 1:
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

                            if DepVar == 1:
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
                                grouped2 = group1.groupby(catFields_edit)
                                for name2, group2 in grouped2:
                                    dataList = []
                                    x = []
                                    y = []
                                    if DepVar == 1:
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

                                    if DepVar == 1:
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
                                if DepVar == 1:
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

                                if DepVar == 1:
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
                if DepVar == 1:
                    yTitle['text'] = 'Relative Abundance'
                elif DepVar == 2:
                    yTitle['text'] = 'Species Richness'
                elif DepVar == 3:
                    yTitle['text'] = 'Shannon Diversity'
                elif DepVar == 4:
                    yTitle['text'] = 'Total Abundance (rRNA gene copies)'
                yAxisDict['title'] = yTitle

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

                base[RID] = 'Step 4 of 4: Formatting graph data for display...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDict['text'] = result
                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                removeRIDANOVA(RID)
                return HttpResponse(res, content_type='application/json')

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "Error with ANcOVA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = simplejson.dumps(myDict)
            removeRIDANOVA(RID)
            return HttpResponse(res, content_type='application/json')


def getTaxaDF(selectAll, taxaDict, savedDF, metaDF, allFields, DepVar, RID, stops, PID):
    global base
    try:
        base[RID] = 'Step 2 of 4: Selecting your chosen taxa or KEGG level...'
        taxaDF = pd.DataFrame(columns=['sampleid', 'rank', 'rank_id', 'rank_name', 'rel_abund', 'abund_16S', 'rich', 'diversity'])
        if selectAll == 0:
            for key in taxaDict:
                taxaList = taxaDict[key]
                if isinstance(taxaList, unicode):
                    if key == 'Kingdom':
                        tempDF = savedDF.loc[savedDF['kingdomid'] == taxaList]
                        tempDF = tempDF[['sampleid', 'kingdomid', 'kingdomName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'kingdomid': 'rank_id', 'kingdomName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Kingdom'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Phyla':
                        tempDF = savedDF.loc[savedDF['phylaid'] == taxaList]
                        tempDF = tempDF[['sampleid', 'phylaid', 'phylaName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'phylaid': 'rank_id', 'phylaName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Phyla'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Class':
                        tempDF = savedDF.loc[savedDF['classid'] == taxaList]
                        tempDF = tempDF[['sampleid', 'classid', 'className', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'classid': 'rank_id', 'className': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Class'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Order':
                        tempDF = savedDF.loc[savedDF['orderid'] == taxaList]
                        tempDF = tempDF[['sampleid', 'orderid', 'orderName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'orderid': 'rank_id', 'orderName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Order'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Family':
                        tempDF = savedDF.loc[savedDF['familyid'] == taxaList]
                        tempDF = tempDF[['sampleid', 'familyid', 'familyName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'familyid': 'rank_id', 'familyName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Family'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Genus':
                        tempDF = savedDF.loc[savedDF['genusid'] == taxaList]
                        tempDF = tempDF[['sampleid', 'genusid', 'genusName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'genusid': 'rank_id', 'genusName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Genus'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Species':
                        tempDF = savedDF.loc[savedDF['speciesid'] == taxaList]
                        tempDF = tempDF[['sampleid', 'speciesid', 'speciesName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'speciesid': 'rank_id', 'speciesName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Species'
                        taxaDF = taxaDF.append(tempDF)
                else:
                    if key == 'Kingdom':
                        tempDF = savedDF.loc[savedDF['kingdomid'].isin(taxaList)]
                        tempDF = tempDF[['sampleid', 'kingdomid', 'kingdomName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'kingdomid': 'rank_id', 'kingdomName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Kingdom'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Phyla':
                        tempDF = savedDF.loc[savedDF['phylaid'].isin(taxaList)]
                        tempDF = tempDF[['sampleid', 'phylaid', 'phylaName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'phylaid': 'rank_id', 'phylaName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Phyla'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Class':
                        tempDF = savedDF.loc[savedDF['classid'].isin(taxaList)]
                        tempDF = tempDF[['sampleid', 'classid', 'className', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'classid': 'rank_id', 'className': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Class'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Order':
                        tempDF = savedDF.loc[savedDF['orderid'].isin(taxaList)]
                        tempDF = tempDF[['sampleid', 'orderid', 'orderName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'orderid': 'rank_id', 'orderName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Order'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Family':
                        tempDF = savedDF.loc[savedDF['familyid'].isin(taxaList)]
                        tempDF = tempDF[['sampleid', 'familyid', 'familyName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'familyid': 'rank_id', 'familyName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Family'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Genus':
                        tempDF = savedDF.loc[savedDF['genusid'].isin(taxaList)]
                        tempDF = tempDF[['sampleid', 'genusid', 'genusName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'genusid': 'rank_id', 'genusName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Genus'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Species':
                        tempDF = savedDF.loc[savedDF['speciesid'].isin(taxaList)]
                        tempDF = tempDF[['sampleid', 'speciesid', 'speciesName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'speciesid': 'rank_id', 'speciesName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Species'
                        taxaDF = taxaDF.append(tempDF)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif selectAll == 1:
            taxaDF = savedDF.loc[:, ['sampleid', 'kingdomid', 'kingdomName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'kingdomid': 'rank_id', 'kingdomName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Kingdom'
        elif selectAll == 2:
            taxaDF = savedDF.loc[:, ['sampleid', 'phylaid', 'phylaName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'phylaid': 'rank_id', 'phylaName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Phyla'
        elif selectAll == 3:
            taxaDF = savedDF.loc[:, ['sampleid', 'classid', 'className', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'classid': 'rank_id', 'className': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Class'
        elif selectAll == 4:
            taxaDF = savedDF.loc[:, ['sampleid', 'orderid', 'orderName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'orderid': 'rank_id', 'orderName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Order'
        elif selectAll == 5:
            taxaDF = savedDF.loc[:, ['sampleid', 'familyid', 'familyName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'familyid': 'rank_id', 'familyName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Family'
        elif selectAll == 6:
            taxaDF = savedDF.loc[:, ['sampleid', 'genusid', 'genusName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'genusid': 'rank_id', 'genusName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Genus'
        elif selectAll == 7:
            taxaDF = savedDF.loc[:, ['sampleid', 'speciesid', 'speciesName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'speciesid': 'rank_id', 'speciesName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Species'

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')

        wantedList = allFields + ['sampleid', 'rank', 'rank_name', 'rank_id']
        if DepVar == 1:
            finalDF = finalDF.groupby(wantedList)[['rel_abund']].sum()
        elif DepVar == 4:
            finalDF = finalDF.groupby(wantedList)[['abund_16S']].sum()
        elif DepVar == 2:
            finalDF = finalDF.groupby(wantedList)[['rich']].sum()
        elif DepVar == 3:
            finalDF = finalDF.groupby(wantedList)[['diversity']].sum()

        finalDF.reset_index(drop=False, inplace=True)

        return finalDF

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "Error with ANcOVA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getKeggDF(keggAll, keggDict, savedDF, tempDF, allFields, DepVar, RID, stops, PID):
    global base
    try:
        base[RID] = 'Step 2 of 4: Selecting your chosen taxa or KEGG level...'
        koDict = {}
        if keggAll == 0:
            for key in keggDict:
                keggList = keggDict[key]
                if isinstance(keggList, unicode):
                    if key == 'Level1':
                        koList = ko_entry.objects.using('picrust').filter(ko_lvl1_id_id=keggList).values_list('ko_orthology', flat=True)
                        if koList:
                            koDict[keggList] = koList
                    elif key == 'Level2':
                        koList = ko_entry.objects.using('picrust').filter(ko_lvl2_id_id=keggList).values_list('ko_orthology', flat=True)
                        if koList:
                            koDict[keggList] = koList
                    elif key == 'Level3':
                        koList = ko_entry.objects.using('picrust').filter(ko_lvl3_id_id=keggList).values_list('ko_orthology', flat=True)
                        if koList:
                            koDict[keggList] = koList
                    elif key == 'Level4':
                        koList = ko_entry.objects.using('picrust').filter(ko_lvl4_id=keggList).values_list('ko_orthology', flat=True)
                        if koList:
                            koDict[keggList] = koList
                else:
                    for i in keggList:
                        if key == 'Level1':
                            koList = ko_entry.objects.using('picrust').filter(ko_lvl1_id_id=i).values_list('ko_orthology', flat=True)
                            if koList:
                                koDict[i] = koList
                        elif key == 'Level2':
                            koList = ko_entry.objects.using('picrust').filter(ko_lvl2_id_id=i).values_list('ko_orthology', flat=True)
                            if koList:
                                koDict[i] = koList
                        elif key == 'Level3':
                            koList = ko_entry.objects.using('picrust').filter(ko_lvl3_id_id=i).values_list('ko_orthology', flat=True)
                            if koList:
                                koDict[i] = koList
                        elif key == 'Level4':
                            koList = ko_entry.objects.using('picrust').filter(ko_lvl4_id=i).values_list('ko_orthology', flat=True)
                            if koList:
                                koDict[i] = koList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif keggAll == 1:
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

        path = 'myPhyloDB/media/temp/anova/' + str(RID)
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
            processes = [threading.Thread(target=sumStuff, args=(listDF[x], koDict, RID, x, PID)) for x in xrange(numcore)]
        else:
            numcore = min(local_cfg.usr_numcore, int(maxCPU))
            listDF = np.array_split(picrustDF, numcore)
            processes = [threading.Thread(target=sumStuff, args=(listDF[x], koDict, RID, x, PID)) for x in xrange(numcore)]

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
            path = 'myPhyloDB/media/temp/anova/'+str(RID)+'/file%d.temp' % i
            frame = pd.read_csv(path)
            picrustDF = picrustDF.append(frame, ignore_index=True)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        shutil.rmtree('myPhyloDB/media/temp/anova/'+str(RID))
        picrustDF.set_index('speciesid', inplace=True)
        picrustDF[picrustDF > 0] = 1

        # merge to get final gene counts for all selected samples
        taxaDF = pd.merge(profileDF, picrustDF, left_index=True, right_index=True, how='inner')

        for level in levelList:
            if DepVar == 1:
                taxaDF[level] = taxaDF['rel_abund'] * taxaDF[level]
            elif DepVar == 2:
                taxaDF[level] = np.where(taxaDF['rel_abund'] * taxaDF[level] > 0, 1, 0)
            elif DepVar == 3:
                taxaDF[level] = taxaDF['rel_abund'] * taxaDF[level]
                taxaDF[level] = taxaDF[level].div(taxaDF[level].sum(), axis=0)
                taxaDF[level] = taxaDF[level].apply(lambda x: -1 * x * math.log(x) if x > 0 else 0)
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
        elif DepVar == 2:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rich')
        elif DepVar == 3:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='diversity')
        elif DepVar == 4:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund_16S')

        wanted = allFields + ['sampleid']
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
            myDict = {'error': "Error with ANcOVA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getNZDF(nzAll, myDict, savedDF, tempDF, allFields, DepVar, RID, stops, PID):
    global base
    try:
        nzDict = {}
        if nzAll == 0:
            for key in myDict:
                keggList = myDict[key]
                if isinstance(keggList, unicode):
                    if key == 'Level1':
                        nzList = nz_entry.objects.using('picrust').filter(nz_lvl1_id_id=keggList).values_list('nz_orthology', flat=True)
                        if nzList:
                            nzDict[keggList] = nzList
                    elif key == 'Level2':
                        nzList = nz_entry.objects.using('picrust').filter(nz_lvl2_id_id=keggList).values_list('nz_orthology', flat=True)
                        if nzList:
                            nzDict[keggList] = nzList
                    elif key == 'Level3':
                        nzList = nz_entry.objects.using('picrust').filter(nz_lvl3_id_id=keggList).values_list('nz_orthology', flat=True)
                        if nzList:
                            nzDict[keggList] = nzList
                    elif key == 'Level4':
                        nzList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=keggList).values_list('nz_orthology', flat=True)
                        if nzList:
                            nzDict[keggList] = nzList
                    elif key == 'Level5':
                        nzList = nz_entry.objects.using('picrust').filter(nz_lvl5_id=keggList).values_list('nz_orthology', flat=True)
                        if nzList:
                            nzDict[keggList] = nzList
                else:
                    for i in keggList:
                        if key == 'Level1':
                            nzList = nz_entry.objects.using('picrust').filter(nz_lvl1_id_id=i).values_list('nz_orthology', flat=True)
                            if nzList:
                                nzDict[i] = nzList
                        elif key == 'Level2':
                            nzList = nz_entry.objects.using('picrust').filter(nz_lvl2_id_id=i).values_list('nz_orthology', flat=True)
                            if nzList:
                                nzDict[i] = nzList
                        elif key == 'Level3':
                            nzList = nz_entry.objects.using('picrust').filter(nz_lvl3_id_id=i).values_list('nz_orthology', flat=True)
                            if nzList:
                                nzDict[i] = nzList
                        elif key == 'Level4':
                            nzList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=i).values_list('nz_orthology', flat=True)
                            if nzList:
                                nzDict[i] = nzList
                        elif key == 'Level5':
                            nzList = nz_entry.objects.using('picrust').filter(nz_lvl5_id=i).values_list('nz_orthology', flat=True)
                            if nzList:
                                nzDict[i] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 1:
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

        path = 'myPhyloDB/media/temp/anova/' + str(RID)
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
            processes = [threading.Thread(target=sumStuff, args=(listDF[x], nzDict, RID, x, PID, stops)) for x in xrange(numcore)]
        else:
            numcore = min(local_cfg.usr_numcore, int(maxCPU))
            listDF = np.array_split(picrustDF, numcore)
            processes = [threading.Thread(target=sumStuff, args=(listDF[x], nzDict, RID, x, PID, stops)) for x in xrange(numcore)]

        for p in processes:
            p.start()

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
            path = 'myPhyloDB/media/temp/anova/'+str(RID)+'/file%d.temp' % i
            frame = pd.read_csv(path)
            picrustDF = picrustDF.append(frame, ignore_index=True)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        shutil.rmtree('myPhyloDB/media/temp/anova/'+str(RID))
        picrustDF.set_index('speciesid', inplace=True)
        picrustDF[picrustDF > 0] = 1

        # merge to get final gene counts for all selected samples
        taxaDF = pd.merge(profileDF, picrustDF, left_index=True, right_index=True, how='inner')
        for level in levelList:
            if DepVar == 1:
                taxaDF[level] = taxaDF['rel_abund'] * taxaDF[level]
            elif DepVar == 2:
                taxaDF[level] = np.where(taxaDF['rel_abund'] * taxaDF[level] > 0, 1, 0)
            elif DepVar == 3:
                taxaDF[level] = taxaDF['rel_abund'] * taxaDF[level]
                taxaDF[level] = taxaDF[level].div(taxaDF[level].sum(), axis=0)
                taxaDF[level] = taxaDF[level].apply(lambda x: -1 * x * math.log(x) if x > 0 else 0)
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
        elif DepVar == 2:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rich')
        elif DepVar == 3:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='diversity')
        elif DepVar == 4:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund_16S')

        wanted = allFields + ['sampleid']
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
            myDict = {'error': "Error with ANcOVA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def sumStuff(slice, koDict, RID, num, PID, stops):
    db.close_old_connections()

    f = open('myPhyloDB/media/temp/anova/'+str(RID)+'/file'+str(num)+".temp", 'w')

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


def getTabANOVA(request):
    if request.is_ajax():
        RID = request.GET["all"]

        myDir = 'myPhyloDB/media/temp/anova/'
        fileName = str(myDir) + str(RID) + '.pkl'
        savedDF = pd.read_pickle(fileName)

        myDir = 'myPhyloDB/media/temp/anova/'
        fileName = str(myDir) + str(RID) + '.csv'
        savedDF.to_csv(fileName)

        myDict = {}
        myDir = 'temp/anova/'
        fileName = str(myDir) + str(RID) + '.csv'
        myDict['name'] = str(fileName)
        res = simplejson.dumps(myDict)

        return HttpResponse(res, content_type='application/json')


def removeANOVAFiles(request):
    if request.is_ajax():
        RID = request.GET["all"]

        file = "myPhyloDB/media/temp/anova/" + str(RID) + ".pkl"
        if os.path.exists(file):
            os.remove(file)

        file = "myPhyloDB/media/temp/anova/" + str(RID) + ".csv"
        if os.path.exists(file):
            os.remove(file)

        return HttpResponse()
