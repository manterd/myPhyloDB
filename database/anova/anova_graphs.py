import datetime
from django.http import HttpResponse
from django.db.models import Sum
import logging
import numpy as np
import pandas as pd
import pickle
from pyper import *
from scipy import stats
import simplejson

from database.anova.anova_DF import UnivMetaDF, normalizeUniv
from database.models import Sample, Profile
from database.utils import multidict, taxaProfileDF, stoppableThread


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}
stops = {}
LOG_FILENAME = 'error_log.txt'


def statusANOVA(request):
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
            stage[RID] = str(base[RID]) + '<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
        except:
            stage[RID] = '<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
        myDict['stage'] = stage[RID]
        json_data = simplejson.dumps(myDict, encoding="Latin-1")
        return HttpResponse(json_data, content_type='application/json')


def removeRIDANOVA(request):
    global base, stage, time1, time2, TimeDiff
    try:
        if request.is_ajax():
            RID = request.GET["all"]
            base.pop(RID, None)
            stage.pop(RID, None)
            time1.pop(RID, None)
            time2.pop(RID, None)
            TimeDiff.pop(RID, None)
            stops.pop(RID, None)
            return True
        else:
            return False
    except:
        return False


thread1 = stoppableThread()
stop1 = False


def stopANOVA(request):
    global stops
    if request.is_ajax():
        RID = request.GET["all"]
        stops[RID] = True
        myDict = {}
        myDict['error'] = 'Your analysis has been stopped!'
        res = simplejson.dumps(myDict)
        return HttpResponse(res, content_type='application/json')


def getCatUnivData(request):
    global res, thread1, stop1
    if request.is_ajax():
        thread1 = stoppableThread(target=loopCat, args=(request,))
        thread1.start()
        thread1.join()
        stop1 = False
        return HttpResponse(res, content_type='application/json')


def getQuantUnivData(request):
    global res, thread1, stop1
    if request.is_ajax():
        thread1 = stoppableThread(target=loopQuant, args=(request,))
        thread1.start()
        thread1.join()
        stop1 = False
        return HttpResponse(res, content_type='application/json')


def loopCat(request):
    global res, base, stage, time1, TimeDiff, stop1, stops
    try:
        while True:
            if request.is_ajax():
                # Get variables from web page
                allJson = request.GET["all"]
                all = simplejson.loads(allJson)

                RID = str(all["RID"])
                stops[RID] = False
                time1[RID] = time.time()  # Moved these down here so RID is available
                base[RID] = 'Step 1 of 4: Selecting your chosen meta-variables...'

                # Get normalized data from cookie
                #savedDF = pickle.loads(request.session['savedDF'])

                path = pickle.loads(request.session['savedDF'])
                savedDF = pd.read_pickle(path)

                selectAll = int(all["selectAll"])
                DepVar = int(all["DepVar"])
                sig_only = int(all["sig_only"])

                # Select samples and meta-variables from savedDF
                catSampleIDs = list(set(all["metaIDsCat"]))
                catFields = list(set(all["metaFieldsCat"]))
                catValues = list(set(all["metaValsCat"]))
                quantSampleIDs = list(set(all["metaIDsQuant"]))
                quantFields = list(set(all["metaFieldsQuant"]))
                quantValues = list(set(all["metaValsQuant"]))

                allSampleIDs = list(set(catSampleIDs + quantSampleIDs))
                allFields = list(set(catFields + quantFields))
                allValues = list(set(catValues + quantValues))

                # Removes samples (rows) that are not in our samplelist
                metaDF = savedDF.loc[savedDF['sampleid'].isin(allSampleIDs)]
                metaDF = metaDF[allFields]

                result = ''
                result += 'Categorical variables selected: ' + ", ".join(catFields) + '\n'
                result += 'Quantitative variables selected: ' + ", ".join(quantFields) + '\n'
                result += '===============================================\n'

                base[RID] = 'Step 1 of 4: Selecting your chosen meta-variables...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                if stops[RID]:
                    print "Received stop code"
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                base[RID] = 'Step 2 of 4: Selecting your chosen taxa...'

                # Select only the taxa of interest if user used the taxa tree
                taxaString = all["taxa"]
                taxaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(taxaString)

                # get selected taxa fro each rank selected in the tree
                taxaDF = pd.DataFrame(columns=['sampleid', 'rank', 'taxa_id', 'taxa_name', 'abund', 'abund_16S', 'rich', 'diversity'])
                if selectAll == 0:
                    for key in taxaDict:
                        taxaList = taxaDict[key]
                        if isinstance(taxaList, unicode):
                            if key == 'Kingdom':
                                tempDF = savedDF.loc[savedDF['kingdomid'] == taxaList]
                                tempDF = tempDF[['sampleid', 'kingdomid', 'kingdomName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'kingdomid': 'taxa_id', 'kingdomName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Kingdom'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Phyla':
                                tempDF = savedDF.loc[savedDF['phylaid'] == taxaList]
                                tempDF = tempDF[['sampleid', 'phylaid', 'phylaName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'phylaid': 'taxa_id', 'phylaName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Phyla'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Class':
                                tempDF = savedDF.loc[savedDF['classid'] == taxaList]
                                tempDF = tempDF[['sampleid', 'classid', 'className', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'classid': 'taxa_id', 'className': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Class'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Order':
                                tempDF = savedDF.loc[savedDF['orderid'] == taxaList]
                                tempDF = tempDF[['sampleid', 'orderid', 'orderName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'orderid': 'taxa_id', 'orderName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Order'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Family':
                                tempDF = savedDF.loc[savedDF['familyid'] == taxaList]
                                tempDF = tempDF[['sampleid', 'familyid', 'familyName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'familyid': 'taxa_id', 'familyName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Family'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Genus':
                                tempDF = savedDF.loc[savedDF['genusid'] == taxaList]
                                tempDF = tempDF[['sampleid', 'genusid', 'genusName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'genusid': 'taxa_id', 'genusName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Genus'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Species':
                                tempDF = savedDF.loc[savedDF['speciesid'] == taxaList]
                                tempDF = tempDF[['sampleid', 'speciesid', 'speciesName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'speciesid': 'taxa_id', 'speciesName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Species'
                                taxaDF = taxaDF.append(tempDF)
                        else:
                            if key == 'Kingdom':
                                tempDF = savedDF.loc[savedDF['kingdomid'].isin(taxaList)]
                                tempDF = tempDF[['sampleid', 'kingdomid', 'kingdomName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'kingdomid': 'taxa_id', 'kingdomName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Kingdom'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Phyla':
                                tempDF = savedDF.loc[savedDF['phylaid'].isin(taxaList)]
                                tempDF = tempDF[['sampleid', 'phylaid', 'phylaName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'phylaid': 'taxa_id', 'phylaName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Phyla'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Class':
                                tempDF = savedDF.loc[savedDF['classid'].isin(taxaList)]
                                tempDF = tempDF[['sampleid', 'classid', 'className', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'classid': 'taxa_id', 'className': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Class'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Order':
                                tempDF = savedDF.loc[savedDF['orderid'].isin(taxaList)]
                                tempDF = tempDF[['sampleid', 'orderid', 'orderName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'orderid': 'taxa_id', 'orderName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Order'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Family':
                                tempDF = savedDF.loc[savedDF['familyid'].isin(taxaList)]
                                tempDF = tempDF[['sampleid', 'familyid', 'familyName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'familyid': 'taxa_id', 'familyName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Family'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Genus':
                                tempDF = savedDF.loc[savedDF['genusid'].isin(taxaList)]
                                tempDF = tempDF[['sampleid', 'genusid', 'genusName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'genusid': 'taxa_id', 'genusName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Genus'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Species':
                                tempDF = savedDF.loc[savedDF['speciesid'].isin(taxaList)]
                                tempDF = tempDF[['sampleid', 'speciesid', 'speciesName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'speciesid': 'taxa_id', 'speciesName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Species'
                                taxaDF = taxaDF.append(tempDF)
                elif selectAll == 1:
                    taxaDF = savedDF[['sampleid', 'kingdomid', 'kingdomName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'kingdomid': 'taxa_id', 'kingdomName': 'taxa_name'}, inplace=True)
                    taxaDF['rank'] = 'Kingdom'
                elif selectAll == 2:
                    taxaDF = savedDF[['sampleid', 'phylaid', 'phylaName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'phylaid': 'taxa_id', 'phylaName': 'taxa_name'}, inplace=True)
                    taxaDF['rank'] = 'Phyla'
                elif selectAll == 3:
                    taxaDF = savedDF[['sampleid', 'classid', 'className', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'classid': 'taxa_id', 'className': 'taxa_name'}, inplace=True)
                    taxaDF['rank'] = 'Class'
                elif selectAll == 4:
                    taxaDF = savedDF[['sampleid', 'orderid', 'orderName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'orderid': 'taxa_id', 'orderName': 'taxa_name'}, inplace=True)
                    taxaDF['rank'] = 'Order'
                elif selectAll == 5:
                    taxaDF = savedDF[['sampleid', 'familyid', 'familyName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'familyid': 'taxa_id', 'familyName': 'taxa_name'}, inplace=True)
                    taxaDF['rank'] = 'Family'
                elif selectAll == 6:
                    taxaDF = savedDF[['sampleid', 'genusid', 'genusName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'genusid': 'taxa_id', 'genusName': 'taxa_name'}, inplace=True)
                    taxaDF['rank'] = 'Genus'
                elif selectAll == 7:
                    taxaDF = savedDF[['sampleid', 'speciesid', 'speciesName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'speciesid': 'taxa_id', 'speciesName': 'taxa_name'}, inplace=True)
                    taxaDF['rank'] = 'Species'

                finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')
                wantedList = allFields + ['sampleid', 'rank', 'taxa_name', 'taxa_id']
                finalDF = finalDF.groupby(wantedList)[['abund', 'abund_16S', 'rich', 'diversity']].sum()
                finalDF.reset_index(drop=False, inplace=True)

                base[RID] = 'Step 2 of 4: Selecting your chosen taxa...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                if stops[RID]:
                    print "Received stop code"
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

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
                grouped1 = finalDF.groupby(['rank', 'taxa_name', 'taxa_id'])
                pValDict = {}
                counter = 1
                for name1, group1 in grouped1:
                    D = ''
                    p_val = 1.0

                    if os.name == 'nt':
                        r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                    else:
                        r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                    r.assign("df", group1)
                    trtString = " * ".join(allFields)
                    if DepVar == 1:
                        anova_string = "fit <- aov(abund ~ " + str(trtString) + ", data=df)"
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
                        anova_string = "fit <- aov(copies ~ " + str(trtString) + ", data=df)"
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
                    else:
                        p_val = 1.0
                        D = 'ANOVA cannot be performed, please check that you have more than one treatment level and appropriate replication.\n'

                    pValDict[name1] = p_val

                    result += 'Taxa level: ' + str(name1[0]) + '\n'
                    result += 'Taxa name: ' + str(name1[1]) + '\n'
                    result += 'Taxa ID: ' + str(name1[2]) + '\n'
                    if DepVar == 1:
                        result += 'Dependent Variable: Relative Abundance (proportion)' + '\n'
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
                    if stops[RID]:
                        print "Received stop code"
                        return None
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                base[RID] = 'Step 3 of 4: Performing statistical test...done!'
                base[RID] = 'Step 4 of 4: Formatting graph data for display...'


                grouped1 = finalDF.groupby(['rank', 'taxa_name', 'taxa_id'])
                for name1, group1 in grouped1:
                    dataList = []
                    pValue = pValDict[name1]

                    if sig_only == 0:
                        if DepVar == 1:
                            grouped2 = group1.groupby(catFields)['abund'].mean()
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

                    elif sig_only == 1:
                        if pValue < 0.05:
                            if DepVar == 1:
                                grouped2 = group1.groupby(catFields)['abund'].mean()
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

                    colors_idx += 1
                    if colors_idx >= len(colors):
                        colors_idx = 0

                    grouped2 = group1.groupby(catFields)['abund'].mean()
                    if catFields.__len__() == 1:
                        xAxisDict['categories'] = catValues
                    else:
                        g2indexvals = grouped2.index.values
                        level = g2indexvals[0].__len__()
                        labelTree = recLabels(g2indexvals, level)
                        xAxisDict['categories'] = [labelTree]

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                    if stops[RID]:
                        print "Received stop code"
                        return None
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                yTitle = {}
                if DepVar == 1:
                    yTitle['text'] = 'Abundance'
                elif DepVar == 2:
                    yTitle['text'] = 'Species Richness'
                elif DepVar == 3:
                    yTitle['text'] = 'Species Diversity'
                elif DepVar == 4:
                    yTitle['text'] = 'Abundance (rRNA gene copies)'
                yAxisDict['title'] = yTitle

                finalDict['series'] = seriesList
                finalDict['xAxis'] = xAxisDict
                finalDict['yAxis'] = yAxisDict
                finalDict['text'] = result

                if not seriesList:
                    finalDict['empty'] = 0
                else:
                    finalDict['empty'] = 1

                base[RID] = 'Step 4 of 4: Formatting graph data for display...done!'

                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                return None

    except:
        if not stop1:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with ANcOVA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            raise


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


def loopQuant(request):
    global res, base, time1, TimeDiff, stop1, stops
    try:
        while True:
            if request.is_ajax():
                # Get variables from web page
                allJson = request.GET["all"]
                all = simplejson.loads(allJson)

                RID = str(all["RID"])
                stops[RID] = False
                time1[RID] = time.time()
                base[RID] = 'Step 1 of 4: Querying database...'

                # Get normalized data from cookie
                #savedDF = pickle.loads(request.session['savedDF'])

                path = pickle.loads(request.session['savedDF'])
                savedDF = pd.read_pickle(path)

                selectAll = int(all["selectAll"])
                DepVar = int(all["DepVar"])
                sig_only = int(all["sig_only"])

                # Select samples and meta-variables from savedDF
                catSampleIDs = list(set(all["metaIDsCat"]))
                catFields = list(set(all["metaFieldsCat"]))
                catValues = list(set(all["metaValsCat"]))
                quantSampleIDs = list(set(all["metaIDsQuant"]))
                quantFields = list(set(all["metaFieldsQuant"]))
                quantValues = list(set(all["metaValsQuant"]))

                allSampleIDs = list(set(catSampleIDs + quantSampleIDs))
                allFields = list(set(catFields + quantFields))
                allValues = list(set(catValues + quantValues))

                # Removes samples (rows) that are not in our samplelist
                metaDF = savedDF.loc[savedDF['sampleid'].isin(allSampleIDs)]
                metaDF = metaDF[allFields]

                result = ''
                result += 'Categorical variables selected: ' + ", ".join(catFields) + '\n'
                result += 'Quantitative variables selected: ' + ", ".join(quantFields) + '\n'
                result += '===============================================\n'

                base[RID] = 'Step 1 of 4: Selecting your chosen meta-variables...done'


                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                if stops[RID]:
                    print "Received stop code"
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                base[RID] = 'Step 2 of 4: Selecting your chosen taxa...'

                # Select only the taxa of interest if user used the taxa tree
                taxaString = all["taxa"]
                taxaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(taxaString)

                # get selected taxa fro each rank selected in the tree
                taxaDF = pd.DataFrame(columns=['sampleid', 'rank', 'taxa_id', 'taxa_name', 'abund', 'abund_16S', 'rich', 'diversity'])
                if selectAll == 0:
                    for key in taxaDict:
                        taxaList = taxaDict[key]
                        if isinstance(taxaList, unicode):
                            if key == 'Kingdom':
                                tempDF = savedDF.loc[savedDF['kingdomid'] == taxaList]
                                tempDF = tempDF[['sampleid', 'kingdomid', 'kingdomName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'kingdomid': 'taxa_id', 'kingdomName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Kingdom'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Phyla':
                                tempDF = savedDF.loc[savedDF['phylaid'] == taxaList]
                                tempDF = tempDF[['sampleid', 'phylaid', 'phylaName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'phylaid': 'taxa_id', 'phylaName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Phyla'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Class':
                                tempDF = savedDF.loc[savedDF['classid'] == taxaList]
                                tempDF = tempDF[['sampleid', 'classid', 'className', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'classid': 'taxa_id', 'className': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Class'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Order':
                                tempDF = savedDF.loc[savedDF['orderid'] == taxaList]
                                tempDF = tempDF[['sampleid', 'orderid', 'orderName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'orderid': 'taxa_id', 'orderName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Order'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Family':
                                tempDF = savedDF.loc[savedDF['familyid'] == taxaList]
                                tempDF = tempDF[['sampleid', 'familyid', 'familyName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'familyid': 'taxa_id', 'familyName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Family'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Genus':
                                tempDF = savedDF.loc[savedDF['genusid'] == taxaList]
                                tempDF = tempDF[['sampleid', 'genusid', 'genusName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'genusid': 'taxa_id', 'genusName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Genus'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Species':
                                tempDF = savedDF.loc[savedDF['speciesid'] == taxaList]
                                tempDF = tempDF[['sampleid', 'speciesid', 'speciesName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'speciesid': 'taxa_id', 'speciesName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Species'
                                taxaDF = taxaDF.append(tempDF)
                        else:
                            if key == 'Kingdom':
                                tempDF = savedDF.loc[savedDF['kingdomid'].isin(taxaList)]
                                tempDF = tempDF[['sampleid', 'kingdomid', 'kingdomName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'kingdomid': 'taxa_id', 'kingdomName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Kingdom'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Phyla':
                                tempDF = savedDF.loc[savedDF['phylaid'].isin(taxaList)]
                                tempDF = tempDF[['sampleid', 'phylaid', 'phylaName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'phylaid': 'taxa_id', 'phylaName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Phyla'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Class':
                                tempDF = savedDF.loc[savedDF['classid'].isin(taxaList)]
                                tempDF = tempDF[['sampleid', 'classid', 'className', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'classid': 'taxa_id', 'className': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Class'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Order':
                                tempDF = savedDF.loc[savedDF['orderid'].isin(taxaList)]
                                tempDF = tempDF[['sampleid', 'orderid', 'orderName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'orderid': 'taxa_id', 'orderName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Order'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Family':
                                tempDF = savedDF.loc[savedDF['familyid'].isin(taxaList)]
                                tempDF = tempDF[['sampleid', 'familyid', 'familyName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'familyid': 'taxa_id', 'familyName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Family'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Genus':
                                tempDF = savedDF.loc[savedDF['genusid'].isin(taxaList)]
                                tempDF = tempDF[['sampleid', 'genusid', 'genusName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'genusid': 'taxa_id', 'genusName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Genus'
                                taxaDF = taxaDF.append(tempDF)
                            elif key == 'Species':
                                tempDF = savedDF.loc[savedDF['speciesid'].isin(taxaList)]
                                tempDF = tempDF[['sampleid', 'speciesid', 'speciesName', 'abund', 'abund_16S', 'rich', 'diversity']]
                                tempDF.rename(columns={'speciesid': 'taxa_id', 'speciesName': 'taxa_name'}, inplace=True)
                                tempDF['rank'] = 'Species'
                                taxaDF = taxaDF.append(tempDF)
                elif selectAll == 1:
                    taxaDF = savedDF[['sampleid', 'kingdomid', 'kingdomName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'kingdomid': 'taxa_id', 'kingdomName': 'taxa_name'}, inplace=True)
                    taxaDF['rank'] = 'Kingdom'
                elif selectAll == 2:
                    taxaDF = savedDF[['sampleid', 'phylaid', 'phylaName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'phylaid': 'taxa_id', 'phylaName': 'taxa_name'}, inplace=True)
                    taxaDF['rank'] = 'Phyla'
                elif selectAll == 3:
                    taxaDF = savedDF[['sampleid', 'classid', 'className', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'classid': 'taxa_id', 'className': 'taxa_name'}, inplace=True)
                    taxaDF['rank'] = 'Class'
                elif selectAll == 4:
                    taxaDF = savedDF[['sampleid', 'orderid', 'orderName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'orderid': 'taxa_id', 'orderName': 'taxa_name'}, inplace=True)
                    taxaDF['rank'] = 'Order'
                elif selectAll == 5:
                    taxaDF = savedDF[['sampleid', 'familyid', 'familyName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'familyid': 'taxa_id', 'familyName': 'taxa_name'}, inplace=True)
                    taxaDF['rank'] = 'Family'
                elif selectAll == 6:
                    taxaDF = savedDF[['sampleid', 'genusid', 'genusName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'genusid': 'taxa_id', 'genusName': 'taxa_name'}, inplace=True)
                    taxaDF['rank'] = 'Genus'
                elif selectAll == 7:
                    taxaDF = savedDF[['sampleid', 'speciesid', 'speciesName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'speciesid': 'taxa_id', 'speciesName': 'taxa_name'}, inplace=True)
                    taxaDF['rank'] = 'Species'

                finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')
                wantedList = allFields + ['sampleid', 'rank', 'taxa_name', 'taxa_id']
                finalDF = finalDF.groupby(wantedList)[['abund', 'abund_16S', 'rich', 'diversity']].sum()
                finalDF.reset_index(drop=False, inplace=True)

                base[RID] = 'Step 2 of 4: Selecting your chosen taxa...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
                if stops[RID]:
                    print "Received stop code"
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                base[RID] = 'Step 3 of 4: Performing statistical test...!'
                finalDict = {}
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
                shapes = ['circle', 'square', 'triangle', 'triangle-down', 'diamond']
                grouped1 = finalDF.groupby(['rank', 'taxa_name', 'taxa_id'])
                pValDict = {}
                p_vals = ''
                counter = 1
                for name1, group1 in grouped1:
                    D = ''
                    p_val = 1.0

                    if os.name == 'nt':
                        r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                    else:
                        r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                    r.assign("df", group1)
                    trtString = " * ".join(allFields)

                    if DepVar == 1:
                        anova_string = "fit <- lm(abund ~ " + str(trtString) + ", data=df)"
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
                        anova_string = "fit <- lm(copies ~ " + str(trtString) + ", data=df)"
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

                    else:
                        p_value = 1.0
                        D = 'ANOVA cannot be performed, please check that you have more than one treatment level and appropriate replication.\n'

                    pValDict[name1] = p_val

                    result += 'Taxa level: ' + str(name1[0]) + '\n'
                    result += 'Taxa name: ' + str(name1[1]) + '\n'
                    result += 'Taxa ID: ' + str(name1[2]) + '\n'
                    if DepVar == 1:
                        result += 'Dependent Variable: Relative Abundance (proportion)' + '\n'
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
                    if stops[RID]:
                        print "Received stop code"
                        return None
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                base[RID] = 'Step 3 of 4: Performing statistical test...done!'
                base[RID] = 'Step 4 of 4: Formatting graph data for display...'

                seriesList = []
                grouped1 = finalDF.groupby(['rank', 'taxa_name', 'taxa_id'])
                for name1, group1 in grouped1:
                    pValue = pValDict[name1]
                    shapes_idx = 0

                    if len(catFields) > 0:
                        grouped2 = group1.groupby(catFields)
                        for name2, group2 in grouped2:
                            dataList = []
                            x = []
                            y = []
                            if sig_only == 0:
                                if DepVar == 1:
                                    dataList = group2[[quantFields[0], 'abund']].values.astype(float).tolist()
                                    x = group2[quantFields[0]].values.astype(float).tolist()
                                    y = group2['abund'].values.astype(float).tolist()
                                elif DepVar == 2:
                                    dataList = group1.groupby(catFields)[[quantFields[0], 'rich']].astype(float).tolist()
                                    x = group1.groupby(catFields)[quantFields[0]].values.astype(float).tolist()
                                    y = group1.groupby(catFields)['rich'].values.astype(float).tolist()
                                elif DepVar == 3:
                                    dataList = group1.groupby(catFields)[[quantFields[0], 'diversity']].values.astype(float).tolist()
                                    x = group1.groupby(catFields)[quantFields[0]].values.astype(float).tolist()
                                    y = group1.groupby(catFields)['diversity'].values.astype(float).tolist()
                                elif DepVar == 4:
                                    dataList = group1.groupby(catFields)[[quantFields[0], 'abund_16S']].values.astype(float).tolist()
                                    x = group1.groupby(catFields)[quantFields[0]].values.astype(float).tolist()
                                    y = group1.groupby(catFields)['abund_16S'].values.astype(float).tolist()

                            elif sig_only == 1:
                                if pValue < 0.05:
                                    if DepVar == 1:
                                        dataList = group1.groupby(catFields)[[quantFields[0], 'abund']].values.astype(float).tolist()
                                        x = group1.groupby(catFields)[quantFields[0]].values.astype(float).tolist()
                                        y = group1.groupby(catFields)['abund'].values.astype(float).tolist()
                                    elif DepVar == 2:
                                        dataList = group1.groupby(catFields)[[quantFields[0], 'rich']].values.astype(float).tolist()
                                        x = group1.groupby(catFields)[quantFields[0]].values.astype(float).tolist()
                                        y = group1.groupby(catFields)['rich'].values.astype(float).tolist()
                                    elif DepVar == 3:
                                        dataList = group1.groupby(catFields)[[quantFields[0], 'diversity']].values.astype(float).tolist()
                                        x = group1.groupby(catFields)[quantFields[0]].values.astype(float).tolist()
                                        y = group1.groupby(catFields)['diversity'].values.astype(float).tolist()
                                    elif DepVar == 4:
                                        dataList = group1.groupby(catFields)[[quantFields[0], 'abund_16S']].values.astype(float).tolist()
                                        x = group1.groupby(catFields)[quantFields[0]].values.astype(float).tolist()
                                        y = group1.groupby(catFields)['abund_16S'].values.astype(float).tolist()

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

                            if not np.isnan(p_vals).any():
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

                    # start loop if catFields == 0
                    else:
                        if sig_only == 0:
                            if DepVar == 1:
                                dataList = group1[[quantFields[0], 'abund']].values.astype(float).tolist()
                                x = group1[quantFields[0]].values.astype(float).tolist()
                                y = group1['abund'].values.astype(float).tolist()
                            elif DepVar == 2:
                                dataList = group1[[quantFields[0], 'rich']].values.astype(float).tolist()
                                x = group1[quantFields[0]].values.astype(float).tolist()
                                y = group1['rich'].values.astype(float).tolist()
                            elif DepVar == 3:
                                dataList = group1[[quantFields[0], 'diversity']].values.astype(float).tolist()
                                x = group1[quantFields[0]].values.astype(float).tolist()
                                y = group1['diversity'].values.astype(float).tolist()
                            elif DepVar == 4:
                                dataList = group1[[quantFields[0], 'abund_16S']].values.astype(float).tolist()
                                x = group1[quantFields[0]].values.astype(float).tolist()
                                y = group1['abund_16S'].values.astype(float).tolist()

                        elif sig_only == 1:
                            if pValue < 0.05:
                                if DepVar == 1:
                                    dataList = group1[[quantFields[0], 'abund']].values.astype(float).tolist()
                                    x = group1[quantFields[0]].values.astype(float).tolist()
                                    y = group1['abund'].values.astype(float).tolist()
                                elif DepVar == 2:
                                    dataList = group1[[quantFields[0], 'rich']].values.astype(float).tolist()
                                    x = group1[quantFields[0]].values.astype(float).tolist()
                                    y = group1['rich'].values.astype(float).tolist()
                                elif DepVar == 3:
                                    dataList = group1[[quantFields[0], 'diversity']].values.astype(float).tolist()
                                    x = group1[quantFields[0]].values.astype(float).tolist()
                                    y = group1['diversity'].values.astype(float).tolist()
                                elif DepVar == 4:
                                    dataList = group1[[quantFields[0], 'abund_16S']].values.astype(float).tolist()
                                    x = group1[quantFields[0]].values.astype(float).tolist()
                                    y = group1['abund_16S'].values.astype(float).tolist()

                        if not np.isnan(p_vals).any():
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
                    if stops[RID]:
                        print "Received stop code"
                        return None
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

                xAxisDict = {}
                xTitle = {}
                xTitle['text'] = quantFields[0]
                xAxisDict['title'] = xTitle

                yAxisDict = {}
                yTitle = {}
                if DepVar == 1:
                    yTitle['text'] = 'Abundance'
                elif DepVar == 2:
                    yTitle['text'] = 'Species Richness'
                elif DepVar == 3:
                    yTitle['text'] = 'Shannon Diversity'
                elif DepVar == 4:
                    yTitle['text'] = 'Abundance (rRNA gene copies)'
                yAxisDict['title'] = yTitle

                finalDict['series'] = seriesList
                finalDict['xAxis'] = xAxisDict
                finalDict['yAxis'] = yAxisDict
                if not seriesList:
                    finalDict['empty'] = 0
                else:
                    finalDict['empty'] = 1

                base[RID] = 'Step 4 of 4: Formatting graph data for display...done!'
                finalDict['text'] = result
                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                return None

    except:
        if not stop1:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with ANcOVA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            raise


