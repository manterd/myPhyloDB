import ast
import datetime
from django import db
from django.http import HttpResponse
from django_pandas.io import read_frame
import logging
import math
import multiprocessing as mp
import numpy as np
import pandas as pd
from pyper import *
import simplejson
import shutil
import threading

from database.models import PICRUSt
from database.models import ko_lvl1, ko_lvl2, ko_lvl3, ko_entry
from database.models import nz_lvl1, nz_lvl2, nz_lvl3, nz_lvl4, nz_entry
from database.utils import multidict, wOdum
import database.queue
from config import local_cfg


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}

LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def statusPCoA(request):
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


def removeRIDPCoA(RID):
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


def getPCoA(request, stops, RID, PID):
    global base, stage, time1, TimeDiff
    try:
        while True:
            if request.is_ajax():
                allJson = request.GET["all"]
                all = simplejson.loads(allJson)

                time1[RID] = time.time()
                base[RID] = 'Step 1 of 8: Selecting your chosen meta-variables...'

                myDir = 'media/usr_temp/' + str(request.user) + '/'
                path = str(myDir) + 'usr_norm_data.csv'

                with open(path, 'rb') as f:
                    savedDF = pd.read_csv(f, index_col=0, sep='\t')

                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])

                distance = int(all["distance"])
                PC1 = int(all["PC1"])
                PC2 = int(all["PC2"])
                test = int(all["test"])
                alpha = float(all["alpha"])
                perms = int(all["perms"])

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

                if distance == 1:
                    result += 'Distance score: Manhattan' + '\n'
                elif distance == 2:
                    result += 'Distance score: Euclidean' + '\n'
                elif distance == 3:
                    result += 'Distance score: Canberra' + '\n'
                elif distance == 4:
                    result += 'Distance score: Bray-Curtis' + '\n'
                elif distance == 5:
                    result += 'Distance score: Kulczynski' + '\n'
                elif distance == 6:
                    result += 'Distance score: Jaccard' + '\n'
                elif distance == 7:
                    result += 'Distance score: Gower' + '\n'
                elif distance == 8:
                    result += 'Distance score: altGower' + '\n'
                elif distance == 9:
                    result += 'Distance score: Morisita' + '\n'
                elif distance == 10:
                    result += 'Distance score: Horn' + '\n'
                elif distance == 11:
                    result += 'Distance score: Mountford' + '\n'
                elif distance == 12:
                    result += 'Distance score: Binomial' + '\n'
                elif distance == 13:
                    result += 'Distance score: Chao' + '\n'
                elif distance == 14:
                    result += 'Distance score: Cao' + '\n'
                elif distance == 15:
                    result += 'Distance score: wOdum' + '\n'
                    result += 'alpha: ' + str(alpha) + '\n'

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

                if not catFields_edit:
                    myDict = {}
                    myDict['error'] = "Selected meta data only has one level.\nPlease select a different variable."
                    res = simplejson.dumps(myDict)
                    return None

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
                tempDF = savedDF.loc[savedDF['sampleid'].isin(allSampleIDs)]

                if metaDictCat:
                    for key in metaDictCat:
                        tempDF = tempDF.loc[tempDF[key].isin(metaDictCat[key])]

                if metaDictQuant:
                    for key in metaDictQuant:
                        valueList = [float(x) for x in metaDictQuant[key]]
                        tempDF = tempDF.loc[tempDF[key].isin(valueList)]

                wantedList = allFields + ['sampleid']
                metaDF = tempDF[wantedList]

                result += 'Categorical variables selected by user: ' + ", ".join(catFields) + '\n'
                result += 'Categorical variables removed from analysis (contains only 1 level): ' + ", ".join(removed) + '\n'
                result += 'Quantitative variables selected by user: ' + ", ".join(quantFields) + '\n'
                result += '===============================================\n'

                base[RID] = 'Step 1 of 8: Selecting your chosen meta-variables...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 2 of 8: Selecting your chosen taxa or KEGG level...'

                DepVar = 1
                finalDF = pd.DataFrame()
                if button3 == 1:
                    DepVar = int(all["DepVar_taxa"])
                    finalDF = getTaxaDF(selectAll, savedDF, metaDF, allFields, DepVar, RID, stops, PID)

                if button3 == 2:
                    DepVar = int(all["DepVar_kegg"])
                    finalDF = getKeggDF(keggAll, savedDF, tempDF, allFields, DepVar, RID, stops, PID)

                if button3 == 3:
                    DepVar = int(all["DepVar_nz"])
                    finalDF = getNZDF(nzAll, savedDF, tempDF, allFields, DepVar, RID, stops, PID)

                # save location info to session
                myDir = 'media/temp/pcoa/'
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
                meta_rDF = meta_rDF.loc[meta_rDF['sampleid'].isin(catSampleIDs)]

                if metaDictCat:
                    for key in metaDictCat:
                        meta_rDF = meta_rDF.loc[meta_rDF[key].isin(metaDictCat[key])]

                if metaDictQuant:
                    for key in metaDictQuant:
                        valueList = [float(x) for x in metaDictQuant[key]]
                        meta_rDF = meta_rDF.loc[meta_rDF[key].isin(valueList)]

                wantedList = allFields + ['sampleid', 'sample_name']
                meta_rDF = meta_rDF[wantedList]
                meta_rDF.set_index('sampleid', drop=True, inplace=True)

                base[RID] = 'Step 2 of 8: Selecting your chosen taxa...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 3 of 8: Calculating distance matrix...'

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                r.assign("data", count_rDF)
                r("library(vegan)")

                if distance == 1:
                    r("dist <- vegdist(data, method='manhattan')")
                elif distance == 2:
                    r("dist <- vegdist(data, method='euclidean')")
                elif distance == 3:
                    r("dist <- vegdist(data, method='canberra')")
                elif distance == 4:
                    r("dist <- vegdist(data, method='bray')")
                elif distance == 5:
                    r("dist <- vegdist(data, method='kulczynski')")
                elif distance == 6:
                    r("dist <- vegdist(data, method='jaccard')")
                elif distance == 7:
                    r("dist <- vegdist(data, method='gower')")
                elif distance == 8:
                    r("dist <- vegdist(data, method='altGower')")
                elif distance == 9:
                    r("dist <- vegdist(data, method='morisita')")
                elif distance == 10:
                    r("dist <- vegdist(data, method='horn')")
                elif distance == 11:
                    r("dist <- vegdist(data, method='mountford')")
                elif distance == 12:
                    r("dist <- vegdist(data, method='binomial')")
                elif distance == 13:
                    r("dist <- vegdist(data, method='chao')")
                elif distance == 14:
                    r("dist <- vegdist(data, method='cao')")
                elif distance == 15:
                    datamtx = np.asarray(count_rDF)
                    dists = wOdum(datamtx, alpha)
                    r.assign("dist", dists)
                    r("dist <- as.dist(dist)")

                r("mat <- as.matrix(dist, diag=TRUE, upper=TRUE)")
                mat = r.get("mat")

                rowList = meta_rDF.sample_name.values.tolist()
                distDF = pd.DataFrame(mat, columns=[rowList], index=rowList)

                base[RID] = 'Step 3 of 8: Calculating distance matrix...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 4 of 8: Principal coordinates analysis...'

                anovaFields = []
                for i in metaDictCat:
                    levels = len(set(metaDictCat[i]))
                    if levels > 1:
                        anovaFields.append(i)

                trtLength = len(set(catValues))
                trtString = " + ".join(anovaFields)

                pcoaDF = pd.DataFrame()
                eigDF = pd.DataFrame()
                bigf = ''
                envFit = ''
                if trtLength > 0:
                    r.assign("meta", meta_rDF)
                    pcoa_string = "ord <- capscale(dist ~ " + str(trtString) + ", meta)"
                    r.assign("cmd", pcoa_string)
                    r("eval(parse(text=cmd))")

                    path = "media/temp/pcoa/Rplots/" + str(RID) + ".pcoa.pdf"
                    if os.path.exists(path):
                        os.remove(path)

                    if not os.path.exists('media/temp/pcoa/Rplots'):
                        os.makedirs('media/temp/pcoa/Rplots')

                    file = "pdf('media/temp/pcoa/Rplots/" + str(RID) + ".pcoa.pdf', height=4, width=6)"
                    r.assign("cmd", file)
                    r("eval(parse(text=cmd))")

                    r("par(mfrow=c(1,2))")
                    r("layout(matrix(c(1,2), 1), widths=c(4,2), heights=c(4,4))")
                    r("ordiplot(ord, type='n')")

                    ellipseVal = all['ellipseVal']
                    if ellipseVal != 'None':
                        myStr = "cat <- factor(meta$" + str(ellipseVal) + ")"
                        r.assign("cmd", myStr)
                        r("eval(parse(text=cmd))")
                        r("pl <- ordiellipse(ord, cat, kind='sd', conf=0.95, draw='polygon', border='black')")


                    surfVal = all['surfVal']
                    if surfVal != 'None':
                        myStr = "quant <- meta$" + str(surfVal)
                        r.assign("cmd", myStr)
                        r("eval(parse(text=cmd))")
                        r("ordisurf(ord, quant, cex=1, labcex=0.6, add=TRUE)")

                    colorVal = all['colorVal']
                    if colorVal != 'None':
                        myStr = "cat <- factor(meta$" + str(colorVal) + ")"
                        r.assign("cmd", myStr)
                        r("eval(parse(text=cmd))")
                        r("points(ord, display='sites', pch=15, col=cat)")
                        r("par(mar=c(0,0,4,0))")
                        r("plot(0, type='n', xlab='', ylab='', axes=FALSE)")
                        r("legend('topleft', legend=levels(cat), bty='n', pch=15, col=1:length(cat))")
                    else:
                        r("points(ord, display='sites', pch=15)")
                        r("plot(0, type='n', xlab='', ylab='', axes=FALSE)")

                    r("dev.off()")

                    if len(quantFields) > 0:
                        trtString = " + ".join(quantFields)
                        envfit_str = "fit <- envfit(ord ~ " + str(trtString) + ", meta, choices=c(" + str(PC1) + "," + str(PC2) + "))"
                        r.assign("cmd", envfit_str)
                        r("eval(parse(text=cmd))")

                        fit_out = r("fit")
                        fit_out = fit_out.replace('try({fit})', '')
                        fit_out = fit_out.replace('***VECTORS', '')
                        tempStuff = fit_out.split('\n')
                        for part in tempStuff:
                            if part > tempStuff[1]:
                                envFit += part + '\n'

                    r("res <- summary(ord)")
                    r("id <- rownames(meta)")
                    r("pcoa <- data.frame(id, meta, res$sites)")
                    pcoaDF = r.get("pcoa")

                    pcoaDF.rename(columns={'id': 'Sample ID'}, inplace=True)
                    pcoaDF.rename(columns={'sample_name': 'Sample Name'}, inplace=True)

                    r("Stat <- c('Eigenvalue', 'Proportion Explained', 'Cumulative Proportion')")
                    r("eig <- data.frame(Stat, res$cont$importance)")
                    eigDF = r.get("eig")

                    base[RID] = 'Step 4 of 8: Principal coordinates analysis...done!'

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    base[RID] = 'Step 5 of 8: Performing perMANOVA...'

                    if perms < 10:
                        bigf = 'A minimum of 10 permutations is required...'
                    elif len(catFields_edit) == 0:
                        bigf = 'No categorical variables are available for perMANOVA/betaDisper analysis'
                    elif perms >= 10 and len(catFields_edit) > 0:
                        if test == 1:
                            for i in catFields_edit:
                                factor_string = str(i) + " <- factor(meta$" + str(i) + ")"
                                r.assign("cmd", factor_string)
                                r("eval(parse(text=cmd))")

                                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                                if stops[PID] == RID:
                                    res = ''
                                    return HttpResponse(res, content_type='application/json')
                                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                            r.assign("perms", perms)
                            trtString = " * ".join(catFields_edit)
                            amova_string = "res <- adonis(dist ~ " + str(trtString) + ", perms=perms)"
                            r.assign("cmd", amova_string)
                            r("eval(parse(text=cmd))")

                            res_aov = r("res$aov.tab")

                            tempStuff = res_aov.split('\n')
                            for part in tempStuff:
                                if part != tempStuff[0]:
                                    bigf += part + '\n'
                            base[RID] = 'Step 5 of 8: Performing perMANOVA...done!'

                        elif test == 2:
                            base[RID] = 'Step 4 of 8: Principal coordinates analysis...done!'
                            base[RID] = 'Step 5 of 8: Performing BetaDisper...'

                            for i in catFields_edit:
                                factor_string = str(i) + " <- factor(meta$" + str(i) + ")"
                                r.assign("cmd", factor_string)
                                r("eval(parse(text=cmd))")

                                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                                if stops[PID] == RID:
                                    res = ''
                                    return HttpResponse(res, content_type='application/json')
                                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                            r.assign("perms", perms)
                            for i in catFields_edit:
                                beta_string = "res <- betadisper(dist, " + str(i) + ")"
                                r.assign("cmd", beta_string)
                                r("eval(parse(text=cmd))")

                                r("something <- anova(res)")
                                beta = r("something")
                                tempStuff = beta.split('\n')
                                bigf += 'group: ' + str(i) + '\n'
                                for part in tempStuff:
                                    if part != tempStuff[0]:
                                        bigf += part + '\n'
                                base[RID] = 'Step 5 of 8: Performing BetaDisper...done!'

                                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                                if stops[PID] == RID:
                                    res = ''
                                    return HttpResponse(res, content_type='application/json')
                                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                else:
                    state = "Your selected variable(s) only have one treatment level, please select additional data!"
                    myDict = {}
                    myDict['error'] = state
                    res = simplejson.dumps(myDict)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 6 of 8: Formatting graph data for display...'
                finalDict = {}
                seriesList = []
                xAxisDict = {}
                yAxisDict = {}

                CAP1 = PC1 + len(catFields_edit) + len(quantFields) + 1
                CAP2 = PC2 + len(catFields_edit) + len(quantFields) + 1

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

                colors_idx = 0
                if catFields_edit:
                    grouped = pcoaDF.groupby(catFields_edit)
                    for name, group in grouped:
                        if len(catFields_edit) > 1:
                            trt = "; ".join(name)
                        else:
                            trt = name

                        dataList = []
                        for index, row in group.iterrows():
                            dataDict = {}
                            dataDict['name'] = row['Sample Name']
                            dataDict['x'] = float(row[CAP1])
                            dataDict['y'] = float(row[CAP2])
                            dataList.append(dataDict)

                        seriesDict = {}
                        seriesDict['name'] = str(trt)
                        seriesDict['data'] = dataList
                        seriesDict['color'] = colors[colors_idx]
                        seriesList.append(seriesDict)
                        colors_idx += 1

                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                        if stops[PID] == RID:
                            res = ''
                            return HttpResponse(res, content_type='application/json')
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                xTitle = {}
                xTitle['text'] = str(eigDF.columns.values.tolist()[PC1]) + " (" + str(eigDF.iloc[1][PC1] * 100) + "%)"
                xTitle['style'] = {'color': 'black', 'fontSize': '18px', 'fontWeight': 'bold'}
                xAxisDict['title'] = xTitle

                yTitle = {}
                yTitle['text'] = str(eigDF.columns.values.tolist()[PC2]) + " (" + str(eigDF.iloc[1][PC2] * 100) + "%)"
                yTitle['style'] = {'color': 'black', 'fontSize': '18px', 'fontWeight': 'bold'}
                yAxisDict['title'] = yTitle

                styleDict = {'style': {'color': 'black', 'fontSize': '14px'}}
                xAxisDict['labels'] = styleDict
                yAxisDict['labels'] = styleDict

                finalDict['series'] = seriesList
                finalDict['xAxis'] = xAxisDict
                finalDict['yAxis'] = yAxisDict

                if test == 1:
                    result += 'perMANOVA results:' + '\n'
                if test == 2:
                    result += 'betaDisper results:' + '\n'

                if len(catFields_edit) == 0:
                    result += 'test cannot be run...' + '\n'
                else:
                    bigf = bigf.decode('utf-8')
                    result += bigf + '\n'
                result += '===============================================\n'

                if len(catFields_edit) > 0:
                    result += '\nenvfit results:\n'
                    envFit = envFit.decode('utf-8')
                    result += envFit
                result += '===============================================\n'

                result += '\nEigenvalues\n'
                eigStr = eigDF.to_string()
                result += str(eigStr) + '\n'
                result += '===============================================\n\n\n\n'

                finalDict['text'] = result

                base[RID] = 'Step 6 of 8: Formatting graph data for display...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 7 of 8: Formatting PCoA table...'

                res_table = pcoaDF.to_html(classes="table display")
                res_table = res_table.replace('border="1"', 'border="0"')
                finalDict['res_table'] = str(res_table)

                base[RID] = 'Step 7 of 8: Formatting PCoA table...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 8 of 8: Formatting distance score table...'

                distDF.sort_index(axis=1, inplace=True)
                distDF.sort_index(axis=0, inplace=True)
                dist_table = distDF.to_html(classes="table display")
                dist_table = dist_table.replace('border="1"', 'border="0"')
                finalDict['dist_table'] = str(dist_table)

                base[RID] = 'Step 8 of 8: Formatting distance score table...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                removeRIDPCoA(RID)
                return HttpResponse(res, content_type='application/json')

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "Error with PCoA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = simplejson.dumps(myDict)
            removeRIDPCoA(RID)
            return HttpResponse(res, content_type='application/json')


def getTaxaDF(selectAll, savedDF, metaDF, allFields, DepVar, RID, stops, PID):
    global base
    try:
        base[RID] = 'Step 2 of 4: Selecting your chosen taxa or KEGG level...'
        taxaDF = pd.DataFrame(columns=['sampleid', 'rank', 'rank_id', 'rank_name', 'rel_abund', 'abund_16S', 'rich', 'diversity'])

        if selectAll == 2:
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

        taxaDF.drop('sampleid', axis=1, inplace=True)
        finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')

        wantedList = allFields + ['sampleid', 'rank', 'rank_name', 'rank_id']
        if DepVar == 1:
            finalDF = finalDF.groupby(wantedList)[['rel_abund']].sum()
        elif DepVar == 2:
            finalDF = finalDF.groupby(wantedList)[['rich']].sum()
        elif DepVar == 3:
            finalDF = finalDF.groupby(wantedList)[['diversity']].sum()
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
            myDict['error'] = "Error with PCoA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getKeggDF(keggAll, savedDF, tempDF, allFields, DepVar, RID, stops, PID):
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
        wanted = ['sampleid', 'speciesid', 'rel_abund', 'abund_16S', 'rich', 'diversity']
        profileDF = tempDF.loc[:, wanted]
        profileDF.set_index('speciesid', inplace=True)

        # get PICRUSt data for species
        speciesList = pd.unique(profileDF.index.ravel().tolist())
        qs = PICRUSt.objects.using('picrust').filter(speciesid__in=speciesList)
        picrustDF = read_frame(qs, fieldnames=['speciesid__speciesid', 'geneCount'])
        picrustDF.set_index('speciesid__speciesid', inplace=True)

        path = 'media/temp/pcoa/' + str(RID)
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
            path = 'media/temp/pcoa/'+str(RID)+'/file%d.temp' % i
            frame = pd.read_csv(path)
            picrustDF = picrustDF.append(frame, ignore_index=True)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        shutil.rmtree('media/temp/pcoa/'+str(RID))
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
        elif DepVar == 4:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund_16S')
        elif DepVar == 2:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rich')
        elif DepVar == 3:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='diversity')

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
            myDict = {}
            myDict['error'] = "Error with PCoA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getNZDF(nzAll, savedDF, tempDF, allFields, DepVar, RID, stops, PID):
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
        wanted = ['sampleid', 'speciesid', 'rel_abund', 'abund_16S', 'rich', 'diversity']
        profileDF = tempDF.loc[:, wanted]
        profileDF.set_index('speciesid', inplace=True)

        # get PICRUSt data for species
        speciesList = pd.unique(profileDF.index.ravel().tolist())
        qs = PICRUSt.objects.using('picrust').filter(speciesid__in=speciesList)
        picrustDF = read_frame(qs, fieldnames=['speciesid__speciesid', 'geneCount'])
        picrustDF.set_index('speciesid__speciesid', inplace=True)

        path = 'media/temp/pcoa/' + str(RID)
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
            path = 'media/temp/pcoa/'+str(RID)+'/file%d.temp' % i
            frame = pd.read_csv(path)
            picrustDF = picrustDF.append(frame, ignore_index=True)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        shutil.rmtree('media/temp/pcoa/'+str(RID))
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
        elif DepVar == 4:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund_16S')
        elif DepVar == 2:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rich')
        elif DepVar == 3:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='diversity')

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
            myDict = {}
            myDict['error'] = "Error with PCoA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def sumStuff(slice, koDict, RID, num, stops, PID):
    global base
    db.close_old_connections()

    f = open('media/temp/pcoa/'+str(RID)+'/file'+str(num)+".temp", 'w')

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


def removePCoAFiles(request):
    if request.is_ajax():
        RID = request.GET["all"]

        file = "media/temp/pcoa/Rplots/" + str(RID) + ".pcoa.pdf"
        if os.path.exists(file):
            os.remove(file)

        file = "media/temp/pcoa/" + str(RID) + ".pkl"
        if os.path.exists(file):
            os.remove(file)

        file = "media/temp/pcoa/" + str(RID) + ".csv"
        if os.path.exists(file):
            os.remove(file)

        return HttpResponse()


def getTabPCoA(request):
    if request.is_ajax():
        RID = request.GET["all"]
        myDir = 'media/temp/pcoa/'
        fileName = str(myDir) + str(RID) + '.pkl'
        savedDF = pd.read_pickle(fileName)

        myDir = 'media/temp/pcoa/'
        fileName = str(myDir) + str(RID) + '.csv'
        savedDF.to_csv(fileName)

        myDict = {}
        myDir = 'temp/pcoa/'
        fileName = str(myDir) + str(RID) + '.csv'
        myDict['name'] = str(fileName)
        res = simplejson.dumps(myDict)

        return HttpResponse(res, content_type='application/json')
