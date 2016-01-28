import datetime
from django.http import HttpResponse
import logging
import numpy as np
import pandas as pd
import pickle
from pyper import *
from scipy.spatial.distance import *
import simplejson

from database.utils import multidict, stoppableThread, wOdum


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}
stop4 = False
stops = {}
thread4 = stoppableThread()
res = ''
LOG_FILENAME = 'error_log.txt'


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
            stage[RID] = str(base[RID]) + '<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
        except:
            stage[RID] = '<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
        myDict['stage'] = stage[RID]
        json_data = simplejson.dumps(myDict, encoding="Latin-1")
        return HttpResponse(json_data, content_type='application/json')


def removeRIDPCoA(request):
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


def stopPCoA(request):
    global thread4, stops, stop4, res
    if request.is_ajax():
        RID = request.GET["all"]
        stops[RID] = True
        stop4 = True
        thread4.terminate()
        thread4.join()
        removeRIDPCoA(request)

        res = ''
        myDict = {}
        myDict['error'] = 'none'
        myDict['message'] = 'Your analysis has been stopped!'
        stop = simplejson.dumps(myDict)
        return HttpResponse(stop, content_type='application/json')


def getPCoA(request):
    global res, thread4, stop4
    if request.is_ajax():
        stop4 = False
        thread4 = stoppableThread(target=loopCat, args=(request,))
        thread4.start()
        thread4.join()
        removeRIDPCoA(request)
        return HttpResponse(res, content_type='application/json')


def loopCat(request):
    global res, base, stage, time1, TimeDiff, stops, stop4
    try:
        while True:
            if request.is_ajax():
                allJson = request.GET["all"]
                all = simplejson.loads(allJson)

                RID = str(all["RID"])
                stops[RID] = False

                time1[RID] = time.time()
                base[RID] = 'Step 1 of 8: Selecting your chosen meta-variables...'

                path = pickle.loads(request.session['savedDF'])
                savedDF = pd.read_pickle(path)

                DepVar = int(all["DepVar"])
                selectAll = int(all["taxa"])
                distance = int(all["distance"])
                PC1 = int(all["PC1"])
                PC2 = int(all["PC2"])
                test = int(all["test"])
                alpha = float(all["alpha"])
                perms = int(all["perms"])

                result = ''
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
                    myDict['error'] = "Selected meta data only has one level.\nPlease select a different variable."
                    res = simplejson.dumps(myDict)
                    return None

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
                metaDF = savedDF.loc[savedDF['sampleid'].isin(allSampleIDs)]

                if metaDictCat:
                    for key in metaDictCat:
                        metaDF = metaDF.loc[metaDF[key].isin(metaDictCat[key])]

                if metaDictQuant:
                    for key in metaDictQuant:
                        valueList = [float(x) for x in metaDictQuant[key]]
                        metaDF = metaDF.loc[metaDF[key].isin(valueList)]

                wantedList = allFields + ['sample_name']
                metaDF = metaDF[wantedList]

                result += 'Categorical variables selected by user: ' + ", ".join(catFields) + '\n'
                result += 'Categorical variables removed from analysis (contains only 1 level): ' + ", ".join(removed) + '\n'
                result += 'Quantitative variables selected by user: ' + ", ".join(quantFields) + '\n'
                result += '===============================================\n'

                base[RID] = 'Step 1 of 4: Selecting your chosen meta-variables...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 2 of 4: Selecting your chosen taxa...'

                # get selected taxa fro each rank selected in the tree
                taxaDF = pd.DataFrame(columns=['sampleid', 'rank', 'taxa_id', 'taxa_name', 'abund', 'abund_16S', 'rich', 'diversity'])

                if selectAll == 1:
                    taxaDF = savedDF.loc[:, ['sampleid', 'kingdomid', 'kingdomName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'kingdomid': 'taxa_id', 'kingdomName': 'taxa_name'}, inplace=True)
                    taxaDF.loc[:, 'rank'] = 'Kingdom'
                elif selectAll == 2:
                    taxaDF = savedDF.loc[:, ['sampleid', 'phylaid', 'phylaName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'phylaid': 'taxa_id', 'phylaName': 'taxa_name'}, inplace=True)
                    taxaDF.loc[:, 'rank'] = 'Phyla'
                elif selectAll == 3:
                    taxaDF = savedDF.loc[:, ['sampleid', 'classid', 'className', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'classid': 'taxa_id', 'className': 'taxa_name'}, inplace=True)
                    taxaDF.loc[:, 'rank'] = 'Class'
                elif selectAll == 4:
                    taxaDF = savedDF.loc[:, ['sampleid', 'orderid', 'orderName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'orderid': 'taxa_id', 'orderName': 'taxa_name'}, inplace=True)
                    taxaDF.loc[:, 'rank'] = 'Order'
                elif selectAll == 5:
                    taxaDF = savedDF.loc[:, ['sampleid', 'familyid', 'familyName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'familyid': 'taxa_id', 'familyName': 'taxa_name'}, inplace=True)
                    taxaDF.loc[:, 'rank'] = 'Family'
                elif selectAll == 6:
                    taxaDF = savedDF.loc[:, ['sampleid', 'genusid', 'genusName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'genusid': 'taxa_id', 'genusName': 'taxa_name'}, inplace=True)
                    taxaDF.loc[:, 'rank'] = 'Genus'
                elif selectAll == 7:
                    taxaDF = savedDF.loc[:, ['sampleid', 'speciesid', 'speciesName', 'abund', 'abund_16S', 'rich', 'diversity']]
                    taxaDF.rename(columns={'speciesid': 'taxa_id', 'speciesName': 'taxa_name'}, inplace=True)
                    taxaDF.loc[:, 'rank'] = 'Species'

                finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')

                wantedList = allFields + ['sampleid', 'sample_name', 'rank', 'taxa_name', 'taxa_id']
                finalDF = finalDF.groupby(wantedList)[['abund', 'abund_16S', 'rich', 'diversity']].sum()

                finalDF.reset_index(drop=False, inplace=True)

                taxaSums = finalDF.groupby('taxa_id').sum()
                goodList = taxaSums[taxaSums['abund'] > 0].index.values.tolist()
                finalDF = finalDF.loc[finalDF['taxa_id'].isin(goodList)]

                count_rDF = pd.DataFrame()
                if DepVar == 1:
                    count_rDF = finalDF.pivot(index='sampleid', columns='taxa_id', values='abund')
                elif DepVar == 2:
                    count_rDF = finalDF.pivot(index='sampleid', columns='taxa_id', values='abund_16S')
                elif DepVar == 3:
                    count_rDF = finalDF.pivot(index='sampleid', columns='taxa_id', values='rich')
                elif DepVar == 4:
                    count_rDF = finalDF.pivot(index='sampleid', columns='taxa_id', values='diversity')

                meta_rDF = finalDF.drop_duplicates(subset='sampleid', keep='last')
                wantedList = allFields + ['sampleid', 'sample_name']
                meta_rDF = meta_rDF[wantedList]
                meta_rDF.set_index('sampleid', drop=True, inplace=True)

                base[RID] = 'Step 2 of 5: Selecting your chosen taxa...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
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

                rowList = meta_rDF['sample_name'].tolist()

                distDF = pd.DataFrame(mat, columns=[rowList], index=rowList)

                base[RID] = 'Step 3 of 8: Calculating distance matrix...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
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
                    pcoa_string = "ord <- capscale(mat ~ " + str(trtString) + ", meta)"
                    r.assign("cmd", pcoa_string)
                    r("eval(parse(text=cmd))")

                    ellipseVal = all['ellipseVal']
                    myStr = "cat <- factor(meta$" + str(ellipseVal) + ")"
                    r.assign("cmd", myStr)
                    r("eval(parse(text=cmd))")

                    surfVal = all['surfVal']
                    if surfVal:
                        myStr = "quant <- meta$" + str(surfVal)
                        r.assign("cmd", myStr)
                        r("eval(parse(text=cmd))")

                    name = request.user
                    ip = request.META.get('REMOTE_ADDR')
                    user = str(name) + "." + str(ip)

                    path = "media/Rplots/" + str(user) + ".pcoa.jpg"
                    if os.path.exists(path):
                        os.remove(path)

                    if not os.path.exists('media/Rplots'):
                        os.makedirs('media/Rplots')

                    file = "jpeg('media/Rplots/" + str(user) + ".pcoa.jpg', height=400, width=400)"
                    r.assign("cmd", file)
                    r("eval(parse(text=cmd))")

                    r("ordiplot(ord, type='n')")
                    r("points(ord, display='sites', pch=15, col=cat)")
                    r("legend('topright', legend=levels(cat), pch=15, col=1:length(cat))")

                    addEllipse = all['addEllipse']
                    if addEllipse == 'yes':
                        r("pl <- ordiellipse(ord, cat, kind='sd', conf=0.95, draw='polygon', border='black')")

                    addSurf = all['addSurf']
                    if addSurf == 'yes':
                        r("ordisurf(ord, quant, add=TRUE)")
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
                    if stops[RID]:
                        res = ''
                        return None
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
                                if stops[RID]:
                                    res = ''
                                    return None
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
                                if stops[RID]:
                                    res = ''
                                    return None
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
                                if stops[RID]:
                                    res = ''
                                    return None
                                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                else:
                    state = "Your selected variable(s) only have one treatment level, please select additional data!"
                    myDict = {}
                    myDict['error'] = state
                    res = simplejson.dumps(myDict)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
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
                        dataList = group[[CAP1, CAP2]].values.astype(float).tolist()
                        if len(catFields_edit) > 1:
                            trt = "; ".join(name)
                        else:
                            trt = name
                        seriesDict = {}
                        seriesDict['name'] = str(trt)
                        seriesDict['data'] = dataList
                        seriesDict['color'] = colors[colors_idx]
                        seriesList.append(seriesDict)
                        colors_idx += 1

                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                        if stops[RID]:
                            res = ''
                            return None
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                xTitle = {}
                xTitle['text'] = str(eigDF.columns.values.tolist()[PC1]) + " (" + str(eigDF.iloc[1][PC1] * 100) + "%)"
                xAxisDict['title'] = xTitle

                yTitle = {}
                yTitle['text'] = str(eigDF.columns.values.tolist()[PC2]) + " (" + str(eigDF.iloc[1][PC2] * 100) + "%)"
                yAxisDict['title'] = yTitle

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
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 7 of 8: Formatting PCoA table...'

                res_table = pcoaDF.to_html(classes="table display")
                res_table = res_table.replace('border="1"', 'border="0"')
                finalDict['res_table'] = str(res_table)

                base[RID] = 'Step 7 of 8: Formatting PCoA table...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 8 of 8: Formatting distance score table...'

                distDF.sort_index(axis=1, inplace=True)
                distDF.sort_index(axis=0, inplace=True)
                dist_table = distDF.to_html(classes="table display")
                dist_table = dist_table.replace('border="1"', 'border="0"')
                finalDict['dist_table'] = str(dist_table)

                base[RID] = 'Step 8 of 8: Formatting distance score table...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                return None

    except:
        if not stop4:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with PCoA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
        return None


def removegraphPCoA(request):
    name = request.user
    ip = request.META.get('REMOTE_ADDR')
    user = str(name) + "." + str(ip)

    file = "media/Rplots/" + str(user) + ".pcoa.jpg"
    if os.path.exists(file):
        os.remove(file)
    return HttpResponse()
