import datetime
from django.http import HttpResponse
import logging
import numpy as np
import pandas as pd
import pickle
from pyper import *
from scipy import stats
from scipy.spatial.distance import *
import simplejson

from database.models import Kingdom, Phyla, Class, Order, Family, Genus, Species
from database.utils import multidict, stoppableThread


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}
stops = {}
stop5 = False
thread5 = stoppableThread()
res = ''
LOG_FILENAME = 'error_log.txt'


def statusSPLS(request):
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


def removeRIDSPLS(request):
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


def stopSPLS(request):
    global thread5, stops, stop5, res
    if request.is_ajax():
        RID = request.GET["all"]
        stops[RID] = True
        stop5 = True
        thread5.terminate()
        thread5.join()
        removeRIDSPLS(request)

        res = ''
        myDict = {}
        myDict['error'] = 'none'
        myDict['message'] = 'Your analysis has been stopped!'
        stop = simplejson.dumps(myDict)
        return HttpResponse(stop, content_type='application/json')


def getSPLS(request):
    global res, thread5, stop5
    if request.is_ajax():
        stop5 = False
        thread5 = stoppableThread(target=loopCat, args=(request,))
        thread5.start()
        thread5.join()
        removeRIDSPLS(request)
        return HttpResponse(res, content_type='application/json')


def loopCat(request):
    global res, base, stage, time1, TimeDiff, stops, stop5
    try:
        while True:
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

                # Select samples and meta-variables from savedDF
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

                # Removes samples (rows) that are not in our samplelist
                metaDF = savedDF.loc[savedDF['sampleid'].isin(quantSampleIDs)]

                if metaDictQuant:
                    for key in metaDictQuant:
                        valueList = [float(x) for x in metaDictQuant[key]]
                        metaDF = metaDF.loc[metaDF[key].isin(valueList)]

                wantedList = quantFields + ['sample_name']
                metaDF = metaDF[wantedList]

                result = ''
                if DepVar == 1:
                    result += 'Dependent variable: Abundance\n'
                if DepVar == 2:
                    result += 'Dependent variable: Species Richness\n'
                if DepVar == 3:
                    result += 'Dependent variable: Species Diversity\n'
                if DepVar == 4:
                    result += 'Dependent variable: Abundance (rRNA gene copies)\n'

                result += 'Quantitative variables selected: ' + ", ".join(quantFields) + '\n'
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

                wantedList = quantFields + ['sampleid', 'sample_name', 'rank', 'taxa_name', 'taxa_id']
                finalDF = finalDF.groupby(wantedList)[['abund', 'abund_16S', 'rich', 'diversity']].sum()

                finalDF.reset_index(drop=False, inplace=True)

                taxaSums = finalDF.groupby('taxa_id').sum()
                goodList = taxaSums[taxaSums['abund'] > 0].index.values.tolist()
                finalDF = finalDF.loc[finalDF['taxa_id'].isin(goodList)]

                count_rDF = pd.DataFrame()
                if DepVar == 1:
                    count_rDF = finalDF.pivot(index='sampleid', columns='taxa_id', values='abund')
                elif DepVar == 2:
                    count_rDF = finalDF.pivot(index='sampleid', columns='taxa_id', values='rich')
                elif DepVar == 3:
                    count_rDF = finalDF.pivot(index='sampleid', columns='taxa_id', values='diversity')
                elif DepVar == 4:
                    count_rDF = finalDF.pivot(index='sampleid', columns='taxa_id', values='abund_16S')

                meta_rDF = finalDF.drop_duplicates(subset='sampleid', keep='last')

                wantedList = quantFields + ['sampleid']
                meta_rDF = meta_rDF[wantedList]
                meta_rDF.set_index('sampleid', drop=True, inplace=True)

                base[RID] = 'Step 2 of 5: Selecting your chosen taxa...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 3 of 5: Calculating sPLS...'

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                r.assign("X", count_rDF)
                r.assign("Y", meta_rDF)
                r.assign("names", count_rDF.columns.values)
                r("colnames(X) <- names")

                r("library(mixOmics)")
                r("ZeroVar <- nearZeroVar(X, freqCut=90/10, uniqueCut=25)")
                r("List <- row.names(ZeroVar$Metrics)")
                r("X_new <- X[,-which(names(X) %in% List)]")
                r("X_scaled <- scale(X_new, center=TRUE, scale=TRUE)")
                r("Y_scaled <- scale(Y, center=TRUE, scale=TRUE)")
                r("detach('package:mixOmics', unload=TRUE)")
                r("library(spls)")
                r("set.seed(1)")

                spls_string = "cv <- cv.spls(X_scaled, Y_scaled, scale.x=FALSE, scale.y=FALSE, eta=seq(0.1, 0.9, 0.1), K=c(1:" + str(len(quantFields)-1) + "), plot.it=FALSE)"
                r.assign("cmd", spls_string)
                r("eval(parse(text=cmd))")

                r("f <- spls(X_scaled, Y_scaled, scale.x=FALSE, scale.y=FALSE, eta=cv$eta.opt, K=cv$K.opt)")

                r("out <- capture.output(print(f))")
                fout = r.get("out")

                if fout is not None:
                    for i in fout:
                        result += str(i) + '\n'
                else:
                    result += 'No significant variables were found\n'
                result += '===============================================\n\n'

                r("set.seed(1)")
                r("ci.f <- ci.spls(f, plot.it=FALSE, plot.fix='y')")
                r("cis <- ci.f$cibeta")
                r("cf <- correct.spls(ci.f, plot.it=FALSE)")
                r("out <- capture.output(cis)")
                fout = r.get("out")

                if fout is not None:
                    result += '\n\nBootstrapped confidence intervals of coefficients:\n'
                    for i in fout:
                        result += str(i) + '\n'
                    result += '===============================================\n\n'

                r("coef.f <- coef(f)")
                r("sum <- sum(coef.f != 0)")
                total = r.get("sum")

                base[RID] = 'Step 3 of 5: Calculating sPLS...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 4 of 5: Formatting graph data for display...'

                finalDict = {}
                if total is not None:
                    r("pred.f <- predict(f, type='fit')")
                    r("library(DMwR)")
                    r("pred.ns <- unscale(pred.f, Y_scaled)")
                    r("pred.ns.rows <- row.names(pred.ns)")
                    pred = r.get("pred.ns")
                    rows = r.get("pred.ns.rows")

                    predList = ['pred_' + s for s in quantFields]
                    predDF = pd.DataFrame(pred,  columns=[predList], index=rows)

                    resultDF = pd.merge(meta_rDF, predDF, left_index=True, right_index=True)
                    result += 'sPLS Model Fit (y = mx + b):\n'
                    result += 'y = predicted\n'
                    result += 'x = observed\n\n'

                    for i in xrange(len(quantFields)):
                        x = resultDF[quantFields[i]].astype(float).values.tolist()
                        y = resultDF[predList[i]].astype(float).values.tolist()

                        slp, inter, r_value, p, se = stats.linregress(x, y)
                        r_sq = r_value * r_value

                        result += 'Variable: ' + str(quantFields[i]) + '\n'
                        result += 'Slope (m): ' + str(slp) + '\n'
                        result += 'Intercept (b): ' + str(inter) + '\n'
                        result += 'R2: ' + str(r_sq) + '\n'
                        result += 'Std Error: ' + str(se) + '\n\n\n'

                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                        if stops[RID]:
                            res = ''
                            return None
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    r("coef.f.rows <- row.names(coef.f)")
                    cf = r.get("coef.f")
                    rows = r.get("coef.f.rows")

                    coeffsDF = pd.DataFrame(cf,  columns=[quantFields], index=rows)
                    coeffsDF = coeffsDF.loc[(coeffsDF != 0).any(axis=1)]
                    coeffsDF.sort_index(inplace=True)
                    taxIDList = coeffsDF.index.values.tolist()

                    namesDF = pd.DataFrame()
                    if selectAll == 1:
                        taxNameList = Kingdom.objects.filter(kingdomid__in=taxIDList).values('kingdomid', 'kingdomName')
                        namesDF = pd.DataFrame(list(taxNameList))
                        namesDF.set_index('kingdomid', inplace=True)
                        namesDF.rename(columns={'kingdomName': 'taxa_name'}, inplace=True)
                    elif selectAll == 2:
                        taxNameList = Phyla.objects.filter(phylaid__in=taxIDList).values('phylaid', 'phylaName')
                        namesDF = pd.DataFrame(list(taxNameList))
                        namesDF.set_index('phylaid', inplace=True)
                        namesDF.rename(columns={'phylaName': 'taxa_name'}, inplace=True)
                    elif selectAll == 3:
                        taxNameList = Class.objects.filter(classid__in=taxIDList).values('classid', 'className')
                        namesDF = pd.DataFrame(list(taxNameList))
                        namesDF.set_index('classid', inplace=True)
                        namesDF.rename(columns={'className': 'taxa_name'}, inplace=True)
                    elif selectAll == 4:
                        taxNameList = Order.objects.filter(orderid__in=taxIDList).values('orderid', 'orderName')
                        namesDF = pd.DataFrame(list(taxNameList))
                        namesDF.set_index('orderid', inplace=True)
                        namesDF.rename(columns={'orderName': 'taxa_name'}, inplace=True)
                    elif selectAll == 5:
                        taxNameList = Family.objects.filter(familyid__in=taxIDList).values('familyid', 'familyName')
                        namesDF = pd.DataFrame(list(taxNameList))
                        namesDF.set_index('familyid', inplace=True)
                        namesDF.rename(columns={'familyName': 'taxa_name'}, inplace=True)
                    elif selectAll == 6:
                        taxNameList = Genus.objects.filter(genusid__in=taxIDList).values('genusid', 'genusName')
                        namesDF = pd.DataFrame(list(taxNameList))
                        namesDF.set_index('genusid', inplace=True)
                        namesDF.rename(columns={'genusName': 'taxa_name'}, inplace=True)
                    elif selectAll == 7:
                        taxNameList = Species.objects.filter(speciesid__in=taxIDList).values('speciesid', 'speciesName')
                        namesDF = pd.DataFrame(list(taxNameList))
                        namesDF.set_index('speciesid', inplace=True)
                        namesDF.rename(columns={'speciesName': 'taxa_name'}, inplace=True)

                    namesDF.sort_index(inplace=True)
                    taxaNameList = namesDF['taxa_name'].values.tolist()

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[RID]:
                        res = ''
                        return None
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    coeffsDF = pd.merge(namesDF, coeffsDF, left_index=True, right_index=True, how='inner')
                    coeffsDF.reset_index(inplace=True)
                    res_table = coeffsDF.to_html(classes="table display")
                    res_table = res_table.replace('border="1"', 'border="0"')
                    finalDict['res_table'] = str(res_table)

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[RID]:
                        res = ''
                        return None
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    if selectAll == 2:
                        clustDF = coeffsDF.drop('phylaid', axis=1)
                        clustDF.set_index('taxa_name', inplace=True)
                    elif selectAll == 3:
                        clustDF = coeffsDF.drop('classid', axis=1)
                        clustDF.set_index('taxa_name', inplace=True)
                    elif selectAll == 4:
                        clustDF = coeffsDF.drop('orderid', axis=1)
                        clustDF.set_index('taxa_name', inplace=True)
                    elif selectAll == 5:
                        clustDF = coeffsDF.drop('familyid', axis=1)
                        clustDF.set_index('taxa_name', inplace=True)
                    elif selectAll == 6:
                        clustDF = coeffsDF.drop('genusid', axis=1)
                        clustDF.set_index('taxa_name', inplace=True)
                    elif selectAll == 7:
                        clustDF = coeffsDF.drop('speciesid', axis=1)
                        clustDF.set_index('taxa_name', inplace=True)

                    resultDF.reset_index(inplace=True)
                    resultDF.rename(columns={'index': 'sampleid'}, inplace=True)
                    pred_table = resultDF.to_html(classes="table display")
                    pred_table = pred_table.replace('border="1"', 'border="0"')
                    finalDict['pred_table'] = str(pred_table)

                    xAxisDict = {}
                    xAxisDict['categories'] = taxaNameList
                    labelsDict = {}
                    labelsDict['rotation'] = 270
                    labelsDict['enabled'] = True
                    xAxisDict['labels'] = labelsDict
                    xAxisDict['title'] = {'text': None}
                    xAxisDict['tickLength'] = 0

                    yAxisDict = {}
                    yAxisDict['categories'] = quantFields
                    yAxisDict['title'] = {'text': None}

                    seriesList = []
                    seriesDict = {}
                    seriesDict['borderWidth'] = '1'

                    row, col = coeffsDF.shape
                    dataList = []
                    for i in xrange(row):
                        for j in xrange(len(quantFields)):
                            val = round(coeffsDF[quantFields[j]].iloc[i], 5)
                            tup = (i, j, val)
                            obsList = list(tup)
                            dataList.append(obsList)

                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                            if stops[RID]:
                                res = ''
                                return None
                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[RID]:
                        res = ''
                        return None
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    seriesDict['data'] = dataList
                    labelDict = {}
                    labelDict['enabled'] = True
                    labelDict['color'] = 'black',
                    labelDict['syle'] = {'textShadow': 'none'}
                    seriesList.append(seriesDict)

                    finalDict['xAxis'] = xAxisDict
                    finalDict['yAxis'] = yAxisDict
                    finalDict['series'] = seriesList

                    # R clustered heatmap
                    clustDF = pd.DataFrame()
                    if selectAll == 2:
                        clustDF = coeffsDF.drop('phylaid', axis=1)
                    elif selectAll == 3:
                        clustDF = coeffsDF.drop('classid', axis=1)
                    elif selectAll == 4:
                        clustDF = coeffsDF.drop('orderid', axis=1)
                    elif selectAll == 5:
                        clustDF = coeffsDF.drop('familyid', axis=1)
                    elif selectAll == 6:
                        clustDF = coeffsDF.drop('genusid', axis=1)
                    elif selectAll == 7:
                        clustDF = coeffsDF.drop('speciesid', axis=1)
                    row, col = clustDF.shape

                    method = all['methodVal']
                    metric = all['metricVal']

                    name = request.user
                    ip = request.META.get('REMOTE_ADDR')
                    user = str(name) + "." + str(ip)

                    path = "media/Rplots/" + str(user) + ".spls.jpg"
                    if os.path.exists(path):
                        os.remove(path)

                    if not os.path.exists('media/Rplots'):
                        os.makedirs('media/Rplots')

                    height = 250 + 15*row
                    width = 250 + 20*(col-1)
                    file = "jpeg('media/Rplots/" + str(user) + ".spls.jpg', height=" + str(height) + ", width=" + str(width) + ")"
                    r.assign("cmd", file)
                    r("eval(parse(text=cmd))")

                    r.assign("df", clustDF[quantFields])
                    r("df <- as.matrix(df)")
                    r.assign("rows", clustDF.taxa_name.values)
                    r("rownames(df) <- rows")
                    r("library(pheatmap)")
                    r("library(RColorBrewer)")
                    r("col.pal <- brewer.pal(9,'RdBu')")

                    if row > 2 and col > 3:
                        hmap_str = "pheatmap(df, fontsize=12, color=col.pal, clustering_method='" + str(method) + "', clustering_distance_rows='" + str(metric) + "', clustering_distance_cols='" + str(metric) + "')"
                        r.assign("cmd", hmap_str)
                        r("eval(parse(text=cmd))")
                        r("dev.off()")

                    if row > 2 and col <= 3:
                        hmap_str = "pheatmap(df, color=col.pal, cluster_col=FALSE, clustering_method='" + str(method) + "', clustering_distance_rows='" + str(metric) + "')"
                        r.assign("cmd", hmap_str)
                        r("eval(parse(text=cmd))")
                        r("dev.off()")

                    if row <= 2 and col > 3:
                        hmap_str = "pheatmap(df, color=col.pal, cluster_row=FALSE, clustering_method='" + str(method) + "', clustering_distance_cols='" + str(metric) + "')"
                        r.assign("cmd", hmap_str)
                        r("eval(parse(text=cmd))")
                        r("dev.off()")

                    if row <= 2 and col <= 3:
                        hmap_str = "pheatmap(df, color=col.pal, cluster_col=FALSE, cluster_row=FALSE)"
                        r.assign("cmd", hmap_str)
                        r("eval(parse(text=cmd))")
                        r("dev.off()")

                finalDict['text'] = result

                base[RID] = 'Step 4 of 5: Formatting graph data for display...done!'
                base[RID] = 'Step 5 of 5: Formatting sPLS coefficient table...'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                return None

    except:
        if not stop5:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with Differential Abundance!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
        return None


def removegraphSPLS(request):
    name = request.user
    ip = request.META.get('REMOTE_ADDR')
    user = str(name) + "." + str(ip)

    file = "media/Rplots/" + str(user) + ".spls.jpg"
    if os.path.exists(file):
        os.remove(file)
    return HttpResponse()
