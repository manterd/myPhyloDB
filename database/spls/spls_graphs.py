import ast
import datetime
from django import db
from django.http import HttpResponse
from django_pandas.io import read_frame
import logging
import math
import numpy as np
import pandas as pd
from pyper import *
from scipy import stats
import shutil
import simplejson
import threading

from database.models import Kingdom, Phyla, Class, Order, Family, Genus, Species
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


def removeRIDSPLS(RID):
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


def getSPLS(request, stops, RID, PID):
    global base, stage, time1, TimeDiff
    try:
        while True:
                allJson = request.GET["all"]
                all = simplejson.loads(allJson)

                time1[RID] = time.time()
                base[RID] = 'Step 1 of 5: Selecting your chosen meta-variables...'

                myDir = 'media/usr_temp/' + str(request.user) + '/'
                path = str(myDir) + 'usr_norm_data.csv'

                with open(path, 'rb') as f:
                    savedDF = pd.read_csv(f, index_col=0, sep='\t')

                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])

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
                tempDF = savedDF.loc[savedDF['sampleid'].isin(quantSampleIDs)]

                if metaDictQuant:
                    for key in metaDictQuant:
                        valueList = [float(x) for x in metaDictQuant[key]]
                        tempDF = tempDF.loc[tempDF[key].isin(valueList)]

                wantedList = quantFields + ['sampleid']
                metaDF = tempDF[wantedList]

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

                result += 'Quantitative variables selected: ' + ", ".join(quantFields) + '\n'
                result += '\n===============================================\n'

                base[RID] = 'Step 1 of 5: Selecting your chosen meta-variables...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 2 of 5: Selecting your chosen taxa or KEGG level...'

                DepVar = 1
                finalDF = pd.DataFrame()
                if button3 == 1:
                    DepVar = int(all["DepVar_taxa"])
                    finalDF = getTaxaDF(selectAll, savedDF, metaDF, quantFields, DepVar, RID, stops, PID)

                if button3 == 2:
                    DepVar = int(all["DepVar_kegg"])
                    finalDF = getKeggDF(keggAll, savedDF, tempDF, quantFields, DepVar, RID, stops, PID)

                if button3 == 3:
                    DepVar = int(all["DepVar_nz"])
                    finalDF = getNZDF(nzAll, savedDF, tempDF, quantFields, DepVar, RID, stops, PID)

                # save location info to session
                myDir = 'media/temp/spls/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                base[RID] = 'Step 2 of 5: Selecting your chosen taxa or KEGG level...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 3 of 5: Calculating sPLS...'

                if DepVar == 1:
                    result += 'Dependent Variable: Relative Abundance' + '\n'
                elif DepVar == 2:
                    result += 'Dependent Variable: Species Richness' + '\n'
                elif DepVar == 3:
                    result += 'Dependent Variable: Species Diversity' + '\n'
                elif DepVar == 4:
                    result += 'Dependent Variable: Abundance (rRNA gene copies)' + '\n'
                result += '\n===============================================\n'

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

                if metaDictQuant:
                    for key in metaDictQuant:
                        valueList = [float(x) for x in metaDictQuant[key]]
                        meta_rDF = meta_rDF.loc[meta_rDF[key].isin(valueList)]

                wantedList = quantFields + ['sampleid']
                meta_rDF = meta_rDF[wantedList]
                meta_rDF.set_index('sampleid', drop=True, inplace=True)

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                r.assign("X", count_rDF)
                r.assign("Y", meta_rDF)
                r.assign("names", count_rDF.columns.values)
                r("colnames(X) <- names")

                r("library(mixOmics)")

                freqCut = all["freqCut"]
                num = int(freqCut.split('/')[0])
                den = int(freqCut.split('/')[1])
                r.assign("num", num)
                r.assign("den", den)

                uniqueCut = int(all["uniqueCut"])
                r.assign("uniqueCut", uniqueCut)

                r("ZeroVar <- nearZeroVar(X, freqCut=num/den, uniqueCut=uniqueCut)")
                r("List <- row.names(ZeroVar$Metrics)")
                r("X_new <- X[,-which(names(X) %in% List)]")
                r("if (length(X_new) == 0) {X_new <- X}")
                columns = r.get("ncol(X_new)")
                if columns == 0:
                    myDict = {}
                    myDict['error'] = "All predictor variables have zero variance.\nsPLS-Regr was aborted!"
                    res = simplejson.dumps(myDict)
                    return HttpResponse(res, content_type='application/json')

                r("maxK <- length(Y)")
                r("X_scaled <- scale(X_new, center=TRUE, scale=TRUE)")
                r("Y_scaled <- scale(Y, center=TRUE, scale=TRUE)")
                r("detach('package:mixOmics', unload=TRUE)")
                r("library(spls)")
                r("set.seed(1)")

                spls_string = "cv <- cv.spls(X_scaled, Y_scaled, scale.x=FALSE, scale.y=FALSE, eta=seq(0.1, 0.9, 0.1), K=c(1:maxK), plot.it=FALSE)"
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
                    result += '\n===============================================\n'

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
                    result += '\n===============================================\n'

                r("coef.f <- coef(f)")
                r("sum <- sum(coef.f != 0)")
                total = r.get("sum")

                base[RID] = 'Step 3 of 5: Calculating sPLS...done!'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
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
                        if stops[PID] == RID:
                            res = ''
                            return HttpResponse(res, content_type='application/json')
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    r("coef.f.rows <- row.names(coef.f)")
                    cf = r.get("coef.f")
                    rows = r.get("coef.f.rows")

                    coeffsDF = pd.DataFrame(cf,  columns=[quantFields], index=rows)
                    coeffsDF = coeffsDF.loc[(coeffsDF != 0).any(axis=1)]
                    coeffsDF.sort_index(inplace=True)
                    taxIDList = coeffsDF.index.values.tolist()

                    namesDF = pd.DataFrame()
                    if button3 == 1:
                        if selectAll == 1:
                            taxNameList = Kingdom.objects.filter(kingdomid__in=taxIDList).values('kingdomid', 'kingdomName')
                            namesDF = pd.DataFrame(list(taxNameList))
                            namesDF.rename(columns={'kingdomName': 'rank_name', 'kingdomid': 'rank_id'}, inplace=True)
                            namesDF.set_index('rank_id', inplace=True)
                        elif selectAll == 2:
                            taxNameList = Phyla.objects.filter(phylaid__in=taxIDList).values('phylaid', 'phylaName')
                            namesDF = pd.DataFrame(list(taxNameList))
                            namesDF.rename(columns={'phylaName': 'rank_name', 'phylaid': 'rank_id'}, inplace=True)
                            namesDF.set_index('rank_id', inplace=True)
                        elif selectAll == 3:
                            taxNameList = Class.objects.filter(classid__in=taxIDList).values('classid', 'className')
                            namesDF = pd.DataFrame(list(taxNameList))
                            namesDF.rename(columns={'className': 'rank_name', 'classid': 'rank_id'}, inplace=True)
                            namesDF.set_index('rank_id', inplace=True)
                        elif selectAll == 4:
                            taxNameList = Order.objects.filter(orderid__in=taxIDList).values('orderid', 'orderName')
                            namesDF = pd.DataFrame(list(taxNameList))
                            namesDF.rename(columns={'orderName': 'rank_name', 'orderid': 'rank_id'}, inplace=True)
                            namesDF.set_index('rank_id', inplace=True)
                        elif selectAll == 5:
                            taxNameList = Family.objects.filter(familyid__in=taxIDList).values('familyid', 'familyName')
                            namesDF = pd.DataFrame(list(taxNameList))
                            namesDF.rename(columns={'familyName': 'rank_name', 'familyid': 'rank_id'}, inplace=True)
                            namesDF.set_index('rank_id', inplace=True)
                        elif selectAll == 6:
                            taxNameList = Genus.objects.filter(genusid__in=taxIDList).values('genusid', 'genusName')
                            namesDF = pd.DataFrame(list(taxNameList))
                            namesDF.rename(columns={'genusName': 'rank_name', 'genusid': 'rank_id'}, inplace=True)
                            namesDF.set_index('rank_id', inplace=True)
                        elif selectAll == 7:
                            taxNameList = Species.objects.filter(speciesid__in=taxIDList).values('speciesid', 'speciesName')
                            namesDF = pd.DataFrame(list(taxNameList))
                            namesDF.rename(columns={'speciesName': 'rank_name', 'speciesid': 'rank_id'}, inplace=True)
                            namesDF.set_index('rank_id', inplace=True)

                    elif button3 == 2:
                        if keggAll == 1:
                            taxNameList = ko_lvl1.objects.using('picrust').filter(ko_lvl1_id__in=taxIDList).values('ko_lvl1_id', 'ko_lvl1_name')
                            namesDF = pd.DataFrame(list(taxNameList))
                            namesDF.rename(columns={'ko_lvl1_name': 'rank_name', 'ko_lvl1_id': 'rank_id'}, inplace=True)
                            namesDF.set_index('rank_id', inplace=True)
                        elif keggAll == 2:
                            taxNameList = ko_lvl2.objects.using('picrust').filter(ko_lvl2_id__in=taxIDList).values('ko_lvl2_id', 'ko_lvl2_name')
                            namesDF = pd.DataFrame(list(taxNameList))
                            namesDF.rename(columns={'ko_lvl2_name': 'rank_name', 'ko_lvl2_id': 'rank_id'}, inplace=True)
                            namesDF.set_index('rank_id', inplace=True)
                        elif keggAll == 3:
                            taxNameList = ko_lvl3.objects.using('picrust').filter(ko_lvl3_id__in=taxIDList).values('ko_lvl3_id', 'ko_lvl3_name')
                            namesDF = pd.DataFrame(list(taxNameList))
                            namesDF.rename(columns={'ko_lvl3_name': 'rank_name', 'ko_lvl3_id': 'rank_id'}, inplace=True)
                            namesDF.set_index('rank_id', inplace=True)

                    elif button3 == 3:
                        if nzAll == 1:
                            taxNameList = nz_lvl1.objects.using('picrust').filter(nz_lvl1_id__in=taxIDList).values('nz_lvl1_id', 'nz_lvl1_name')
                            namesDF = pd.DataFrame(list(taxNameList))
                            namesDF.rename(columns={'nz_lvl1_name': 'rank_name', 'nz_lvl1_id': 'rank_id'}, inplace=True)
                            namesDF.set_index('rank_id', inplace=True)
                        elif nzAll == 2:
                            taxNameList = nz_lvl2.objects.using('picrust').filter(nz_lvl2_id__in=taxIDList).values('nz_lvl2_id', 'nz_lvl2_name')
                            namesDF = pd.DataFrame(list(taxNameList))
                            namesDF.rename(columns={'nz_lvl2_name': 'rank_name', 'nz_lvl2_id': 'rank_id'}, inplace=True)
                            namesDF.set_index('rank_id', inplace=True)
                        elif nzAll == 3:
                            taxNameList = nz_lvl3.objects.using('picrust').filter(nz_lvl3_id__in=taxIDList).values('nz_lvl3_id', 'nz_lvl3_name')
                            namesDF = pd.DataFrame(list(taxNameList))
                            namesDF.rename(columns={'nz_lvl3_name': 'rank_name', 'nz_lvl3_id': 'rank_id'}, inplace=True)
                            namesDF.set_index('rank_id', inplace=True)
                        elif nzAll == 4:
                            taxNameList = nz_lvl4.objects.using('picrust').filter(nz_lvl4_id__in=taxIDList).values('nz_lvl4_id', 'nz_lvl4_name')
                            namesDF = pd.DataFrame(list(taxNameList))
                            namesDF.rename(columns={'nz_lvl4_name': 'rank_name', 'nz_lvl4_id': 'rank_id'}, inplace=True)
                            namesDF.set_index('rank_id', inplace=True)
                        elif nzAll == 5:
                            taxNameList = nz_lvl4.objects.using('picrust').filter(nz_lvl4_id__in=taxIDList).values('nz_lvl4_id', 'nz_lvl4_name')
                            df1 = pd.DataFrame(list(taxNameList))
                            df1.rename(columns={'nz_lvl4_name': 'rank_name', 'nz_lvl4_id': 'rank_id'}, inplace=True)
                            df1.set_index('rank_id', inplace=True)
                            taxNameList = nz_entry.objects.using('picrust').filter(nz_lvl5_id__in=taxIDList).values('nz_lvl5_id', 'nz_name')
                            df2 = pd.DataFrame(list(taxNameList))
                            df2.rename(columns={'nz_name': 'rank_name', 'nz_lvl5_id': 'rank_id'}, inplace=True)
                            df2.set_index('rank_id', inplace=True)
                            namesDF = pd.concat([df1, df2])
                        elif nzAll == 6:
                            taxNameList = nz_lvl4.objects.using('picrust').filter(nz_lvl4_id__in=taxIDList).values('nz_lvl4_id', 'nz_lvl4_name')
                            namesDF = pd.DataFrame(list(taxNameList))
                            namesDF.rename(columns={'nz_lvl1_name': 'rank_name', 'nz_lvl1_id': 'rank_id'}, inplace=True)
                            namesDF.set_index('rank_id', inplace=True)

                    namesDF.sort_index(inplace=True)
                    taxNameList = namesDF['rank_name'].values.tolist()

                    if button3 == 2:
                        if keggAll > 1:
                            taxNameList[:] = (item[:20] + '...' if len(item) > 20 else item for item in taxNameList)

                    elif button3 == 3:
                        if nzAll > 1:
                            taxNameList[:] = (item.split()[0] for item in taxNameList)

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    coeffsDF = pd.merge(namesDF, coeffsDF, left_index=True, right_index=True, how='inner')
                    coeffsDF.reset_index(inplace=True)
                    res_table = coeffsDF.to_html(classes="table display")
                    res_table = res_table.replace('border="1"', 'border="0"')
                    finalDict['res_table'] = str(res_table)

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    resultDF.reset_index(inplace=True)
                    resultDF.rename(columns={'index': 'sampleid'}, inplace=True)
                    pred_table = resultDF.to_html(classes="table display")
                    pred_table = pred_table.replace('border="1"', 'border="0"')
                    finalDict['pred_table'] = str(pred_table)

                    xAxisDict = {}
                    xAxisDict['categories'] = taxNameList
                    labelsDict = {}
                    labelsDict['rotation'] = 270
                    labelsDict['enabled'] = True
                    labelsDict['style'] = {'color': 'black', 'fontSize': '14px'}
                    xAxisDict['labels'] = labelsDict
                    xAxisDict['title'] = {'text': None}
                    xAxisDict['tickLength'] = 0

                    yAxisDict = {}
                    yAxisDict['categories'] = quantFields
                    yAxisDict['labels'] = {'style': {'color': 'black', 'fontSize': '14px'}}
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
                            if stops[PID] == RID:
                                res = ''
                                return HttpResponse(res, content_type='application/json')
                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
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
                    clustDF = coeffsDF.drop('rank_id', axis=1)

                    row, col = clustDF.shape

                    method = all['methodVal']
                    metric = all['metricVal']

                    path = "media/temp/spls/Rplots/" + str(RID) + ".spls.pdf"
                    if os.path.exists(path):
                        os.remove(path)

                    if not os.path.exists('media/temp/spls/Rplots'):
                        os.makedirs('media/temp/spls/Rplots')

                    height = 2.5 + 0.2*row
                    width = 3.5 + 0.2*(col-1)
                    file = "pdf('media/temp/spls/Rplots/" + str(RID) + ".spls.pdf', height=" + str(height) + ", width=" + str(width) + ", onefile=FALSE)"
                    r.assign("cmd", file)
                    r("eval(parse(text=cmd))")

                    r.assign("df", clustDF[quantFields])
                    r("df <- as.matrix(df)")

                    r.assign("rows", taxNameList)
                    r("rownames(df) <- rows")
                    r("library(pheatmap)")
                    r("library(RColorBrewer)")
                    r("col.pal <- brewer.pal(9,'RdBu')")

                    if row > 2 and col > 3:
                        hmap_str = "pheatmap(df, fontsize=12, color=col.pal, clustering_method='" + str(method) + "', clustering_distance_rows='" + str(metric) + "', clustering_distance_cols='" + str(metric) + "')"
                        r.assign("cmd", hmap_str)
                        r("eval(parse(text=cmd))")

                    if row > 2 and col <= 3:
                        hmap_str = "pheatmap(df, color=col.pal, cluster_col=FALSE, clustering_method='" + str(method) + "', clustering_distance_rows='" + str(metric) + "')"
                        r.assign("cmd", hmap_str)
                        r("eval(parse(text=cmd))")

                    if row <= 2 and col > 3:
                        hmap_str = "pheatmap(df, color=col.pal, cluster_row=FALSE, clustering_method='" + str(method) + "', clustering_distance_cols='" + str(metric) + "')"
                        r.assign("cmd", hmap_str)
                        r("eval(parse(text=cmd))")

                    if row <= 2 and col <= 3:
                        hmap_str = "pheatmap(df, color=col.pal, cluster_col=FALSE, cluster_row=FALSE)"
                        r.assign("cmd", hmap_str)
                        r("eval(parse(text=cmd))")

                    r("dev.off()")

                finalDict['text'] = result

                base[RID] = 'Step 4 of 5: Formatting graph data for display...done!'
                base[RID] = 'Step 5 of 5: Formatting sPLS coefficient table...'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                removeRIDSPLS(RID)
                return HttpResponse(res, content_type='application/json')

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "Error with sPLS-Regr!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = simplejson.dumps(myDict)
            removeRIDSPLS(RID)
            return HttpResponse(res, content_type='application/json')


def getTaxaDF(selectAll, savedDF, metaDF, allFields, DepVar, RID, stops, PID):
    global base
    try:
        base[RID] = 'Step 2 of 5: Selecting your chosen taxa or KEGG level...'
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
            myDict = {}
            myDict['error'] = "Error with sPLS-Regr!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getKeggDF(keggAll, savedDF, tempDF, allFields, DepVar, RID, stops, PID):
    global base
    try:
        base[RID] = 'Step 2 of 5: Selecting your chosen taxa or KEGG level...'

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

        path = 'media/temp/spls/' + str(RID)
        if not os.path.exists(path):
            os.makedirs(path)

        if os.name == 'nt':
            numcore = 1
            listDF = np.array_split(picrustDF, numcore)
            processes = [threading.Thread(target=sumStuff, args=(listDF[x], koDict, RID, x, stops, PID)) for x in xrange(numcore)]
        else:
            numcore = local_cfg.usr_numcore
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
            path = 'media/temp/spls/'+str(RID)+'/file%d.temp' % i
            frame = pd.read_csv(path)
            picrustDF = picrustDF.append(frame, ignore_index=True)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        shutil.rmtree('media/temp/spls/'+str(RID))
        picrustDF.set_index('speciesid', inplace=True)

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
            myDict['error'] = "Error with sPLS-Regr!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
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

        path = 'media/temp/spls/' + str(RID)
        if not os.path.exists(path):
            os.makedirs(path)

        if os.name == 'nt':
            numcore = 1
            listDF = np.array_split(picrustDF, numcore)
            processes = [threading.Thread(target=sumStuff, args=(listDF[x], nzDict, RID, x, stops, PID)) for x in xrange(numcore)]
        else:
            numcore = local_cfg.usr_numcore
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
            path = 'media/temp/spls/'+str(RID)+'/file%d.temp' % i
            frame = pd.read_csv(path)
            picrustDF = picrustDF.append(frame, ignore_index=True)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        shutil.rmtree('media/temp/spls/'+str(RID))
        picrustDF.set_index('speciesid', inplace=True)

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
            myDict['error'] = "Error with sPLS-Regr!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def sumStuff(slice, koDict, RID, num, stops, PID):
    global base
    db.close_old_connections()

    f = open('media/temp/spls/'+str(RID)+'/file'+str(num)+".temp", 'w')

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


def removeSPLSFiles(request):
    if request.is_ajax():
        RID = request.GET["all"]

        file = "media/temp/spls/Rplots/" + str(RID) + ".spls.pdf"
        if os.path.exists(file):
            os.remove(file)

        file = "media/temp/spls/" + str(RID) + ".pkl"
        if os.path.exists(file):
            os.remove(file)

        file = "media/temp/spls/" + str(RID) + ".csv"
        if os.path.exists(file):
            os.remove(file)

        return HttpResponse()


def getTabSPLS(request):
    if request.is_ajax():
        RID = request.GET["all"]
        myDir = 'media/temp/spls/'
        fileName = str(myDir) + str(RID) + '.pkl'
        savedDF = pd.read_pickle(fileName)

        myDir = 'media/temp/spls/'
        fileName = str(myDir) + str(RID) + '.csv'
        savedDF.to_csv(fileName)

        myDict = {}
        myDir = 'temp/spls/'
        fileName = str(myDir) + str(RID) + '.csv'
        myDict['name'] = str(fileName)
        res = simplejson.dumps(myDict)

        return HttpResponse(res, content_type='application/json')
