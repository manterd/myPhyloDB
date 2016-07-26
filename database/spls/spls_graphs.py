import datetime
from django.http import HttpResponse
import logging
import pandas as pd
from pyper import *
from scipy import stats
import simplejson

from database.models import Kingdom, Phyla, Class, Order, Family, Genus, Species
from database.models import ko_lvl1, ko_lvl2, ko_lvl3
from database.models import nz_lvl1, nz_lvl2, nz_lvl3, nz_lvl4
from database.utils import multidict
from database.utils_kegg import getTaxaDF, getKeggDF, getNZDF
import database.queue


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getSPLS(request, stops, RID, PID):
    try:
        while True:
                allJson = request.body.split('&')[0]
                all = simplejson.loads(allJson)
                database.queue.base(RID, 'Step 1 of 5: Selecting your chosen meta-variables...')
                myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
                path = str(myDir) + 'usr_norm_data.csv'

                with open(path, 'rb') as f:
                    savedDF = pd.read_csv(f, index_col=0, sep=',')

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

                database.queue.base(RID, 'Step 1 of 5: Selecting your chosen meta-variables...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.base(RID, 'Step 2 of 5: Selecting your chosen taxa or KEGG level...')

                DepVar = 1
                finalDF = pd.DataFrame()
                if button3 == 1:
                    DepVar = int(all["DepVar_taxa"])
                    finalDF, missingList = getTaxaDF('rel_abund', selectAll, '', savedDF, metaDF, quantFields, DepVar, RID, stops, PID)
                    result += 'The following PGPRs were not detected: ' + ", ".join(missingList) + '\n'
                    result += '===============================================\n'

                if button3 == 2:
                    DepVar = int(all["DepVar_kegg"])
                    finalDF = getKeggDF('rel_abund', keggAll, '', savedDF, tempDF, quantFields, DepVar, RID, stops, PID)

                if button3 == 3:
                    DepVar = int(all["DepVar_nz"])
                    finalDF = getNZDF('rel_abund', nzAll, '', savedDF, tempDF, quantFields, DepVar, RID, stops, PID)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/spls/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                database.queue.base(RID, 'Step 2 of 5: Selecting your chosen taxa or KEGG level...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.base(RID, 'Step 3 of 5: Calculating sPLS...')

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
                    myDict = {'error': "All predictor variables have zero variance.\nsPLS-Regr was aborted!"}
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
                    myDict = {'error': "Analysis did not converge.\nsPLS-Regr was aborted!"}
                    res = simplejson.dumps(myDict)
                    return HttpResponse(res, content_type='application/json')


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

                database.queue.base(RID, 'Step 3 of 5: Calculating sPLS...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.base(RID, 'Step 4 of 5: Formatting graph data for display...')

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
                            print taxNameList
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
                            namesDF = pd.DataFrame(list(taxNameList))
                            namesDF.rename(columns={'nz_lvl4_name': 'rank_name', 'nz_lvl4_id': 'rank_id'}, inplace=True)
                            namesDF.set_index('rank_id', inplace=True)
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

                    path = "myPhyloDB/media/temp/spls/Rplots/" + str(RID) + ".spls.pdf"
                    if os.path.exists(path):
                        os.remove(path)

                    if not os.path.exists('myPhyloDB/media/temp/spls/Rplots'):
                        os.makedirs('myPhyloDB/media/temp/spls/Rplots')

                    height = 2.5 + 0.2*row
                    width = 3.5 + 0.2*(col-1)
                    file = "pdf('myPhyloDB/media/temp/spls/Rplots/" + str(RID) + ".spls.pdf', height=" + str(height) + ", width=" + str(width) + ", onefile=FALSE)"
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

                database.queue.base(RID, 'Step 4 of 5: Formatting graph data for display...done!')
                database.queue.base(RID, 'Step 5 of 5: Formatting sPLS coefficient table...')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                return HttpResponse(res, content_type='application/json')

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "Error with sPLS-Regr!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')
