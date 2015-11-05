from django.http import HttpResponse
from django.db.models import Sum
import numpy as np
import pandas as pd
import pickle
from pyper import *
from scipy import stats
from scipy.spatial.distance import *
import simplejson

from database.spls.spls_DF import metaData, normalizeData
from database.models import Sample, Profile, Kingdom, Phyla, Class, Order, Family, Genus, Species
from database.utils import multidict, taxaProfileDF


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}


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
            return True
        else:
            return False
    except:
        return False


def getSPLS(request):
    try:
        global base, time1, TimeDiff
        samples = Sample.objects.all()
        samples.query = pickle.loads(request.session['selected_samples'])
        selected = samples.values_list('sampleid')
        qs1 = Sample.objects.all().filter(sampleid__in=selected)

        if request.is_ajax():
            allJson = request.GET["all"]
            all = simplejson.loads(allJson)

            RID = str(all["RID"])
            time1[RID] = time.time()
            base[RID] = 'Step 1 of 8: Querying database...'

            DepVar = int(all["DepVar"])
            taxaLevel = int(all["taxa"])
            NormMeth = int(all["NormMeth"])
            Iters = int(all["Iters"])
            NormVal = all["NormVal"]
            size = int(all["MinSize"])

            countList = []
            for sample in qs1:
                total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                if total['count__sum'] is not None:
                    countList.append(total['count__sum'])

            minSize = int(min(countList))
            medianSize = int(np.median(np.array(countList)))
            maxSize = int(max(countList))

            if NormVal == "min":
                NormReads = minSize
            elif NormVal == "median":
                NormReads = medianSize
            elif NormVal == "max":
                NormReads = maxSize
            elif NormVal == "none":
                NormReads = -1
            else:
                NormReads = int(all["NormVal"])

            # Remove samples if below the sequence threshold set by user (rarefaction)
            newList = []
            result = ''
            if taxaLevel == 1:
                result = result + 'Taxa level: Kingdom' + '\n'
            elif taxaLevel == 2:
                result = result + 'Taxa level: Phyla' + '\n'
            elif taxaLevel == 3:
                result = result + 'Taxa level: Class' + '\n'
            elif taxaLevel == 4:
                result = result + 'Taxa level: Order' + '\n'
            elif taxaLevel == 5:
                result = result + 'Taxa level: Family' + '\n'
            elif taxaLevel == 6:
                result = result + 'Taxa level: Genus' + '\n'
            elif taxaLevel == 7:
                result = result + 'Taxa level: Species' + '\n'

            metaStr = all["metaQuant"]
            fieldList = []
            valueList = []
            metaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStr)
            for key in sorted(metaDict):
                fieldList.append(key)
                valueList.extend(metaDict[key])

            idStr = all["metaIDs"]
            idDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(idStr)

            if DepVar == 4:
                idList = []
                for i in xrange(len(selected)):
                    idList.append(selected[i][0])
                idDict['rRNA_copies'] = idList

            result += 'Quantitative variables selected: ' + ", ".join(fieldList) + '\n'
            result += '===============================================\n'
            result += '\nData Normalization:\n'

            # Limit reads to max value
            if NormMeth == 1:
                for sample in qs1:
                    total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                    if total['count__sum'] is not None:
                        id = sample.sampleid
                        newList.append(id)

            elif NormMeth == 2 or NormMeth == 3:
                if NormReads > maxSize:
                    NormReads = medianSize
                    result += 'The subsample size was too high and automatically reset to the median value...\n'

                for sample in qs1:
                    total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                    if NormMeth == 2:
                        if total['count__sum'] is not None and int(total['count__sum']) >= NormReads:
                            id = sample.sampleid
                            newList.append(id)
                    else:
                        if total['count__sum'] is not None:
                            id = sample.sampleid
                            newList.append(id)

                # If user set reads too high sample list will be blank
                if not newList:
                    NormReads = medianSize
                    for sample in qs1:
                        total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                        if total['count__sum'] is not None and int(total['count__sum']) >= NormReads:
                            id = sample.sampleid
                            newList.append(id)

            elif NormMeth == 4 or NormMeth == 5 or NormMeth == 6:
                if size > maxSize:
                    size = medianSize
                    result += 'The minimum sample size was too high and automatically reset to the median value...\n'
                for sample in qs1:
                    total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                    if total['count__sum'] is not None and int(total['count__sum']) >= size:
                        id = sample.sampleid
                        newList.append(id)

                # If user set reads too high sample list will be blank
                if not newList:
                    size = medianSize
                    for sample in qs1:
                        total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                        if total['count__sum'] is not None and int(total['count__sum']) >= size:
                            id = sample.sampleid
                            newList.append(id)

            metaDF = metaData(idDict)

            lenA, col = metaDF.shape

            metaDF = metaDF.ix[newList]
            metaDF.dropna(inplace=True)
            lenB, col = metaDF.shape

            selectRem = len(selected) - lenA
            normRem = lenA - lenB

            result += str(lenB) + ' selected samples were included in the final analysis.\n'
            if normRem > 0:
                result += str(normRem) + ' samples did not met the desired normalization criteria.\n'
            if selectRem:
                result += str(selectRem) + ' samples were deselected by the user.\n'

            # Create unique list of samples in meta dataframe (may be different than selected samples)
            myList = metaDF.index.values.tolist()

            taxaDF = taxaProfileDF(myList)

            base[RID] = 'Step 1 of 5: Querying database...done!'

            try:
                base[RID] = 'Step 2 of 5: Normalizing data...'

                # Select only the taxa of interest if user used the selectAll button
                taxaDict = {}
                if taxaLevel == 1:
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('kingdomid', flat='True').distinct()
                    taxaDict['Kingdom'] = qs3
                elif taxaLevel == 2:
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('phylaid', flat='True').distinct()
                    taxaDict['Phyla'] = qs3
                elif taxaLevel == 3:
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('classid', flat='True').distinct()
                    taxaDict['Class'] = qs3
                elif taxaLevel == 4:
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('orderid', flat='True').distinct()
                    taxaDict['Order'] = qs3
                elif taxaLevel == 5:
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('familyid', flat='True').distinct()
                    taxaDict['Family'] = qs3
                elif taxaLevel == 6:
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('genusid', flat='True').distinct()
                    taxaDict['Genus'] = qs3
                elif taxaLevel == 7:
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('speciesid', flat='True').distinct()
                    taxaDict['Species'] = qs3

                normDF, DESeq_error = normalizeData(taxaDF, taxaDict, myList, NormMeth, NormReads, metaDF, Iters)

                finalDict = {}
                if NormMeth == 1:
                    result += 'No normalization was performed...\n'
                elif NormMeth == 2 or NormMeth == 3:
                    result = result + 'Data were rarefied to ' + str(NormReads) + ' sequence reads...\n'
                elif NormMeth == 4:
                    result += 'Data were normalized by the total number of sequence reads...\n'
                elif NormMeth == 5 and DESeq_error == 'no':
                    result += 'Data were normalized by DESeq2...\n'
                elif NormMeth == 5 and DESeq_error == 'yes':
                    result += 'DESeq2 cannot run estimateSizeFactors...\n'
                    result += 'Analysis was run without normalization...\n'
                    result += 'To try again, please select fewer samples or another normalization method...\n'
                elif NormMeth == 6 and DESeq_error == 'no':
                    result += 'Data were normalized by DESeq2 with variance stabilization...\n'
                elif NormMeth == 6 and DESeq_error == 'yes':
                    result += 'DESeq2 cannot run estimateSizeFactors...\n'
                    result += 'Analysis was run without normalization...\n'
                    result += 'To try again, please select fewer samples or another normalization method...\n'
                result += '===============================================\n\n'

                normDF.set_index('sampleid', inplace=True)

                finalDF = pd.merge(metaDF, normDF, left_index=True, right_index=True)

                if DepVar == 4:
                    finalDF['copies'] = finalDF.abund / NormReads * finalDF.rRNA_copies
                    finalDF[['abund', 'copies', 'rich', 'diversity']] = finalDF[['abund', 'copies', 'rich', 'diversity']].astype(float)
                else:
                    finalDF[['abund', 'rich', 'diversity']] = finalDF[['abund', 'rich', 'diversity']].astype(float)

                finalDF.reset_index(inplace=True)

                if DepVar == 1:
                    normDF = finalDF.pivot(index='sampleid', columns='taxa_id', values='abund')
                if DepVar == 2:
                    normDF = finalDF.pivot(index='sampleid', columns='taxa_id', values='rich')
                if DepVar == 3:
                    normDF = finalDF.pivot(index='sampleid', columns='taxa_id', values='diversity')
                if DepVar == 4:
                    normDF = finalDF.pivot(index='sampleid', columns='taxa_id', values='copies')

                base[RID] = 'Step 2 of 5: Normalizing data...done!'

            except:
                myDict = {}
                myDict['error'] = 'Your selections resulted in no valid observations'
                res = simplejson.dumps(myDict)
                return HttpResponse(res, content_type='application/json')

            base[RID] = 'Step 3 of 5: Calculating sPLS...'

            if os.name == 'nt':
                r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
            else:
                r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

            r.assign("X", normDF)
            r.assign("Y", metaDF[fieldList])
            r.assign("names", normDF.columns.values)
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
            r("cv <- cv.spls(X_scaled, Y_scaled, scale.x=FALSE, scale.y=FALSE, eta=seq(0.1, 0.9, 0.1), K=c(1:5), plot.it=FALSE)")
            r("f <- spls(X_scaled, Y_scaled, scale.x=FALSE, scale.y=FALSE, eta=cv$eta.opt, K=cv$K.opt)")

            r("out <- capture.output(print(f))")
            fout = r.get("out")

            for i in fout:
                result += str(i) + '\n'
            result += '===============================================\n\n'

            r("set.seed(1)")
            r("ci.f <- ci.spls(f, plot.it=FALSE, plot.fix='y')")
            r("cis <- ci.f$cibeta")
            r("cf <- correct.spls(ci.f, plot.it=FALSE)")
            r("out <- capture.output(cis)")
            fout = r.get("out")
            result += '\n\nBootstrapped confidence intervals of coefficients:\n'
            for i in fout:
                result += str(i) + '\n'
            result += '===============================================\n\n'


            r("coef.f <- coef(f)")
            r("sum <- sum(coef.f != 0)")
            total = r.get("sum")

            base[RID] = 'Step 3 of 5: Calculating sPLS...done!'
            base[RID] = 'Step 4 of 5: Formatting graph data for display...'

            if total > 0:
                r("pred.f <- predict(f, type='fit')")

                r("library(DMwR)")
                r("pred.ns <- unscale(pred.f, Y_scaled)")
                r("pred.ns.rows <- row.names(pred.ns)")
                pred = r.get("pred.ns")
                rows = r.get("pred.ns.rows")

                predList = ['pred_' + s for s in fieldList]
                predDF = pd.DataFrame(pred,  columns=[predList], index=rows)
                predDF = predDF.applymap(round)
                predDF.sort_index(inplace=True)
                finalDF = pd.merge(metaDF, predDF, left_index=True, right_index=True)
                result += 'sPLS Model Fit (y = mx + b):\n'
                result += 'y = predicted\n'
                result += 'x = observed\n\n'

                for i in xrange(len(fieldList)):
                    x = finalDF[fieldList[i]].astype(float).values.tolist()
                    y = finalDF[predList[i]].astype(float).values.tolist()

                    slp, inter, r_value, p, se = stats.linregress(x, y)
                    r_sq = r_value * r_value

                    result += 'Variable: ' + str(fieldList[i]) + '\n'
                    result += 'Slope (m): ' + str(slp) + '\n'
                    result += 'Intercept (b): ' + str(inter) + '\n'
                    result += 'R2: ' + str(r_sq) + '\n'
                    result += 'Std Error: ' + str(se) + '\n\n\n'

                r("coef.f.rows <- row.names(coef.f)")
                cf = r.get("coef.f")
                rows = r.get("coef.f.rows")

                coeffsDF = pd.DataFrame(cf,  columns=[fieldList], index=rows)
                coeffsDF = coeffsDF.loc[(coeffsDF != 0).any(axis=1)]
                coeffsDF.sort_index(inplace=True)
                taxIDList = coeffsDF.index.values.tolist()

                namesDF = pd.DataFrame()
                if taxaLevel == 1:
                    taxNameList = Kingdom.objects.filter(kingdomid__in=taxIDList).values('kingdomid', 'kingdomName')
                    namesDF = pd.DataFrame(list(taxNameList))
                    namesDF.set_index('kingdomid', inplace=True)
                    namesDF.rename(columns={'kingdomName': 'taxa_name'}, inplace=True)
                elif taxaLevel == 2:
                    taxNameList = Phyla.objects.filter(phylaid__in=taxIDList).values('phylaid', 'phylaName')
                    namesDF = pd.DataFrame(list(taxNameList))
                    namesDF.set_index('phylaid', inplace=True)
                    namesDF.rename(columns={'phylaName': 'taxa_name'}, inplace=True)
                elif taxaLevel == 3:
                    taxNameList = Class.objects.filter(classid__in=taxIDList).values('classid', 'className')
                    namesDF = pd.DataFrame(list(taxNameList))
                    namesDF.set_index('classid', inplace=True)
                    namesDF.rename(columns={'className': 'taxa_name'}, inplace=True)
                elif taxaLevel == 4:
                    taxNameList = Order.objects.filter(orderid__in=taxIDList).values('orderid', 'orderName')
                    namesDF = pd.DataFrame(list(taxNameList))
                    namesDF.set_index('orderid', inplace=True)
                    namesDF.rename(columns={'orderName': 'taxa_name'}, inplace=True)
                elif taxaLevel == 5:
                    taxNameList = Family.objects.filter(familyid__in=taxIDList).values('familyid', 'familyName')
                    namesDF = pd.DataFrame(list(taxNameList))
                    namesDF.set_index('familyid', inplace=True)
                    namesDF.rename(columns={'familyName': 'taxa_name'}, inplace=True)
                elif taxaLevel == 6:
                    taxNameList = Genus.objects.filter(genusid__in=taxIDList).values('genusid', 'genusName')
                    namesDF = pd.DataFrame(list(taxNameList))
                    namesDF.set_index('genusid', inplace=True)
                    namesDF.rename(columns={'genusName': 'taxa_name'}, inplace=True)
                elif taxaLevel == 7:
                    taxNameList = Species.objects.filter(speciesid__in=taxIDList).values('speciesid', 'speciesName')
                    namesDF = pd.DataFrame(list(taxNameList))
                    namesDF.set_index('speciesid', inplace=True)
                    namesDF.rename(columns={'speciesName': 'taxa_name'}, inplace=True)

                namesDF.sort_index(inplace=True)
                taxaNameList = namesDF['taxa_name'].values.tolist()

                coeffsDF = pd.merge(namesDF, coeffsDF, left_index=True, right_index=True)
                coeffsDF.reset_index(inplace=True)
                res_table = coeffsDF.to_html(classes="table display")
                res_table = res_table.replace('border="1"', 'border="0"')
                finalDict['res_table'] = str(res_table)

                if taxaLevel == 2:
                    clustDF = coeffsDF.drop('phylaid', axis=1)
                    clustDF.set_index('taxa_name', inplace=True)
                elif taxaLevel == 3:
                    clustDF = coeffsDF.drop('classid', axis=1)
                    clustDF.set_index('taxa_name', inplace=True)
                elif taxaLevel == 4:
                    clustDF = coeffsDF.drop('orderid', axis=1)
                    clustDF.set_index('taxa_name', inplace=True)
                elif taxaLevel == 5:
                    clustDF = coeffsDF.drop('familyid', axis=1)
                    clustDF.set_index('taxa_name', inplace=True)
                elif taxaLevel == 6:
                    clustDF = coeffsDF.drop('genusid', axis=1)
                    clustDF.set_index('taxa_name', inplace=True)
                elif taxaLevel == 7:
                    clustDF = coeffsDF.drop('speciesid', axis=1)
                    clustDF.set_index('taxa_name', inplace=True)

                finalDF.reset_index(inplace=True)
                finalDF.rename(columns={'index': 'sampleid'}, inplace=True)
                pred_table = finalDF.to_html(classes="table display")
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
                catList = fieldList
                yAxisDict['categories'] = catList
                yAxisDict['title'] = {'text': None}

                seriesList = []
                seriesDict = {}
                seriesDict['borderWidth'] = '1'

                row, col = coeffsDF.shape
                dataList = []
                for i in xrange(row):
                    for j in xrange(len(fieldList)):
                        val = round(coeffsDF[fieldList[j]].iloc[i], 5)
                        tup = (i, j, val)
                        obsList = list(tup)
                        dataList.append(obsList)

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
                if taxaLevel == 2:
                    clustDF = coeffsDF.drop('phylaid', axis=1)
                elif taxaLevel == 3:
                    clustDF = coeffsDF.drop('classid', axis=1)
                elif taxaLevel == 4:
                    clustDF = coeffsDF.drop('orderid', axis=1)
                elif taxaLevel == 5:
                    clustDF = coeffsDF.drop('familyid', axis=1)
                elif taxaLevel == 6:
                    clustDF = coeffsDF.drop('genusid', axis=1)
                elif taxaLevel == 7:
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

                height = 200 + 15*row
                width = 200 + 20*(col-1)
                file = "jpeg('media/Rplots/" + str(user) + ".spls.jpg', height=" + str(height) + ", width=" + str(width) + ")"
                r.assign("cmd", file)
                r("eval(parse(text=cmd))")

                r.assign("df", clustDF[fieldList])
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

            res = simplejson.dumps(finalDict)
            return HttpResponse(res, content_type='application/json')

    except Exception as e:
        print "Error with SPLS CAT: ", e
        state = "Error with SPLS CAT: " + str(e)

        myDict = {}
        myDict['error'] = state
        res = simplejson.dumps(myDict)
        return HttpResponse(res, content_type='application/json')


def removegraphSPLS(request):
    name = request.user
    ip = request.META.get('REMOTE_ADDR')
    user = str(name) + "." + str(ip)

    file = "media/Rplots/" + str(user) + ".spls.jpg"
    if os.path.exists(file):
        os.remove(file)