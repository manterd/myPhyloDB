from pcoa_DF import catPCoAMetaDF, normalizePCoA
from django.http import HttpResponse
from database.models import Sample, Profile
from django.db.models import Sum
import numpy as np
from numpy import *
import pandas as pd
import pickle
from stats.distance import wOdum
from scipy import stats
from scipy.spatial.distance import *
import simplejson
from database.utils import multidict, ordered_set, taxaProfileDF
from pyper import *


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}


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


def removeRIDPCOA(request):
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


def getCatPCoAData(request):
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

            taxaLevel = int(all["taxa"])
            distance = int(all["distance"])
            PC1 = int(all["PC1"])
            PC2 = int(all["PC2"])
            test = int(all["test"])
            alpha = float(all["alpha"])
            perms = int(all["perms"])
            NormMeth = int(all["NormMeth"])
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
            if taxaLevel == 0:
                result = result + 'Taxa level: Kingdom' + '\n'
            elif taxaLevel == 1:
                result = result + 'Taxa level: Phyla' + '\n'
            elif taxaLevel == 2:
                result = result + 'Taxa level: Class' + '\n'
            elif taxaLevel == 3:
                result = result + 'Taxa level: Order' + '\n'
            elif taxaLevel == 4:
                result = result + 'Taxa level: Family' + '\n'
            elif taxaLevel == 5:
                result = result + 'Taxa level: Genus' + '\n'
            elif taxaLevel == 6:
                result = result + 'Taxa level: Species' + '\n'

            if distance == 1:
                result = result + 'Distance score: Manhattan' + '\n'
            elif distance == 2:
                result = result + 'Distance score: Euclidean' + '\n'
            elif distance == 3:
                result = result + 'Distance score: Canberra' + '\n'
            elif distance == 4:
                result = result + 'Distance score: Bray-Curtis' + '\n'
            elif distance == 5:
                result = result + 'Distance score: Kulczynski' + '\n'
            elif distance == 6:
                result = result + 'Distance score: Jaccard' + '\n'
            elif distance == 7:
                result = result + 'Distance score: Gower' + '\n'
            elif distance == 8:
                result = result + 'Distance score: altGower' + '\n'
            elif distance == 9:
                result = result + 'Distance score: Morisita' + '\n'
            elif distance == 10:
                result = result + 'Distance score: Horn' + '\n'
            elif distance == 11:
                result = result + 'Distance score: Mountford' + '\n'
            elif distance == 12:
                result = result + 'Distance score: Binomial' + '\n'
            elif distance == 13:
                result = result + 'Distance score: Chao' + '\n'
            elif distance == 14:
                result = result + 'Distance score: Cao' + '\n'
            elif distance == 15:
                result = result + 'Distance score: wOdum' + '\n'

            metaStrCat = all["metaCat"]
            fieldListCat = []
            if metaStrCat:
                metaDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStrCat)
                for key in metaDictCat:
                    fieldListCat.append(key)

            metaStrQuant = all["metaQuant"]
            fieldListQuant = []
            if metaStrQuant:
                metaDictQuant = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStrQuant)
                for key in metaDictQuant:
                    fieldListQuant.append(key)

            metaStr = all["meta"]
            fieldList = []
            valueList = []
            metaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStr)
            for key in metaDict:
                fieldList.append(key)
                valueList.append(metaDict[key])

            result = result + 'Categorical variables selected: ' + ", ".join(fieldListCat) + '\n'
            result = result + 'Quantitative variables selected: ' + ", ".join(fieldListQuant) + '\n'
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

            qs2 = Sample.objects.all().filter(sampleid__in=newList)
            metaDF = catPCoAMetaDF(qs2, metaDict)
            metaDF.dropna(subset=fieldList, inplace=True)
            metaDF.sort(columns='sample_name', inplace=True)
            totalSamp, cols = metaDF.shape

            normRem = len(countList) - len(newList)
            selectRem = len(newList) - totalSamp

            result += str(totalSamp) + ' selected samples were included in the final analysis.\n'
            if normRem > 0:
                result += str(normRem) + ' samples did not met the desired normalization criteria.\n'
            if selectRem:
                result += str(selectRem) + ' samples were deselected by the user.\n'

            # Create combined metadata column
            if len(fieldListCat) > 1:
                for index, row in metaDF.iterrows():
                    metaDF.ix[index, 'merge'] = " & ".join(row[fieldListCat])
            else:
                metaDF['merge'] = metaDF[fieldListCat[0]]

            # Create unique list of samples in meta dataframe (may be different than selected samples)
            myList = metaDF['sampleid'].tolist()
            mySet = list(ordered_set(myList))

            taxaDF = taxaProfileDF(mySet)

            base[RID] = 'Step 1 of 8: Querying database...done!'
            base[RID] = 'Step 2 of 8: Normalizing data...'

            # Sum by taxa level
            taxaDF = taxaDF.groupby(level=taxaLevel).sum()

            normDF, DESeq_error = normalizePCoA(taxaDF, taxaLevel, mySet, NormMeth, NormReads, metaDF)

            normDF.sort_index(inplace=True)

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

            base[RID] = 'Step 2 of 8: Normalizing data...done!'
            base[RID] = 'Step 3 of 8: Calculating distance matrix...'

            metaDF.set_index('sampleid', inplace=True)
            metaDF.sort_index(inplace=True)

            if os.name == 'nt':
                r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
            else:
                r = R(RCMD="R/R-Linux/bin/R")

            r.assign("data", normDF)
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
                datamtx = asarray(normDF)
                dists = wOdum(datamtx, alpha)
                r.assign("dist", dists)

            r("mat <- as.matrix(dist, diag=TRUE, upper=TRUE)")
            mat = r.get("mat")
            rowList = metaDF['sample_name'].tolist()
            distDF = pd.DataFrame(mat, columns=[rowList], index=rowList)

            base[RID] = 'Step 3 of 8: Calculating distance matrix...done!'
            base[RID] = 'Step 4 of 8: Principal coordinates analysis...'

            trtLength = 0
            for i in valueList:
                if len(i) > trtLength:
                    trtLength = len(i)

            if trtLength > 1:
                r.assign("meta", metaDF)
                trtString = " + ".join(fieldListCat)
                pcoa_string = "ord <- capscale(mat ~ " + str(trtString) + ", meta)"
                r.assign("cmd", pcoa_string)
                r("eval(parse(text=cmd))")

                envFit = ''
                if len(fieldListQuant) > 0:
                    trtString = " + ".join(fieldListQuant)
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

                bigf = ''
                if perms <= 10:
                    bigf = 'Please increase the number of permutations...'
                elif len(fieldListCat) == 0:
                    bigf = 'No categorical variables are available for perMANOVA/betaDisper analysis'
                elif perms > 10 and len(fieldListCat) > 0:
                    if test == 1:

                        base[RID] = 'Step 4 of 8: Principal coordinates analysis...done!'
                        base[RID] = 'Step 5 of 8: Performing perMANOVA...'

                        for i in fieldListCat:
                            factor_string = str(i) + " <- factor(meta$" + str(i) + ")"
                            r.assign("cmd", factor_string)
                            r("eval(parse(text=cmd))")

                        r.assign("perms", perms)
                        trtString = " + ".join(fieldListCat)
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

                        for i in fieldList:
                            factor_string = str(i) + " <- factor(meta$" + str(i) + ")"
                            r.assign("cmd", factor_string)
                            r("eval(parse(text=cmd))")

                        r.assign("perms", perms)
                        for i in fieldListCat:
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
            else:
                state = "Your selected variable(s) only have one treatment level, please select additional data!"
                myDict = {}
                myDict['error'] = state
                res = simplejson.dumps(myDict)
                return HttpResponse(res, content_type='application/json')

            base[RID] = 'Step 6 of 8: Formatting graph data for display...'

            seriesList = []
            xAxisDict = {}
            yAxisDict = {}

            CAP1 = PC1 + int(len(fieldList)) + 2
            CAP2 = PC2 + int(len(fieldList)) + 2

            if fieldListCat:
                grouped = pcoaDF.groupby(fieldListCat)
                for name, group in grouped:
                    dataList = group.icol([CAP1, CAP2]).values.astype(float).tolist()
                    trt = name
                    seriesDict = {}
                    seriesDict['name'] = str(trt)
                    seriesDict['data'] = dataList
                    seriesList.append(seriesDict)
            else:
                dataList = pcoaDF.icol([CAP1, CAP2]).values.astype(float).tolist()
                trt = fieldListQuant[0]
                seriesDict = {}
                seriesDict['name'] = str(trt)
                seriesDict['data'] = dataList
                seriesList.append(seriesDict)

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
                result = result + 'perMANOVA results:' + '\n'
            if test == 2:
                result = result + 'betaDisper results:' + '\n'

            if len(fieldListCat) == 0:
                result = result + 'test cannot be run...' + '\n'
            else:
                result = result + str(bigf) + '\n'
            result += '===============================================\n'

            if len(fieldListCat) > 0:
                result = result + '\nenvfit results:\n'
                result = result + str(envFit)
            result += '===============================================\n'

            result = result + '\nEigenvalues\n'
            eigStr = eigDF.to_string()
            result = result + str(eigStr) + '\n'
            result += '===============================================\n\n\n\n'

            finalDict['text'] = result

            base[RID] = 'Step 6 of 8: Formatting graph data for display...done!'
            base[RID] = 'Step 7 of 8: Formatting PCoA table...'

            pcoaDF.reset_index(drop=True, inplace=True)
            res_table = pcoaDF.to_html(classes="table display")
            res_table = res_table.replace('border="1"', 'border="0"')
            finalDict['res_table'] = str(res_table)

            base[RID] = 'Step 7 of 8: Formatting PCoA table...done!'
            base[RID] = 'Step 8 of 8: Formatting distance score table...'

            dist_table = distDF.to_html(classes="table display")
            dist_table = dist_table.replace('border="1"', 'border="0"')
            finalDict['dist_table'] = str(dist_table)

            base[RID] = 'Step 8 of 8: Formatting distance score table...done!'
            finalDict['error'] = 'none'

            res = simplejson.dumps(finalDict)
            return HttpResponse(res, content_type='application/json')

    except Exception as e:
        print "Error with PCoA CAT: ", e
        state = "Error with PCoA CAT: " + str(e)

        myDict = {}
        myDict['error'] = state
        res = simplejson.dumps(myDict)
        return HttpResponse(res, content_type='application/json')
