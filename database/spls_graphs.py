from spls_DF import catSPLSMetaDF, normalizeSPLS
from django.http import HttpResponse
from database.models import Sample, Profile
from django.db.models import Sum
import numpy as np
from numpy import *
import pandas as pd
import pickle
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


def getCatSPLSAData(request):
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

            metaStr = all["metaQuant"]
            fieldList = []
            valueList = []
            metaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStr)
            for key in metaDict:
                fieldList.append(key)
                valueList.append(metaDict[key])

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
            metaDF = catSPLSMetaDF(qs2, metaDict)
            print 'metaDF\n', metaDF

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

            # Create unique list of samples in meta dataframe (may be different than selected samples)
            myList = metaDF['sampleid'].tolist()
            mySet = list(ordered_set(myList))

            taxaDF = taxaProfileDF(mySet)

            # Sum by taxa level
            taxaDF = taxaDF.groupby(level=taxaLevel).sum()

            base[RID] = 'Step 1 of 8: Querying database...done!'
            base[RID] = 'Step 2 of 8: Normalizing data...'

            normDF, DESeq_error = normalizeSPLS(taxaDF, taxaLevel, mySet, NormMeth, NormReads, metaDF)
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
            base[RID] = 'Step 3 of 8: Calculating sPLS...'

            metaDF.sort_index(inplace=True)

            if os.name == 'nt':
                r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
            else:
                r = R(RCMD="R/R-Linux/bin/R")

            metaDF.set_index('sampleid', inplace=True)
            #finalDF = metaDF.join(normDF)

            r.assign("X", normDF)   ### taxa abundances
            r.assign("Y", metaDF[fieldList])    ### quantitative data

            #print r("X")
            #print r("Y")

            taxaIDs = normDF.columns.values
            r.assign("taxaIDs", taxaIDs)
            r("taxaIDs")
            r("row.names(X)")
            r("colnames(X) <- taxaIDs")

            # remove predictors with zero variance
            r("library(mixOmics)")
            r("ZeroVar <- nearZeroVar(X)")
            r("List <- row.names(ZeroVar$Metrics)")
            r("X_new <- X[,-which(names(X) %in% List)]")

            r("library(spls)")
            r("set.seed(1)")
            r("cv <- cv.spls(X_new, Y, eta=seq(0.1, 0.9, 0.1), K=c(1:10))")
            r("f <- spls(X_new, Y, eta=cv$eta.opt, K=cv$K.opt)")

            r("set.seed(1)")
            r("ci.f <- ci.spls(f)")
            r("cf <- correct.spls(ci.f, plot.it=FALSE)")

            ### DELETE 'Rplots.pdf'
            r("cf.rows <- row.names(cf)")
            cf = r.get("cf")
            rows = r.get("cf.rows")

            coeffsDF = pd.DataFrame(cf,  columns=[fieldList], index=rows)
           # print 'coeffsDF\n', coeffsDF        ## index = taxaids, cols = variables
            coeffsDF = coeffsDF[(coeffsDF.T != 0).any()]

            xAxisDict = {}
            catList = coeffsDF.index.values.tolist()
            xAxisDict['categories'] = catList
            labelsDict = {}
            labelsDict['rotation'] = 90
            xAxisDict['labels'] = labelsDict
            print 'xAxis:', xAxisDict

            yAxisDict = {}
            catList = fieldList
            yAxisDict['categories'] = catList
            print 'yAxis:', yAxisDict

            seriesList = []
            seriesDict = {}
            seriesDict['name'] = 'test'
            dataList = coeffsDF[fieldList[0]].values.astype(np.float).tolist()
            seriesDict['data'] = dataList
            seriesList.append(seriesDict)
            print 'seriesList:', seriesList

            finalDict['xAxis'] = xAxisDict
            finalDict['yAxis'] = yAxisDict
            finalDict['series'] = seriesList

            finalDict['text'] = result

            res = simplejson.dumps(finalDict)
            return HttpResponse(res, content_type='application/json')

    except Exception as e:
        print "Error with SPLS CAT: ", e
        state = "Error with SPLS CAT: " + str(e)

        myDict = {}
        myDict['error'] = state
        res = simplejson.dumps(myDict)
        return HttpResponse(res, content_type='application/json')
