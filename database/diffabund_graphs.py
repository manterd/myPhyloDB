from diffabund_DF import catDiffAbundDF
from django.http import HttpResponse
from database.models import Sample, Profile
from django.db.models import Sum
import pandas as pd
import pickle
from scipy import stats
import simplejson
from database.utils import multidict, ordered_set, taxaProfileDF
from numpy import *
import numpy as np
import datetime
from pyper import *
import math
from models import Kingdom, Phyla, Class, Order, Family, Genus, Species


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}


def updateDiffAbund(request):
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


def removeRIDDIFF(request):
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


def getDiffAbund(request):
    try:
        global base, time1, TimeDiff
        # Get selected samples from cookie and query database for sample info
        samples = Sample.objects.all()
        samples.query = pickle.loads(request.session['selected_samples'])
        selected = samples.values_list('sampleid')
        qs1 = Sample.objects.all().filter(sampleid__in=selected)

        if request.is_ajax():
            # Get variables from web page
            allJson = request.GET["all"]
            all = simplejson.loads(allJson)

            RID = str(all["RID"])
            time1[RID] = time.time()
            base[RID] = 'Step 1 of 6: Querying database...'

            taxaLevel = int(all["taxaLevel"])
            NormMeth = int(all["NormMeth"])
            FdrVal = float(all["FdrVal"])
            size = int(all["NormVal"])
            theta = float(all["Theta"])

            # Generate a list of sequence reads per sample
            countList = []
            for sample in qs1:
                total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                if total['count__sum'] is not None:
                    countList.append(total['count__sum'])

            # Remove blank samples if below the sequence threshold set by user (rarefaction)
            newList = []
            result = 'Data Normalization:\n'
            for sample in qs1:
                total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                if total['count__sum'] is not None and total['count__sum'] >= size:
                    id = sample.sampleid
                    newList.append(id)
            qs2 = Sample.objects.all().filter(sampleid__in=newList)

            # Get dict of selected meta variables
            metaString = all["meta"]
            metaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaString)

            # Convert dict to list
            fieldList = []
            for key in metaDict:
                fieldList.append(key)
            metaDF = catDiffAbundDF(qs2, metaDict)
            metaDF.dropna(subset=fieldList, inplace=True)
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

            base[RID] = 'Step 1 of 6: Querying database...done!'
            base[RID] = 'Step 2 of 6: Normalizing data...'

            # Create combined metadata column
            if len(fieldList) > 1:
                for index, row in metaDF.iterrows():
                    metaDF.ix[index, 'merge'] = " & ".join(row[fieldList])
            else:
                metaDF['merge'] = metaDF[fieldList[0]]

            # Sum by taxa level
            taxaDF = taxaDF.groupby(level=taxaLevel).sum()

            # Normalization
            finalDict = {}
            if os.name == 'nt':
                r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
            else:
                r = R(RCMD="R/R-Linux/bin/R")
            r.assign("metaDF", metaDF)
            r("trt <- factor(metaDF$merge)")
            r.assign("count", taxaDF)

            DESeq_error = ''
            if NormMeth == 1:
                r("countFilt <- count")
                r("trtFilt <- trt")

                r("library(DESeq2)")
                r("colData <- data.frame(row.names=colnames(countFilt), trt=trtFilt)")
                r("dds <- DESeqDataSetFromMatrix(countData=countFilt, colData=colData, design= ~ trt)")

                r("dds <- estimateSizeFactors(dds)")
                pycds = r.get("sizeFactors(dds)")

                if pycds is not None:
                    DESeq_error = 'no'
                    r("dds <- estimateDispersions(dds)")
                    r("dds <- nbinomWaldTest(dds)")

                elif pycds is None:
                    DESeq_error = 'yes'
                    r("sizeFactor <- rep(1, length(trtFilt))")
                    r("dds$sizeFactor <- sizeFactor")
                    r("dds <- estimateDispersions(dds)")
                    r("dds <- nbinomWaldTest(dds)")

            elif NormMeth == 2:
                r("rs = rowSums(count)")
                r.assign("theta", theta)
                r("row <- (rs > quantile(rs, probs=theta))")
                r("countFilt <- count[row,]")
                r("trtFilt <- trt")

                r("library(DESeq2)")
                r("colData <- data.frame(row.names=colnames(countFilt), trt=trtFilt)")
                r("dds <- DESeqDataSetFromMatrix(countData=countFilt, colData=colData, design= ~ trt)")

                r("dds <- estimateSizeFactors(dds)")
                pycds = r.get("sizeFactors(dds)")

                if pycds is not None:
                    DESeq_error = 'no'
                    r("dds <- estimateDispersions(dds)")
                    r("dds <- nbinomWaldTest(dds)")

                elif pycds is None:
                    DESeq_error = 'yes'
                    r("sizeFactor <- rep(1, length(trtFilt))")
                    r("dds$sizeFactor <- sizeFactor")
                    r("dds <- estimateDispersions(dds)")
                    r("dds <- nbinomWaldTest(dds)")

            if NormMeth == 1 and DESeq_error == 'no':
                result += 'Data were normalized by DESeq2...\n'
            elif NormMeth == 1 and DESeq_error == 'yes':
                result += 'DESeq2 cannot run estimateSizeFactors...\n'
                result += 'Analysis was run without size normalization...\n'
                result += 'To try again, select a different sample combination or increase the minimum sample size...\n'
            elif NormMeth == 2 and DESeq_error == 'no':
                result += 'Data were normalized by DESeq2+Filter...\n'
            elif NormMeth == 2 and DESeq_error == 'yes':
                result += 'DESeq2 cannot run estimateSizeFactors...\n'
                result += 'Analysis was run without size normalization...\n'
                result += 'To try again, select a different sample combination or  increase the minimum sample size...\n'
            result += '===============================================\n\n\n'

            base[RID] = 'Step 2 of 6: Normalizing data...done!'
            base[RID] = 'Step 3 of 6: Performing statistical test...'
            try:
                mergeList = metaDF['merge'].tolist()
                mergeSet = list(set(mergeList))

                finalDF = pd.DataFrame()
                for i, val in enumerate(mergeSet):
                    start = i + 1
                    stop = int(len(mergeSet))
                    for j in range(start, stop):
                        if i != j:
                            r.assign("trt1", mergeSet[i])
                            r.assign("trt2", mergeSet[j])
                            r("res <- results(dds, contrast=c('trt', trt1, trt2))")
                            r("baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,dds$trt==trt1, drop=FALSE])")
                            r("baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,dds$trt==trt2, drop=FALSE])")
                            r("df <- data.frame(id=rownames(res), baseMean=res$baseMean, baseMeanA=baseMeanA, baseMeanB=baseMeanB, log2FoldChange=res$log2FoldChange, stderr=res$lfcSE, stat=res$stat, pval=res$pvalue, padj=res$padj)")
                            nbinom_res = r.get("df")

                            names = []
                            try:
                                for item in nbinom_res["id"]:
                                    names.append(findTaxa(item))
                            except Exception as e:
                                print("Failed at loop, "+str(e.args))

                            nbinom_res['Taxa Name'] = names
                            nbinom_res.rename(columns={'id': 'Taxa ID'}, inplace=True)
                            stuff = ['Taxa ID', 'Taxa Name', ' baseMean ', ' baseMeanA ', ' baseMeanB ', ' log2FoldChange ', ' stderr ', ' stat ', ' pval ', ' padj ']
                            nbinom_res = nbinom_res.reindex(columns=stuff)

                            iterationName = str(mergeSet[i]) + ' vs ' + str(mergeSet[j])

                            nbinom_res['Comparison'] = iterationName


                            nbinom_res.rename(columns={' baseMean ': 'baseMean'}, inplace=True)

                            nbinom_res.rename(columns={' baseMeanA ': 'baseMeanA'}, inplace=True)

                            nbinom_res.rename(columns={' baseMeanB ': 'baseMeanB'}, inplace=True)

                            nbinom_res.rename(columns={' log2FoldChange ': 'log2FoldChange'}, inplace=True)

                            nbinom_res.rename(columns={' stderr ': 'StdErr'}, inplace=True)

                            nbinom_res.rename(columns={' stat ': 'Stat'}, inplace=True)

                            nbinom_res.rename(columns={' pval ': 'p-value'}, inplace=True)

                            nbinom_res.rename(columns={' padj ': 'p-adjusted'}, inplace=True)

                            nbinom_res[['p-value', 'p-adjusted']].astype(float)


                            finalDF = pd.concat([finalDF, nbinom_res])

                            base[RID] = 'Step 3 of 6: Performing statistical test...' + str(iterationName) + ' is done!'


            except Exception as e:
                print ("Exception! "+str(e.args))

            base[RID] = 'Step 3 of 6: Performing statistical test...done!'
            base[RID] = 'Step 4 of 6: Formatting graph data for display...'

            seriesList = []
            xAxisDict = {}
            yAxisDict = {}

            grouped = finalDF.groupby('Comparison')

            listOfShapes = ['circle', 'square', 'triangle', 'triangle-down', 'diamond',]
            shapeIterator = 0

            for name, group in grouped:
                nosigDF = group[group["p-adjusted"] > FdrVal]
                allData = nosigDF[["baseMean", "log2FoldChange"]].values.astype(np.float).tolist()
                sigDF = group[group["p-adjusted"] <= FdrVal]
                sigData = sigDF[["baseMean", "log2FoldChange"]].values.astype(np.float).tolist()

                seriesDict = {}
                seriesDict['name'] = "NotSig: " + str(name)
                seriesDict['data'] = allData
                markerDict = {}
                markerDict['symbol'] = listOfShapes[shapeIterator]
                seriesDict['marker'] = markerDict
                seriesList.append(seriesDict)

                seriesDict = {}
                seriesDict['name'] = "Sig: " + str(name)
                seriesDict['data'] = sigData
                markerDict = {}
                markerDict['symbol'] = listOfShapes[shapeIterator]
                seriesDict['marker'] = markerDict
                seriesList.append(seriesDict)

                shapeIterator += 1
                if shapeIterator >= len(listOfShapes):
                    shapeIterator = 0

                base[RID] = 'Step 4 of 6: Formatting graph data for display...' + str(name) + ' is done!'

            xTitle = {}
            xTitle['text'] = "baseMean"
            xAxisDict['title'] = xTitle
            xAxisDict['type'] = 'logarithmic'

            yTitle = {}
            yTitle['text'] = "log2FoldChange"
            yAxisDict['title'] = yTitle
            yAxisDict['type'] = 'linear'

            finalDict['series'] = seriesList
            finalDict['xAxis'] = xAxisDict
            finalDict['yAxis'] = yAxisDict

            base[RID] = 'Step 4 of 6: Formatting graph data for display...done!'
            base[RID] = 'Step 5 of 6:  Formatting nbinomTest results for display...'

            finalDF = finalDF[['Comparison', 'Taxa ID', 'Taxa Name', 'baseMean', 'baseMeanA', 'baseMeanB', 'log2FoldChange', 'StdErr', 'Stat', 'p-value', 'p-adjusted']]
            res_table = finalDF.to_html(classes="table display")
            res_table = res_table.replace('border="1"', 'border="0"')
            finalDict['res_table'] = str(res_table)
            finalDict['text'] = result

            base[RID] = 'Step 6 of 6: Formatting results for display...done!'
            finalDict['error'] = 'none'

            res = simplejson.dumps(finalDict)

            return HttpResponse(res, content_type='application/json')

    except Exception as e:
        print "Error with DIFFABUND: ", e
        print "RID: ", RID
        state = "Error with DIFFABUND: " + str(e)

        myDict = {}
        myDict['error'] = state
        res = simplejson.dumps(myDict)
        return HttpResponse(res, content_type='application/json')


def findTaxa(id):
    taxa = ""
    try:
        temp = Kingdom.objects.filter(kingdomid=id)
        taxa += temp[0].kingdomName
    except:
        # not kingdom, try next one
        try:
            temp = Phyla.objects.filter(phylaid=id)
            taxa += temp[0].phylaName
        except:
            try:
                temp = Class.objects.filter(classid=id)
                taxa += temp[0].className
            except:
                try:
                    temp = Order.objects.filter(orderid=id)
                    taxa += temp[0].orderName
                except:
                    try:
                        temp = Family.objects.filter(familyid=id)
                        taxa += temp[0].familyName
                    except:
                        try:
                            temp = Genus.objects.filter(genusid=id)
                            taxa += temp[0].genusName
                        except:
                            try:
                                temp = Species.objects.filter(speciesid=id)
                                taxa += temp[0].speciesName
                            except:
                                # not found, error!
                                print("Could not find taxa for "+str(id))
    return taxa