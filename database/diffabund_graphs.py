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


stage = ''


def updateDiffAbund(request):
    global stage
    if request.is_ajax():
        myDict = {}
        myDict['stage'] = stage
        json_data = simplejson.dumps(myDict, encoding="Latin-1")
        return HttpResponse(json_data, content_type='application/json')


def getDiffAbund(request):
    global stage
    stage = 'Step 1 of 4: Querying database...'
    # Get selected samples from cookie and query database for sample info
    samples = Sample.objects.all()
    samples.query = pickle.loads(request.session['selected_samples'])
    selected = samples.values_list('sampleid')
    qs1 = Sample.objects.all().filter(sampleid__in=selected)

    if request.is_ajax():
        # Get variables from web page
        allJson = request.GET["all"]
        all = simplejson.loads(allJson)

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
            if total['count__sum'] is not None:
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

        # Create unique list of samples in meta dataframe (may be different than selected samples)
        myList = metaDF['sampleid'].tolist()
        mySet = list(ordered_set(myList))

        taxaDF = taxaProfileDF(mySet)

        stage = 'Step 1 of 4: Querying database...complete'
        stage = 'Step 2 of 4: Normalizing data...'

        # Create combined metadata column
        if len(fieldList) > 1:
            metaDF['merge'] = reduce(lambda x, y: metaDF[x] + ' & ' + metaDF[y], fieldList)
        else:
            metaDF['merge'] = metaDF[fieldList[0]]

        # Sum by taxa level
        taxaDF = taxaDF.groupby(level=taxaLevel).sum()

        # Normalization
        finalDict = {}
        r = R(RCMD="R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
        r.assign("metaDF", metaDF)
        r("trt <- factor(metaDF$merge)")
        r.assign("count", taxaDF)

        #Create filter for samples below threshold
        r.assign("size", size)
        r("total <- colSums(count)")
        r("col <- (total > size)")

        DESeq_error = ''
        if NormMeth == 1:
            r("countFilt <- count[,col]")
            r("trtFilt <- trt[col]")

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
            r("countFilt <- count[row,col]")
            r("trtFilt <- trt[col]")

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
            result += 'Analysis was run without normalization...\n'
            result += 'To try again, please increase the minimum sample size...\n'
        elif NormMeth == 2 and DESeq_error == 'no':
            result += 'Data were normalized by DESeq+Filter...\n'
        elif NormMeth == 2 and DESeq_error == 'yes':
            result += 'DESeq2 cannot run estimateSizeFactors...\n'
            result += 'Analysis was run without normalization...\n'
            result += 'To try again, please increase the minimum sample size...\n'
        result += '===============================================\n\n\n'

        stage = 'Step 2 of 4: Normalizing data...complete'
        stage = 'Step 3 of 4: Performing statistical test...'

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
                    r("df <- data.frame(id=rownames(res), baseMean=res$baseMean, log2FoldChange=res$log2FoldChange, stderr=res$lfcSE, stat=res$stat, pval=res$pvalue, padj=res$padj)")
                    nbinom_res = r.get("df")

                    names = []
                    for item in nbinom_res["id"]:
                        names.append(findTaxa(item))

                    nbinom_res['Taxa Name'] = names
                    nbinom_res.rename(columns={'id': 'Taxa ID'}, inplace=True)
                    stuff = ['Taxa ID', 'Taxa Name', ' baseMean ', ' log2FoldChange ', ' stderr ', ' stat ', ' pval ', ' padj ']
                    nbinom_res = nbinom_res.reindex(columns=stuff)
                    iterationName = str(mergeSet[i]) + ' vs ' + str(mergeSet[j])
                    nbinom_res['Comparison'] = iterationName

                    nbinom_res.rename(columns={' baseMean ': 'baseMean'}, inplace=True)
                    nbinom_res.rename(columns={' log2FoldChange ': 'log2FoldChange'}, inplace=True)
                    nbinom_res.rename(columns={' stderr ': 'StdErr'}, inplace=True)
                    nbinom_res.rename(columns={' stat ': 'Stat'}, inplace=True)
                    nbinom_res.rename(columns={' pval ': 'p-value'}, inplace=True)
                    nbinom_res.rename(columns={' padj ': 'p-adjusted'}, inplace=True)
                    nbinom_res[['p-value', 'p-adjusted']].astype(float)

                    if finalDF .empty:
                        finalDF = nbinom_res
                    else:
                        finalDF = pd.concat([finalDF, nbinom_res])

        stage = 'Step 3 of 4: Performing statistical test...completed'

        stage = 'Step 4 of 4: Preparing graph data...'
        seriesList = []
        xAxisDict = {}
        yAxisDict = {}

        grouped = finalDF.groupby('Comparison')
        for name, group in grouped:
            nosigDF = group[group["p-adjusted"] > FdrVal]
            allData = nosigDF[["baseMean", "log2FoldChange"]].values.astype(np.float).tolist()
            sigDF = group[group["p-adjusted"] <= FdrVal]
            sigData = sigDF[["baseMean", "log2FoldChange"]].values.astype(np.float).tolist()

            seriesDict = {}
            seriesDict['name'] = "NotSig: " + str(name)
            seriesDict['data'] = allData
            seriesList.append(seriesDict)

            seriesDict = {}
            seriesDict['name'] = "Sig: " + str(name)
            seriesDict['data'] = sigData
            seriesList.append(seriesDict)

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

        finalDF = finalDF[['Comparison', 'Taxa ID', 'Taxa Name', 'baseMean', 'log2FoldChange', 'StdErr', 'Stat', 'p-value', 'p-adjusted']]
        res_table = finalDF.to_html(classes="table display")
        res_table = res_table.replace('border="1"', 'border="0"')
        finalDict['res_table'] = str(res_table)
        finalDict['text'] = result

        res = simplejson.dumps(finalDict)
        stage = 'Step 4 of 4: Preparing graph data...completed'

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