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
        FDR = float(all["FdrVal"])
        StatTest = int(all["StatTest"])
        sig_only = int(all["sig_only"])

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
        metaDF.sort(columns='sample_name', inplace=True)

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
        r("condition <- factor(metaDF$merge)")
        r("library(DESeq)")
        r.assign("countTable", taxaDF)
        r("cds <- newCountDataSet(countTable, condition)")
        r("cds <- estimateSizeFactors(cds)")
        r("cds <- estimateDispersions(cds, method='blind', fitType='local')")

        pycds = r.get("sizeFactors(cds)")
        found = 0
        for thing in pycds:
            if str(thing) == "None":
                found += 1

        DESeq_error = 'no'
        if found > 0:
            DESeq_error = 'yes'
            sizeFactor = []
            for i in myList:
                sizeFactor.append(1)
            r.assign("sizeFactor", sizeFactor)
            r("cds$sizeFactor <- sizeFactor")
            r("cds <- estimateDispersions(cds, method='blind', fitType='local')")

        if NormMeth == 1 and DESeq_error == 'no':
            result += 'Data were normalized by DESeq...\n'
        elif NormMeth == 1 and DESeq_error == 'yes':
            result += 'DESeq cannot run estimateSizeFactors...\n'
            result += 'Analysis was run without normalization...\n'
            result += 'To try again, please select fewer samples or another normalization method...\n'
        result += '===============================================\n\n\n'

        stage = 'Step 2 of 4: Normalizing data...complete'
        stage = 'Step 3 of 4: Performing statistical test...'

        # nbinomTest
        if StatTest == 1:
            mergeList = metaDF['merge'].tolist()
            mergeSet = list(set(mergeList))
            for i, val in enumerate(mergeSet):
                start = i + 1
                stop = int(len(mergeSet))
                for j in range(start, stop):
                    if i != j:
                        result += '===============================================\n'
                        result = result + 'Comparison ' + str(mergeSet[i]) + ' vs ' + str(mergeSet[j]) + '\n'
                        r.assign("trt1", mergeSet[i])
                        r.assign("trt2", mergeSet[j])
                        r("res <- nbinomTest(cds, trt1, trt2)")

                        nbinom_res = r.get("res")
                        result += str(nbinom_res)
                result += '===============================================\n\n\n'

                # TODO modify output sent to result string
                    # write nbinom_res to datafame
                    # add taxa name and rank columns
                    # remove baseMean
                    # output to string

                # TODO output means to datatable
                    # write



        stage = 'Step 3 of 4: Performing statistical test...completed'
        stage = 'Step 4 of 4: Preparing graph data...'

        # TODO add fold change graph
            # create highcharts scatter plot
            # x = taxa name
            # y = log2 fold change
            # significant points in red
            # OR see if we can a

        stage = 'Step 4 of 4: Preparing graph data...completed'
        finalDict['text'] = result
        res = simplejson.dumps(finalDict)
        return HttpResponse(res, content_type='application/json')
