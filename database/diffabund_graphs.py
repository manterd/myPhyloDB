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
        DESeq_error = ''
        r = R(RCMD="R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
        r.assign("metaDF", metaDF)
        r("condition <- factor(metaDF$merge)")
        r("library(DESeq)")

        if NormMeth == 1:
            r.assign("countTable", taxaDF)
            r("cds <- newCountDataSet(countTable, condition)")
            sizeFactor = []
            for i in mySet:
                sizeFactor.append(1)
            r.assign("sizeFactor", sizeFactor)
            r("cds$sizeFactor <- sizeFactor")
            r("cds <- estimateDispersions(cds, method='blind', fitType='local')")

        #FIXME proportion data not working (remove from options?)
        #NormMeth 1 & 3 do work
        elif NormMeth == 2:
            sizeFactor = []
            newDF = taxaDF.drop(mySet, axis=1)
            for i in mySet:
                newDF[i] = taxaDF[i].div(taxaDF[i].sum(), axis=0)
                sizeFactor.append(1)
            r.assign("relabundTable", newDF)
            r("cds <- newCountDataSet(relabundTable, condition)")
            r.assign("sizeFactor", sizeFactor)
            r("cds$sizeFactor <- sizeFactor")
            print r("counts(cds)")
            #r("cds <- estimateDispersions(cds, method='blind', fitType='local')")

        elif NormMeth == 3:
            r.assign("countTable", taxaDF)
            r("cds <- newCountDataSet(countTable, condition)")

            #FIXME size factors will fail if too many counts in dataset (INT16 problem in R)
            r("cds <- estimateSizeFactors(cds)")
            pycds = r.get("sizeFactors(cds)")

            found = 0
            for thing in pycds:
                if str(thing) == "None":
                    found += 1

            if found == 0:
                DESeq_error = 'no'
                r("cds <- estimateDispersions(cds, method='blind', fitType='local')")

            else:
                DESeq_error = 'yes'
                sizeFactor = []
                for i in mySet:
                    sizeFactor.append(1)
                r.assign("sizeFactor", sizeFactor)
                r("cds$sizeFactor <- sizeFactor")
                r("cds <- estimateDispersions(cds, method='blind', fitType='local')")

        #elif NormMeth == 4:
            #TODO add DESeq with independent filtering

        if NormMeth == 1:
            result += 'No normalization was performed...\n'
        elif NormMeth == 2:
            result += 'Data were normalized by the total number of sequence reads...\n'
        elif NormMeth == 3 and DESeq_error == 'no':
            result += 'Data were normalized by DESeq variance stabilization ...\n'
        elif NormMeth == 3 and DESeq_error == 'yes':
            result += 'DESeq cannot run estimateSizeFactors...\n'
            result += 'Analysis was run with size factors set to 1)...\n'
            result += 'To try again, please select fewer samples or another normalization method...\n'
        result += '===============================================\n\n\n'

        stage = 'Step 2 of 4: Normalizing data...complete'
        stage = 'Step 3 of 4: Performing statistical test...'

        # nbinomTest
        if StatTest == 1:
            mergeList = metaDF['merge'].tolist()
            mergeSet = list(set(mergeList))

            #if len(mergeSet) == 1:
                #TODO need to send an error back (No statistical test performed only 1 level in data)

            if len(mergeSet) == 2:
                r.assign("trt1", mergeSet[0])
                r.assign("trt2", mergeSet[1])
                r("res <- nbinomTest(cds, trt1, trt2)")
                print r("res")

                # TODO add res to result string

            #if len(mergeSet) == 3:
                #TODO need to add loop to deal with each pairwise comparison
                #TODO or potentially switch to fitNbinomGLMs

        #elif StatTest == 2:
        #TODO add multifactor (fitNbinomGLMs)

        stage = 'Step 3 of 4: Performing statistical test...completed'
        stage = 'Step 4 of 4: Preparing graph data...'

        #TODO finish result string and add foldchange graph
        #TODO add significant only capability

        stage = 'Step 4 of 4: Preparing graph data...completed'
        finalDict['text'] = result
        res = simplejson.dumps(finalDict)
        return HttpResponse(res, content_type='application/json')
