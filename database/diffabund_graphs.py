from diffabund_DF import catDiffAbundDF
from django.http import HttpResponse
from database.models import Sample, Profile
from django.db.models import Sum
import pandas as pd
import pickle
from scipy import stats
import simplejson
from database.utils import multidict, ordered_set, taxaProfileDF
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

        if NormMeth == 1:
            result = result + 'No normalization was performed...\n'
        elif NormMeth == 2:
            result = result + 'Data were normalized by the total number of sequence reads...\n'
        elif NormMeth == 3:
            result = result + 'Data were normalized by DESeqVS...\n'

        result = result + '===============================================\n\n\n'

        #countDF = pd.DataFrame()
        #if NormMeth == 1:
            #TODO add no normalization

        #if NormMeth == 2:
            #TODO add independent filtering

        if NormMeth == 3:
            # Create CountDataSet for DESeq
            r = R(RCMD="R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
            r.assign("countTable", taxaDF)
            r.assign("metaDF", metaDF)
            r("condition <- factor(metaDF$merge)")
            r("library(DESeq)")
            r("cds <- newCountDataSet(countTable, condition)")
            r("cds <- estimateSizeFactors(cds)")
            pycds = r.get("sizeFactors(cds)")
            r("cds <- estimateDispersions(cds, method='blind', fitType='local')")
            r("vsd <- varianceStabilizingTransformation(cds)")
            print r("varianceStabilizingTransformation(cds)")

            ### Create a heatmap of genes vs sample.
            r("pdf('test.pdf')")
            r("library(RColorBrewer)")
            r("library(gplots)")
            r("hmcol=colorRampPalette(brewer.pal(9,'GnBu'))(100)")
            r("heatmap.2(exprs(vsd), col=hmcol, trace='none', margin=c(15,15), cexRow=0.8, cexCol=0.8)")

            found = False
            for thing in pycds:
                if str(thing) == "None":
                    found = True

            stage = 'Step 2 of 4: Normalizing data...complete'
            stage = 'Step 3 of 4: Performing statistical test...'

            # nbinomTest
            mergeList = metaDF['merge'].tolist()
            mergeSet = list(set(mergeList))

            #if len(mergeSet) == 1:
                #TODO

            if len(mergeSet) == 2:
                r.assign("trt1", mergeSet[0])
                r.assign("trt2", mergeSet[1])

                if found:
                    DESeq_error = 'no'
                    r("res <- nbinomTest(cds, trt1, trt2)")
                    print "cds"
                    print r("res")
                else:
                    DESeq_error = 'yes'
                    r("res <- nbinomTest(exprs(vsd), trt1, trt2)")
                    print "vsd"
                    print r("res")

                # TODO output res to dataframe

            #if len(mergeSet) == 3:
                #TODO
                #not sure how to deal with this?
                #multiple nbinomTests?

            #TODO add multifactor (fitNbinomGLMs)

        #TODO finish result string and add foldchange graph
        #TODO add significant only capability

        result = ''
        return HttpResponse(result, content_type='application/json')
