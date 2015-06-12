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
        r("cds <- estimateDispersions(cds, method='blind', fitType='local', sharingMode='fit-only')")

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
            r("cds <- estimateDispersions(cds, method='blind', fitType='local', sharingMode='fit-only')")

        if NormMeth == 1 and DESeq_error == 'no':
            result += 'Data were normalized by DESeq...\n'
        elif NormMeth == 1 and DESeq_error == 'yes':
            result += 'DESeq cannot run estimateSizeFactors...\n'
            result += 'Analysis was run without normalization...\n'
            result += 'To try again, please select fewer samples or another normalization method...\n'
        result += '===============================================\n\n\n'

        stage = 'Step 2 of 4: Normalizing data...complete'
        stage = 'Step 3 of 4: Performing statistical test...'
        AxisArray = {}
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

                        names = []
                        for item in nbinom_res["id"]:
                            names.append(findTaxa(item))

                        try:
                            nbinom_res['Taxa name'] = names
                            nbinom_res.rename(columns={'id': 'Taxa ID'}, inplace=True)
                            stuff = ['Taxa ID', 'Taxa name', ' baseMean ', ' baseMeanA ', ' baseMeanB ', ' foldChange ', ' log2FoldChange ', ' pval ', ' padj ']
                            nbinom_res = nbinom_res.reindex(columns=stuff)
                            nbinom_res.rename(columns={' log2FoldChange ': 'log2FoldChange'}, inplace=True)
                            nbinom_res.rename(columns={' baseMeanA ': 'baseMeanA'}, inplace=True)
                            nbinom_res.rename(columns={' baseMeanB ': 'baseMeanB'}, inplace=True)
                            nbinom_res.rename(columns={' pval ': 'pval'}, inplace=True)
                            nbinom_res.rename(columns={' padj ': 'padj'}, inplace=True)
                        except:
                            print ("Join failed")

                        try:
                            iterationName = str(mergeSet[i]) + ' vs ' + str(mergeSet[j])  # + taxa
                            print ("Iteration: "+iterationName)
                            dataSet = [None]
                            print "DataSet: "+str(dataSet)
                            for thing in nbinom_res["log2FoldChange"]:
                                print("Element!")
                                print("thing: "+str(thing))
                                dataSet.append(thing)
                            AxisArray[iterationName] = dataSet
                        except:
                            print("Failed to add xAxis data")

                        try:
                            nbinom_res.drop(' baseMean ', axis=1, inplace=True)

                        except:
                            print"Deletion failed"
                            print sys.exc_info()[0]

                        result += nbinom_res.to_string()
                        result += '\n===============================================\n\n\n'

                # TODO output means to datatable ???
                    # write



        stage = 'Step 3 of 4: Performing statistical test...completed'
        stage = 'Step 4 of 4: Preparing graph data...'

        # TODO add fold change graph
            # create highcharts scatter plot
            # x = taxa name
            # y = log2 fold change
            # significant points in red

        '''regrDict = {}
        regrDict['type'] = 'line'
        regrDict['name'] = 'R2: ' + str(r_square) + '; p-value: ' + str(p_value) + '<br>' + '(y = ' + str(slope) + 'x' + ' + ' + str(intercept) + ')'
        regrDict['data'] = regrList
        seriesList.append(regrDict)'''

        xAxisDict = {}
        yAxisDict = {}

        xTitle = {}
        xTitle['text'] = taxaDF  # taxalist here?
        xAxisDict['title'] = xTitle

        yTitle = {}
        yTitle['text'] = fieldList[0]
        yAxisDict['title'] = yTitle

        finalDict['series'] = seriesList
        finalDict['xAxis'] = xAxisDict
        finalDict['yAxis'] = yAxisDict

        stage = 'Step 4 of 4: Preparing graph data...completed'
        finalDict['text'] = result
        res = simplejson.dumps(finalDict)
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

    '''qs1 = Species.objects.values('kingdomid__kingdomName', 'kingdomid', 'phylaid__phylaName', 'phylaid', 'classid__className', 'classid', 'orderid__orderName', 'orderid', 'familyid__familyName', 'familyid', 'genusid__genusName', 'genusid', 'speciesName', 'speciesid')
    speciesDF = pd.DataFrame.from_records(qs1)
    speciesDF.rename(columns={'kingdomid__kingdomName': 'Kingdom Name', 'kingdomid': 'Kingdom ID'}, inplace=True)
    speciesDF.rename(columns={'phylaid__phylaName': 'Phylum Name', 'phylaid': 'Phylum ID'}, inplace=True)
    speciesDF.rename(columns={'classid__className': 'Class Name', 'classid': 'Class ID'}, inplace=True)
    speciesDF.rename(columns={'orderid__orderName': 'Order Name', 'orderid': 'Order ID'}, inplace=True)
    speciesDF.rename(columns={'familyid__familyName': 'Family Name', 'familyid': 'Family ID'}, inplace=True)
    speciesDF.rename(columns={'genusid__genusName': 'Genus Name', 'genusid': 'Genus ID'}, inplace=True)
    speciesDF.rename(columns={'speciesName': 'Species Name', 'speciesid': 'Species ID'}, inplace=True)
    table = speciesDF.to_html(classes="table display", columns=['Kingdom Name', 'Kingdom ID', 'Phylum Name', 'Phylum ID', 'Class Name', 'Class ID', 'Order Name', 'Order ID', 'Family Name', 'Family ID', 'Genus Name', 'Genus ID', 'Species Name', 'Species ID'])
    table = table.replace('border="1"', 'border="0"')'''
