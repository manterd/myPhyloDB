from diffabund_DF import catDiffAbundDF, normalizeDiffAbundDA
from django.http import HttpResponse
from database.models import Sample, Profile
from django.db.models import Sum
import pandas as pd
import pickle
from pyvttbl.Anova1way import Anova1way
from scipy import stats
import simplejson
from database.utils import multidict, ordered_set, taxaProfileDF
import numpy as np
import datetime
from pyper import *
import math


stage = ''


def statusDiffAbund(request):
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
        #taxaLevel = int(all["taxa"])    ### need to add selectbox to html
        allJson = request.GET["all"]
        all = simplejson.loads(allJson)
        selectAll = int(all["selectAll"])
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

        '''
        # Calculate min/median/max of sequence reads for rarefaction
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
        '''

        # Remove samples if below the sequence threshold set by user (rarefaction)
        #newList = []
        result = 'Data Normalization:\n'

        '''
        # Limit reads to max value
        if NormReads > maxSize:
            NormReads = medianSize
            result = result + 'The desired sample size was too high and automatically reset to the median value...\n'
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
            result = result + 'The desired sample size was too high and automatically reset to the median value...\n'
            for sample in qs1:
                total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                if total['count__sum'] is not None and int(total['count__sum']) >= NormReads:
                    id = sample.sampleid
                    newList.append(id)
        qs2 = Sample.objects.all().filter(sampleid__in=newList)
        '''

        # Get dict of selected meta variables
        metaString = all["meta"]
        metaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaString)

        # Convert dict to list
        fieldList = []
        for key in metaDict:
            fieldList.append(key)

        # Create dataframe of meta variables by sample
        metaDF = catDiffAbundDF(qs1, metaDict)
        metaDF.dropna(subset=fieldList, inplace=True)
        metaDF.sort(columns='sample_name', inplace=True)

        # Create unique list of samples in meta dataframe (may be different than selected samples due to null values)
        myList = metaDF['sampleid'].tolist()
        mySet = list(ordered_set(myList))

        # Create dataframe with all taxa/count data by sample
        taxaDF = taxaProfileDF(mySet)

        stage = 'Step 1 of 4: Querying database...complete'
        # Normalize data
        stage = 'Step 2 of 4: Normalizing data...'

        # Create combined metadata column
        # Only need this for DESeq
        if len(fieldList) > 1:
            metaDF['merge'] = reduce(lambda x, y: metaDF[x] + ' & ' + metaDF[y], fieldList)
        else:
            metaDF['merge'] = metaDF[fieldList[0]]

        # Sum by taxa level
        taxaLevel = 0
        taxaDF = taxaDF.groupby(level=taxaLevel).sum()
        print taxaDF

        normDF = normalizeDiffAbundDA(taxaDF, taxaLevel, mySet, NormMeth, metaDF)
        finalDF = metaDF.merge(normDF, on='sampleid', how='outer')
        finalDF[['abund', 'rich', 'diversity']] = finalDF[['abund', 'rich', 'diversity']].astype(float)
        pd.set_option('display.max_rows', finalDF.shape[0], 'display.max_columns', finalDF.shape[1], 'display.width', 1000)

        finalDict = {}
        if NormMeth == 1:
            result = result + 'No normalization was performed...\n'
        if NormMeth == 2:
            result = result + 'Data were normalized by the total number of sequence reads...\n'
        if NormMeth == 3 :
            result = result + 'Data were normalized by DESeq...\n'
        result = result + '===============================================\n\n\n'

        stage = 'Step 2 of 4: Normalizing data...complete'
        stage = 'Step 3 of 4: Performing statistical test...'

        seriesList = []
        xAxisDict = {}
        yAxisDict = {}

        ### group DataFrame by each taxa level selected
        grouped1 = finalDF.groupby(['rank', 'taxa_name', 'taxa_id'])

        ### group DataFrame by each meta variable selected
        for name1, group1 in grouped1:
            trtList = []
            valList = []
            grouped2 = group1.groupby(fieldList)['abund']

            for name2, group2 in grouped2:
                if isinstance(name2, unicode):
                    trt = name2
                else:
                    trt = ' & '.join(list(name2))
                trtList.append(trt)
                valList.append(list(group2.T))

            D = ""
            p_val = 1.0
            if StatTest == 1:
                # transform data to help equalize variances
                for row in valList:
                    if NormMeth == 4:
                        row[:] = [math.asin(x) for x in row]
                    else:
                        row[:] = [math.log(x + 1, 2) for x in row]

                D = Anova1way()
                try:
                    D.run(valList, conditions_list=trtList)
                    p_val = float(D['p'])
                except:
                    p_val = 1.0
                    D = 'Samples sizes are either too unequal or n=1.\n'

            elif StatTest == 2:
                rows_list = []
                m = 0.0
                for i, val1 in enumerate(trtList):
                    smp1 = valList[i]
                    start = i + 1
                    stop = int(len(trtList))
                    for j in range(start, stop):
                        if i != j:
                            val2 = trtList[j]
                            smp2 = valList[j]
                            if len(smp1) > 1 and len(smp2) > 1:
                                (t_val, p_val) = stats.ttest_ind(smp1, smp2, equal_var=False)
                                m += 1.0
                            else:
                                p_val = "Nan"
                            dict1 = {'taxa_level': name1[0], 'taxa_name': name1[1], 'taxa_id': name1[2], 'sample1': val1, 'sample2': val2, 'p_value': p_val}
                            rows_list.append(dict1)

                pvalDF = pd.DataFrame(rows_list)
                pvalDF.sort(columns='p_value', inplace=True)
                pvalDF.reset_index(drop=True, inplace=True)
                pvalDF.index += 1
                pvalDF['BH'] = pvalDF.index / m * FDR
                pvalDF['Sig'] = pvalDF.p_value <= pvalDF.BH

                pvalDF[['p_value', 'BH']] = pvalDF[['p_value', 'BH']].astype(float)
                D = pvalDF.to_string(columns=['sample1', 'sample2', 'p_value', 'BH', 'Sig'])
                if True in pvalDF['Sig'].values:
                    p_val = 0.0
                else:
                    p_val = 1.0

            stage = 'Step 3 of 4: Performing statistical test...complete'
            stage = 'Step 4 of 4: Preparing graph data...'
            if sig_only == 1:
                if p_val <= 0.05:
                    result = result + '===============================================\n'
                    result = result + 'Taxa level: ' + str(name1[0]) + '\n'
                    result = result + 'Taxa name: ' + str(name1[1]) + '\n'
                    result = result + 'Taxa ID: ' + str(name1[2]) + '\n'
                    result = result + 'Dependent Variable: Abundance' + '\n'

                    indVar = ' x '.join(fieldList)
                    result = result + 'Independent Variable: ' + str(indVar) + '\n\n'

                    result = result + str(D) + '\n'
                    result = result + '===============================================\n'
                    result = result + '\n\n\n\n'

                    dataList = []
                    grouped2 = group1.groupby(fieldList).mean()

                    dataList.extend(list(grouped2['abund'].T))

                    seriesDict = {}
                    seriesDict['name'] = name1
                    seriesDict['data'] = dataList
                    seriesList.append(seriesDict)

                    xTitle = {}
                    xTitle['text'] = indVar
                    xAxisDict['title'] = xTitle
                    xAxisDict['categories'] = trtList

                    yTitle = {}
                    yTitle['text'] = 'Abundance'
                    yAxisDict['title'] = yTitle

            if sig_only == 0:
                result = result + '===============================================\n'
                result = result + 'Taxa level: ' + str(name1[0]) + '\n'
                result = result + 'Taxa name: ' + str(name1[1]) + '\n'
                result = result + 'Taxa ID: ' + str(name1[2]) + '\n'
                result = result + 'Dependent Variable: Abundance' + '\n'

                indVar = ' x '.join(fieldList)
                result = result + 'Independent Variable: ' + str(indVar) + '\n\n'

                result = result + str(D) + '\n'
                result = result + '===============================================\n'
                result = result + '\n\n\n\n'

                dataList = []
                grouped2 = group1.groupby(fieldList).mean()

                dataList.extend(list(grouped2['abund'].T))

                seriesDict = {}
                seriesDict['name'] = name1
                seriesDict['data'] = dataList
                seriesList.append(seriesDict)

                xTitle = {}
                xTitle['text'] = indVar
                xAxisDict['title'] = xTitle
                xAxisDict['categories'] = trtList

                yTitle = {}
                yTitle['text'] = 'Abundance'
                yAxisDict['title'] = yTitle

        finalDict['series'] = seriesList
        finalDict['xAxis'] = xAxisDict
        finalDict['yAxis'] = yAxisDict
        finalDict['text'] = result
        if not seriesList:
            finalDict['empty'] = 0
        else:
            finalDict['empty'] = 1

        finalDF.reset_index(drop=True, inplace=True)

        biome = {}
        newList = ['sampleid', 'sample_name']
        newList.extend(fieldList)
        grouped = finalDF.groupby(newList, sort=False)
        nameList = []
        for name, group in grouped:
            metaDict = {}
            for i in xrange(1, len(newList)):
                metaDict[str(newList[i])] = str(name[i])
            nameList.append({"id": str(name[0]), "metadata": metaDict})

        grouped = finalDF.groupby(['rank', 'taxa_name', 'taxa_id'], sort=False)
        taxaList = []
        dataList = []
        for name, group in grouped:
            metaDict ={}
            taxonList = []
            taxonList.append(str(name[1]))
            metaDict['taxonomy'] = taxonList
            taxaList.append({"id": str(name[2]), "metadata": metaDict})
            dataList.append(group['abund'].tolist())

        biome['format'] = 'Biological Observation Matrix 0.9.1-dev'
        biome['format_url'] = 'http://biom-format.org/documentation/format_versions/biom-1.0.html'
        biome['type'] = 'OTU table'
        biome['generated_by'] = 'myPhyloDB'
        biome['date'] = str(datetime.datetime.now())
        biome['matrix_type'] = 'dense'
        biome['matrix_element_type'] = 'float'
        biome['rows'] = taxaList
        biome['columns'] = nameList
        biome['data'] = dataList

        res_table = finalDF.to_html(classes="table display")
        res_table = res_table.replace('border="1"', 'border="0"')
        finalDict['res_table'] = str(res_table)
        stage = 'Step 4 of 4: Preparing graph data...complete'

        biome_json = simplejson.dumps(biome, ensure_ascii=True, indent=4, sort_keys=True)
        finalDict['biome'] = str(biome_json)

        res = simplejson.dumps(finalDict)
        return HttpResponse(res, content_type='application/json')
