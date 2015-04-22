from anova_DF import catUnivMetaDF, quantUnivMetaDF, normalizeUniv
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
import json


stage = ''


def statusANOVA(request):
    global stage
    if request.is_ajax():
        myDict = {}
        myDict['stage'] = stage
        json_data = simplejson.dumps(myDict, encoding="Latin-1")
        return HttpResponse(json_data, content_type='application/json')


def getCatUnivData(request):
    global stage
    stage = 'Step 1 of 4: Querying database...'
    samples = Sample.objects.all()
    samples.query = pickle.loads(request.session['selected_samples'])
    selected = samples.values_list('sampleid')
    qs1 = Sample.objects.all().filter(sampleid__in=selected)

    if request.is_ajax():
        allJson = request.GET["all"]
        all = simplejson.loads(allJson)
        button = int(all["button"])
        sig_only = int(all["sig_only"])
        factor = all["normalize"]
        selectAll = int(all["selectAll"])
        remove = int(all["remove"])
        stat = int(all["statistic"])

        countList = []
        for sample in qs1:
            total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
            if total['count__sum'] is not None:
                countList.append(total['count__sum'])

        minSize = int(min(countList))
        medianSize = int(np.median(np.array(countList)))
        maxSize = int(max(countList))

        if factor == "min":
            norm = minSize
        elif factor == "median":
            norm = medianSize
        elif factor == "max":
            norm = maxSize
        elif factor == "none":
            norm = -1
        else:
            norm = int(all["normalize"])

        newList = []
        for sample in qs1:
            total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
            if remove == 1:
                if total['count__sum'] is not None and int(total['count__sum']) >= norm:
                    id = sample.sampleid
                    newList.append(id)
            else:
                if total['count__sum'] is not None:
                    id = sample.sampleid
                    newList.append(id)

        qs2 = Sample.objects.all().filter(sampleid__in=newList)

        metaString = all["meta"]
        metaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaString)

        fieldList = []
        for key in metaDict:
            fieldList.append(key)

        metaDF = catUnivMetaDF(qs2, metaDict)
        metaDF.dropna(subset=fieldList, inplace=True)
        metaDF.sort(columns='sample_name', inplace=True)

        myList = metaDF['sampleid'].tolist()
        mySet = list(ordered_set(myList))

        taxaDF = taxaProfileDF(mySet)

        taxaString = all["taxa"]

        taxaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(taxaString)
        if selectAll == 1:
            taxaDict = {}
            qs3 = Profile.objects.all().filter(sampleid__in=mySet).values_list('kingdomid', flat='True').distinct()
            taxaDict['Kingdom'] = qs3
        elif selectAll == 2:
            taxaDict = {}
            qs3 = Profile.objects.all().filter(sampleid__in=mySet).values_list('phylaid', flat='True').distinct()
            taxaDict['Phyla'] = qs3
        elif selectAll == 3:
            taxaDict = {}
            qs3 = Profile.objects.all().filter(sampleid__in=mySet).values_list('classid', flat='True').distinct()
            taxaDict['Class'] = qs3
        elif selectAll == 4:
            taxaDict = {}
            qs3 = Profile.objects.all().filter(sampleid__in=mySet).values_list('orderid', flat='True').distinct()
            taxaDict['Order'] = qs3
        elif selectAll == 5:
            taxaDict = {}
            qs3 = Profile.objects.all().filter(sampleid__in=mySet).values_list('familyid', flat='True').distinct()
            taxaDict['Family'] = qs3
        elif selectAll == 6:
            taxaDict = {}
            qs3 = Profile.objects.all().filter(sampleid__in=mySet).values_list('genusid', flat='True').distinct()
            taxaDict['Genus'] = qs3
        elif selectAll == 7:
            taxaDict = {}
            qs3 = Profile.objects.all().filter(sampleid__in=mySet).values_list('speciesid', flat='True').distinct()
            taxaDict['Species'] = qs3

        stage = 'Step 1 of 4: Querying database...complete'
        stage = 'Step 2 of 4: Normalizing data...'
        normDF = normalizeUniv(taxaDF, taxaDict, mySet, norm, button)
        stage = 'Step 2 of 4: Normalizing data...complete'
        stage = 'Step 3 of 4: Performing ANOVA...'
        finalDF = metaDF.merge(normDF, on='sampleid', how='outer')
        finalDF[['count', 'rel_abund', 'rich', 'diversity']] = finalDF[['count', 'rel_abund', 'rich', 'diversity']].astype(float)
        pd.set_option('display.max_rows', finalDF.shape[0], 'display.max_columns', finalDF.shape[1], 'display.width', 1000)

        finalDict = {}
        result = 'Data Normalization:\n'
        if norm == -1:
            result = result + 'Data were not normalized...\n'
        else:
            result = result + 'Data normalized to ' + str(norm) + ' sequence reads...\n'
        if remove == 1:
            result = result + 'Samples below this threshold were removed\n\n'
        if remove == 0:
            result = result + 'All selected samples were included in the analysis\n\n'

        seriesList = []
        xAxisDict = {}
        yAxisDict = {}

        ### group DataFrame by each taxa level selected
        grouped1 = finalDF.groupby(['rank', 'taxa_name', 'taxa_id'])

        ### group DataFrame by each meta variable selected
        for name1, group1 in grouped1:
            trtList = []
            valList = []
            grouped2 = pd.DataFrame()
            if button == 1:
                grouped2 = group1.groupby(fieldList)['count']
            elif button == 2:
                grouped2 = group1.groupby(fieldList)['rel_abund']
            elif button == 3:
                grouped2 = group1.groupby(fieldList)['rich']
            elif button == 4:
                grouped2 = group1.groupby(fieldList)['diversity']

            for name2, group2 in grouped2:
                if isinstance(name2, unicode):
                    trt = name2
                else:
                    trt = ' & '.join(list(name2))
                trtList.append(trt)
                valList.append(list(group2.T))

            D = Anova1way()
            try:
                D.run(valList, conditions_list=trtList)
                anova_error = 'no'
            except:
                D = 'a or b too big, or ITMAX too small in Betacf.\n' + 'Samples sizes are either too unequal or n=1.\n'
                anova_error = 'yes'

            stage = 'Step 3 of 4: Performing ANOVA...complete'
            stage = 'Step 4 of 4: Preparing graph data...'
            if sig_only == 1:
                if float(D['p']) <= 0.05:
                    result = result + '===============================================\n'
                    result = result + 'Taxa level: ' + str(name1[0]) + '\n'
                    result = result + 'Taxa name: ' + str(name1[1]) + '\n'
                    result = result + 'Taxa ID: ' + str(name1[2]) + '\n'
                    if button == 1:
                        result = result + 'Dependent Variable: Sequence Reads' + '\n'
                    elif button == 2:
                        result = result + 'Dependent Variable: Relative Abundance' + '\n'
                    elif button == 3:
                        result = result + 'Dependent Variable: Species Richness' + '\n'
                    elif button == 4:
                        result = result + 'Dependent Variable: Shannon Diversity' + '\n'

                    indVar = ' x '.join(fieldList)
                    result = result + 'Independent Variable: ' + str(indVar) + '\n'

                    if anova_error == 'yes':
                        result = result + '\nANOVA cannot be performed...' + '\n' + D + '\n'
                    else:
                        result = result + str(D) + '\n'
                    result = result + '===============================================\n'
                    result = result + '\n\n\n\n'

                    dataList = []
                    if stat == 1:
                        grouped2 = group1.groupby(fieldList).mean()
                    else:
                        grouped2 = group1.groupby(fieldList).sum()

                    if button == 1:
                        dataList.extend(list(grouped2['count'].T))
                    elif button == 2:
                        dataList.extend(list(grouped2['rel_abund'].T))
                    elif button == 3:
                        dataList.extend(list(grouped2['rich'].T))
                    elif button == 4:
                        dataList.extend(list(grouped2['diversity'].T))

                    seriesDict = {}
                    seriesDict['name'] = name1
                    seriesDict['data'] = dataList
                    seriesList.append(seriesDict)

                    xTitle = {}
                    xTitle['text'] = indVar
                    xAxisDict['title'] = xTitle
                    xAxisDict['categories'] = trtList

                    yTitle = {}
                    if button == 1:
                        yTitle['text'] = 'Sequence Reads'
                    elif button == 2:
                        yTitle['text'] = 'Relative Abundance'
                    elif button == 3:
                        yTitle['text'] = 'Species Richness'
                    elif button == 4:
                        yTitle['text'] = 'Shannon Diversity'
                    yAxisDict['title'] = yTitle

            if sig_only == 0:
                result = result + '===============================================\n'
                result = result + 'Taxa level: ' + str(name1[0]) + '\n'
                result = result + 'Taxa name: ' + str(name1[1]) + '\n'
                result = result + 'Taxa ID: ' + str(name1[2]) + '\n'
                if button == 1:
                    result = result + 'Dependent Variable: Sequence Reads' + '\n'
                elif button == 2:
                    result = result + 'Dependent Variable: Relative Abundance' + '\n'
                elif button == 3:
                    result = result + 'Dependent Variable: Species Richness' + '\n'
                elif button == 4:
                    result = result + 'Dependent Variable: Shannon Diversity' + '\n'

                indVar = ' x '.join(fieldList)
                result = result + 'Independent Variable: ' + str(indVar) + '\n'

                if anova_error == 'yes':
                        result = result + '\nANOVA cannot be performed...' + '\n' + D + '\n'
                else:
                    result = result + str(D) + '\n'
                result = result + '===============================================\n'
                result = result + '\n\n\n\n'

                dataList = []
                if stat == 1:
                    grouped2 = group1.groupby(fieldList).mean()
                else:
                    grouped2 = group1.groupby(fieldList).sum()
                if button == 1:
                    dataList.extend(list(grouped2['count'].T))
                elif button == 2:
                    dataList.extend(list(grouped2['rel_abund'].T))
                elif button == 3:
                    dataList.extend(list(grouped2['rich'].T))
                elif button == 4:
                    dataList.extend(list(grouped2['diversity'].T))

                seriesDict = {}
                seriesDict['name'] = name1
                seriesDict['data'] = dataList
                seriesList.append(seriesDict)

                xTitle = {}
                xTitle['text'] = indVar
                xAxisDict['title'] = xTitle
                xAxisDict['categories'] = trtList

                yTitle = {}
                if button == 1:
                    yTitle['text'] = 'Sequence Reads'
                elif button == 2:
                    yTitle['text'] = 'Relative Abundance'
                elif button == 3:
                    yTitle['text'] = 'Species Richness'
                elif button == 4:
                    yTitle['text'] = 'Shannon Diversity'
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
        if button == 1:
            finalDF.drop(['rel_abund', 'rich', 'diversity'], axis=1, inplace=True)
        elif button == 2:
            finalDF.drop(['count', 'rich', 'diversity'], axis=1, inplace=True)
        elif button == 3:
            finalDF.drop(['count', 'rel_abund', 'diversity'], axis=1, inplace=True)
        elif button == 4:
            finalDF.drop(['count', 'rel_abund', 'rich'], axis=1, inplace=True)

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
            if button == 1:
                dataList.append(group['count'].tolist())
            if button == 2:
                dataList.append(group['rel_abund'].tolist())
            if button == 3:
                dataList.append(group['rich'].tolist())
            if button == 4:
                dataList.append(group['diversity'].tolist())

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


def getQuantUnivData(request):
    global stage
    stage = 'Step 1 of 4: Querying database...'
    samples = Sample.objects.all()
    samples.query = pickle.loads(request.session['selected_samples'])
    selected = samples.values_list('sampleid')
    qs1 = Sample.objects.all().filter(sampleid__in=selected)

    if request.is_ajax():
        allJson = request.GET["all"]
        all = simplejson.loads(allJson)
        button = int(all["button"])
        sig_only = int(all["sig_only"])
        factor = all["normalize"]
        selectAll = int(all["selectAll"])
        remove = int(all["remove"])

        countList = []
        for sample in qs1:
            total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
            if total['count__sum'] is not None:
                countList.append(total['count__sum'])

        minSize = int(min(countList))
        medianSize = int(np.median(np.array(countList)))
        maxSize = int(max(countList))

        if factor == "min":
            norm = minSize
        elif factor == "median":
            norm = medianSize
        elif factor == "max":
            norm = maxSize
        elif factor == "none":
            norm = -1
        else:
            norm = int(all["normalize"])

        newList = []
        for sample in qs1:
            total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
            if remove == 1:
                if total['count__sum'] is not None and int(total['count__sum']) >= norm:
                    id = sample.sampleid
                    newList.append(id)
            else:
                if total['count__sum'] is not None:
                    id = sample.sampleid
                    newList.append(id)

        qs2 = Sample.objects.all().filter(sampleid__in=newList)

        metaString = all["meta"]
        metaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaString)

        fieldList = []
        for key in metaDict:
            fieldList.append(metaDict[key])

        metaDF = quantUnivMetaDF(qs2, metaDict)
        metaDF.dropna(subset=fieldList, inplace=True)
        metaDF.sort(columns='sample_name', inplace=True)

        myList = metaDF['sampleid'].tolist()
        mySet = list(ordered_set(myList))
        taxaDF = taxaProfileDF(mySet)

        taxaString = all["taxa"]

        taxaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(taxaString)
        if selectAll == 1:
            taxaDict = {}
            qs3 = Profile.objects.all().filter(sampleid__in=mySet).values_list('kingdomid', flat='True').distinct()
            taxaDict['Kingdom'] = qs3
        elif selectAll == 2:
            taxaDict = {}
            qs3 = Profile.objects.all().filter(sampleid__in=mySet).values_list('phylaid', flat='True').distinct()
            taxaDict['Phyla'] = qs3
        elif selectAll == 3:
            taxaDict = {}
            qs3 = Profile.objects.all().filter(sampleid__in=mySet).values_list('classid', flat='True').distinct()
            taxaDict['Class'] = qs3
        elif selectAll == 4:
            taxaDict = {}
            qs3 = Profile.objects.all().filter(sampleid__in=mySet).values_list('orderid', flat='True').distinct()
            taxaDict['Order'] = qs3
        elif selectAll == 5:
            taxaDict = {}
            qs3 = Profile.objects.all().filter(sampleid__in=mySet).values_list('familyid', flat='True').distinct()
            taxaDict['Family'] = qs3
        elif selectAll == 6:
            taxaDict = {}
            qs3 = Profile.objects.all().filter(sampleid__in=mySet).values_list('genusid', flat='True').distinct()
            taxaDict['Genus'] = qs3
        elif selectAll == 7:
            taxaDict = {}
            qs3 = Profile.objects.all().filter(sampleid__in=mySet).values_list('speciesid', flat='True').distinct()
            taxaDict['Species'] = qs3

        stage = 'Step 1 of 4: Querying database...complete'
        stage = 'Step 2 of 4: Normalizing data...'
        normDF = normalizeUniv(taxaDF, taxaDict, mySet, norm, button)
        stage = 'Step 2 of 4: Normalizing data...complete'
        stage = 'Step 3 of 4: Performing linear regression...'

        finalDF = metaDF.merge(normDF, on='sampleid', how='outer')
        finalDF[[fieldList[0], 'count', 'rel_abund', 'rich', 'diversity']] = finalDF[[fieldList[0], 'count', 'rel_abund', 'rich', 'diversity']].astype(float)
        pd.set_option('display.max_rows', finalDF.shape[0], 'display.max_columns', finalDF.shape[1], 'display.width', 1000)

        finalDict = {}
        seriesList = []
        xAxisDict = {}
        yAxisDict = {}
        grouped1 = finalDF.groupby(['rank', 'taxa_name', 'taxa_id'])
        for name1, group1 in grouped1:
            dataList = []
            x = []
            y = []
            if button == 1:
                dataList = group1[[fieldList[0], 'count']].values.tolist()
                x = group1[fieldList[0]].values.tolist()
                y = group1['count'].values.tolist()
            elif button == 2:
                dataList = group1[[fieldList[0], 'rel_abund']].values.tolist()
                x = group1[fieldList[0]].values.tolist()
                y = group1['rel_abund'].values.tolist()
            elif button == 3:
                dataList = group1[[fieldList[0], 'rich']].values.tolist()
                x = group1[fieldList[0]].values.tolist()
                y = group1['rich'].values.tolist()
            elif button == 4:
                dataList = group1[[fieldList[0], 'diversity']].values.tolist()
                x = group1[fieldList[0]].values.tolist()
                y = group1['diversity'].values.tolist()

            if max(x) == min(x):
                stop = 0
            else:
                stop = 1
                slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
                p_value = "%0.3f" % p_value
                r_square = r_value * r_value
                r_square = "%0.4f" % r_square
                min_y = slope*min(x) + intercept
                max_y = slope*max(x) + intercept
                slope = "%.3E" % slope
                intercept = "%.3E" % intercept

                regrList = []
                regrList.append([min(x), min_y])
                regrList.append([max(x), max_y])

            stage = 'Step 3 of 4: Performing linear regression...complete'
            stage = 'Step 4 of 4: Preparing graph data...'
            if sig_only == 0:
                seriesDict = {}
                seriesDict['type'] = 'scatter'
                seriesDict['name'] = name1
                seriesDict['data'] = dataList
                seriesList.append(seriesDict)
                if stop == 0:
                    regDict = {}
                elif stop == 1:
                    regrDict = {}
                    regrDict['type'] = 'line'
                    name2 = list(name1)
                    temp = 'R2: ' + str(r_square) + '; p-value: ' + str(p_value) + '<br>' + '(y = ' + str(slope) + 'x' + ' + ' + str(intercept) + ')'
                    name2.append(temp)
                    regrDict['name'] = name2
                    regrDict['data'] = regrList
                    seriesList.append(regrDict)

            if sig_only == 1:
                if p_value <= 0.05:
                    seriesDict = {}
                    seriesDict['type'] = 'scatter'
                    name2 = list(name1)
                    temp = 'R2: ' + str(r_square) + '; p-value: ' + str(p_value) + '<br>' + '(y = ' + str(slope) + 'x' + ' + ' + str(intercept) + ')'
                    name2.append(temp)
                    seriesDict['name'] = name2
                    seriesDict['data'] = dataList
                    seriesList.append(seriesDict)

                    regrDict = {}
                    regrDict['type'] = 'line'
                    regrDict['name'] = name1
                    regrDict['data'] = regrList
                    seriesList.append(regrDict)

            xTitle = {}
            xTitle['text'] = fieldList[0]
            xAxisDict['title'] = xTitle

            yTitle = {}
            if button == 1:
                yTitle['text'] = 'Sequence Reads'
            elif button == 2:
                yTitle['text'] = 'Relative Abundance'
            elif button == 3:
                yTitle['text'] = 'Species Richness'
            elif button == 4:
                yTitle['text'] = 'Shannon Diversity'
            yAxisDict['title'] = yTitle

        finalDict['series'] = seriesList
        finalDict['xAxis'] = xAxisDict
        finalDict['yAxis'] = yAxisDict
        if not seriesList:
            finalDict['empty'] = 0
        else:
            finalDict['empty'] = 1

        finalDF.reset_index(drop=True, inplace=True)
        if button == 1:
            finalDF.drop(['rel_abund', 'rich', 'diversity'], axis=1, inplace=True)
        elif button == 2:
            finalDF.drop(['count', 'rich', 'diversity'], axis=1, inplace=True)
        elif button == 3:
            finalDF.drop(['count', 'rel_abund', 'diversity'], axis=1, inplace=True)
        elif button == 4:
            finalDF.drop(['count', 'rel_abund', 'rich'], axis=1, inplace=True)

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
            if button == 1:
                dataList.append(group['count'].tolist())
            if button == 2:
                dataList.append(group['rel_abund'].tolist())
            if button == 3:
                dataList.append(group['rich'].tolist())
            if button == 4:
                dataList.append(group['diversity'].tolist())

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

        result = 'Data Normalization:\n'
        if norm == -1:
            result = result + 'Data were not normalized...\n'
        else:
            result = result + 'Data normalized to ' + str(norm) + ' sequence reads...\n'
        if remove == 1:
            result = result + 'Samples below this threshold were removed\n\n'
        if remove == 0:
            result = result + 'All selected samples were included in the analysis\n\n'
        finalDict['text'] = result

        stage = 'Step 4 of 4: Preparing graph data...complete'

        biome_json = simplejson.dumps(biome, ensure_ascii=True, indent=4, sort_keys=True)
        finalDict['biome'] = str(biome_json)

        res = simplejson.dumps(finalDict)
        return HttpResponse(res, content_type='application/json')


