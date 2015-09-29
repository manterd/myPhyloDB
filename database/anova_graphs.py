from anova_DF import catUnivMetaDF, quantUnivMetaDF, normalizeUniv
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


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}


def statusANOVA(request):
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


def removeRIDANOVA(request):
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


def getCatUnivData(request):
    try:
        try:
            global base, stage, time1, TimeDiff

            # Get selected samples from cookie and query database for sample info
            samples = Sample.objects.all()
            samples.query = pickle.loads(request.session['selected_samples'])
            selected = samples.values_list('sampleid')
            qs1 = Sample.objects.all().filter(sampleid__in=selected)
        except Exception as e:
            print("Error starting ANOVA: ", e)

        if request.is_ajax():
            try:
                # Get variables from web page
                allJson = request.GET["all"]
                all = simplejson.loads(allJson)

                RID = str(all["RID"])
                time1[RID] = time.time()  # Moved these down here so RID is available
                base[RID] = 'Step 1 of 6: Querying database...'

                selectAll = int(all["selectAll"])
                DepVar = int(all["DepVar"])
                NormMeth = int(all["NormMeth"])
                NormVal = all["NormVal"]
                FDR = float(all["FdrVal"])
                StatTest = int(all["StatTest"])
                sig_only = int(all["sig_only"])
                size = int(all["MinSize"])

                # Generate a list of sequence reads per sample
                countList = []
                for sample in qs1:
                    total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                    if total['count__sum'] is not None:
                        countList.append(total['count__sum'])

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

                # Remove samples if below the sequence threshold set by user (rarefaction)
                newList = []
                result = 'Data Normalization:\n'

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

                elif NormMeth == 4 or NormMeth == 5:
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
                numberRem = len(countList) - len(newList)
                if numberRem > 0:
                    result += str(numberRem) + ' samples did not met the desired normalization criteria; and were not included in the analysis...\n'
                    result += str(len(newList)) + ' samples met the desired normalization criteria; and were included in the analysis...\n'
                else:
                    result += 'All ' + str(len(countList)) + ' selected samples were included in the analysis...\n'

                # Get dict of selected meta variables
                metaString = all["meta"]
                metaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaString)

                # Convert dict to list
                fieldList = []
                for key in metaDict:
                    fieldList.append(key)

                # Create dataframe of meta variables by sample
                metaDF = catUnivMetaDF(qs2, metaDict)
                metaDF.dropna(subset=fieldList, inplace=True)
                metaDF.sort(columns='sampleid', inplace=True)

                # Create unique list of samples in meta dataframe (may be different than selected samples)
                myList = metaDF['sampleid'].tolist()
                mySet = list(ordered_set(myList))

                # Create dataframe with all taxa/count data by sample
                taxaDF = taxaProfileDF(mySet)

                # Select only the taxa of interest if user used the taxa tree
                taxaString = all["taxa"]
                taxaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(taxaString)

                # Select only the taxa of interest if user used the selectAll button
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

                base[RID] = 'Step 1 of 6: Querying database...done!'
            except Exception as e:
                print("Error querying database (ANOVA CAT): ", e)
            # Normalize data
            try:
                base[RID] = 'Step 2 of 6: Normalizing data...'

                normDF, DESeq_error = normalizeUniv(taxaDF, taxaDict, mySet, NormMeth, NormReads, metaDF)
                normDF.sort('sampleid')

                finalDF = metaDF.merge(normDF, on='sampleid', how='outer')
                finalDF[['abund', 'rich', 'diversity']] = finalDF[['abund', 'rich', 'diversity']].astype(float)
                pd.set_option('display.max_rows', finalDF.shape[0], 'display.max_columns', finalDF.shape[1], 'display.width', 1000)

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
                result += '===============================================\n\n\n'

                base[RID] = 'Step 2 of 6: Normalizing data...done!'
            except Exception as e:
                print("Error with normalization (ANOVA CAT): ", e)
            base[RID] = 'Step 3 of 6: Performing statistical test...'

            seriesList = []
            xAxisDict = {}
            yAxisDict = {}

            # group DataFrame by each taxa level selected
            grouped1 = finalDF.groupby(['rank', 'taxa_name', 'taxa_id'])

            for name1, group1 in grouped1:
                try:
                    trtList = []
                    valList = []
                    grouped2 = pd.DataFrame()
                    if DepVar == 1:
                        grouped2 = group1.groupby(fieldList)['abund']
                    elif DepVar == 2:
                        grouped2 = group1.groupby(fieldList)['rich']
                    elif DepVar == 3:
                        grouped2 = group1.groupby(fieldList)['diversity']
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
                        if os.name == 'nt':
                            r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                        else:
                            r = R(RCMD="R/R-Linux/bin/R")
                        r.assign("df", group1)
                        trtString = " * ".join(fieldList)

                        if DepVar == 1:
                            anova_string = "fit <- aov(abund ~ " + str(trtString) + ", data=df)"
                            r.assign("cmd", anova_string)
                            r("eval(parse(text=cmd))")
                        elif DepVar == 2:
                            anova_string = "fit <- aov(rich ~ " + str(trtString) + ", data=df)"
                            r.assign("cmd", anova_string)
                            r("eval(parse(text=cmd))")
                        elif DepVar == 3:
                            anova_string = "fit <- aov(diversity ~ " + str(trtString) + ", data=df)"
                            r.assign("cmd", anova_string)
                            r("eval(parse(text=cmd))")

                        aov = r("summary(fit)")
                        pString = r("summary(fit)[[1]][['Pr(>F)']]")
                        tempStuff = pString.split(' ')
                        pList = []
                        for part in tempStuff:
                            try:
                                pList.append(float(part))
                            except:
                                placeholder = ''
                        p_val = min(pList)

                        tempStuff = aov.split('\n')
                        for part in tempStuff:
                            if part != tempStuff[0]:
                                D += part + '\n'

                        fList = []
                        for part in tempStuff:
                            if part != tempStuff[0] and part != tempStuff[1]:
                                part = part.replace('\r', '')
                                part1 = part.split(' ')
                                if part1[0] == 'Residuals':
                                    break
                                fList.append(part1[0])

                        D += '\nTukey Honest Significant Differences:\n'
                        r("tukey <- TukeyHSD(fit, ordered=TRUE)")
                        for i in fList:
                            D += str(i) + '\n'
                            hsd_string = "table <- tukey$" + str(i)
                            r.assign("cmd", hsd_string)
                            r("eval(parse(text=cmd))")
                            table = r("table")
                            tempStuff = table.split('\n')
                            for part in tempStuff:
                                if part != tempStuff[0]:
                                    D += part + '\n'

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
                                        p_val = np.nan
                                    dict1 = {'taxa_level': name1[0], 'taxa_name': name1[1], 'taxa_id': name1[2], 'sample1': val1, 'sample2': val2, 'mean1': np.mean(smp1), 'mean2': np.mean(smp2), 'stdev1': np.std(smp1), 'stdev2': np.std(smp2), 'p_value': p_val}
                                    rows_list.append(dict1)

                        pvalDF = pd.DataFrame(rows_list)
                        pvalDF.sort(columns='p_value', inplace=True)
                        pvalDF.reset_index(drop=True, inplace=True)
                        pvalDF.index += 1
                        pvalDF['p_adj'] = pvalDF.p_value * m / pvalDF.index
                        pvalDF['Sig'] = pvalDF.p_adj <= FDR
                        pvalDF[['mean1', 'mean2', 'stdev1', 'stdev2', 'p_value', 'p_adj']] = pvalDF[['mean1', 'mean2', 'stdev1', 'stdev2', 'p_value', 'p_adj']].astype(float)
                        D = pvalDF.to_string(columns=['sample1', 'sample2', 'mean1', 'mean2', 'stdev1', 'stdev2', 'p_value', 'p_adj', 'Sig'])
                        if True in pvalDF['Sig'].values:
                            p_val = 0.0
                        else:
                            p_val = 1.0
                    base[RID] = 'Step 3 of 6: Performing statistical test...done!'
                except Exception as e:
                    print("Error performing statistical test (ANOVA CAT): ", e)
                try:
                    base[RID] = 'Step 4 of 6: Formatting graph data for display...'

                    if sig_only == 1:
                        if p_val <= 0.05:
                            result += '===============================================\n'
                            result = result + 'Taxa level: ' + str(name1[0]) + '\n'
                            result = result + 'Taxa name: ' + str(name1[1]) + '\n'
                            result = result + 'Taxa ID: ' + str(name1[2]) + '\n'
                            if DepVar == 1:
                                result = result + 'Dependent Variable: Abundance' + '\n'
                            elif DepVar == 2:
                                result = result + 'Dependent Variable: Species Richness' + '\n'
                            elif DepVar == 3:
                                result = result + 'Dependent Variable: Species Diversity' + '\n'

                            result = result + '\nANOVA table:\n'

                            result = result + str(D) + '\n'
                            result += '===============================================\n'
                            result += '\n\n\n\n'

                            dataList = []
                            grouped2 = group1.groupby(fieldList).mean()

                            if DepVar == 1:
                                dataList.extend(list(grouped2['abund'].T))
                            elif DepVar == 2:
                                dataList.extend(list(grouped2['rich'].T))
                            elif DepVar == 3:
                                dataList.extend(list(grouped2['diversity'].T))

                            seriesDict = {}
                            seriesDict['name'] = name1
                            seriesDict['data'] = dataList
                            seriesList.append(seriesDict)

                            #xTitle = {}
                            #xTitle['text'] = indVar
                            #xAxisDict['title'] = xTitle
                            xAxisDict['categories'] = trtList

                            yTitle = {}
                            if DepVar == 1:
                                yTitle['text'] = 'Abundance'
                            elif DepVar == 2:
                                yTitle['text'] = 'Species Richness'
                            elif DepVar == 3:
                                yTitle['text'] = 'Species Diversity'
                            yAxisDict['title'] = yTitle

                    if sig_only == 0:
                        result += '===============================================\n'
                        result = result + 'Taxa level: ' + str(name1[0]) + '\n'
                        result = result + 'Taxa name: ' + str(name1[1]) + '\n'
                        result = result + 'Taxa ID: ' + str(name1[2]) + '\n'
                        if DepVar == 1:
                            result = result + 'Dependent Variable: Abundance' + '\n'
                        elif DepVar == 2:
                            result = result + 'Dependent Variable: Species Richness' + '\n'
                        elif DepVar == 3:
                            result = result + 'Dependent Variable: Species Diversity' + '\n'

                        result = result + '\nANOVA table:\n'

                        result = result + str(D) + '\n'
                        result += '===============================================\n'
                        result += '\n\n\n\n'

                        dataList = []
                        grouped2 = group1.groupby(fieldList).mean()

                        if DepVar == 1:
                            dataList.extend(list(grouped2['abund'].T))
                        elif DepVar == 2:
                            dataList.extend(list(grouped2['rich'].T))
                        elif DepVar == 3:
                            dataList.extend(list(grouped2['diversity'].T))

                        seriesDict = {}
                        seriesDict['name'] = name1
                        seriesDict['data'] = dataList
                        seriesList.append(seriesDict)

                        #xTitle = {}
                        #xTitle['text'] = indVar
                        #xAxisDict['title'] = xTitle
                        xAxisDict['categories'] = trtList

                        yTitle = {}
                        if DepVar == 1:
                            yTitle['text'] = 'Abundance'
                        elif DepVar == 2:
                            yTitle['text'] = 'Species Richness'
                        elif DepVar == 3:
                            yTitle['text'] = 'Species Diversity'
                        yAxisDict['title'] = yTitle
                except Exception as e:
                    print("Error with formatting (ANOVA CAT): ", e)
            try:
                finalDict['series'] = seriesList
                finalDict['xAxis'] = xAxisDict
                finalDict['yAxis'] = yAxisDict
                finalDict['text'] = result
                if not seriesList:
                    finalDict['empty'] = 0
                else:
                    finalDict['empty'] = 1

                finalDF.reset_index(drop=True, inplace=True)
                if DepVar == 1:
                    finalDF.drop(['rich', 'diversity'], axis=1, inplace=True)
                elif DepVar == 2:
                    finalDF.drop(['abund', 'diversity'], axis=1, inplace=True)
                elif DepVar == 3:
                    finalDF.drop(['abund', 'rich'], axis=1, inplace=True)

                base[RID] = 'Step 4 of 6: Formatting graph data for display...done!'
                base[RID] = 'Step 5 of 6: Formatting biome data...'

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
                    if DepVar == 1:
                        dataList.append(group['abund'].tolist())
                    if DepVar == 2:
                        dataList.append(group['rich'].tolist())
                    if DepVar == 3:
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

                base[RID] = 'Step 5 of 6: Formatting biome data...done!'
            except Exception as e:
                print("Error finishing formatting (ANOVA CAT): ", e)
            base[RID] = 'Step 6 of 6: Formatting result table...'

            res_table = finalDF.to_html(classes="table display")
            res_table = res_table.replace('border="1"', 'border="0"')
            finalDict['res_table'] = str(res_table)


            biome_json = simplejson.dumps(biome, ensure_ascii=True, indent=4, sort_keys=True)
            finalDict['error'] = 'none'
            finalDict['biome'] = str(biome_json)
            res = simplejson.dumps(finalDict)

            base[RID] = 'Step 6 of 6: Formatting result table...done!'

            return HttpResponse(res, content_type='application/json')

    except Exception as e:
        print "Error with ANOVA: ", e
        state = "Error with ANOVA: " + str(e)

        myDict = {}
        myDict['error'] = state
        res = simplejson.dumps(myDict)
        return HttpResponse(res, content_type='application/json')


def getQuantUnivData(request):
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
            base[RID] = 'Step 1 of 6: Querying database...'

            selectAll = int(all["selectAll"])
            DepVar = int(all["DepVar"])
            NormMeth = int(all["NormMeth"])
            NormVal = all["NormVal"]
            sig_only = int(all["sig_only"])
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
            result = 'Data Normalization:\n'

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

            elif NormMeth == 4 or NormMeth == 5:
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
            numberRem = len(countList) - len(newList)
            if numberRem > 0:
                result += str(numberRem) + ' samples did not met the desired normalization criteria; and were not included in the analysis...\n'
                result += str(len(newList)) + ' samples met the desired normalization criteria; and were included in the analysis...\n'
            else:
                result += 'All ' + str(len(countList)) + ' selected samples were included in the analysis...\n'

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
            metaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStr)
            for key in metaDict:
                fieldList.append(key)

            metaDF = catUnivMetaDF(qs2, metaDict)
            metaDF.dropna(subset=fieldList, inplace=True)
            metaDF.sort(columns='sampleid', inplace=True)

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

            base[RID] = 'Step 1 of 6: Querying database...done!'
            base[RID] = 'Step 2 of 6: Normalizing data...'

            normDF, DESeq_error = normalizeUniv(taxaDF, taxaDict, mySet, NormMeth, NormReads, metaDF)
            normDF.sort('sampleid')

            finalDF = metaDF.merge(normDF, on='sampleid', how='outer')
            finalDF[['abund', 'rich', 'diversity']] = finalDF[['abund', 'rich', 'diversity']].astype(float)
            pd.set_option('display.max_rows', finalDF.shape[0], 'display.max_columns', finalDF.shape[1], 'display.width', 1000)

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
            result += '===============================================\n\n\n'

            base[RID] = 'Step 2 of 6: Normalizing data...done!'
            base[RID] = 'Step 3 of 6: Performing linear regression...'

            seriesList = []
            xAxisDict = {}
            yAxisDict = {}
            colors = ['#000000', '#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#cccc00', '#a65628', '#f781bf']
            shapes = ['circle', 'square', 'triangle', 'triangle-down', 'diamond']
            colors_idx = 0
            shapes_idx = 0
            grouped1 = finalDF.groupby(['rank', 'taxa_name', 'taxa_id'])
            for name1, group1 in grouped1:
                dataList = []
                x = []

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R")

                r.assign("df", group1)
                trtString = "*".join(fieldList)

                D = ""
                if DepVar == 1:
                    anova_string = "fit <- lm(abund ~ " + str(trtString) + ", data=df)"
                    r.assign("cmd", anova_string)
                    r("eval(parse(text=cmd))")

                elif DepVar == 2:
                    anova_string = "fit <- lm(rich ~ " + str(trtString) + ", data=df)"
                    r.assign("cmd", anova_string)
                    r("eval(parse(text=cmd))")

                elif DepVar == 3:
                    anova_string = "fit <- lm(diversity ~ " + str(trtString) + ", data=df)"
                    r.assign("cmd", anova_string)
                    r("eval(parse(text=cmd))")

                anova_string = "df$pred <- predict(fit, df)"
                r.assign("cmd", anova_string)
                r("eval(parse(text=cmd))")
                aov = r("summary(fit)")

                tempStuff = aov.split('\n')
                for part in tempStuff:
                    if part != tempStuff[0]:
                        D += part + '\n'

                r("p_vals <- summary(fit)$coefficients[,4]")
                p_vals = r.get("p_vals")
                p_value = min(p_vals)

                resultDF = r.get("df")
                resultDF = resultDF.rename(columns=lambda x: x.strip())

                if fieldListCat:
                    grouped2 = resultDF.groupby(fieldListCat)
                    for name2, group2 in grouped2:
                        if DepVar == 1:
                            dataList = group2[[fieldListQuant[0], 'abund']].values.astype(float).tolist()
                            x = group2[fieldListQuant[0]].values.tolist()
                        elif DepVar == 2:
                            dataList = group2[[fieldListQuant[0], 'rich']].values.astype(float).tolist()
                            x = group2[fieldListQuant[0]].values.tolist()
                        elif DepVar == 3:
                            dataList = group2[[fieldListQuant[0], 'diversity']].values.astype(float).tolist()
                            x = group2[fieldListQuant[0]].values.tolist()

                        if max(x) == min(x):
                            stop = 0
                        else:
                            stop = 1
                            seriesDict = {}
                            seriesDict['turboThreshold'] = 0
                            seriesDict['regression'] = 'true'
                            regDict = {}
                            regDict['color'] = colors[colors_idx]
                            seriesDict['regressionSettings'] = regDict
                            name = str(name1[1]) + ": " + str(name2)
                            seriesDict['name'] = name
                            seriesDict['color'] = colors[colors_idx]
                            markerDict = {}
                            markerDict['symbol'] = shapes[shapes_idx]
                            seriesDict['marker'] = markerDict
                            seriesDict['data'] = dataList
                            seriesList.append(seriesDict)

                        shapes_idx += 1
                        if shapes_idx >= len(shapes):
                            shapes_idx = 0

                if not fieldListCat:
                    if DepVar == 1:
                        dataList = resultDF[[str(fieldListQuant[0]), 'abund']].values.astype(float).tolist()
                        x = resultDF[fieldListQuant[0]].values.tolist()
                    elif DepVar == 2:
                        dataList = resultDF[[fieldListQuant[0], 'rich']].values.astype(float).tolist()
                        x = resultDF[fieldListQuant[0]].values.tolist()
                    elif DepVar == 3:
                        dataList = resultDF[[fieldListQuant[0], 'diversity']].values.astype(float).tolist()
                        x = resultDF[fieldListQuant[0]].values.tolist()
                    print 'dataList:', dataList
                    if max(x) == min(x):
                        stop = 0
                    else:
                        stop = 1
                        seriesDict = {}
                        seriesDict['turboThreshold'] = 0
                        seriesDict['regression'] = 'true'
                        regDict = {}
                        regDict['color'] = colors[colors_idx]
                        seriesDict['regressionSettings'] = regDict
                        name = name1[1]
                        seriesDict['name'] = name
                        seriesDict['color'] = colors[colors_idx]
                        markerDict = {}
                        markerDict['symbol'] = shapes[shapes_idx]
                        seriesDict['marker'] = markerDict
                        seriesDict['data'] = dataList
                        seriesList.append(seriesDict)

                base[RID] = 'Step 3 of 6: Performing linear regression...done!'
                base[RID] = 'Step 4 of 6: Formatting graph data for display...'

                if sig_only == 0:
                    result = result + '\nANCOVA table:\n'

                    result = result + str(D) + '\n'
                    result += '===============================================\n'
                    result += '\n\n\n\n'

                if sig_only == 1:
                    if p_value <= 0.05:
                        result = result + '\nANCOVA table:\n'

                        result = result + str(D) + '\n'
                        result += '===============================================\n'
                        result += '\n\n\n\n'

                xTitle = {}
                xTitle['text'] = fieldList[0]
                xAxisDict['title'] = xTitle

                yTitle = {}
                if DepVar == 1:
                    yTitle['text'] = 'Abundance'
                elif DepVar == 2:
                    yTitle['text'] = 'Species Richness'
                elif DepVar == 3:
                    yTitle['text'] = 'Shannon Diversity'
                yAxisDict['title'] = yTitle

                colors_idx += 1
                if colors_idx >= len(colors):
                    colors_idx = 0

            finalDict['series'] = seriesList
            finalDict['xAxis'] = xAxisDict
            finalDict['yAxis'] = yAxisDict
            if not seriesList:
                finalDict['empty'] = 0
            else:
                finalDict['empty'] = 1

            finalDF.reset_index(drop=True, inplace=True)
            if DepVar == 1:
                finalDF.drop(['rich', 'diversity'], axis=1, inplace=True)
            elif DepVar == 2:
                finalDF.drop(['abund', 'diversity'], axis=1, inplace=True)
            elif DepVar == 3:
                finalDF.drop(['abund', 'rich'], axis=1, inplace=True)

            base[RID] = 'Step 4 of 6: Formatting graph data for display...'
            base[RID] = 'Step 5 of 6: Formatting biome data...'

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
                if DepVar == 1:
                    dataList.append(group['abund'].tolist())
                if DepVar == 2:
                    dataList.append(group['rich'].tolist())
                if DepVar == 3:
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

            base[RID] = 'Step 5 of 6: Formatting biome data...done!'
            base[RID] = 'Step 6 of 6: Formatting result table...'

            res_table = finalDF.to_html(classes="table display")
            res_table = res_table.replace('border="1"', 'border="0"')
            finalDict['res_table'] = str(res_table)

            finalDict['text'] = result
            base[RID] = 'Step 4 of 4: Preparing graph data...completed'

            biome_json = simplejson.dumps(biome, ensure_ascii=True, indent=4, sort_keys=True)
            finalDict['biome'] = str(biome_json)

            base[RID] = 'Step 6 of 6: Formatting result table...done!'

            finalDict['error'] = 'none'
            res = simplejson.dumps(finalDict)

            return HttpResponse(res, content_type='application/json')

    except Exception as e:
        print "Error with Linear Regression: ", e
        state = "Error with Linear Regression: " + str(e)

        myDict = {}
        myDict['error'] = state
        res = simplejson.dumps(myDict)
        return HttpResponse(res, content_type='application/json')


