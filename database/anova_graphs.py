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
                totalSamp, columns = metaDF.shape

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
            colors = [
                "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
                "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
                "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
                "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58",
                "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
                "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
                "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5",
                "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
                "#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01",
                "#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966",
                "#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
                "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", "#C6DC99", "#203B3C",
                "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
                "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
                "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
                "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F",
                "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
                "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
                "#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00",
                "#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
                "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25",
                "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"
            ]

            colors_idx = 0
            for name1, group1 in grouped1:
                try:
                    trtList = []
                    valList = []

                    grouped2 = group1.groupby(fieldList)

                    xAxisList = []
                    for name2, group2 in grouped2:
                        if isinstance(name2, unicode):
                            trt = name2
                        else:
                            trt = ' & '.join(list(name2))
                        trtList.append(trt)
                        valList.append(list(group2.T))

                    grouped3 = group1.groupby(fieldList).mean()

                    catList = []
                    for i in xrange(len(fieldList)):
                        catList.append(grouped3.index.get_level_values(i).unique().tolist())

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
                        if p_val >= 0.05:
                            result += '\n===============================================\n\n'
                            result += 'Taxa Rank: ' + str(name1[0]) + '; Taxa Name: ' + str(name1[1]) + '; Taxa ID: ' + str(name1[2]) + ' is not significant,\n'
                            result += 'it was removed from your analysis.\n'
                            result += '\n===============================================\n\n'

                        if p_val < 0.05:
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
                            seriesDict['color'] = colors[colors_idx]
                            seriesDict['data'] = dataList
                            seriesList.append(seriesDict)

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
                        seriesDict['color'] = colors[colors_idx]
                        seriesDict['data'] = dataList
                        seriesList.append(seriesDict)

                        if catList.__len__() == 1:
                            xAxisDict['categories'] = catList[0]
                        else:
                            finalList = []
                            finalList.append(recLabel("", catList))
                            xAxisDict['categories'] = finalList

                        yTitle = {}
                        if DepVar == 1:
                            yTitle['text'] = 'Abundance'
                        elif DepVar == 2:
                            yTitle['text'] = 'Species Richness'
                        elif DepVar == 3:
                            yTitle['text'] = 'Species Diversity'
                        yAxisDict['title'] = yTitle

                    colors_idx += 1
                    if colors_idx >= len(colors):
                        colors_idx = 0

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


def recLabel(name, list):
    retDict = {}
    if list.__len__() == 1:
        retDict['name'] = name
        retDict['categories'] = list[0]
        return retDict
    subList = []
    for stuff in list[0]:
        subList.append(recLabel(stuff, list[1:]))
    retDict['name'] = name
    retDict['categories'] = subList
    return retDict


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

            metaStrCat = all["metaCat"]
            fieldListCat = []
            fieldListVals = []
            if metaStrCat:
                metaDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStrCat)
                for key in metaDictCat:
                    fieldListCat.append(key)
                    fieldListVals.append(metaDictCat[key])

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
            metaDF = catUnivMetaDF(qs2, metaDict)
            metaDF.dropna(subset=fieldList, inplace=True)
            metaDF.sort(columns='sampleid', inplace=True)
            totalSamp, cols = metaDF.shape

            normRem = len(countList) - len(newList)
            selectRem = len(newList) - totalSamp

            result += str(totalSamp) + ' selected samples were included in the final analysis.\n'
            if normRem > 0:
                result += str(normRem) + ' samples did not met the desired normalization criteria.\n'
            if selectRem:
                result += str(selectRem) + ' samples were deselected by the user.\n'

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
            colors = [
                "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
                "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
                "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
                "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58",
                "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
                "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
                "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5",
                "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
                "#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01",
                "#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966",
                "#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
                "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", "#C6DC99", "#203B3C",
                "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
                "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
                "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
                "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F",
                "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
                "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
                "#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00",
                "#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
                "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25",
                "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"
            ]
            shapes = ['circle', 'square', 'triangle', 'triangle-down', 'diamond']
            colors_idx = 0
            shapes_idx = 0
            grouped1 = finalDF.groupby(['rank', 'taxa_name', 'taxa_id'])
            for name1, group1 in grouped1:
                dataList = []
                x = []
                y = []
                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R")

                r.assign("df", group1)
                if len(fieldListVals) > 1:
                    trtString = "*".join(fieldList)
                else:
                    trtString = fieldListQuant[0]

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

                if p_vals is not None:
                    p_value = min(p_vals)  # TODO this one
                else:
                    p_value = 1.0

                resultDF = r.get("df")
                resultDF = resultDF.rename(columns=lambda x: x.strip())

                if len(fieldListVals) > 1:
                    grouped2 = resultDF.groupby(fieldListCat)
                    for name2, group2 in grouped2:
                        obs = group2['abund'].loc[group2['abund'] > 0].count()
                        if obs < 3:
                            result += '\n===============================================\n\n'
                            result += 'Taxa Rank: ' + str(name1[0]) + '; Taxa Name: ' + str(name1[1]) + '; Taxa ID: ' + str(name1[2]) + ' has ' + str(obs) + ' observation(s),\n'
                            result += 'it was removed from your analysis.\n'
                            result += '\n===============================================\n\n'
                            break

                        else:
                            if DepVar == 1:
                                dataList = group2[[fieldListQuant[0], 'abund']].values.astype(float).tolist()
                                x = group2[fieldListQuant[0]].astype(float).values.tolist()
                                y = group2['abund'].values.astype(float).tolist()
                            elif DepVar == 2:
                                dataList = group2[[fieldListQuant[0], 'rich']].values.astype(float).tolist()
                                x = group2[fieldListQuant[0]].values.astype(float).tolist()
                                y = group2['rich'].values.astype(float).tolist()
                            elif DepVar == 3:
                                dataList = group2[[fieldListQuant[0], 'diversity']].values.astype(float).tolist()
                                x = group2[fieldListQuant[0]].astype(float).values.tolist()
                                y = group2['diversity'].astype(float).values.tolist()

                            slp, inter, r_value, p, std_err = stats.linregress(x, y)
                            min_y = float(slp*min(x) + inter)
                            max_y = float(slp*max(x) + inter)
                            slope = "%0.3f" % slp
                            intercept = "%0.3f" % inter

                            if sig_only == 0:
                                result += 'Taxa level: ' + str(name1[0]) + '\n'
                                result += 'Taxa name: ' + str(name1[1]) + '\n'
                                result += 'Taxa ID: ' + str(name1[2]) + '\n'
                                if DepVar == 1:
                                    result = result + 'Dependent Variable: Abundance' + '\n'
                                elif DepVar == 2:
                                    result = result + 'Dependent Variable: Species Richness' + '\n'
                                elif DepVar == 3:
                                    result = result + 'Dependent Variable: Species Diversity' + '\n'
                                result += '\nANCOVA table:\n'
                                result += str(D) + '\n'
                                result += '===============================================\n'
                                result += '\n'

                                regrList = []
                                regrList.append([float(min(x)), min_y])
                                regrList.append([float(max(x)), max_y])

                                seriesDict = {}
                                seriesDict['turboThreshold'] = 0
                                seriesDict['type'] = 'scatter'
                                seriesDict['name'] = str(name1[1]) + ": " + str(name2)
                                seriesDict['data'] = dataList
                                seriesDict['color'] = colors[colors_idx]

                                markerDict = {}
                                markerDict['symbol'] = shapes[shapes_idx]
                                seriesDict['marker'] = markerDict
                                seriesDict['data'] = dataList
                                seriesList.append(seriesDict)

                                regrDict = {}
                                regrDict['type'] = 'line'
                                regrDict['name'] = 'y = ' + str(slope) + 'x' + ' + ' + str(intercept)
                                regrDict['data'] = regrList
                                regrDict['color'] = colors[colors_idx]

                                markerDict = {}
                                markerDict['enabled'] = False
                                regrDict['marker'] = markerDict

                                seriesList.append(regrDict)

                            if sig_only == 1:
                                if p_value >= 0.05:
                                    result += '\n===============================================\n\n'
                                    result += 'Treatment Level: ' + str(name2) + '\n'
                                    result += 'Taxa Rank: ' + str(name1[0]) + '; Taxa Name: ' + str(name1[1]) + '; Taxa ID: ' + str(name1[2]) + ' is not significant,\n'
                                    result += 'it was removed from your analysis.\n'
                                    result += '\n===============================================\n\n'
                                if p_value <= 0.05:
                                    result += 'Taxa level: ' + str(name1[0]) + '\n'
                                    result += 'Taxa name: ' + str(name1[1]) + '\n'
                                    result += 'Taxa ID: ' + str(name1[2]) + '\n'
                                    if DepVar == 1:
                                        result = result + 'Dependent Variable: Abundance' + '\n'
                                    elif DepVar == 2:
                                        result = result + 'Dependent Variable: Species Richness' + '\n'
                                    elif DepVar == 3:
                                        result = result + 'Dependent Variable: Species Diversity' + '\n'
                                    result += '\nANCOVA table:\n'
                                    result += str(D) + '\n'
                                    result += '===============================================\n'
                                    result += '\n'

                                    regrList = []
                                    regrList.append([float(min(x)), min_y])
                                    regrList.append([float(max(x)), max_y])

                                    seriesDict = {}
                                    seriesDict['turboThreshold'] = 0
                                    seriesDict['type'] = 'scatter'
                                    seriesDict['name'] = str(name1[1]) + ": " + str(name2)
                                    seriesDict['data'] = dataList
                                    seriesDict['color'] = colors[colors_idx]

                                    markerDict = {}
                                    markerDict['symbol'] = shapes[shapes_idx]
                                    seriesDict['marker'] = markerDict
                                    seriesDict['data'] = dataList
                                    seriesList.append(seriesDict)

                                    regrDict = {}
                                    regrDict['type'] = 'line'
                                    regrDict['name'] = 'y = ' + str(slope) + 'x' + ' + ' + str(intercept)
                                    regrDict['data'] = regrList
                                    regrDict['color'] = colors[colors_idx]

                                    markerDict = {}
                                    markerDict['enabled'] = False
                                    regrDict['marker'] = markerDict

                                    seriesList.append(regrDict)

                            shapes_idx += 1
                            if shapes_idx >= len(shapes):
                                shapes_idx = 0

                else:
                    obs = resultDF['abund'].loc[resultDF['abund'] > 0].count()
                    if obs < 3:
                        result += '\n===============================================\n\n'
                        result += 'Taxa Rank: ' + str(name1[0]) + '; Taxa Name: ' + str(name1[1]) + '; Taxa ID: ' + str(name1[2]) + ' has ' + str(obs) + ' observation(s),\n'
                        result += 'it was removed from your analysis.\n'
                        result += '\n===============================================\n\n'
                    else:
                        if DepVar == 1:
                            dataList = resultDF[[str(fieldListQuant[0]), 'abund']].values.astype(float).tolist()
                            x = resultDF[fieldListQuant[0]].astype(float).values.tolist()
                            y = resultDF['abund'].astype(float).values.tolist()
                        elif DepVar == 2:
                            dataList = resultDF[[fieldListQuant[0], 'rich']].values.astype(float).tolist()
                            x = resultDF[fieldListQuant[0]].astype(float).values.tolist()
                            y = resultDF['rich'].astype(float).values.tolist()
                        elif DepVar == 3:
                            dataList = resultDF[[fieldListQuant[0], 'diversity']].values.astype(float).tolist()
                            x = resultDF[fieldListQuant[0]].astype(float).values.tolist()
                            y = resultDF['diversity'].astype(float).values.tolist()

                        slp, inter, r_value, p, std_err = stats.linregress(x, y)
                        min_y = float(slp*min(x) + inter)
                        max_y = float(slp*max(x) + inter)
                        slope = "%0.3f" % slp
                        intercept = "%0.3f" % inter

                        if sig_only == 0:
                            result += 'Taxa level: ' + str(name1[0]) + '\n'
                            result += 'Taxa name: ' + str(name1[1]) + '\n'
                            result += 'Taxa ID: ' + str(name1[2]) + '\n'
                            if DepVar == 1:
                                result = result + 'Dependent Variable: Abundance' + '\n'
                            elif DepVar == 2:
                                result = result + 'Dependent Variable: Species Richness' + '\n'
                            elif DepVar == 3:
                                result = result + 'Dependent Variable: Species Diversity' + '\n'
                            result += '\nANCOVA table:\n'
                            result += str(D) + '\n'
                            result += '===============================================\n'
                            result += '\n'

                            regrList = []
                            regrList.append([float(min(x)), min_y])
                            regrList.append([float(max(x)), max_y])

                            seriesDict = {}
                            seriesDict['turboThreshold'] = 0
                            seriesDict['type'] = 'scatter'
                            seriesDict['name'] = str(name1[1])
                            seriesDict['data'] = dataList
                            seriesDict['color'] = colors[colors_idx]

                            markerDict = {}
                            markerDict['symbol'] = shapes[shapes_idx]
                            seriesDict['marker'] = markerDict
                            seriesDict['data'] = dataList
                            seriesList.append(seriesDict)

                            regrDict = {}
                            regrDict['type'] = 'line'
                            regrDict['name'] = 'y = ' + str(slope) + 'x' + ' + ' + str(intercept)
                            regrDict['data'] = regrList
                            regrDict['color'] = colors[colors_idx]

                            markerDict = {}
                            markerDict['enabled'] = False
                            regrDict['marker'] = markerDict

                            seriesList.append(regrDict)

                        if sig_only == 1:
                            if p_value >= 0.05:
                                result += '\n===============================================\n\n'
                                result += 'Taxa Rank: ' + str(name1[0]) + '; Taxa Name: ' + str(name1[1]) + '; Taxa ID: ' + str(name1[2]) + ' is not significant,\n'
                                result += 'it was removed from your analysis.\n'
                                result += '\n===============================================\n\n'
                            if p_value < 0.05:
                                result += 'Taxa level: ' + str(name1[0]) + '\n'
                                result += 'Taxa name: ' + str(name1[1]) + '\n'
                                result += 'Taxa ID: ' + str(name1[2]) + '\n'
                                if DepVar == 1:
                                    result = result + 'Dependent Variable: Abundance' + '\n'
                                elif DepVar == 2:
                                    result = result + 'Dependent Variable: Species Richness' + '\n'
                                elif DepVar == 3:
                                    result = result + 'Dependent Variable: Species Diversity' + '\n'
                                result += '\nANCOVA table:\n'
                                result += str(D) + '\n'
                                result += '===============================================\n'
                                result += '\n'

                                regrList = []
                                regrList.append([float(min(x)), min_y])
                                regrList.append([float(max(x)), max_y])

                                seriesDict = {}
                                seriesDict['turboThreshold'] = 0
                                seriesDict['type'] = 'scatter'
                                seriesDict['name'] = str(name1[1])
                                seriesDict['data'] = dataList
                                seriesDict['color'] = colors[colors_idx]

                                markerDict = {}
                                markerDict['symbol'] = shapes[shapes_idx]
                                seriesDict['marker'] = markerDict
                                seriesDict['data'] = dataList
                                seriesList.append(seriesDict)

                                regrDict = {}
                                regrDict['type'] = 'line'
                                regrDict['name'] = 'y = ' + str(slope) + 'x' + ' + ' + str(intercept)
                                regrDict['data'] = regrList
                                regrDict['color'] = colors[colors_idx]

                                markerDict = {}
                                markerDict['enabled'] = False
                                regrDict['marker'] = markerDict

                                seriesList.append(regrDict)

                base[RID] = 'Step 3 of 6: Performing linear regression...done!'
                base[RID] = 'Step 4 of 6: Formatting graph data for display...'


                xTitle = {}
                xTitle['text'] = fieldListQuant[0]
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

