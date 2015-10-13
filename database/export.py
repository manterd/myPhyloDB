__author__ = 'KorsaM'
import datetime
from django.http import HttpResponse
from django.db.models import Sum
import numpy as np
import pandas as pd
import pickle
from pyper import *
from scipy import stats
import simplejson

from database.anova.anova_DF import UnivMetaDF, normalizeUniv
from database.models import Sample, Profile
from database.utils import multidict, taxaProfileDF


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}


def statusExport(request):
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


def removeRIDExport(request):
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


def getExCatData(request):
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
                Iters = int(all["Iters"])
                NormVal = all["NormVal"]
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
                result = ''
                metaStrCat = all["metaValsCat"]
                fieldListCat = []
                valueListCat = []
                idDictCat = {}
                try:
                    metaDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStrCat)
                    for key in sorted(metaDictCat):
                        fieldListCat.append(key)
                        valueListCat.append(metaDictCat[key])

                    idStrCat = all["metaIDsCat"]
                    idDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(idStrCat)
                except:
                    placeholder = ''

                metaStrQuant = all["metaValsQuant"]
                fieldListQuant = []
                valueListQuant = []
                idDictQuant = {}
                try:
                    metaDictQuant = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStrQuant)
                    for key in sorted(metaDictQuant):
                        fieldListQuant.append(key)
                        valueListQuant.extend(metaDictQuant[key])

                    idStrQuant = all["metaIDsQuant"]
                    idDictQuant = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(idStrQuant)
                except:
                    placeholder = ''

                metaStr = all["metaVals"]
                fieldList = []
                valueList = []
                idDict = {}
                try:
                    metaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStr)
                    for key in sorted(metaDict):
                        fieldList.append(key)
                        valueList.append(metaDict[key])

                    idStr = all["metaIDs"]
                    idDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(idStr)

                except:
                    placeholder = ''

                result = result + 'Categorical variables selected: ' + ", ".join(fieldListCat) + '\n'
                result = result + 'Quantitative variables selected: ' + ", ".join(fieldListQuant) + '\n'
                result += '===============================================\n'
                result += '\nData Normalization:\n'

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

                metaDF = UnivMetaDF(idDict)

                lenA, col = metaDF.shape

                metaDF = metaDF.ix[newList]
                metaDF.dropna(inplace=True)
                lenB, col = metaDF.shape

                selectRem = len(selected) - lenA
                normRem = lenA - lenB

                result += str(lenB) + ' selected samples were included in the final analysis.\n'
                if normRem > 0:
                    result += str(normRem) + ' samples did not met the desired normalization criteria.\n'
                if selectRem:
                    result += str(selectRem) + ' samples were deselected by the user.\n'

                # Create unique list of samples in meta dataframe (may be different than selected samples)
                myList = metaDF.index.values.tolist()

                # Create dataframe with all taxa/count data by sample
                taxaDF = taxaProfileDF(myList)

                # Select only the taxa of interest if user used the taxa tree
                taxaString = all["taxa"]
                taxaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(taxaString)

                # Select only the taxa of interest if user used the selectAll button
                if selectAll == 1:
                    taxaDict = {}
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('kingdomid', flat='True').distinct()
                    taxaDict['Kingdom'] = qs3
                elif selectAll == 2:
                    taxaDict = {}
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('phylaid', flat='True').distinct()
                    taxaDict['Phyla'] = qs3
                elif selectAll == 3:
                    taxaDict = {}
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('classid', flat='True').distinct()
                    taxaDict['Class'] = qs3
                elif selectAll == 4:
                    taxaDict = {}
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('orderid', flat='True').distinct()
                    taxaDict['Order'] = qs3
                elif selectAll == 5:
                    taxaDict = {}
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('familyid', flat='True').distinct()
                    taxaDict['Family'] = qs3
                elif selectAll == 6:
                    taxaDict = {}
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('genusid', flat='True').distinct()
                    taxaDict['Genus'] = qs3
                elif selectAll == 7:
                    taxaDict = {}
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('speciesid', flat='True').distinct()
                    taxaDict['Species'] = qs3

                base[RID] = 'Step 1 of 6: Querying database...done!'
            except Exception as e:
                print("Error querying database (ANOVA CAT): ", e)

            # Normalize data
            try:
                base[RID] = 'Step 2 of 6: Normalizing data...'

                normDF, DESeq_error = normalizeUniv(taxaDF, taxaDict, myList, NormMeth, NormReads, metaDF, Iters)

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

                normDF.set_index('sampleid', inplace=True)

                finalDF = pd.merge(metaDF, normDF, left_index=True, right_index=True)
                finalDF[['abund', 'rich', 'diversity']] = finalDF[['abund', 'rich', 'diversity']].astype(float)

                base[RID] = 'Step 2 of 6: Normalizing data...done!'

            except Exception as e:
                print("Error with normalization (ANOVA CAT): ", e)


            # Done with normalization, send data back

    except Exception as e:
        print "Error with ANOVA: ", e
        state = "Error with ANOVA: " + str(e)

        myDict = {}
        myDict['error'] = state
        res = simplejson.dumps(myDict)
        return HttpResponse(res, content_type='application/json')


def getExQuantData(request):
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
            Iters = int(all["Iters"])
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
            result = ''
            metaStrCat = all["metaValsCat"]
            fieldListCat = []
            valueListCat = []
            idDictCat = {}
            try:
                metaDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStrCat)
                for key in sorted(metaDictCat):
                    fieldListCat.append(key)
                    valueListCat.append(metaDictCat[key])

                idStrCat = all["metaIDsCat"]
                idDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(idStrCat)
            except:
                placeholder = ''

            metaStrQuant = all["metaValsQuant"]
            fieldListQuant = []
            valueListQuant = []
            idDictQuant = {}
            try:
                metaDictQuant = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStrQuant)
                for key in sorted(metaDictQuant):
                    fieldListQuant.append(key)
                    valueListQuant.extend(metaDictQuant[key])

                idStrQuant = all["metaIDsQuant"]
                idDictQuant = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(idStrQuant)
            except:
                placeholder = ''

            metaStr = all["metaVals"]
            fieldList = []
            valueList = []
            idDict = {}
            try:
                metaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStr)
                for key in sorted(metaDict):
                    fieldList.append(key)
                    valueList.append(metaDict[key])

                idStr = all["metaIDs"]
                idDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(idStr)

            except:
                placeholder = ''

            result = result + 'Categorical variables selected: ' + ", ".join(fieldListCat) + '\n'
            result = result + 'Quantitative variables selected: ' + ", ".join(fieldListQuant) + '\n'
            result += '===============================================\n'
            result += '\nData Normalization:\n'

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

            metaDF = UnivMetaDF(idDict)

            lenA, col = metaDF.shape

            metaDF = metaDF.ix[newList]
            metaDF.dropna(inplace=True)
            lenB, col = metaDF.shape

            selectRem = len(selected) - lenA
            normRem = lenA - lenB

            result += str(lenB) + ' selected samples were included in the final analysis.\n'
            if normRem > 0:
                result += str(normRem) + ' samples did not met the desired normalization criteria.\n'
            if selectRem:
                result += str(selectRem) + ' samples were deselected by the user.\n'

            # Create unique list of samples in meta dataframe (may be different than selected samples)
            myList = metaDF.index.values.tolist()

            # Create dataframe with all taxa/count data by sample
            taxaDF = taxaProfileDF(myList)

            # Select only the taxa of interest if user used the taxa tree
            taxaString = all["taxa"]
            taxaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(taxaString)

            # Select only the taxa of interest if user used the selectAll button
            if selectAll == 1:
                taxaDict = {}
                qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('kingdomid', flat='True').distinct()
                taxaDict['Kingdom'] = qs3
            elif selectAll == 2:
                taxaDict = {}
                qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('phylaid', flat='True').distinct()
                taxaDict['Phyla'] = qs3
            elif selectAll == 3:
                taxaDict = {}
                qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('classid', flat='True').distinct()
                taxaDict['Class'] = qs3
            elif selectAll == 4:
                taxaDict = {}
                qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('orderid', flat='True').distinct()
                taxaDict['Order'] = qs3
            elif selectAll == 5:
                taxaDict = {}
                qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('familyid', flat='True').distinct()
                taxaDict['Family'] = qs3
            elif selectAll == 6:
                taxaDict = {}
                qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('genusid', flat='True').distinct()
                taxaDict['Genus'] = qs3
            elif selectAll == 7:
                taxaDict = {}
                qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('speciesid', flat='True').distinct()
                taxaDict['Species'] = qs3

            base[RID] = 'Step 1 of 6: Querying database...done!'
            base[RID] = 'Step 2 of 6: Normalizing data...'

            normDF, DESeq_error = normalizeUniv(taxaDF, taxaDict, myList, NormMeth, NormReads, metaDF, Iters)

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

            normDF.set_index('sampleid', inplace=True)

            finalDF = pd.merge(metaDF, normDF, left_index=True, right_index=True)
            finalDF[['abund', 'rich', 'diversity']] = finalDF[['abund', 'rich', 'diversity']].astype(float)

            base[RID] = 'Step 2 of 6: Normalizing data...done!'
            # Done normalizing, export!

    except Exception as e:
        print "Error with Linear Regression: ", e
        state = "Error with Linear Regression: " + str(e)

        myDict = {}
        myDict['error'] = state
        res = simplejson.dumps(myDict)
        return HttpResponse(res, content_type='application/json')



