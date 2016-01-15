import datetime
from django.http import HttpResponse
from django.db.models import Sum
import logging
import numpy as np
import pandas as pd
import pickle
from pyper import *
import simplejson

from database.export.export_DF import UnivMetaDF, normalizeUniv
from database.models import Sample, Profile
from database.utils import multidict, taxaProfileDF, stoppableThread


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}
LOG_FILENAME = 'error_log.txt'


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


thread3 = stoppableThread()
stop3 = False


def stopExport(request):
    global res, thread3, stop3
    if request.is_ajax():
        stop3 = True
        try:
            thread3.terminate()
            thread3.join()
            myDict = {}
            myDict['error'] = 'Your analysis has been stopped!'
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')
        except:
            pass


def getExCatData(request):
    global res, thread3, stop3
    if request.is_ajax():
        thread3 = stoppableThread(target=loopCat, args=(request,))
        thread3.start()
        thread3.join()
        stop3 = False
        return HttpResponse(res, content_type='application/json')


def loopCat(request):
    global res, base, stage, time1, TimeDiff, stop3
    try:
        while True:
            # Get selected samples from cookie and query database for sample info
            samples = Sample.objects.all()
            samples.query = pickle.loads(request.session['selected_samples'])
            selected = samples.values_list('sampleid')
            qs1 = Sample.objects.all().filter(sampleid__in=selected)

            if request.is_ajax():
                # Get variables from web page
                allJson = request.GET["all"]
                all = simplejson.loads(allJson)

                RID = str(all["RID"])
                time1[RID] = time.time()  # Moved these down here so RID is available
                base[RID] = 'Step 1 of 4: Querying database...'

                selectAll = int(all["selectAll"])
                NormMeth = int(all["NormMeth"])
                Iters = int(all["Iters"])
                NormVal = all["NormVal"]
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
                metaStrCat = all["metaValsCat"]
                fieldListCat = []
                valueListCat = []
                try:
                    metaDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStrCat)
                    for key in sorted(metaDictCat):
                        fieldListCat.append(key)
                        valueListCat.append(metaDictCat[key])

                    idStrCat = all["metaIDsCat"]
                except:
                    placeholder = ''

                metaStrQuant = all["metaValsQuant"]
                fieldListQuant = []
                valueListQuant = []
                try:
                    metaDictQuant = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStrQuant)
                    for key in sorted(metaDictQuant):
                        fieldListQuant.append(key)
                        valueListQuant.extend(metaDictQuant[key])

                    idStrQuant = all["metaIDsQuant"]
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
                metaDF = metaDF.ix[newList]
                metaDF.dropna(inplace=True)

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

                base[RID] = 'Step 1 of 4: Querying database...done!'
                base[RID] = 'Step 2 of 4: Normalizing data...'

                normDF, DESeq_error = normalizeUniv(taxaDF, taxaDict, myList, NormMeth, NormReads, metaDF, Iters)

                finalDict = {}
                normDF.set_index('sampleid', inplace=True)

                finalDF = pd.merge(metaDF, normDF, left_index=True, right_index=True)
                #finalDF['abund'] = finalDF['abund'].div(finalDF['abund'].groupby(finalDF.index).sum())

                finalDF[['abund', 'rich', 'diversity']] = finalDF[['abund', 'rich', 'diversity']].astype(float)
                finalDF.reset_index(drop=False, inplace=True)
                finalDF.rename(columns={'index': 'sampleid'}, inplace=True)

                base[RID] = 'Step 2 of 4: Normalizing data...done!'
                base[RID] = 'Step 3 of 4: Formatting biome data...'

                biome = {}
                if not u'sample_name' in fieldList:
                    newList = ['sampleid', 'sample_name']
                    newList.extend(fieldList)
                else:
                    newList = ['sampleid']
                    newList.extend(fieldList)

                grouped = finalDF.groupby(newList, sort=False)
                nameList = []
                for name, group in grouped:
                    metaDict = {}
                    for i in xrange(1, len(newList)):
                        metaDict[str(newList[i])] = str(name[i])
                    nameList.append({"id": str(name[0]), "metadata": metaDict})

                if key == 'Kingdom':
                    grouped = finalDF.groupby(['rank', 'kingdomName', 'kingdomid'], sort=False)
                elif key == 'Phyla':
                    grouped = finalDF.groupby(['rank', 'phylaName', 'phylaid'], sort=False)
                elif key == 'Class':
                    grouped = finalDF.groupby(['rank', 'className', 'classid'], sort=False)
                elif key == 'Order':
                    grouped = finalDF.groupby(['rank', 'orderName', 'orderid'], sort=False)
                elif key == 'Family':
                    grouped = finalDF.groupby(['rank', 'familyName', 'familyid'], sort=False)
                elif key == 'Genus':
                    grouped = finalDF.groupby(['rank', 'genusName', 'genusid'], sort=False)
                elif key == 'Species':
                    grouped = finalDF.groupby(['rank', 'speciesName', 'speciesid'], sort=False)

                taxaList = []
                dataList = []
                for name, group in grouped:
                    metaDict ={}
                    taxonList = []
                    taxonList.append(str(name[1]))
                    metaDict['taxonomy'] = taxonList
                    taxaList.append({"id": str(name[2]), "metadata": metaDict})
                    dataList.append(group['abund'].values.astype(float).tolist())

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

                base[RID] = 'Step 3 of 4: Formatting biome data...done!'
                base[RID] = 'Step 4 of 4: Formatting result table...'

                res_table = finalDF.to_html(classes="table display")
                res_table = res_table.replace('border="1"', 'border="0"')
                finalDict['res_table'] = str(res_table)


                biome_json = simplejson.dumps(biome, ensure_ascii=True, indent=4, sort_keys=True)
                finalDict['error'] = 'none'
                finalDict['biome'] = str(biome_json)
                res = simplejson.dumps(finalDict)
                base[RID] = 'Step 4 of 4: Formatting result table...done!'
                return None
    except:
        if not stop3:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with Differential Abundance!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            raise

    #finally:
    #    return None



