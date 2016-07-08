import ast
import datetime
from django import db
from django.http import HttpResponse
from django_pandas.io import read_frame
import logging
import math
import multiprocessing as mp
from natsort import natsorted
import numpy as np
import pandas as pd
from pyper import *
from PyPDF2 import PdfFileReader, PdfFileMerger
import simplejson
import shutil
import threading

from database.models import PICRUSt
from database.models import ko_lvl1, ko_lvl2, ko_lvl3, ko_entry
from database.models import nz_lvl1, nz_lvl2, nz_lvl3, nz_lvl4, nz_entry
from database.utils import multidict
import database.queue
from config import local_cfg


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}

LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def statusSpAC(request):
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
            if TimeDiff[RID] == 0:
                stage[RID] = 'Analysis has been placed in queue, there are '+str(database.queue.stat(RID))+' others in front of you.'
            else:
                stage[RID] = str(base[RID]) + '<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
        except:
            if TimeDiff[RID] == 0:
                stage[RID] = 'In queue'
            else:
                stage[RID] = '<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
        myDict['stage'] = stage[RID]
        json_data = simplejson.dumps(myDict, encoding="Latin-1")
        return HttpResponse(json_data, content_type='application/json')


def removeRIDSpAC(RID):
    global base, stage, time1, time2, TimeDiff
    try:
        base.pop(RID, None)
        stage.pop(RID, None)
        time1.pop(RID, None)
        time2.pop(RID, None)
        TimeDiff.pop(RID, None)
        return True
    except:
        return False


def getSpAC(request, stops, RID, PID):
    global base, stage, time1, TimeDiff
    try:
        while True:
            if request.is_ajax():
                allJson = request.GET["all"]
                all = simplejson.loads(allJson)

                time1[RID] = time.time()
                base[RID] = 'Step 1 of 4: Selecting your chosen meta-variables...'

                myDir = 'media/usr_temp/' + str(request.user) + '/'
                path = str(myDir) + 'usr_norm_data.csv'

                with open(path, 'rb') as f:
                    savedDF = pd.read_csv(f, index_col=0, sep='\t')

                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])

                result = ''
                button3 = int(all['button3'])
                if button3 == 1:
                    if selectAll == 1:
                        result += 'Taxa level: Kingdom' + '\n'
                    elif selectAll == 2:
                        result += 'Taxa level: Phyla' + '\n'
                    elif selectAll == 3:
                        result += 'Taxa level: Class' + '\n'
                    elif selectAll == 4:
                        result += 'Taxa level: Order' + '\n'
                    elif selectAll == 5:
                        result += 'Taxa level: Family' + '\n'
                    elif selectAll == 6:
                        result += 'Taxa level: Genus' + '\n'
                    elif selectAll == 7:
                        result += 'Taxa level: Species' + '\n'
                elif button3 == 2:
                    if keggAll == 1:
                        result += 'KEGG Pathway level: 1' + '\n'
                    elif keggAll == 2:
                        result += 'KEGG Pathway level: 2' + '\n'
                    elif keggAll == 3:
                        result += 'KEGG Pathway level: 3' + '\n'
                elif button3 == 3:
                    if nzAll == 1:
                        result += 'KEGG Enzyme level: 1' + '\n'
                    elif nzAll == 2:
                        result += 'KEGG Enzyme level: 2' + '\n'
                    elif nzAll == 3:
                        result += 'KEGG Enzyme level: 3' + '\n'
                    elif nzAll == 4:
                        result += 'KEGG Enzyme level: 4' + '\n'
                    elif keggAll == 5:
                        result += 'KEGG Enzyme level: GIBBs' + '\n'
                    elif keggAll == 6:
                        result += 'KEGG Enzyme level: Nitrogen cycle' + '\n'

                # Select samples and meta-variables from savedDF
                metaValsCat = all['metaValsCat']
                metaIDsCat = all['metaIDsCat']

                metaDictCat = {}
                catFields = []
                catValues = []
                if metaValsCat:
                    metaDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaValsCat)
                    for key in sorted(metaDictCat):
                        catFields.append(key)
                        catValues.extend(metaDictCat[key])

                catFields_edit = []
                removed = []
                for i in metaDictCat:
                    levels = len(set(metaDictCat[i]))
                    if levels > 1:
                        catFields_edit.append(i)
                    else:
                        removed.append(i)

                catSampleIDs = []
                if metaIDsCat:
                    idDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaIDsCat)
                    for key in sorted(idDictCat):
                        catSampleIDs.extend(idDictCat[key])

                if not catFields_edit:
                    catSampleIDs = savedDF['sampleid'].tolist()

                allSampleIDs = catSampleIDs
                allFields = catFields_edit

                # Removes samples (rows) that are not in our samplelist
                tempDF = savedDF.loc[savedDF['sampleid'].isin(allSampleIDs)]

                if metaDictCat:
                    for key in metaDictCat:
                        tempDF = tempDF.loc[tempDF[key].isin(metaDictCat[key])]

                wantedList = allFields + ['sampleid']
                metaDF = tempDF[wantedList]

                result += 'Categorical variables selected by user: ' + ", ".join(catFields) + '\n'
                result += 'Categorical variables removed from analysis (contains only 1 level): ' + ", ".join(removed) + '\n'
                result += '===============================================\n'

                base[RID] = 'Step 1 of 4: Selecting your chosen meta-variables...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 2 of 4: Selecting your chosen taxa or KEGG level...'

                DepVar = 1
                finalDF = pd.DataFrame()
                if button3 == 1:
                    DepVar = int(all["DepVar_taxa"])
                    finalDF = getTaxaDF(selectAll, savedDF, metaDF, allFields, DepVar, RID, stops, PID)

                if button3 == 2:
                    DepVar = int(all["DepVar_kegg"])
                    finalDF = getKeggDF(keggAll, savedDF, tempDF, allFields, DepVar, RID, stops, PID)

                if button3 == 3:
                    DepVar = int(all["DepVar_nz"])
                    finalDF = getNZDF(nzAll, savedDF, tempDF, allFields, DepVar, RID, stops, PID)

                # save location info to session
                myDir = 'media/temp/spac/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                count_rDF = pd.DataFrame()
                if DepVar == 1:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='abund')
                elif DepVar == 2:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='rich')
                elif DepVar == 3:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='diversity')
                elif DepVar == 4:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='abund_16S')

                meta_rDF = savedDF.drop_duplicates(subset='sampleid', take_last=True)

                # Removes samples (rows) that are not in our samplelist
                meta_rDF = meta_rDF.loc[meta_rDF['sampleid'].isin(catSampleIDs)]

                if metaDictCat:
                    for key in metaDictCat:
                        meta_rDF = meta_rDF.loc[meta_rDF[key].isin(metaDictCat[key])]

                wantedList = allFields + ['sampleid', 'sample_name']
                meta_rDF = meta_rDF[wantedList]
                meta_rDF.set_index('sampleid', drop=True, inplace=True)

                base[RID] = 'Step 2 of 4: Selecting your chosen taxa...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 3 of 4: Calculating Species Accumulation Curves...'

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                path = os.path.join('media', 'temp', 'spac', 'Rplots', RID)
                if not os.path.exists(path):
                    os.makedirs(path)

                r.assign("path", path)
                r("setwd(path)")
                r("options('width'=5000)")
                r.assign("RID", RID)
                r("library(vegan)")

                method = all['method']
                r.assign('method', method)

                r("pdf_counter <- 1")
                if allFields:
                    grouped = meta_rDF.groupby(allFields)
                    for name, group in grouped:
                        print name
                        if isinstance(name, tuple):
                            name = "; ".join(name)
                        r.assign('name', name)
                        IDs = group.index.values.tolist()
                        data = count_rDF.ix[IDs]
                        r.assign("data", data)

                        r("pdf(paste('SpAC_temp', pdf_counter, '.pdf', sep=''))")
                        if method == 'estaccumR':
                            r("x <- estaccumR(data)")
                        if method == 'poolaccum':
                            r("x <- poolaccum(data)")

                        values = r("summary(x)")
                        values = values.lstrip('try({summary(x)})')
                        result += str(name) + '\n'
                        values = values.lstrip('try({y})')
                        result += str(values) + '\n'
                        result += '===============================================\n'

                        nSamps, cols = group.shape
                        if nSamps > 2:
                            r("plot(x, main=paste(name))")
                        else:
                            r("plot(0:10, type='n', axes=FALSE, bty='n', xlab='', ylab='', \
                                main=paste(name, '\nInsufficient data to generate plot (n < 3)!', sep='')) \
                            ")
                        r("dev.off()")
                        r("pdf_counter <- pdf_counter + 1")

                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                        if stops[PID] == RID:
                            res = ''
                            return HttpResponse(res, content_type='application/json')
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                else:
                    r.assign("data", count_rDF)

                    r("pdf(paste('SpAC_temp', pdf_counter, '.pdf', sep=''))")
                    if method == 'estaccumR':
                        r("x <- estaccumR(data)")
                    if method == 'poolaccum':
                        r("x <- poolaccum(data)")

                    values = r("summary(x)")
                    values = values.lstrip('try({summary(x)})')
                    result += str(values) + '\n'

                    r("plot(x)")
                    r("dev.off()")
                    r("pdf_counter <- pdf_counter + 1")

                base[RID] = 'Step 3 of 4: Calculating Species Accumulation Curves...done!'

                base[RID] = 'Step 4 of 4: Pooling pdf files for display...'

                # Combining Pdf files
                finalFile = 'media/temp/spac/Rplots/' + str(RID) + '/SpAC_final.pdf'

                pdf_files = [f for f in os.listdir(path) if f.endswith("pdf")]
                pdf_files = natsorted(pdf_files, key=lambda y: y.lower())

                merger = PdfFileMerger()
                for filename in pdf_files:
                    merger.append(PdfFileReader(os.path.join(path, filename), 'rb'))

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                merger.write(finalFile)

                finalDict = {}
                finalDict['text'] = result

                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                removeRIDSpAC(RID)
                return HttpResponse(res, content_type='application/json')

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "Error with SpAC!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = simplejson.dumps(myDict)
            removeRIDSpAC(RID)
            return HttpResponse(res, content_type='application/json')


def getTaxaDF(selectAll, savedDF, metaDF, allFields, DepVar, RID, stops, PID):
    global base
    try:
        base[RID] = 'Step 2 of 4: Selecting your chosen taxa or KEGG level...'
        taxaDF = pd.DataFrame(columns=['sampleid', 'rank', 'rank_id', 'rank_name', 'abund', 'abund_16S', 'rich', 'diversity'])

        if selectAll == 2:
            taxaDF = savedDF.loc[:, ['sampleid', 'phylaid', 'phylaName', 'abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'phylaid': 'rank_id', 'phylaName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Phyla'
        elif selectAll == 3:
            taxaDF = savedDF.loc[:, ['sampleid', 'classid', 'className', 'abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'classid': 'rank_id', 'className': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Class'
        elif selectAll == 4:
            taxaDF = savedDF.loc[:, ['sampleid', 'orderid', 'orderName', 'abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'orderid': 'rank_id', 'orderName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Order'
        elif selectAll == 5:
            taxaDF = savedDF.loc[:, ['sampleid', 'familyid', 'familyName', 'abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'familyid': 'rank_id', 'familyName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Family'
        elif selectAll == 6:
            taxaDF = savedDF.loc[:, ['sampleid', 'genusid', 'genusName', 'abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'genusid': 'rank_id', 'genusName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Genus'
        elif selectAll == 7:
            taxaDF = savedDF.loc[:, ['sampleid', 'speciesid', 'speciesName', 'abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'speciesid': 'rank_id', 'speciesName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Species'

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        taxaDF.drop('sampleid', axis=1, inplace=True)
        finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')

        wantedList = allFields + ['sampleid', 'rank', 'rank_name', 'rank_id']
        if DepVar == 1:
            finalDF = finalDF.groupby(wantedList)[['abund']].sum()
        elif DepVar == 2:
            finalDF = finalDF.groupby(wantedList)[['rich']].sum()
        elif DepVar == 3:
            finalDF = finalDF.groupby(wantedList)[['diversity']].sum()
        elif DepVar == 4:
            finalDF = finalDF.groupby(wantedList)[['abund_16S']].sum()

        finalDF.reset_index(drop=False, inplace=True)
        return finalDF

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with SpAC!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getKeggDF(keggAll, savedDF, tempDF, allFields, DepVar, RID, stops, PID):
    global base
    try:
        base[RID] = 'Step 2 of 8: Selecting your chosen taxa or KEGG level...'
        koDict = {}
        if keggAll == 1:
            keys = ko_lvl1.objects.using('picrust').values_list('ko_lvl1_id', flat=True)
            for key in keys:
                koList = ko_entry.objects.using('picrust').filter(ko_lvl1_id_id=key).values_list('ko_orthology', flat=True)
                if koList:
                    koDict[key] = koList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif keggAll == 2:
            keys = ko_lvl2.objects.using('picrust').values_list('ko_lvl2_id', flat=True)
            for key in keys:
                koList = ko_entry.objects.using('picrust').filter(ko_lvl2_id_id=key).values_list('ko_orthology', flat=True)
                if koList:
                    koDict[key] = koList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif keggAll == 3:
            keys = ko_lvl3.objects.using('picrust').values_list('ko_lvl3_id', flat=True)
            for key in keys:
                koList = ko_entry.objects.using('picrust').filter(ko_lvl3_id_id=key).values_list('ko_orthology', flat=True)
                if koList:
                    koDict[key] = koList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        if stops[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        # create sample and species lists based on meta data selection
        wanted = ['sampleid', 'speciesid', 'abund', 'abund_16S', 'rich', 'diversity']
        profileDF = tempDF.loc[:, wanted]
        profileDF.set_index('speciesid', inplace=True)

        # get PICRUSt data for species
        speciesList = pd.unique(profileDF.index.ravel().tolist())
        qs = PICRUSt.objects.using('picrust').filter(speciesid__in=speciesList)
        picrustDF = read_frame(qs, fieldnames=['speciesid__speciesid', 'geneCount'])
        picrustDF.set_index('speciesid__speciesid', inplace=True)

        path = 'media/temp/spac/' + str(RID)
        if not os.path.exists(path):
            os.makedirs(path)

        maxCPU = mp.cpu_count()
        maxCPU /= 3
        maxCPU = math.trunc(maxCPU)
        if maxCPU < 1:
            maxCPU = 1

        if os.name == 'nt':
            numcore = 1
            listDF = np.array_split(picrustDF, numcore)
            processes = [threading.Thread(target=sumStuff, args=(listDF[x], koDict, RID, x, stops, PID)) for x in xrange(numcore)]
        else:
            numcore = min(local_cfg.usr_numcore, maxCPU)
            listDF = np.array_split(picrustDF, numcore)
            processes = [threading.Thread(target=sumStuff, args=(listDF[x], koDict, RID, x, stops, PID)) for x in xrange(numcore)]

        for p in processes:
            p.start()

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        for p in processes:
            p.join()

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        levelList = []
        for key in koDict:
            levelList.append(key)

        picrustDF = pd.DataFrame()
        for i in xrange(numcore):
            path = 'media/temp/spac/'+str(RID)+'/file%d.temp' % i
            frame = pd.read_csv(path)
            picrustDF = picrustDF.append(frame, ignore_index=True)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        shutil.rmtree('media/temp/spac/'+str(RID))
        picrustDF.set_index('speciesid', inplace=True)
        picrustDF[picrustDF > 0] = 1

        # merge to get final gene counts for all selected samples
        taxaDF = pd.merge(profileDF, picrustDF, left_index=True, right_index=True, how='inner')

        for level in levelList:
            if DepVar == 1:
                taxaDF[level] = taxaDF['abund'] * taxaDF[level]
            elif DepVar == 2:
                taxaDF[level] = np.where(taxaDF['abund'] * taxaDF[level] > 0, 1, 0)
            elif DepVar == 3:
                taxaDF[level] = taxaDF['abund'] * taxaDF[level]
                taxaDF[level] = taxaDF[level].div(taxaDF[level].sum(), axis=0)
                taxaDF[level] = taxaDF[level].apply(lambda x: -1 * x * math.log(x) if x > 0 else 0)
            elif DepVar == 4:
                taxaDF[level] = taxaDF['abund_16S'] * taxaDF[level]

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        taxaDF = taxaDF.groupby('sampleid')[levelList].agg('sum')
        taxaDF.reset_index(drop=False, inplace=True)

        if DepVar == 1:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund')
        elif DepVar == 4:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund_16S')
        elif DepVar == 2:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rich')
        elif DepVar == 3:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='diversity')

        wanted = allFields + ['sampleid']
        metaDF = savedDF.loc[:, wanted]
        metaDF.set_index('sampleid', drop=True, inplace=True)
        grouped = metaDF.groupby(level=0)
        metaDF = grouped.last()

        taxaDF.set_index('sampleid', drop=True, inplace=True)
        finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')
        finalDF.reset_index(drop=False, inplace=True)

        finalDF['rank'] = ''
        finalDF['rank_name'] = ''
        for index, row in finalDF.iterrows():
            if ko_lvl1.objects.using('picrust').filter(ko_lvl1_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl1'
                finalDF.loc[index, 'rank_name'] = ko_lvl1.objects.using('picrust').get(ko_lvl1_id=row['rank_id']).ko_lvl1_name
            elif ko_lvl2.objects.using('picrust').filter(ko_lvl2_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl2'
                finalDF.loc[index, 'rank_name'] = ko_lvl2.objects.using('picrust').get(ko_lvl2_id=row['rank_id']).ko_lvl2_name
            elif ko_lvl3.objects.using('picrust').filter(ko_lvl3_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl3'
                finalDF.loc[index, 'rank_name'] = ko_lvl3.objects.using('picrust').get(ko_lvl3_id=row['rank_id']).ko_lvl3_name
            elif ko_entry.objects.using('picrust').filter(ko_lvl4_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl4'
                finalDF.loc[index, 'rank_name'] = ko_entry.objects.using('picrust').get(ko_lvl4_id=row['rank_id']).ko_name

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        return finalDF

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with SpAC!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getNZDF(nzAll, savedDF, tempDF, allFields, DepVar, RID, stops, PID):
    global base
    try:
        nzDict = {}
        if nzAll == 1:
            keys = nz_lvl1.objects.using('picrust').values_list('nz_lvl1_id', flat=True)
            for key in keys:
                nzList = nz_entry.objects.using('picrust').filter(nz_lvl1_id_id=key).values_list('nz_orthology', flat=True)
                if nzList:
                    nzDict[key] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 2:
            keys = nz_lvl2.objects.using('picrust').values_list('nz_lvl2_id', flat=True)
            for key in keys:
                nzList = nz_entry.objects.using('picrust').filter(nz_lvl2_id_id=key).values_list('nz_orthology', flat=True)
                if nzList:
                    nzDict[key] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 3:
            keys = nz_lvl3.objects.using('picrust').values_list('nz_lvl3_id', flat=True)
            for key in keys:
                nzList = nz_entry.objects.using('picrust').filter(nz_lvl3_id_id=key).values_list('nz_orthology', flat=True)
                if nzList:
                    nzDict[key] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 4:
            keys = nz_lvl4.objects.using('picrust').values_list('nz_lvl4_id', flat=True)
            for key in keys:
                nzList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=key).values_list('nz_orthology', flat=True)
                if nzList:
                    nzDict[key] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 5:
            # 1.18.6.1  nitrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.18.6.1  nitrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.3.3.11  pyrroloquinoline-quinone synthase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.3.3.11  pyrroloquinoline-quinone synthase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.4.99.5  glycine dehydrogenase (cyanide-forming)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.4.99.5  glycine dehydrogenase (cyanide-forming)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.1.1.76  (S,S)-butanediol dehydrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.1.1.76  (S,S)-butanediol dehydrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.14  chitinase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='3.2.1.14  chitinase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 4.1.1.74  indolepyruvate decarboxylase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='4.1.1.74  indolepyruvate decarboxylase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 6.3.2.39  aerobactin synthase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='6.3.2.39  aerobactin synthase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.4  cellulase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.4  cellulase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.91  cellulose 1,4-beta-cellobiosidase (non-reducing end)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.91  cellulose 1,4-beta-cellobiosidase (non-reducing end)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.21  beta-glucosidase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.21  beta-glucosidase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.8  endo-1,4-beta-xylanase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.8  endo-1,4-beta-xylanase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.37  xylan 1,4-beta-xylosidase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.37  xylan 1,4-beta-xylosidase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.5.1.4  amidase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.5.1.4  amidase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.5.1.5  urease
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.5.1.5  urease').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.1.3.1  alkaline phosphatase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.1.3.1  alkaline phosphatase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.1.3.2  acid phosphatase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.1.3.2  acid phosphatase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.1.6.1  arylsulfatase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.1.6.1  arylsulfatase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

        elif nzAll == 6:
            # 1.18.6.1  nitrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.18.6.1  nitrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.5.1.5  urease
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='3.5.1.5  urease').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.14.99.39  ammonia monooxygenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.14.99.39  ammonia monooxygenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.6  hydroxylamine dehydrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.6  hydroxylamine dehydrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.99.4  nitrate reductase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.99.4  nitrate reductase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.1  nitrite reductase (NO-forming)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.1  nitrite reductase (NO-forming)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.5  nitric oxide reductase (cytochrome c)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.5  nitric oxide reductase (cytochrome c)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.4  nitrous-oxide reductase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.4  nitrous-oxide reductase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

        # create sample and species lists based on meta data selection
        wanted = ['sampleid', 'speciesid', 'abund', 'abund_16S', 'rich', 'diversity']
        profileDF = tempDF.loc[:, wanted]
        profileDF.set_index('speciesid', inplace=True)

        # get PICRUSt data for species
        speciesList = pd.unique(profileDF.index.ravel().tolist())
        qs = PICRUSt.objects.using('picrust').filter(speciesid__in=speciesList)
        picrustDF = read_frame(qs, fieldnames=['speciesid__speciesid', 'geneCount'])
        picrustDF.set_index('speciesid__speciesid', inplace=True)

        path = 'media/temp/spac/' + str(RID)
        if not os.path.exists(path):
            os.makedirs(path)

        maxCPU = mp.cpu_count()
        maxCPU /= 3
        maxCPU = math.trunc(maxCPU)
        if maxCPU < 1:
            maxCPU = 1

        if os.name == 'nt':
            numcore = 1
            listDF = np.array_split(picrustDF, numcore)
            processes = [threading.Thread(target=sumStuff, args=(listDF[x], nzDict, RID, x, stops, PID)) for x in xrange(numcore)]
        else:
            numcore = min(local_cfg.usr_numcore, maxCPU)
            listDF = np.array_split(picrustDF, numcore)
            processes = [threading.Thread(target=sumStuff, args=(listDF[x], nzDict, RID, x, stops, PID)) for x in xrange(numcore)]

        for p in processes:
            p.start()

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        for p in processes:
            p.join()

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        levelList = []
        for key in nzDict:
            levelList.append(key)

        picrustDF = pd.DataFrame()
        for i in xrange(numcore):
            path = 'media/temp/spac/'+str(RID)+'/file%d.temp' % i
            frame = pd.read_csv(path)
            picrustDF = picrustDF.append(frame, ignore_index=True)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        shutil.rmtree('media/temp/spac/'+str(RID))
        picrustDF.set_index('speciesid', inplace=True)
        picrustDF[picrustDF > 0] = 1

        # merge to get final gene counts for all selected samples
        taxaDF = pd.merge(profileDF, picrustDF, left_index=True, right_index=True, how='inner')
        for level in levelList:
            if DepVar == 1:
                taxaDF[level] = taxaDF['abund'] * taxaDF[level]
            elif DepVar == 2:
                taxaDF[level] = np.where(taxaDF['abund'] * taxaDF[level] > 0, 1, 0)
            elif DepVar == 3:
                taxaDF[level] = taxaDF['abund'] * taxaDF[level]
                taxaDF[level] = taxaDF[level].div(taxaDF[level].sum(), axis=0)
                taxaDF[level] = taxaDF[level].apply(lambda x: -1 * x * math.log(x) if x > 0 else 0)
            elif DepVar == 4:
                taxaDF[level] = taxaDF['abund_16S'] * taxaDF[level]

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        taxaDF = taxaDF.groupby('sampleid')[levelList].agg('sum')
        taxaDF.reset_index(drop=False, inplace=True)

        if DepVar == 1:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund')
        elif DepVar == 4:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund_16S')
        elif DepVar == 2:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rich')
        elif DepVar == 3:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='diversity')

        wanted = allFields + ['sampleid']
        metaDF = savedDF.loc[:, wanted]
        metaDF.set_index('sampleid', drop=True, inplace=True)
        grouped = metaDF.groupby(level=0)
        metaDF = grouped.last()

        taxaDF.set_index('sampleid', drop=True, inplace=True)
        finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')
        finalDF.reset_index(drop=False, inplace=True)

        finalDF['rank'] = ''
        finalDF['rank_name'] = ''
        for index, row in finalDF.iterrows():
            if nz_lvl1.objects.using('picrust').filter(nz_lvl1_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl1'
                finalDF.loc[index, 'rank_name'] = nz_lvl1.objects.using('picrust').get(nz_lvl1_id=row['rank_id']).nz_lvl1_name
            elif nz_lvl2.objects.using('picrust').filter(nz_lvl2_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl2'
                finalDF.loc[index, 'rank_name'] = nz_lvl2.objects.using('picrust').get(nz_lvl2_id=row['rank_id']).nz_lvl2_name
            elif nz_lvl3.objects.using('picrust').filter(nz_lvl3_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl3'
                finalDF.loc[index, 'rank_name'] = nz_lvl3.objects.using('picrust').get(nz_lvl3_id=row['rank_id']).nz_lvl3_name
            elif nz_lvl4.objects.using('picrust').filter(nz_lvl4_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl4'
                finalDF.loc[index, 'rank_name'] = nz_lvl4.objects.using('picrust').get(nz_lvl4_id=row['rank_id']).nz_lvl4_name
            elif nz_entry.objects.using('picrust').filter(nz_lvl5_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl5'
                finalDF.loc[index, 'rank_name'] = nz_entry.objects.using('picrust').get(nz_lvl5_id=row['rank_id']).nz_name

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        return finalDF

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with spac!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def sumStuff(slice, koDict, RID, num, stops, PID):
    global base
    db.close_old_connections()

    f = open('media/temp/spac/'+str(RID)+'/file'+str(num)+".temp", 'w')

    keyList = []
    for key in koDict:
        keyList.append(key)

    f.write('speciesid,'+",".join(keyList)+'\n')

    for index, row in slice.iterrows():
        d = ast.literal_eval(row['geneCount'])

        f.write(str(index)+',')
        sumList = []
        for key in koDict:
            sum = 0.0
            myList = koDict[key]
            for k in myList:
                if k in d:
                    sum += d[k]
            sumList.append(sum)

        f.write(','.join(map(str, sumList)))
        f.write('\n')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    f.close()


def removeSpACFiles(request):
    if request.is_ajax():
        RID = request.GET["all"]

        file = "media/temp/spac/Rplots/" + str(RID) + ".spac.pdf"
        if os.path.exists(file):
            os.remove(file)

        file = "media/temp/spac/" + str(RID) + ".pkl"
        if os.path.exists(file):
            os.remove(file)

        file = "media/temp/spac/" + str(RID) + ".csv"
        if os.path.exists(file):
            os.remove(file)

        return HttpResponse()


def getTabSpAC(request):
    if request.is_ajax():
        RID = request.GET["all"]
        myDir = 'media/temp/spac/'
        fileName = str(myDir) + str(RID) + '.pkl'
        savedDF = pd.read_pickle(fileName)

        myDir = 'media/temp/spac/'
        fileName = str(myDir) + str(RID) + '.csv'
        savedDF.to_csv(fileName)

        myDict = {}
        myDir = 'temp/spac/'
        fileName = str(myDir) + str(RID) + '.csv'
        myDict['name'] = str(fileName)
        res = simplejson.dumps(myDict)

        return HttpResponse(res, content_type='application/json')
