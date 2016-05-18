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
import pickle
from pyper import *
from PyPDF2 import PdfFileReader, PdfFileMerger
import simplejson
import shutil
import threading

from database.models import PICRUSt
from database.models import nz_lvl1, nz_lvl2, nz_lvl3, nz_lvl4, nz_entry
from database.utils import multidict, stoppableThread


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}
stop_soil_index = False
stops = {}
thread_soil_index = stoppableThread()
res = ''
LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def statussoil_index(request):
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


def removeRIDsoil_index(request):
    global base, stage, time1, time2, TimeDiff
    try:
        if request.is_ajax():
            RID = request.GET["all"]
            base.pop(RID, None)
            stage.pop(RID, None)
            time1.pop(RID, None)
            time2.pop(RID, None)
            TimeDiff.pop(RID, None)
            stops.pop(RID, None)
            return True
        else:
            return False
    except:
        return False


def stopsoil_index(request):
    global thread_soil_index, stops, stop_soil_index, res
    if request.is_ajax():
        RID = request.GET["all"]
        stops[RID] = True
        stop_soil_index = True
        thread_soil_index.terminate()
        thread_soil_index.join()
        removeRIDsoil_index(request)

        res = ''
        myDict = {}
        myDict['error'] = 'none'
        myDict['message'] = 'Your analysis has been stopped!'
        stop = simplejson.dumps(myDict)
        return HttpResponse(stop, content_type='application/json')


def getsoil_index(request):
    global res, thread_soil_index, stop_soil_index
    if request.is_ajax():
        stop_soil_index = False
        thread_soil_index = stoppableThread(target=loopCat, args=(request,))
        thread_soil_index.start()
        thread_soil_index.join()
        removeRIDsoil_index(request)
        return HttpResponse(res, content_type='application/json')


def loopCat(request):
    global res, base, stage, time1, TimeDiff, stops, stop_soil_index
    try:
        while True:
            if request.is_ajax():
                allJson = request.GET["all"]
                all = simplejson.loads(allJson)
                RID = str(all["RID"])
                stops[RID] = False

                time1[RID] = time.time()
                base[RID] = 'Step 1 of X: Selecting your chosen meta-variables...'

                path = pickle.loads(request.session['savedDF'])
                savedDF = pd.read_pickle(path)

                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])

                result = ''
                button3 = int(all['button3'])

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

                base[RID] = 'Step 1 of X: Selecting your chosen meta-variables...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 2 of X: Selecting your chosen taxa or KEGG level...'

                finalDF, keggDF, levelList = getTaxaDF(savedDF, metaDF, RID)
                finalDF = finalDF[finalDF.rel_abund != 0]

                # save location info to session
                myDir = 'media/temp/soil_index/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                meta_rDF = savedDF.drop_duplicates(subset='sampleid', take_last=True)

                # Removes samples (rows) that are not in our samplelist
                meta_rDF = meta_rDF.loc[meta_rDF['sampleid'].isin(catSampleIDs)]

                if metaDictCat:
                    for key in metaDictCat:
                        meta_rDF = meta_rDF.loc[meta_rDF[key].isin(metaDictCat[key])]

                wantedList = allFields + ['sampleid', 'sample_name']
                meta_rDF = meta_rDF[wantedList]
                meta_rDF.set_index('sampleid', drop=True, inplace=True)

                rnaDF = finalDF.replace(r'\s+', np.nan, regex=True)
                rnaDF.fillna(0.0, axis=1, inplace=True)
                rnaDF['bins'] = 'bin1'
                rnaDF.ix[rnaDF.rRNACount == 1, 'bins'] = 'bin1'
                rnaDF.ix[(rnaDF.rRNACount > 1) & (rnaDF.rRNACount <= 2), 'bins'] = 'bin2'
                rnaDF.ix[(rnaDF.rRNACount > 2) & (rnaDF.rRNACount <= 3), 'bins'] = 'bin3'
                rnaDF.ix[(rnaDF.rRNACount > 3) & (rnaDF.rRNACount <= 4), 'bins'] = 'bin4'
                rnaDF.ix[(rnaDF.rRNACount > 4) & (rnaDF.rRNACount <= 5), 'bins'] = 'bin5'
                rnaDF.ix[(rnaDF.rRNACount > 5) & (rnaDF.rRNACount <= 6), 'bins'] = 'bin6'
                rnaDF.ix[(rnaDF.rRNACount > 6) & (rnaDF.rRNACount <= 7), 'bins'] = 'bin7'
                rnaDF.ix[(rnaDF.rRNACount > 7) & (rnaDF.rRNACount <= 8), 'bins'] = 'bin8'
                rnaDF.ix[(rnaDF.rRNACount > 8) & (rnaDF.rRNACount <= 9), 'bins'] = 'bin9'
                rnaDF.ix[(rnaDF.rRNACount > 9) & (rnaDF.rRNACount <= 10), 'bins'] = 'bin10'
                rnaDF.ix[rnaDF.rRNACount > 10, 'bins'] = 'bin11'
                myBins = list(set(rnaDF['bins'].tolist()))

                binDF = rnaDF.groupby(['sampleid', 'bins'])['rel_abund'].sum()
                binDF = binDF.unstack(level=-1)
                binDF.fillna(0.0, axis=1, inplace=True)
                binDF.reset_index(drop=False, inplace=True)
                sumDF1 = binDF.groupby('sampleid')[myBins].sum()

                sumDF2 = finalDF.groupby('sampleid')['abund_16S', 'rich', 'diversity'].sum()
                allDF = pd.merge(sumDF2, sumDF1, left_index=True, right_index=True, how='outer')

                kegg_rDF = keggDF.pivot(index='sampleid', columns='rank_name', values='rel_abund')
                allDF = pd.merge(allDF, kegg_rDF, left_index=True, right_index=True, how='outer')

                allDF = pd.merge(allDF, meta_rDF, left_index=True, right_index=True, how='outer')
                wantList = ['abund_16S', 'rich', 'diversity'] + myBins
                bytrt1 = allDF.groupby(catFields_edit)[wantList]
                #print bytrt1.mean()  # means by other means

                nzDF = keggDF.groupby(['sampleid', 'rank_name'])['rel_abund'].sum()
                nzDF = nzDF.unstack(level=-1)
                starDF = pd.merge(nzDF, meta_rDF, left_index=True, right_index=True, how='outer')
                bytrt2 = starDF.groupby(catFields_edit)

                myList = [u'1.1.1.76  (S,S)-butanediol dehydrogenase', u'1.18.6.1  nitrogenase', u'1.3.3.11  pyrroloquinoline-quinone synthase', u'1.4.99.5  glycine dehydrogenase (cyanide-forming)', u'3.1.3.1  alkaline phosphatase', u'3.1.3.2  acid phosphatase', u'3.1.6.1  arylsulfatase', u'3.2.1.14  chitinase', u'3.2.1.21  beta-glucosidase', u'3.2.1.37  xylan 1,4-beta-xylosidase', u'3.2.1.4  cellulase', u'3.2.1.8  endo-1,4-beta-xylanase', u'3.2.1.91  cellulose 1,4-beta-cellobiosidase (non-reducing end)', u'3.5.1.4  amidase', u'3.5.1.5  urease', u'3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase', u'4.1.1.74  indolepyruvate decarboxylase', u'6.3.2.39  aerobactin synthase']
                newList = [u'(S,S)-butanediol dehydrogenase', u'nitrogenase', u'pyrroloquinoline-quinone synthase', u'glycine dehydrogenase (cyanide-forming)', u'alkaline phosphatase', u'3.1.3.2  acid phosphatase', u'3.1.6.1  arylsulfatase', u'3.2.1.14  chitinase', u'3.2.1.21  beta-glucosidase', u'3.2.1.37  xylan 1,4-beta-xylosidase', u'3.2.1.4  cellulase', u'3.2.1.8  endo-1,4-beta-xylanase', u'3.2.1.91  cellulose 1,4-beta-cellobiosidase (non-reducing end)', u'3.5.1.4  amidase', u'3.5.1.5  urease', u'3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase', u'4.1.1.74  indolepyruvate decarboxylase', u'6.3.2.39  aerobactin synthase']
                ### rename and reorder bytrt2

                df = bytrt2.mean()
                print df
                enzymes = df.columns.values.tolist()

                print enzymes

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                r.assign('data', bytrt2.mean())
                r('max <- apply(data, 2, max)')
                r('final <- scale(data, center=F, scale=max)')
                r('trt <- row.names(final)')
                r.assign('enzymes', enzymes)





                base[RID] = 'Step 2 of X: Selecting your chosen taxa...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 3 of X: Calculating Species Accumulation Curves...'


                path = os.path.join('media', 'temp', 'soil_index', 'Rplots', RID)
                if not os.path.exists(path):
                    os.makedirs(path)

                r.assign("path", path)
                r("setwd(path)")
                r("options('width'=5000)")
                r.assign("RID", RID)

                base[RID] = 'Step 3 of X: Calculating Species Accumulation Curves...done!'

                base[RID] = 'Step 4 of X: Pooling pdf files for display...'

                # Combining Pdf files
                finalFile = 'media/temp/soil_index/Rplots/' + str(RID) + '/soil_index_final.pdf'

                pdf_files = [f for f in os.listdir(path) if f.endswith("pdf")]
                pdf_files = natsorted(pdf_files, key=lambda y: y.lower())

                merger = PdfFileMerger()
                for filename in pdf_files:
                    merger.append(PdfFileReader(os.path.join(path, filename), 'rb'))

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[RID]:
                        res = ''
                        return None
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                merger.write(finalFile)

                finalDict = {}
                finalDict['text'] = result

                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                return None

    except:
        if not stop_soil_index:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with soil_index!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
        return None


def getTaxaDF(savedDF, metaDF, RID):
    global base, stops, stop_soil_index, res
    try:
        base[RID] = 'Step 2 of 4: Selecting your chosen taxa or KEGG level...'

        taxaDF = savedDF.loc[:, ['sampleid', 'speciesid', 'speciesName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
        taxaDF.rename(columns={'speciesid': 'rank_id', 'speciesName': 'rank_name'}, inplace=True)

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[RID]:
            res = ''
            return None
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        taxaDF.drop('sampleid', axis=1, inplace=True)
        finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')

        speciesList = taxaDF['rank_id'].tolist()
        speciesList = list(set(speciesList))
        qs = PICRUSt.objects.using('picrust').filter(speciesid__in=speciesList)
        df = read_frame(qs, fieldnames=['speciesid__speciesid', 'rRNACount'])
        df.rename(columns={'speciesid__speciesid': 'rank_id'}, inplace=True)
        finalDF = pd.merge(finalDF, df, left_on='rank_id', right_on='rank_id', how='outer')

        keggDF, levelList = getNZDF(metaDF, finalDF, RID)

        return finalDF, keggDF, levelList

    except:
        if not stop_soil_index:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with soil_index!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
        return None


def getNZDF(metaDF, finalDF, RID):
    global base, stops, stop_soil_index, res

    try:
        nzDict = {}
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

        # create sample and species lists based on meta data selection
        wanted = ['sampleid', 'rank_id', 'rel_abund', 'abund_16S']
        profileDF = finalDF.loc[:, wanted]
        profileDF.set_index('rank_id', inplace=True)

        # get PICRUSt data for species
        speciesList = pd.unique(finalDF['rank_id'].tolist())
        qs = PICRUSt.objects.using('picrust').filter(speciesid__in=speciesList)
        picrustDF = read_frame(qs, fieldnames=['speciesid__speciesid', 'geneCount'])
        picrustDF.set_index('speciesid__speciesid', inplace=True)

        path = 'media/temp/soil_index/' + str(RID)
        if not os.path.exists(path):
            os.makedirs(path)

        if os.name == 'nt':
            numcore = 1
            listDF = np.array_split(picrustDF, numcore)
            processes = [threading.Thread(target=sumStuff, args=(listDF[x], nzDict, RID, x)) for x in xrange(numcore)]
        else:
            numcore = mp.cpu_count()
            listDF = np.array_split(picrustDF, numcore)
            processes = [mp.Process(target=sumStuff, args=(listDF[x], nzDict, RID, x)) for x in xrange(numcore)]

        for p in processes:
            p.start()

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[RID]:
                res = ''
                return None
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        for p in processes:
            p.join()

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[RID]:
                res = ''
                return None
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        levelList = []
        for key in nzDict:
            levelList.append(key)

        picrustDF = pd.DataFrame()
        for i in xrange(numcore):
            path = 'media/temp/soil_index/'+str(RID)+'/file%d.temp' % i
            frame = pd.read_csv(path)
            picrustDF = picrustDF.append(frame, ignore_index=True)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[RID]:
                res = ''
                return None
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        shutil.rmtree('media/temp/soil_index/'+str(RID))
        picrustDF.set_index('speciesid', inplace=True)

        # merge to get final gene counts for all selected samples
        taxaDF = pd.merge(profileDF, picrustDF, left_index=True, right_index=True, how='inner')

        for level in levelList:
            taxaDF[level] = taxaDF['rel_abund'] * taxaDF[level]

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[RID]:
                res = ''
                return None
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        taxaDF = taxaDF.groupby('sampleid')[levelList].agg('sum')
        taxaDF.reset_index(drop=False, inplace=True)

        taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rel_abund')

        metaDF.set_index('sampleid', drop=True, inplace=True)
        grouped = metaDF.groupby(level=0)
        metaDF = grouped.last()

        taxaDF.set_index('sampleid', drop=True, inplace=True)
        keggDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')
        keggDF.reset_index(drop=False, inplace=True)

        keggDF['rank'] = ''
        keggDF['rank_name'] = ''

        for index, row in keggDF.iterrows():
            if nz_lvl1.objects.using('picrust').filter(nz_lvl1_id=row['rank_id']).exists():
                keggDF.loc[index, 'rank'] = 'Lvl1'
                keggDF.loc[index, 'rank_name'] = nz_lvl1.objects.using('picrust').get(nz_lvl1_id=row['rank_id']).nz_lvl1_name
            elif nz_lvl2.objects.using('picrust').filter(nz_lvl2_id=row['rank_id']).exists():
                keggDF.loc[index, 'rank'] = 'Lvl2'
                keggDF.loc[index, 'rank_name'] = nz_lvl2.objects.using('picrust').get(nz_lvl2_id=row['rank_id']).nz_lvl2_name
            elif nz_lvl3.objects.using('picrust').filter(nz_lvl3_id=row['rank_id']).exists():
                keggDF.loc[index, 'rank'] = 'Lvl3'
                keggDF.loc[index, 'rank_name'] = nz_lvl3.objects.using('picrust').get(nz_lvl3_id=row['rank_id']).nz_lvl3_name
            elif nz_lvl4.objects.using('picrust').filter(nz_lvl4_id=row['rank_id']).exists():
                keggDF.loc[index, 'rank'] = 'Lvl4'
                keggDF.loc[index, 'rank_name'] = nz_lvl4.objects.using('picrust').get(nz_lvl4_id=row['rank_id']).nz_lvl4_name
            elif nz_entry.objects.using('picrust').filter(nz_lvl5_id=row['rank_id']).exists():
                keggDF.loc[index, 'rank'] = 'Lvl5'
                keggDF.loc[index, 'rank_name'] = nz_entry.objects.using('picrust').get(nz_lvl5_id=row['rank_id']).nz_name

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[RID]:
                res = ''
                return None
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        return keggDF, levelList

    except:
        if not stop_soil_index:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with soil_index!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
        return None


def sumStuff(slice, koDict, RID, num):
    global base, stops, res
    db.close_old_connections()

    f = open('media/temp/soil_index/'+str(RID)+'/file'+str(num)+".temp", 'w')

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
        if stops[RID]:
            res = ''
            return None
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    f.close()


def removesoil_indexFiles(request):
    if request.is_ajax():
        RID = request.GET["all"]

        file = "media/temp/soil_index/Rplots/" + str(RID) + ".soil_index.pdf"
        if os.path.exists(file):
            os.remove(file)

        file = "media/temp/soil_index/" + str(RID) + ".pkl"
        if os.path.exists(file):
            os.remove(file)

        file = "media/temp/soil_index/" + str(RID) + ".csv"
        if os.path.exists(file):
            os.remove(file)

        return HttpResponse()


def getTabsoil_index(request):
    if request.is_ajax():
        RID = request.GET["all"]
        myDir = 'media/temp/soil_index/'
        fileName = str(myDir) + str(RID) + '.pkl'
        savedDF = pd.read_pickle(fileName)

        myDir = 'media/temp/soil_index/'
        fileName = str(myDir) + str(RID) + '.csv'
        savedDF.to_csv(fileName)

        myDict = {}
        myDir = 'temp/soil_index/'
        fileName = str(myDir) + str(RID) + '.csv'
        myDict['name'] = str(fileName)
        res = simplejson.dumps(myDict)

        return HttpResponse(res, content_type='application/json')
