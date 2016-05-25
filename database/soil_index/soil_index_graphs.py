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
                base[RID] = 'Step 1 of 3: Selecting your chosen meta-variables...'

                path = pickle.loads(request.session['savedDF'])
                savedDF = pd.read_pickle(path)

                result = ''

                # Select samples and meta-variables from savedDF
                metaValsCat = all['metaValsCat']
                metaIDsCat = all['metaIDsCat']

                # Get maxvals for plotting
                maxVals = all['maxVals']

                locMax = all['locMax']


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

                base[RID] = 'Step 1 of 3: Selecting your chosen meta-variables...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 2 of 3: Selecting your chosen taxa or KEGG level...'

                finalDF, keggDF, levelList = getTaxaDF(savedDF, metaDF, RID)
                finalDF = finalDF[finalDF.rel_abund != 0]

                count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='abund')
                count_rDF.fillna(0, inplace=True)
                count_rDF['sum'] = count_rDF[count_rDF.gt(0)].count(axis=1)
                count_rDF['singletons'] = count_rDF[count_rDF.lt(2)].sum(axis=1)
                count_rDF['coverage'] = 1 - count_rDF['singletons'] / count_rDF['sum']
                count_rDF['total'] = count_rDF.sum(axis=1)

                # save location info to session
                myDir = 'media/temp/soil_index/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                meta_rDF = savedDF.drop_duplicates(subset='sampleid', take_last=True)

                # Removes samples (rows) that are not in our sampleeta_rlist
                meta_rDF = meta_rDF.loc[meta_rDF['sampleid'].isin(catSampleIDs)]

                bioDF  = meta_rDF[catFields_edit].copy()
                chemDF = meta_rDF[catFields_edit].copy()
                physDF = meta_rDF[catFields_edit].copy()

                # replace NaN's with 0's
                bioDF = bioDF.fillna(0)
                chemDF = chemDF.fillna(0)
                physDF = physDF.fillna(0)

                # check for null columns, fill with blanks if found (0's, not nulls)
                if 'microbial_biomass_C' in meta_rDF.columns:
                    bioDF['microbial_biomass_C'] = meta_rDF['microbial_biomass_C']
                else:
                    bioDF['microbial_biomass_C'] = 0.0

                if 'microbial_biomass_N' in meta_rDF.columns:
                    bioDF['microbial_biomass_N'] = meta_rDF['microbial_biomass_N']
                else:
                    bioDF['microbial_biomass_N'] = 0.0

                if 'microbial_respiration' in meta_rDF.columns:
                    bioDF['microbial_respiration'] = meta_rDF['microbial_respiration']
                else:
                    bioDF['microbial_respiration'] = 0.0

                if 'soil_pH' in meta_rDF.columns:
                    chemDF['soil_pH'] = meta_rDF['soil_pH']
                else:
                    chemDF['soil_pH'] = 0.0
                if 'soil_C' in meta_rDF.columns:
                    chemDF['soil_C'] = meta_rDF['soil_C']
                else:
                    chemDF['soil_C'] = 0.0
                if 'soil_OM' in meta_rDF.columns:
                    chemDF['soil_OM'] = meta_rDF['soil_OM']
                else:
                    chemDF['soil_OM'] = 0.0
                if 'soil_EC' in meta_rDF.columns:
                    chemDF['soil_EC'] = meta_rDF['soil_EC']
                else:
                    chemDF['soil_EC'] = 0.0
                if 'soil_N' in meta_rDF.columns:
                    chemDF['soil_N'] = meta_rDF['soil_N']
                else:
                    chemDF['soil_N'] = 0.0
                if 'soil_P' in meta_rDF.columns:
                    chemDF['soil_P'] = meta_rDF['soil_P']
                else:
                    chemDF['soil_P'] = 0.0
                if 'soil_K' in meta_rDF.columns:
                    chemDF['soil_K'] = meta_rDF['soil_K']
                else:
                    chemDF['soil_K'] = 0.0

                if 'water_content_soil' in meta_rDF.columns:
                    physDF['water_content_soil'] = meta_rDF['water_content_soil']
                else:
                    physDF['water_content_soil'] = 0.0
                if 'bulk_density' in meta_rDF.columns:
                    physDF['bulk_density'] = meta_rDF['bulk_density']
                else:
                    physDF['bulk_density'] = 0.0
                if 'porosity' in meta_rDF.columns:
                    physDF['porosity'] = meta_rDF['porosity']
                else:
                    physDF['porosity'] = 0.0
                if 'soil_surf_hard' in meta_rDF.columns:
                    physDF['soil_surf_hard'] = meta_rDF['soil_surf_hard']
                else:
                    physDF['soil_surf_hard'] = 0.0
                if 'soil_subsurf_hard' in meta_rDF.columns:
                    physDF['soil_subsurf_hard'] = meta_rDF['soil_subsurf_hard']
                else:
                    physDF['soil_subsurf_hard'] = 0.0

                # get means
                bioDF  = bioDF.groupby(catFields_edit)
                bioDF = bioDF.mean()
                chemDF = chemDF.groupby(catFields_edit)
                chemDF = chemDF.mean()
                physDF = physDF.groupby(catFields_edit)
                physDF = physDF.mean()

                if metaDictCat:
                    for key in metaDictCat:
                        meta_rDF = meta_rDF.loc[meta_rDF[key].isin(metaDictCat[key])]

                wantedList = allFields + ['sampleid', 'sample_name']
                meta_rDF = meta_rDF[wantedList]
                meta_rDF.set_index('sampleid', drop=True, inplace=True)

                rnaDF = finalDF.replace(r'\s+', np.nan, regex=True)
                rnaDF.fillna(0.0, axis=1, inplace=True)
                rnaDF['bins'] = 'bin1'
                rnaDF.ix[(rnaDF.rRNACount > 0) & (rnaDF.rRNACount <= 2), 'bins'] = 'bin1'
                rnaDF.ix[(rnaDF.rRNACount > 2) & (rnaDF.rRNACount < 4), 'bins'] = 'bin2'
                rnaDF.ix[rnaDF.rRNACount >= 4, 'bins'] = 'bin3'
                myBins = list(set(rnaDF['bins'].tolist()))

                binDF = rnaDF.groupby(['sampleid', 'bins'])['rel_abund'].sum()
                binDF = binDF.unstack(level=-1)
                binDF.fillna(0.0, axis=1, inplace=True)
                binDF.reset_index(drop=False, inplace=True)
                sumDF1 = binDF.groupby('sampleid')[myBins].sum()

                sumDF2 = finalDF.groupby('sampleid')['abund_16S', 'rich', 'diversity'].sum()
                allDF = pd.merge(sumDF2, sumDF1, left_index=True, right_index=True, how='outer')

                allDF = pd.merge(allDF, meta_rDF, left_index=True, right_index=True, how='outer')
                allDF = pd.merge(allDF, count_rDF, left_index=True, right_index=True, how='inner')

                wantList = ['abund_16S', 'rich', 'diversity', 'coverage'] + myBins

                bytrt1 = allDF.groupby(catFields_edit)[wantList]
                df1 = bytrt1.mean()

                # merge with other df's before sending to R
                df1 = pd.merge(df1, bioDF, left_index=True, right_index=True, how='outer')
                df1 = pd.merge(df1, chemDF, left_index=True, right_index=True, how='outer')
                df1 = pd.merge(df1, physDF, left_index=True, right_index=True, how='outer')

                df1.reset_index(drop=False, inplace=True)
                if len(catFields_edit) > 1:
                    for index, row in df1.iterrows():
                        df1.loc[index, 'merge'] = "; ".join(row[catFields_edit])
                else:
                    df1.loc[:, 'merge'] = df1.loc[:, catFields_edit[0]]

                df1.set_index('merge', inplace=True)

                nzDF = keggDF.groupby(['sampleid', 'rank_name'])['abund'].sum()
                nzDF = nzDF.unstack(level=-1)
                starDF = pd.merge(nzDF, meta_rDF, left_index=True, right_index=True, how='outer')
                bytrt2 = starDF.groupby(catFields_edit)

                myList = [u'3.2.1.4  cellulase', u'3.2.1.8  endo-1,4-beta-xylanase',
                          u'3.2.1.21  beta-glucosidase', u'3.2.1.37  xylan 1,4-beta-xylosidase',
                          u'3.2.1.91  cellulose 1,4-beta-cellobiosidase (non-reducing end)',
                          u'3.5.1.4  amidase', u'3.5.1.5  urease', u'3.1.3.1  alkaline phosphatase',
                          u'3.1.3.2  acid phosphatase', u'3.1.6.1  arylsulfatase', u'1.18.6.1  nitrogenase',
                          u'1.3.3.11  pyrroloquinoline-quinone synthase', u'6.3.2.39  aerobactin synthase',
                          u'3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase',
                          u'4.1.1.74  indolepyruvate decarboxylase',
                          u'1.1.1.76  (S,S)-butanediol dehydrogenase',
                          u'1.4.99.5  glycine dehydrogenase (cyanide-forming)',
                          u'3.2.1.14  chitinase', u'1.14.99.39  ammonia monooxygenase',
                          u'1.7.2.6  hydroxylamine dehydrogenase', u'1.7.99.4  nitrate reductase',
                          u'1.7.2.1  nitrite reductase (NO-forming)',
                          u'1.7.2.5  nitric oxide reductase (cytochrome c)',
                          u'1.7.2.4  nitrous-oxide reductase']

                newList = [u'cellulase', u'endo-1,4-beta-xylanase', u'beta-glucosidase',
                           u'xylan 1,4-beta-xylosidase',
                           u'cellulose 1,4-beta-cellobiosidase (non-reducing end)', u'amidase', u'urease',
                           u'alkaline phosphatase', u'acid phosphatase', u'arylsulfatase', u'nitrogenase',
                           u'pyrroloquinoline-quinone synthase', u'aerobactin synthase',
                           u'1-aminocyclopropane-1-carboxylate deaminase', u'indolepyruvate decarboxylase',
                           u'(S,S)-butanediol dehydrogenase', u'glycine dehydrogenase (cyanide-forming)',
                           u'chitinase', u'ammonia monooxygenase', u'hydroxylamine dehydrogenase', u'nitrate reductase',
                           u'nitrite reductase (NO-forming)', u'nitric oxide reductase (cytochrome c)',
                           u'nitrous-oxide reductase']

                df2 = bytrt2[myList].mean()
                df2.rename(columns={u'3.2.1.4  cellulase': 'cellulase',
                                    u'3.2.1.8  endo-1,4-beta-xylanase': 'endo-1,4-beta-xylanase',
                                    u'3.2.1.21  beta-glucosidase': 'beta-glucosidase',
                                    u'3.2.1.37  xylan 1,4-beta-xylosidase': 'xylan 1,4-beta-xylosidase',
                                    u'3.2.1.91  cellulose 1,4-beta-cellobiosidase (non-reducing end)': 'cellulose 1,4-beta-cellobiosidase (non-reducing end)',
                                    u'3.5.1.4  amidase': 'amidase',
                                    u'3.5.1.5  urease': 'urease',
                                    u'3.1.3.1  alkaline phosphatase': 'alkaline phosphatase',
                                    u'3.1.3.2  acid phosphatase': 'acid phosphatase',
                                    u'3.1.6.1  arylsulfatase': 'arylsulfatase',
                                    u'1.18.6.1  nitrogenase': 'nitrogenase',
                                    u'1.3.3.11  pyrroloquinoline-quinone synthase': 'pyrroloquinoline-quinone synthase',
                                    u'6.3.2.39  aerobactin synthase': 'aerobactin synthase',
                                    u'3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase': '1-aminocyclopropane-1-carboxylate deaminase',
                                    u'4.1.1.74  indolepyruvate decarboxylase': 'indolepyruvate decarboxylase',
                                    u'1.1.1.76  (S,S)-butanediol dehydrogenase': '(S,S)-butanediol dehydrogenase',
                                    u'1.4.99.5  glycine dehydrogenase (cyanide-forming)': 'glycine dehydrogenase (cyanide-forming)',
                                    u'3.2.1.14  chitinase': 'chitinase',
                                    u'1.14.99.39  ammonia monooxygenase': 'ammonia monooxygenase',
                                    u'1.7.2.6  hydroxylamine dehydrogenase': 'hydroxylamine dehydrogenase',
                                    u'1.7.99.4  nitrate reductase': 'nitrate reductase',
                                    u'1.7.2.1  nitrite reductase (NO-forming)': 'nitrite reductase (NO-forming)',
                                    u'1.7.2.5  nitric oxide reductase (cytochrome c)': 'nitric oxide reductase (cytochrome c)',
                                    u'1.7.2.4  nitrous-oxide reductase': 'nitrous-oxide reductase'}, inplace=True)

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                path = os.path.join('media', 'temp', 'soil_index', 'Rplots', RID)
                if not os.path.exists(path):
                    os.makedirs(path)

                r.assign("path", path)
                r("setwd(path)")
                r("options('width'=5000)")
                r.assign("RID", RID)

                df2.reset_index(drop=False, inplace=True)
                if len(catFields_edit) > 1:
                    for index, row in df2.iterrows():
                        df2.loc[index, 'merge'] = "; ".join(row[catFields_edit])
                else:
                    df2.loc[:, 'merge'] = df2.loc[:, catFields_edit[0]]

                df2.set_index('merge', inplace=True)
                df2 = df2[newList]
                result += str(df2)

                maxDict = {
                    'cellulase': 1e9,
                    'endo-1,4-beta-xylanase': 1e9,
                    'beta-glucosidase': 1e9,
                    'xylan 1,4-beta-xylosidase': 1e9,
                    'cellulose 1,4-beta-cellobiosidase (non-reducing end)': 1e9,
                    'amidase': 1e9,
                    'urease': 1e9,
                    'alkaline phosphatase': 1e9,
                    'acid phosphatase': 1e9,
                    'arylsulfatase': 1e9,
                    'nitrogenase': 1e9,
                    'pyrroloquinoline-quinone synthase': 1e9,
                    'aerobactin synthase': 1e9,
                    '1-aminocyclopropane-1-carboxylate deaminase': 1e9,
                    'indolepyruvate decarboxylase': 1e9,
                    '(S,S)-butanediol dehydrogenase': 1e9,
                    'glycine dehydrogenase (cyanide-forming)': 1e9,
                    'chitinase': 1e9,
                    'ammonia monooxygenase': 1e9,
                    'hydroxylamine dehydrogenase': 1e9,
                    'nitrate reductase': 1e9,
                    'nitrite reductase (NO-forming)': 1e9,
                    'nitric oxide reductase (cytochrome c)': 1e9,
                    'nitrous-oxide reductase': 1e9
                }

                scaleDF = df2.copy()
                for key in maxDict:
                    if locMax:
                        maxVal = max(scaleDF[key])
                    else:
                        maxVal = float(maxVals[key])
                    if maxVal == 0:
                        maxVal = 1
                    scaleDF[key] /= maxVal

                r.assign('data', scaleDF)
                r('trt <- row.names(data)')
                r.assign('enzymes', newList)

                rowcount = len(scaleDF)
                r.assign('rowcount', rowcount)

                r.assign('dat', df1)  # use dat for displaying rounded data
                r('odat <- dat')  # use odat for calculating graphs with more precision

                r('dat$abund_16S <- signif(dat$abund_16S, 3)')
                r('dat$rich <- signif(dat$rich, 3)')
                r('dat$diversity <- signif(dat$diversity, 3)')
                r('dat$coverage <- signif(dat$coverage, 3)')

                # load max vals from maxVals dictionary into R
                if locMax:
                    # bio
                    r('maxAbund <- max(odat$abund_16S)')
                    r('if(maxAbund==0) maxAbund=1')
                    r('maxRich <- max(odat$rich)')
                    r('if(maxRich==0) maxRich=1')
                    r('maxDiversity <- max(odat$diversity)')
                    r('if(maxDiversity==0) maxDiversity=1')

                    # chem
                    r('maxPH <- 14')
                    r('maxC <- max(odat$soil_C)')
                    r('if(maxC==0) maxC=1')
                    r('maxOM <- max(odat$soil_OM)')
                    r('if(maxOM==0) maxOM=1')
                    r('maxP <- max(odat$soil_P)')
                    r('if(maxP==0) maxP=1')
                    r('maxK <- max(odat$soil_K)')
                    r('if(maxK==0) maxK=1')

                    # phys
                    r('maxWCS <- max(odat$water_content_soil)')
                    r('if(maxWCS==0) maxWCS=1')
                    r('maxBD <- max(odat$bulk_density)')
                    r('if(maxBD==0) maxBD=1')
                    r('maxPORO <- max(odat$porosity)')
                    r('if(maxPORO==0) maxPORO=1')
                    r('maxSH <- max(odat$soil_surf_hard)')
                    r('if(maxSH==0) maxSH=1')
                    r('maxSSH <- max(odat$soil_subsurf_hard)')
                    r('if(maxSSH==0) maxSSH=1')

                else:
                    # bio
                    r.assign('maxAbund', float(maxVals['Abundance']))
                    r.assign('maxRich', float(maxVals['Richness']))
                    r.assign('maxDiversity', float(maxVals['Diversity']))
                    # chem
                    r.assign('maxC', float(maxVals['Carbon']))
                    r.assign('maxOM', float(maxVals['Organic Matter']))
                    r.assign('maxPH', float(maxVals['pH']))
                    r.assign('maxP', float(maxVals['Phosphorus']))
                    r.assign('maxK', float(maxVals['Potassium']))
                    # physical
                    r.assign('maxWCS', float(maxVals['Water Content']))
                    r.assign('maxBD', float(maxVals['Bulk Density']))
                    r.assign('maxPORO', float(maxVals['Porosity']))

                r.assign('myBins', myBins)  # load bins

                # bio
                r('odat$abund_16S <- odat$abund_16S / maxAbund')
                r('odat$rich <- odat$rich / maxRich')
                r('odat$diversity <- odat$diversity / maxDiversity')

                # chem
                r('odat$soil_pH <- odat$soil_pH / maxPH')
                r('odat$soil_C <- odat$soil_C / maxC')
                r('odat$soil_OM <- odat$soil_OM / maxOM')
                r('odat$soil_P <- odat$soil_P / maxP')
                r('odat$soil_K <- odat$soil_K / maxK')

                # phys
                r('odat$water_content_soil <- odat$water_content_soil / maxWCS')
                r('odat$bulk_density <- odat$bulk_density / maxBD')
                r('odat$porosity <- odat$porosity / maxPORO')
                r('odat$soil_surf_hard <- odat$odat$soil_surf_hard / maxSH')
                r('odat$soil_subsurf_hard <- odat$soil_subsurf_hard / maxSSH')

                r("colGIBBs <- c('blue', 'blue', 'blue', 'blue', 'blue', \
                    'green', 'green', 'red', 'red', 'turquoise', \
                    'magenta', 'yellow', 'salmon', 'darkgray', 'darkgray', \
                    'brown', 'brown', 'brown')")

                r("colN <- c('blue', 'green', 'red', 'turquoise', \
                    'magenta', 'yellow')")

                r("pdf_counter <- 1")
                for row in range(1, rowcount+1):
                    r("pdf(paste('soil_index_temp', pdf_counter, '.pdf', sep=''), height=4, width=8.5)")
                    r.assign('off', row)
                    r('curName <- row.names(data)[off]')
                    print r('curName')
                    r('layout(matrix(c(1,3,5,6,2,4,4,4), 2, 4, byrow=T), widths=c(2.5,2,2,2))')

                    # GIBBs
                    r('stars(matrix(1, ncol=18, nrow=1), \
                        draw.segments=T, col.segments=rep("white", 18), scale=F, full=T, labels=NA, \
                        cex=0.8, ncol=1, mar=c(1,1,5,1), lty=3)')
                    r('title("GIBBs", line=0)')
                    r('stars(matrix(0.75, ncol=18, nrow=1), \
                        draw.segments=T, col.segments=rep("white", 18), scale=F, full=T, labels=NA, \
                        cex=0.8, ncol=1, mar=c(1,1,5,1), add=T, lty=3)')
                    r('stars(matrix(0.5, ncol=18, nrow=1), \
                        draw.segments=T, col.segments=rep("white", 18), scale=F, full=T, labels=NA, \
                        cex=0.8, ncol=1, mar=c(1,1,5,1), add=T, lty=3)')
                    r('stars(matrix(0.25, ncol=18, nrow=1), \
                        draw.segments=T, col.segments=rep("white", 18), scale=F, full=T, labels=NA, \
                        cex=0.8, ncol=1, mar=c(1,1,5,1), add=T, lty=3)')

                    r("gibbs <- data[,c('cellulase', 'endo.1.4.beta.xylanase', 'beta.glucosidase', \
                        'xylan.1.4.beta.xylosidase', 'cellulose.1.4.beta.cellobiosidase..non.reducing.end.', \
                        'amidase', 'urease','alkaline.phosphatase', 'acid.phosphatase', 'arylsulfatase', \
                        'nitrogenase', 'pyrroloquinoline.quinone.synthase', 'aerobactin.synthase', \
                        'X1.aminocyclopropane.1.carboxylate.deaminase', 'indolepyruvate.decarboxylase', \
                        'X.S.S..butanediol.dehydrogenase', 'glycine.dehydrogenase..cyanide.forming.', \
                        'chitinase')]")

                    r('stars(gibbs[off,], lwd=1, draw.segments=T, \
                        col.segments=colGIBBs, scale=F, full=T, labels=NA, \
                        cex=0.8, ncol=1, mar=c(1,1,5,1), add=T)')

                    # N cycle
                    r('stars(matrix(1, ncol=6, nrow=1), \
                        draw.segments=T, col.segments=rep("white", 6), scale=F, full=T, labels=NA, \
                        cex=0.8, ncol=1, mar=c(1,1,5,1), lty=3)')
                    r('title("N cycle", line=0)')
                    r('stars(matrix(0.75, ncol=6, nrow=1), \
                        draw.segments=T, col.segments=rep("white", 6), scale=F, full=T, labels=NA, \
                        cex=0.8, ncol=1, mar=c(1,1,5,1), add=T, lty=3)')
                    r('stars(matrix(0.5, ncol=6, nrow=1), \
                        draw.segments=T, col.segments=rep("white", 6), scale=F, full=T, labels=NA, \
                        cex=0.8, ncol=1, mar=c(1,1,5,1), add=T, lty=3)')
                    r('stars(matrix(0.25, ncol=6, nrow=1), \
                        draw.segments=T, col.segments=rep("white", 6), scale=F, full=T, labels=NA, \
                        cex=0.8, ncol=1, mar=c(1,1,5,1), add=T, lty=3)')


                    r("ncycle <- data[,c('ammonia.monooxygenase', 'hydroxylamine.dehydrogenase', \
                        'nitrate.reductase', 'nitrite.reductase..NO.forming.', \
                         'nitric.oxide.reductase..cytochrome.c.', 'nitrous.oxide.reductase')]")

                    r('stars(ncycle[off,], lwd=1, draw.segments=T, \
                        col.segments=colN, scale=F, full=T, labels=NA, \
                        cex=0.8, ncol=1, mar=c(1,1,5,1), add=T)')

                    # Biological
                    r('vals <- as.numeric(odat[off, c("microbial_respiration", "coverage", "diversity", "rich", "abund_16S")])')
                    r('par(mar=c(0,5,5,5))')
                    r('bar <- barplot(height=vals, width=0.2, xlim=c(0,1), \
                        ylim=c(0,1), horiz=T, space=0, \
                        names.arg=c("Microbial Respiration", "Coverage", "Diversity", "Richness", "Abundance"),\
                        cex.names=0.8, las=2, axes=F, col=c(2,3,4), beside=T)')
                    r('axis(1, at=c(0, 0.25, 0.5, 0.75, 1), lwd.ticks=0.2, cex.axis=0.8)')
                    r('values <- c(dat[off, c("microbial_respiration", "coverage", "diversity", "rich", "abund_16S")])')
                    r('text(x=vals+0.27, y=bar, labels=values, cex=0.8, xpd=TRUE)')
                    r('title("Biological", line=1)')

                    r('vals <- as.matrix(t(odat[off, myBins]))')
                    r('par(mar=c(5,15,1,18))')
                    r('bar <- barplot(height=vals, width=0.2, xlim=c(0,1), \
                        ylim=c(0,1), horiz=T, space=0, \
                        names.arg=c("rRNA Copy Number"), \
                        cex.names=0.8, las=2, axes=F, col=c("red", "darkgray", "green"), beside=F)')
                    r('axis(1, at=c(0, 0.25, 0.5, 0.75, 1), lwd.ticks=0.2, cex.axis=0.8)')
                    r('text(x=c(0.2, 0.5, 0.8), y=bar+0.2, labels=round(vals, 3), cex=0.8, xpd=TRUE)')
                    r('par(xpd=T)')
                    r('legend(1, 0.25, c("x <= 2", "2 < x < 4", "x >= 4"), cex=0.8, bty="n", fill=c("red", "darkgray", "green"))')

                    # Chemical, needs 'Minor Elements'
                    r('vals <- as.numeric(odat[off, c("soil_K", "soil_P", "soil_pH", "soil_C", "soil_OM")])')
                    r('par(mar=c(0,5,5,5))')
                    r('bar <- barplot(height=vals, width=0.2, xlim=c(0,1), \
                        ylim=c(0,1), horiz=T, space=0, \
                        names.arg=c("Potassium", "Phosphorus", "pH", "Carbon", "Organic Matter"),\
                        cex.names=0.8, las=2, axes=F, col=c(2,3,4), beside=T)')
                    r('axis(1, at=c(0, 0.25, 0.5, 0.75, 1), lwd.ticks=0.2, cex.axis=0.8)')
                    r('values <- dat[off, c("soil_K", "soil_P", "soil_pH", "soil_C", "soil_OM")]')
                    r('text(x=vals+0.2, y=bar, labels=signif(values, 3), cex=0.8, xpd=TRUE)')
                    r('title("Chemical", line=1)')

                    # Physical, needs 'Aggregate Stability'
                    r('vals <- as.numeric(odat[off, c("soil_subsurf_hard", "soil_surf_hard", "water_content_soil", "bulk_density", "porosity")])')
                    r('par(mar=c(0,5,5,5))')
                    r('bar <- barplot(height=vals, width=0.2, xlim=c(0,1), \
                        ylim=c(0,1), horiz=T, space=0, \
                        names.arg=c("Subsurface Hardness", "Surface Hardness", "Water Content", "Bulk Density", "Porosity"),\
                        cex.names=0.8, las=2, axes=F, col=c(2,3,4), beside=T)')
                    r('axis(1, at=c(0, 0.25, 0.5, 0.75, 1), lwd.ticks=0.2, cex.axis=0.8)')
                    r('values <- dat[off, c("soil_subsurf_hard", "soil_surf_hard", "water_content_soil", "bulk_density", "porosity")]')
                    r('text(x=vals+0.2, y=bar, labels=signif(values, 3), cex=0.8, xpd=TRUE)')
                    r('title("Physical", line=1)')
                    r('dev.off()')
                    r("pdf_counter <- pdf_counter + 1")

                base[RID] = 'Step 2 of 3: Selecting your chosen taxa...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[RID]:
                    res = ''
                    return None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 3 of 3: Pooling pdf files for display...'

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
        base[RID] = 'Step 2 of 3: Selecting your chosen taxa or KEGG level...'
        taxaDF = savedDF.loc[:, ['sampleid', 'speciesid', 'speciesName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
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
            taxaDF[level] = taxaDF['abund_16S'] * taxaDF[level]

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[RID]:
                res = ''
                return None
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        taxaDF = taxaDF.groupby('sampleid')[levelList].agg('sum')
        taxaDF.reset_index(drop=False, inplace=True)

        taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund')

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
