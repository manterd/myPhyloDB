import datetime
from django.http import HttpResponse
import logging
from natsort import natsorted
import pandas as pd
from pyper import *
from PyPDF2 import PdfFileReader, PdfFileMerger
import simplejson

from database.utils import multidict
from database.utils_kegg import getTaxaDF, getKeggDF, getNZDF
import database.queue


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getSpAC(request, stops, RID, PID):
    try:
        while True:
            if request.is_ajax():
                allJson = request.body.split('&')[0]
                all = simplejson.loads(allJson)
                database.queue.setBase(RID, 'Step 1 of 4: Selecting your chosen meta-variables...')
                myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
                path = str(myDir) + 'usr_norm_data.csv'

                with open(path, 'rb') as f:
                    savedDF = pd.read_csv(f, index_col=0, sep=',')

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
                        if idDictCat[key] not in catSampleIDs:
                            catSampleIDs.extend(idDictCat[key])

                if not catFields_edit:
                    catSampleIDs = savedDF.sampleid.unique().tolist()

                allSampleIDs = catSampleIDs
                allFields = catFields_edit

                # Removes samples (rows) that are not in our samplelist
                metaDF = savedDF.drop_duplicates(subset='sampleid', take_last=True)
                if allSampleIDs:
                    metaDF = metaDF.loc[metaDF['sampleid'].isin(allSampleIDs)]

                # make sure column types are correct
                metaDF[catFields_edit] = metaDF[catFields_edit].astype(str)

                if metaDictCat:
                    for key in metaDictCat:
                        metaDF = metaDF.loc[metaDF[key].isin(metaDictCat[key])]

                wantedList = allFields + ['sampleid', 'sample_name']
                metaDF = metaDF[wantedList]
                metaDF.set_index('sampleid', drop=True, inplace=True)

                result += 'Categorical variables selected by user: ' + ", ".join(catFields) + '\n'
                result += 'Categorical variables removed from analysis (contains only 1 level): ' + ", ".join(removed) + '\n'
                result += '===============================================\n'

                database.queue.setBase(RID, 'Step 1 of 4: Selecting your chosen meta-variables...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 2 of 4: Selecting your chosen taxa or KEGG level...')

                DepVar = 1
                finalDF = pd.DataFrame()
                if button3 == 1:
                    DepVar = int(all["DepVar_taxa"])
                    finalDF, missingList = getTaxaDF('abund', selectAll, '', savedDF, metaDF, allFields, DepVar, RID, stops, PID)

                if button3 == 2:
                    DepVar = int(all["DepVar_kegg"])
                    finalDF = getKeggDF('abund', keggAll, '', savedDF, metaDF, allFields, DepVar, RID, stops, PID)

                if button3 == 3:
                    DepVar = int(all["DepVar_nz"])
                    finalDF = getNZDF('abund', nzAll, '', savedDF, metaDF, allFields, DepVar, RID, stops, PID)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/spac/'
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

                database.queue.setBase(RID, 'Step 2 of 4: Selecting your chosen taxa...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 3 of 4: Calculating Species Accumulation Curves...')

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                path = os.path.join('myPhyloDB', 'media', 'temp', 'spac', 'Rplots', RID)
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
                    grouped = metaDF.groupby(allFields)
                    for name, group in grouped:
                        # print name
                        if isinstance(name, tuple):
                            name = "; ".join(name)
                        r.assign('name', name)
                        IDs = group.index.values.tolist()
                        data = count_rDF.ix[IDs]
                        r.assign("data", data)

                        r("pdf(paste('SpAC_temp', pdf_counter, '.pdf', sep=''))")
                        if method == 'estaccumR':
                            r("x <- estaccumR(data)")
                        elif method == 'poolaccum':
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
                                main=paste(name, '\nInsufficient data to generate plot (n <= 3)!', sep='')) \
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

                database.queue.setBase(RID, 'Step 3 of 4: Calculating Species Accumulation Curves...done!')

                database.queue.setBase(RID, 'Step 4 of 4: Pooling pdf files for display...')

                # Combining Pdf files
                finalFile = 'myPhyloDB/media/temp/spac/Rplots/' + str(RID) + '/SpAC_final.pdf'

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
                return HttpResponse(res, content_type='application/json')

    except Exception as e:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "There was an error during your analysis:\nError: " + str(e.message) + "\nTimestamp: " + str(datetime.datetime.now())
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')
