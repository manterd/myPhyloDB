import ast
import datetime
from django.http import HttpResponse
from django_pandas.io import read_frame
import logging
from natsort import natsorted
import pandas as pd
from PyPDF2 import PdfFileReader, PdfFileMerger
from pyper import *
import simplejson
import time

from database.models import PICRUSt
from database.models import ko_lvl1, ko_lvl2, ko_lvl3
from database.models import nz_lvl1, nz_lvl2, nz_lvl3, nz_lvl4, nz_entry
from database.utils import multidict
from database.models import Phyla, Class, Order, Family, Genus, Species
import database.queue


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}
done = {}

LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def updateGAGE(request):
    global base, stage, time1, time2, TimeDiff
    if request.is_ajax():
        RID = request.GET["all"]
        time2[RID] = time.time()

        try:
            TimeDiff[RID] = time2[RID] - time1[RID]
        except:
            TimeDiff[RID] = 0

        if done[RID]:
            stage[RID] = 'Analysis is complete, results are loading'
        else:
            try:
                if TimeDiff[RID] == 0:
                    stage[RID] = 'Analysis has been placed in queue, there are ' + str(database.queue.stat(RID)) + ' others in front of you.'
                else:
                    stage[RID] = str(base[RID]) + '<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
            except:
                if TimeDiff[RID] == 0:
                    stage[RID] = 'In queue'
                else:
                    stage[RID] = '<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]

        myDict = {'stage': stage[RID]}
        json_data = simplejson.dumps(myDict, encoding="Latin-1")
        return HttpResponse(json_data, content_type='application/json')


def removeRIDGAGE(RID):
    global base, stage, time1, time2, TimeDiff, done
    try:
        base.pop(RID, None)
        stage.pop(RID, None)
        time1.pop(RID, None)
        time2.pop(RID, None)
        TimeDiff.pop(RID, None)
        done.pop(RID, None)
        return True
    except:
        return False


def getGAGE(request, stops, RID, PID):
    global base, stage, time1, TimeDiff, done
    done[RID] = False
    try:
        while True:
            if request.is_ajax():
                # Get variables from web page
                allJson = request.body.split('&')[0]
                all = simplejson.loads(allJson)

                time1[RID] = time.time()  # Moved these down here so RID is available
                base[RID] = 'Step 1 of 4: Selecting your chosen meta-variables...'

                myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
                path = str(myDir) + 'usr_norm_data.csv'

                with open(path, 'rb') as f:
                    savedDF = pd.read_csv(f, index_col=0, sep='\t')

                result = ''

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

                if not catFields_edit:
                    myDict = {}
                    myDict['error'] = "Selected meta data only has one level.\nPlease select a different variable(s)."
                    res = simplejson.dumps(myDict)
                    return HttpResponse(res, content_type='application/json')

                catSampleIDs = []
                if metaIDsCat:
                    idDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaIDsCat)
                    for key in sorted(idDictCat):
                        catSampleIDs.extend(idDictCat[key])

                # Removes samples (rows) that are not in our samplelist
                tempDF = savedDF.loc[savedDF['sampleid'].isin(catSampleIDs)]

                if metaDictCat:
                    for key in metaDictCat:
                        tempDF = tempDF.loc[tempDF[key].isin(metaDictCat[key])]

                result += 'Categorical variables selected by user: ' + ", ".join(catFields) + '\n'
                result += 'Categorical variables removed from analysis (contains only 1 level): ' + ", ".join(removed) + '\n'
                result += '\n===============================================\n'

                base[RID] = 'Step 1 of 4: Selecting your chosen meta-variables...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 2 of 4: Mapping phylotypes to KEGG pathways...'

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                DepVar = int(all["DepVar"])
                keggString = all["kegg"]
                keggDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(keggString)
                nameList = []
                for value in keggDict.itervalues():
                    if isinstance(value, list):
                        nameList.extend(value)
                    else:
                        nameList.append(value)

                r("library(gage)")
                r("load('myPhyloDB/media/kegg/kegg.gs.RData')")

                keggDict = {}
                r("selPaths <- vector()")
                for i in nameList:
                    pathStr = i.split('[PATH:')[1].split(']')[0]
                    r.assign("pathStr", pathStr)
                    r("selPath <- kegg.gs[grepl(paste(pathStr), names(kegg.gs))]")
                    key = r.get("names(selPath)")
                    value = r.get("selPath$ko")
                    keggDict[key] = value.tolist()
                    r("selPaths <- append(selPaths, names(selPath))")

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDF = getKeggDF(savedDF, tempDF, DepVar, RID, stops, PID)

                base[RID] = 'Step 3 of 4: Performing GAGE analysis...'

                # save location info to session
                myDir = 'myPhyloDB/media/temp/gage/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                count_rDF = pd.DataFrame()
                if DepVar == 1:
                    finalDF['abund'] = finalDF['abund'].round(0).astype(int)
                    count_rDF = finalDF.pivot(index='rank_id', columns='sampleid', values='abund')
                elif DepVar == 4:
                    finalDF['abund_16S'] = finalDF['abund_16S'].round(0).astype(int)
                    count_rDF = finalDF.pivot(index='rank_id', columns='sampleid', values='abund_16S')

                temp_rDF = savedDF.drop_duplicates(subset='sampleid', take_last=True)

                # Removes samples (rows) that are not in our samplelist
                temp_rDF = temp_rDF.loc[temp_rDF['sampleid'].isin(catSampleIDs)]

                if metaDictCat:
                    for key in metaDictCat:
                        temp_rDF = temp_rDF.loc[temp_rDF[key].isin(metaDictCat[key])]

                # Create combined metadata column - GAGE only
                meta_rDF = temp_rDF.copy()
                if len(catFields_edit) > 1:
                    for index, row in temp_rDF.iterrows():
                        meta_rDF.loc[index, 'merge'] = "; ".join(row[catFields_edit])
                else:
                    meta_rDF.loc[:, 'merge'] = temp_rDF.loc[:, catFields_edit[0]]

                wantedList = ['sampleid', 'merge', 'sample_name']
                meta_rDF = meta_rDF.loc[:, wantedList]
                meta_rDF.set_index('sampleid', drop=True, inplace=True)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDict = {}
                r.assign("metaDF", meta_rDF)
                r("trt <- factor(metaDF$merge)")

                r.assign("count", count_rDF)
                r.assign("sampleIDs", count_rDF.columns.values.tolist())
                r("names(count) <- sampleIDs")

                r("library(DESeq2)")
                r("colData <- data.frame(row.names=colnames(count), trt=trt)")
                r("dds <- DESeqDataSetFromMatrix(countData=count, colData=colData, design = ~ trt)")

                r("sizeFactor <- rep(1, length(trt))")
                r("dds$sizeFactor <- sizeFactor")
                r("dds <- estimateDispersions(dds)")
                r("dds <- nbinomWaldTest(dds)")

                if DepVar == 1:
                    result += 'Dependent Variable: Relative Abundance' + '\n'
                elif DepVar == 4:
                    result += 'Dependent Variable: Abundance (rRNA gene copies)' + '\n'
                result += '\n===============================================\n\n\n'

                levels = list(set(meta_rDF['merge'].tolist()))
                levels = natsorted(levels, key=lambda y: y.lower())

                r("library(pathview)")
                r("library(png)")
                r("library(png)")
                r("library(grid)")
                r("pdf_counter <- 1")

                path = os.path.join('myPhyloDB', 'media', 'temp', 'gage', 'Rplots', RID)
                if not os.path.exists(path):
                    os.makedirs(path)

                r.assign("path", path)
                r("setwd(path)")
                r("options(width=5000)")
                r.assign("RID", RID)

                gageDF = pd.DataFrame(columns=['comparison', 'pathway', ' p.geomean ', ' stat.mean ', ' p.val ', ' q.val ', ' set.size '])
                diffDF = pd.DataFrame(columns=['comparison', 'kegg', ' baseMean ', ' baseMeanA ', ' baseMeanB ', ' log2FoldChange ', ' stderr ', ' pval ', ' padj '])
                for i in xrange(len(levels)-1):
                    for j in xrange(i+1, len(levels)):
                        trt1 = levels[i]
                        trt2 = levels[j]
                        r.assign("trt1", trt1)
                        r.assign("trt2", trt2)

                        # get sign based on log2FoldChange
                        r("res <- results(dds, contrast=c('trt', trt1, trt2))")
                        r("change <- -res$log2FoldChange")
                        r("names(change) <- row.names(res)")

                        # output DiffAbund to DataTable
                        r("baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,dds$trt==trt1, drop=FALSE])")
                        r("baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,dds$trt==trt2, drop=FALSE])")
                        r("df <- data.frame(kegg=rownames(res), baseMean=res$baseMean, baseMeanA=baseMeanA, baseMeanB=baseMeanB, log2FoldChange=-res$log2FoldChange, stderr=res$lfcSE, stat=res$stat, pval=res$pvalue, padj=res$padj)")

                        compDF = r.get("df")
                        comparison = str(trt1) + ' vs. ' + str(trt2)
                        compDF.insert(0, 'comparison', comparison)  # nonetype! compDF can be null?
                        diffDF = diffDF.append(compDF, ignore_index=True)
                        all_columns = compDF.columns
                        diffDF = diffDF.ix[:, all_columns]

                        ### GAGE analysis on all pathways...
                        r("gage.res <- gage(change, gsets=kegg.gs, species='ko', same.dir=FALSE)")
                        r("df <- data.frame(pathway=rownames(gage.res$greater), p.geomean=gage.res$greater[, 1], stat.mean=gage.res$greater[, 2], \
                            p.val=gage.res$greater[, 3], q.val=gage.res$greater[, 4], \
                            set.size=gage.res$greater[, 5])")

                        compDF = r.get("df")
                        compDF.insert(0, 'comparison', comparison)
                        compDF.dropna(axis=0, how='any', inplace=True)
                        gageDF = gageDF.append(compDF, ignore_index=True)

                        ### Get data way for pathview
                        # merge sign and sig to get vector (1=sig. postive, 0=not sig., -1=sig. negative)
                        r("binary <- change / abs(change)")
                        r("sig <- (res$pvalue <= 0.05)")
                        r("sig <- sig * 1")

                        r("sig <- sig * binary")
                        r("names(sig) <- row.names(res)")


                        for key in keggDict.iterkeys():
                            r.assign("pathway", key)
                            r("pid <- substr(pathway, start=1, stop=7)")
                            r("pv <- pathview(gene.data=sig, pathway.id=pid, species='ko', kegg.dir='../../../../kegg/pathways', kegg.native=T,  multi.state=F, same.layer=T, low='red', mid='gray', high='green')")

                            # convert to pdf
                            r("pdf(paste('gage_temp', pdf_counter, '.pdf', sep=''))")
                            r("plot.new()")
                            r("pngRaster <- readPNG(paste(pid, 'pathview.png', sep='.'))")
                            r("grid.raster(pngRaster, width=unit(0.8, 'npc'), height=unit(0.8, 'npc'))")
                            r("title(paste(trt1, ' vs ', trt2, sep=''))")
                            r("dev.off()")
                            r("pdf_counter <- pdf_counter + 1")

                        base[RID] = 'Step 3 of 4: Performing GAGE Analysis...\nComparison: ' + str(trt1) + ' vs ' + str(trt2) + ' is done!'

                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                        if stops[PID] == RID:
                            res = ''
                            return HttpResponse(res, content_type='application/json')
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 4 of 4: Pooling pdf files for display...'

                # Combining Pdf files
                finalFile = 'myPhyloDB/media/temp/gage/Rplots/' + str(RID) + '/gage_final.pdf'

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

                #Export tables to html
                gage_table = gageDF.to_html(classes="table display")
                gage_table = gage_table.replace('border="1"', 'border="0"')
                finalDict['gage_table'] = str(gage_table)

                diff_table = diffDF.to_html(classes="table display")
                diff_table = diff_table.replace('border="1"', 'border="0"')
                finalDict['diff_table'] = str(diff_table)

                finalDict['text'] = result
                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                removeRIDGAGE(RID)
                done[RID] = True
                return HttpResponse(res, content_type='application/json')

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with GAGE Analysis!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            removeRIDGAGE(RID)
            return HttpResponse(res, content_type='application/json')


def getKeggDF(savedDF, tempDF, DepVar, RID, stops, PID):
    global base
    try:
        # create sample and species lists based on meta data selection
        wanted = ['sampleid', 'speciesid', 'abund', 'abund_16S']
        profileDF = tempDF.loc[:, wanted]
        profileDF.set_index('speciesid', inplace=True)

        # get PICRUSt data for species
        speciesList = pd.unique(profileDF.index.ravel().tolist())
        qs = PICRUSt.objects.using('picrust').filter(speciesid__in=speciesList)
        picrustDF = read_frame(qs, fieldnames=['speciesid__speciesid', 'geneCount'])
        picrustDF.set_index('speciesid__speciesid', inplace=True)

        finalKeys = []
        total, col = picrustDF.shape
        counter = 0
        for index, row in picrustDF.iterrows():
            d = ast.literal_eval(row['geneCount'])
            for x in d.iterkeys():
                picrustDF.loc[index, x] = d[x]
                if x not in finalKeys:
                    finalKeys.append(x)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            counter += 1
            base[RID] = 'Step 2 of 4: Mapping phylotypes to KEGG pathways...phylotype ' + str(counter) + ' out of ' + str(total) + ' is done!'

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        picrustDF.drop('geneCount', axis=1, inplace=True)
        picrustDF.fillna(0, inplace=True)
        picrustDF[picrustDF > 0] = 1

        # merge to get final gene counts for all selected samples
        taxaDF = pd.merge(profileDF, picrustDF, left_index=True, right_index=True, how='inner')

        for key in finalKeys:
            if DepVar == 1:
                taxaDF[key] = taxaDF['abund'] * taxaDF[key]
            elif DepVar == 4:
                taxaDF[key] = taxaDF['abund_16S'] * taxaDF[key]

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        taxaDF = taxaDF.groupby('sampleid')[finalKeys].agg('sum')
        taxaDF.reset_index(drop=False, inplace=True)

        if DepVar == 1:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund')
        elif DepVar == 4:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund_16S')

        wanted = ['sampleid']
        metaDF = savedDF.loc[:, wanted]
        metaDF.set_index('sampleid', drop=True, inplace=True)
        grouped = metaDF.groupby(level=0)
        metaDF = grouped.last()

        taxaDF.set_index('sampleid', drop=True, inplace=True)
        finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')

        finalDF.reset_index(drop=False, inplace=True)

        return finalDF

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with GAGE Analysis!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getTabGAGE(request):
    if request.is_ajax():
        RID = request.GET["all"]

        myDir = 'myPhyloDB/media/temp/gage/'
        fileName = str(myDir) + str(RID) + '.pkl'
        savedDF = pd.read_pickle(fileName)

        myDir = 'myPhyloDB/media/temp/gage/'
        fileName = str(myDir) + str(RID) + '.csv'
        savedDF.to_csv(fileName)

        myDict = {}
        myDir = '/myPhyloDB/media/temp/gage/'
        fileName = str(myDir) + str(RID) + '.csv'
        myDict['name'] = str(fileName)
        res = simplejson.dumps(myDict)

        return HttpResponse(res, content_type='application/json')


def removeGAGEFiles(request):
    if request.is_ajax():
        RID = request.GET["all"]

        file = "myPhyloDB/media/temp/gage/" + str(RID) + ".pkl"
        if os.path.exists(file):
            os.remove(file)

        file = "myPhyloDB/media/temp/gage/" + str(RID) + ".csv"
        if os.path.exists(file):
            os.remove(file)

        return HttpResponse()


def getFullTaxonomy(level, id):
    record = []

    if level == 2:
        record = Phyla.objects.all().filter(phylaid__in=id).values_list('kingdomid_id__kingdomName', 'phylaName')
    elif level == 3:
        record = Class.objects.all().filter(classid__in=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'className')
    elif level == 4:
        record = Order.objects.all().filter(orderid__in=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderName')
    elif level == 5:
        record = Family.objects.all().filter(familyid__in=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyName')
    elif level == 6:
        record = Genus.objects.all().filter(genusid__in=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyid_id__familyName', 'genusName')
    elif level == 7:
        record = Species.objects.all().filter(speciesid__in=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyid_id__familyName', 'genusid_id__genusName', 'speciesName')

    return record


def getFullKO(level, id):
    record = []

    if level == 1:
        record = ko_lvl1.objects.using('picrust').all().filter(ko_lvl1_id__in=id).values_list('ko_lvl1_name')
    elif level == 2:
        record = ko_lvl2.objects.using('picrust').all().filter(ko_lvl2_id__in=id).values_list('ko_lvl1_id_id__ko_lvl1_name', 'ko_lvl2_name')
    elif level == 3:
        record = ko_lvl3.objects.using('picrust').all().filter(ko_lvl3_id__in=id).values_list('ko_lvl1_id_id__ko_lvl1_name', 'ko_lvl2_id_id__ko_lvl2_name', 'ko_lvl3_name')

    return record


def getFullNZ(level, id):
    record = []

    if level == 1:
        record = nz_lvl1.objects.using('picrust').all().filter(nz_lvl1_id__in=id).values_list('nz_lvl1_name')
    elif level == 2:
        record = nz_lvl2.objects.using('picrust').all().filter(nz_lvl2_id__in=id).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_name')
    elif level == 3:
        record = nz_lvl3.objects.using('picrust').all().filter(nz_lvl3_id__in=id).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_name')
    elif level == 4:
        record = nz_lvl4.objects.using('picrust').all().filter(nz_lvl4_id__in=id).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_id_id__nz_lvl3_name', 'nz_lvl4_name')
    elif level == 5:
        for item in id:
            if nz_lvl3.objects.using('picrust').all().filter(nz_lvl3_id=item).exists():
                qs = nz_lvl3.objects.using('picrust').all().filter(nz_lvl3_id=item).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_name')
                record.extend(qs)
            elif nz_lvl4.objects.using('picrust').all().filter(nz_lvl4_id=item).exists():
                qs = nz_lvl4.objects.using('picrust').all().filter(nz_lvl4_id=item).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_id_id__nz_lvl3_name', 'nz_lvl4_name')
                record.extend(qs)
            elif nz_entry.objects.using('picrust').all().filter(nz_lvl5_id=item).exists():
                qs = nz_entry.objects.using('picrust').all().filter(nz_lvl5_id=item).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_id_id__nz_lvl3_name', 'nz_lvl4_id_id__nz_lvl4_name')
                record.extend(qs)
    elif level == 6:
        record = nz_lvl4.objects.using('picrust').all().filter(nz_lvl4_id__in=id).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_id_id__nz_lvl3_name', 'nz_lvl4_name')

    return record
