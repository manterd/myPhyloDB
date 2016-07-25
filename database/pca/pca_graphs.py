import datetime
from django.http import HttpResponse
from django_pandas.io import read_frame
import logging
import math
import numpy as np
import pandas as pd
from pyper import *
import simplejson

from database.models import PICRUSt
from database.models import Phyla, Class, Order, Family, Genus, Species
from database.models import ko_lvl1, ko_lvl2, ko_lvl3, ko_entry
from database.models import nz_lvl1, nz_lvl2, nz_lvl3, nz_lvl4, nz_entry
from database.utils import multidict, sumKEGG
import database.queue


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}
done = {}

LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def statusPCA(request):
    global base, stage, time1, time2, TimeDiff
    if request.is_ajax():
        RID = request.GET["all"]
        time2[RID] = time.time()

        try:
            TimeDiff[RID] = time2[RID] - time1[RID]
        except:
            TimeDiff[RID] = 0

        try:
            if done[RID]:
                stage[RID] = '<br>Analysis has been running for %.1f seconds<br>Analysis is complete, results are loading' % TimeDiff[RID]
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
                        stage[RID] = str(base[RID]) + '<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
        except:
            stage[RID] = 'Analysis is initializing...'

        myDict = {'stage': stage[RID]}
        json_data = simplejson.dumps(myDict, encoding="Latin-1")
        return HttpResponse(json_data, content_type='application/json')


def removeRIDPCA(RID):
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


def getPCA(request, stops, RID, PID):
    global base, stage, time1, TimeDiff, done
    done[RID] = False
    try:
        while True:
            if request.is_ajax():
                allJson = request.body.split('&')[0]
                all = simplejson.loads(allJson)

                time1[RID] = time.time()
                base[RID] = 'Step 1 of 4 Selecting your chosen meta-variables...'

                myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
                path = str(myDir) + 'usr_norm_data.csv'

                with open(path, 'rb') as f:
                    savedDF = pd.read_csv(f, index_col=0, sep=',')

                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])

                PC1 = int(all["PC1"])
                PC2 = int(all["PC2"])
                vecScale = float(all["vecScale"])

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

                allSampleIDs = catSampleIDs
                allFields = catFields_edit

                # Removes samples (rows) that are not in our samplelist
                if allSampleIDs:
                    tempDF = savedDF.loc[savedDF['sampleid'].isin(allSampleIDs)]
                else:
                    tempDF = savedDF

                if metaDictCat:
                    for key in metaDictCat:
                        tempDF = tempDF.loc[tempDF[key].isin(metaDictCat[key])]

                if allFields:
                    wantedList = allFields + ['sampleid']
                else:
                    wantedList = ['sampleid']
                metaDF = tempDF[wantedList]

                result += 'Categorical variables selected by user: ' + ", ".join(catFields) + '\n'
                result += 'Categorical variables removed from analysis (contains only 1 level): ' + ", ".join(removed) + '\n'
                result += '===============================================\n\n'

                base[RID] = 'Step 1 of 8: Selecting your chosen meta-variables...done'

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
                    finalDF = getTaxaDF(selectAll, savedDF, metaDF, DepVar, RID, stops, PID)

                if button3 == 2:
                    DepVar = int(all["DepVar_kegg"])
                    finalDF = getKeggDF(keggAll, savedDF, tempDF, DepVar, RID, stops, PID)

                if button3 == 3:
                    DepVar = int(all["DepVar_nz"])
                    finalDF = getNZDF(nzAll, savedDF, tempDF, DepVar, RID, stops, PID)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/pca/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

                count_rDF = pd.DataFrame()
                if DepVar == 1:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='rel_abund')
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

                base[RID] = 'Step 2 of 4: Selecting your chosen taxa or KEGG level...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 3 of 4: Performing statistical test...'

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                r.assign("data", count_rDF)
                r.assign("cols", count_rDF.columns.values.tolist())
                r("colnames(data) <- cols")

                r.assign("meta", meta_rDF)
                r.assign("rows", meta_rDF.index.values.tolist())
                r("rownames(meta) <- rows")

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                path = "myPhyloDB/media/temp/pca/Rplots/" + str(RID) + ".pca.pdf"
                if os.path.exists(path):
                    os.remove(path)

                if not os.path.exists('myPhyloDB/media/temp/pca/Rplots'):
                    os.makedirs('myPhyloDB/media/temp/pca/Rplots')

                file = "pdf('myPhyloDB/media/temp/pca/Rplots/" + str(RID) + ".pca.pdf')"
                r.assign("cmd", file)
                r("eval(parse(text=cmd))")

                r("library(FactoMineR)")
                r("library(factoextra)")
                r("library(fpc)")

                r("res.pca <- PCA(data, ncp=(nrow(data)-1), graph=FALSE)")
                r("fviz_screeplot(res.pca)")

                r.assign("PC1", PC1)
                r.assign("PC2", PC2)

                addContrib = all['addContrib']
                if addContrib == 'yes':
                    contrib = int(all["contribVal"])
                    row, taxa = count_rDF.shape
                    samples, col = meta_rDF.shape
                    contrib = min(contrib, taxa, samples)
                    r.assign("contrib", contrib)
                if addContrib == 'no':
                    row, taxa = count_rDF.shape
                    samples, col = meta_rDF.shape
                    contrib = max(taxa, samples)
                    r.assign("contrib", contrib)

                r("fviz_contrib(res.pca, choice='var', axes=PC1, top=contrib) + \
                    theme(axis.text.x=element_text(angle=90, hjust=1))")
                r("fviz_contrib(res.pca, choice='var', axes=PC2, top=contrib) + \
                    theme(axis.text.x=element_text(angle=90, hjust=1))")
                r("fviz_contrib(res.pca, choice='ind', axes=PC1, top=contrib) + \
                    theme(axis.text.x=element_text(angle=90, hjust=1))")
                r("fviz_contrib(res.pca, choice='ind', axes=PC2, top=contrib) + \
                    theme(axis.text.x=element_text(angle=90, hjust=1))")


                ellipseVal = all['ellipseVal']
                if ellipseVal == 'None':
                    r("ellipseTrt <- c('All')")
                if ellipseVal != 'None' and ellipseVal != 'k-means':
                    r.assign("ellipseVal", ellipseVal)
                    r("ellipseTrt <- as.factor(meta[,paste(ellipseVal)])")
                if ellipseVal != 'None' and ellipseVal == 'k-means':
                    r("scores <- res.pca$ind$coord")
                    r("pamk.best <- pamk(scores)")
                    r("km <- kmeans(scores, centers=pamk.best$nc)")
                    r("ellipseTrt <- as.factor(paste('k-cluster: ', km$cluster, sep=''))")

                colorVal = all['colorVal']
                if colorVal == 'None':
                    r("colorTrt <- c('All')")
                if colorVal != 'None' and colorVal != 'k-means':
                    r.assign("colorVal", colorVal)
                    r("colorTrt <- as.factor(meta[,paste(colorVal)])")
                if colorVal != 'None' and colorVal == 'k-means':
                    r("scores <- res.pca$ind$coord")
                    r("pamk.best <- pamk(scores)")
                    r("km <- kmeans(scores, centers=pamk.best$nc)")
                    r("colorTrt <- as.factor(paste('k-cluster: ', km$cluster, sep=''))")

                shapeVal = all['shapeVal']
                if shapeVal == 'None':
                    r("shapeTrt <- 'All'")
                if shapeVal != 'None' and shapeVal != 'k-means':
                    r.assign("shapeVal", shapeVal)
                    r("shapeTrt <- as.factor(meta[,paste(shapeVal)])")
                if shapeVal != 'None' and shapeVal == 'k-means':
                    r("scores <- res.pca$ind$coord")
                    r("pamk.best <- pamk(scores)")
                    r("km <- kmeans(scores, centers=pamk.best$nc)")
                    r("shapeTrt <- as.factor(paste('k-cluster: ', km$cluster, sep=''))")

                r("indDF <- data.frame(x=res.pca$ind$coord[,PC1], \
                    y=res.pca$ind$coord[,PC2], \
                    Color=colorTrt, \
                    Shape=shapeTrt, \
                    Fill=ellipseTrt) \
                ")

                addContrib2 = all['addContrib2']
                if addContrib2 == 'yes':
                    contrib2 = int(all["contribVal2"])
                    r.assign("contrib2", contrib2)
                if addContrib2 == 'no':
                    r("contrib2 <- length(res.pca$var$cos2)")

                r("varDF <- data.frame(varnames=rownames(res.pca$var$coord), \
                    x=res.pca$var$coord[,PC1], \
                    y=res.pca$var$coord[,PC2]) \
                ")

                # rescale variable coordinates (from factoextra)
                r("mult <- min(max(indDF$x)-min(indDF$x)/(max(varDF$x)-min(varDF$x)), max(indDF$y)-min(indDF$y)/(max(varDF$y)-min(varDF$y)))")

                r.assign("vecScale", vecScale)

                r("varDF$v1 <- varDF$x * mult * vecScale")
                r("varDF$v2 <- varDF$y * mult * vecScale")

                # filter data to top contributors (max correlation with selected axes)
                r("x <- as.vector(apply(abs(res.pca$var$cor[,c(PC1, PC2)]), 1, max))")
                r("rank <- rank(-x, ties.method='random')")
                r("rank <- (rank <= contrib2)")
                r("varDF <- varDF[rank,]")

                # Create biplot using ggplot
                r("p <- ggplot(indDF)")

                if not colorVal == 'None':
                    if not shapeVal == 'None':
                        r("p <- p + geom_point(aes(x=x, y=y, color=as.factor(Color), shape=as.factor(Shape)), size=4)")
                        r("p <- p + scale_color_brewer(palette='Set1')")
                        r("p <- p + guides(color=guide_legend('Symbols'), shape=guide_legend('Symbols'))")
                    else:
                        r("p <- p + geom_point(aes(x=x, y=y, color=factor(Color)), size=4)")
                        r("p <- p + scale_color_brewer(palette='Set1')")
                        r("p <- p + guides(color=guide_legend('Symbols'))")
                else:
                    if not shapeVal == 'None':
                        r("p <- p + geom_point(aes(x=x, y=y, shape=factor(Shape)), size=4)")
                        r("p <- p + guides(shape=guide_legend('Symbols'))")
                    else:
                        r("p <- p + geom_point(aes(x=x, y=y), size=4)")

                if not ellipseVal == 'None':
                    r("p <- p + stat_ellipse(aes(x=x, y=y, fill=factor(Fill)), geom='polygon', level=0.95, alpha=0.2)")
                    r("p <- p + scale_fill_brewer(palette='Set1')")
                    r("p <- p + guides(fill=guide_legend('Ellipses'))")

                r("p <- p + geom_hline(aes(yintercept=0), linetype='dashed')")
                r("p <- p + geom_vline(aes(xintercept=0), linetype='dashed')")

                if addContrib2 == 'yes':
                    r("p <- p + geom_segment(data=varDF, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,'cm')), alpha=0.75, color='blue')")
                    r("p <- p + geom_text(data=varDF, aes(x=v1, y=v2, label=varnames, vjust=ifelse(v2 >= 0, -1, 2)), size=2, position=position_jitter(width=0, height=0), color='blue')")

                r("p <- p + ggtitle('Biplot of variables and individuals')")
                r("p <- p + xlab(paste('Dim.', PC1, ' (', round(res.pca$eig[PC1,2], 1), '%)', sep=''))")
                r("p <- p + ylab(paste('Dim.', PC2, ' (', round(res.pca$eig[PC2,2], 1), '%)', sep=''))")
                r("print(p)")

                r("dev.off()")
                base[RID] = 'Step 3 of 4: Performing statistical test...done'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                base[RID] = 'Step 4 of 4: Formatting graph data...'

                finalDict = {}
                r("options(width=5000)")

                result += 'Eigenvalues\n'
                out = r("res.pca$eig")
                temp = out.split('\n')
                for line in temp:
                    if line != temp[0]:
                        result += line + '\n'
                result += '\n===============================================\n'
                finalDict['text'] = result

                ## variables
                nameDF = finalDF[['rank_id']].drop_duplicates(subset='rank_id', take_last=True)
                nameDF.set_index('rank_id', inplace=True)

                r("var <- get_pca_var(res.pca)")

                r("df <- data.frame(var$coord)")
                tempDF = r.get("df")
                tempDF['id'] = count_rDF.columns.values
                tempDF.set_index('id', inplace=True)
                varCoordDF = pd.merge(nameDF, tempDF, left_index=True, right_index=True, how='inner')
                varCoordDF.reset_index(drop=False, inplace=True)
                varCoordDF.rename(columns={'index': 'rank_id'}, inplace=True)

                zipped = []
                if button3 == 1:
                    zipped = getFullTaxonomy(selectAll, varCoordDF['rank_id'])
                elif button3 == 2:
                    zipped = getFullKO(keggAll, varCoordDF['rank_id'])
                elif button3 == 3:
                    zipped = getFullNZ(nzAll, varCoordDF['rank_id'])

                if button3 == 1:
                    if selectAll == 2:
                        k, p = map(None, *zipped)
                        varCoordDF.insert(1, 'Kingdom', k)
                        varCoordDF.insert(2, 'Phyla', p)
                    elif selectAll == 3:
                        k, p, c = map(None, *zipped)
                        varCoordDF.insert(1, 'Kingdom', k)
                        varCoordDF.insert(2, 'Phyla', p)
                        varCoordDF.insert(3, 'Class', c)
                    elif selectAll == 4:
                        k, p, c, o = map(None, *zipped)
                        varCoordDF.insert(1, 'Kingdom', k)
                        varCoordDF.insert(2, 'Phyla', p)
                        varCoordDF.insert(3, 'Class', c)
                        varCoordDF.insert(4, 'Order', o)
                    elif selectAll == 5:
                        k, p, c, o, f = map(None, *zipped)
                        varCoordDF.insert(1, 'Kingdom', k)
                        varCoordDF.insert(2, 'Phyla', p)
                        varCoordDF.insert(3, 'Class', c)
                        varCoordDF.insert(4, 'Order', o)
                        varCoordDF.insert(5, 'Family', f)
                    elif selectAll == 6:
                        k, p, c, o, f, g = map(None, *zipped)
                        varCoordDF.insert(1, 'Kingdom', k)
                        varCoordDF.insert(2, 'Phyla', p)
                        varCoordDF.insert(3, 'Class', c)
                        varCoordDF.insert(4, 'Order', o)
                        varCoordDF.insert(5, 'Family', f)
                        varCoordDF.insert(6, 'Genus', g)
                    elif selectAll == 7:
                        k, p, c, o, f, g, s = map(None, *zipped)
                        varCoordDF.insert(1, 'Kingdom', k)
                        varCoordDF.insert(2, 'Phyla', p)
                        varCoordDF.insert(3, 'Class', c)
                        varCoordDF.insert(4, 'Order', o)
                        varCoordDF.insert(5, 'Family', f)
                        varCoordDF.insert(6, 'Genus', g)
                        varCoordDF.insert(7, 'Species', s)
                if button3 == 2:
                    if keggAll == 1:
                        L1 = [i[0] for i in zipped]
                        varCoordDF.insert(1, 'Level_1', L1)
                    if keggAll == 2:
                        L1, L2 = map(None, *zipped)
                        varCoordDF.insert(1, 'Level_1', L1)
                        varCoordDF.insert(2, 'Level_2', L2)
                    if keggAll == 3:
                        L1, L2, L3 = map(None, *zipped)
                        varCoordDF.insert(1, 'Level_1', L1)
                        varCoordDF.insert(2, 'Level_2', L2)
                        varCoordDF.insert(3, 'Level_3', L3)
                if button3 == 3:
                    if nzAll == 1:
                        L1 = [i[0] for i in zipped]
                        varCoordDF.insert(1, 'Level_1', L1)
                    if nzAll == 2:
                        L1, L2 = map(None, *zipped)
                        varCoordDF.insert(1, 'Level_1', L1)
                        varCoordDF.insert(2, 'Level_2', L2)
                    if nzAll == 3:
                        L1, L2, L3 = map(None, *zipped)
                        varCoordDF.insert(1, 'Level_1', L1)
                        varCoordDF.insert(2, 'Level_2', L2)
                        varCoordDF.insert(3, 'Level_3', L3)
                    if nzAll == 4:
                        L1, L2, L3, L4 = map(None, *zipped)
                        varCoordDF.insert(1, 'Level_1', L1)
                        varCoordDF.insert(2, 'Level_2', L2)
                        varCoordDF.insert(3, 'Level_3', L3)
                        varCoordDF.insert(4, 'Level_4', L4)
                    if nzAll == 5:
                        L1, L2, L3, L4 = map(None, *zipped)
                        varCoordDF.insert(1, 'Level_1', L1)
                        varCoordDF.insert(2, 'Level_2', L2)
                        varCoordDF.insert(3, 'Level_3', L3)
                        varCoordDF.insert(4, 'Level_4', L4)
                    if nzAll == 6:
                        L1, L2, L3, L4 = map(None, *zipped)
                        varCoordDF.insert(1, 'Level_1', L1)
                        varCoordDF.insert(2, 'Level_2', L2)
                        varCoordDF.insert(3, 'Level_3', L3)
                        varCoordDF.insert(4, 'Level_4', L4)

                    varCoordDF.fillna(value='N/A', inplace=True)

                table = varCoordDF.to_html(classes="table display")
                table = table.replace('border="1"', 'border="0"')
                finalDict['varCoordDF'] = str(table)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                r("df <- data.frame(var$contrib)")
                tempDF = r.get("df")
                tempDF['id'] = count_rDF.columns.values
                tempDF.set_index('id', inplace=True)
                varContribDF = pd.merge(nameDF, tempDF, left_index=True, right_index=True, how='inner')
                varContribDF.reset_index(drop=False, inplace=True)
                varContribDF.rename(columns={'index': 'rank_id'}, inplace=True)

                zipped = []
                if button3 == 1:
                    zipped = getFullTaxonomy(selectAll, varCoordDF['rank_id'])
                elif button3 == 2:
                    zipped = getFullKO(keggAll, varCoordDF['rank_id'])
                elif button3 == 3:
                    zipped = getFullNZ(nzAll, varCoordDF['rank_id'])

                if button3 == 1:
                    if selectAll == 2:
                        k, p = map(None, *zipped)
                        varContribDF.insert(1, 'Kingdom', k)
                        varContribDF.insert(2, 'Phyla', p)
                    elif selectAll == 3:
                        k, p, c = map(None, *zipped)
                        varContribDF.insert(1, 'Kingdom', k)
                        varContribDF.insert(2, 'Phyla', p)
                        varContribDF.insert(3, 'Class', c)
                    elif selectAll == 4:
                        k, p, c, o = map(None, *zipped)
                        varContribDF.insert(1, 'Kingdom', k)
                        varContribDF.insert(2, 'Phyla', p)
                        varContribDF.insert(3, 'Class', c)
                        varContribDF.insert(4, 'Order', o)
                    elif selectAll == 5:
                        k, p, c, o, f = map(None, *zipped)
                        varContribDF.insert(1, 'Kingdom', k)
                        varContribDF.insert(2, 'Phyla', p)
                        varContribDF.insert(3, 'Class', c)
                        varContribDF.insert(4, 'Order', o)
                        varContribDF.insert(5, 'Family', f)
                    elif selectAll == 6:
                        k, p, c, o, f, g = map(None, *zipped)
                        varContribDF.insert(1, 'Kingdom', k)
                        varContribDF.insert(2, 'Phyla', p)
                        varContribDF.insert(3, 'Class', c)
                        varContribDF.insert(4, 'Order', o)
                        varContribDF.insert(5, 'Family', f)
                        varContribDF.insert(6, 'Genus', g)
                    elif selectAll == 7:
                        k, p, c, o, f, g, s = map(None, *zipped)
                        varContribDF.insert(1, 'Kingdom', k)
                        varContribDF.insert(2, 'Phyla', p)
                        varContribDF.insert(3, 'Class', c)
                        varContribDF.insert(4, 'Order', o)
                        varContribDF.insert(5, 'Family', f)
                        varContribDF.insert(6, 'Genus', g)
                        varContribDF.insert(7, 'Species', s)
                if button3 == 2:
                    if keggAll == 1:
                        L1 = [i[0] for i in zipped]
                        varContribDF.insert(1, 'Level_1', L1)
                    if keggAll == 2:
                        L1, L2 = map(None, *zipped)
                        varContribDF.insert(1, 'Level_1', L1)
                        varContribDF.insert(2, 'Level_2', L2)
                    if keggAll == 3:
                        L1, L2, L3 = map(None, *zipped)
                        varContribDF.insert(1, 'Level_1', L1)
                        varContribDF.insert(2, 'Level_2', L2)
                        varContribDF.insert(3, 'Level_3', L3)
                if button3 == 3:
                    if nzAll == 1:
                        L1 = [i[0] for i in zipped]
                        varContribDF.insert(1, 'Level_1', L1)
                    if nzAll == 2:
                        L1, L2 = map(None, *zipped)
                        varContribDF.insert(1, 'Level_1', L1)
                        varContribDF.insert(2, 'Level_2', L2)
                    if nzAll == 3:
                        L1, L2, L3 = map(None, *zipped)
                        varContribDF.insert(1, 'Level_1', L1)
                        varContribDF.insert(2, 'Level_2', L2)
                        varContribDF.insert(3, 'Level_3', L3)
                    if nzAll == 4:
                        L1, L2, L3, L4 = map(None, *zipped)
                        varContribDF.insert(1, 'Level_1', L1)
                        varContribDF.insert(2, 'Level_2', L2)
                        varContribDF.insert(3, 'Level_3', L3)
                        varContribDF.insert(4, 'Level_4', L4)
                    if nzAll == 5:
                        L1, L2, L3, L4 = map(None, *zipped)
                        varContribDF.insert(1, 'Level_1', L1)
                        varContribDF.insert(2, 'Level_2', L2)
                        varContribDF.insert(3, 'Level_3', L3)
                        varContribDF.insert(4, 'Level_4', L4)
                    if nzAll == 6:
                        L1, L2, L3, L4 = map(None, *zipped)
                        varContribDF.insert(1, 'Level_1', L1)
                        varContribDF.insert(2, 'Level_2', L2)
                        varContribDF.insert(3, 'Level_3', L3)
                        varContribDF.insert(4, 'Level_4', L4)

                    varContribDF.fillna(value='N/A', inplace=True)

                table = varContribDF.to_html(classes="table display")
                table = table.replace('border="1"', 'border="0"')
                finalDict['varContribDF'] = str(table)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                r("ind <- get_pca_ind(res.pca)")

                if ellipseVal == 'k-means' or colorVal == 'k-means' or shapeVal == 'k-means':
                    r("df <- data.frame(km$cluster, ind$coord)")
                else:
                    r("df <- data.frame(ind$coord)")

                tempDF = r.get("df")
                if not meta_rDF.empty:
                    tempDF['id'] = meta_rDF.index.values.tolist()
                    tempDF.set_index('id', inplace=True)
                    indCoordDF = pd.merge(meta_rDF, tempDF, left_index=True, right_index=True, how='inner')
                    indCoordDF.reset_index(drop=False, inplace=True)
                    indCoordDF.rename(columns={'index': 'rank_id', ' km.cluster ': 'k-means cluster'}, inplace=True)
                else:
                    indCoordDF = tempDF.copy()
                    indCoordDF.rename(columns={' km.cluster ': 'k-means cluster'}, inplace=True)
                table = indCoordDF.to_html(classes="table display")
                table = table.replace('border="1"', 'border="0"')
                finalDict['indCoordDF'] = str(table)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                r("df <- data.frame(ind$contrib)")
                tempDF = r.get("df")
                if not meta_rDF.empty:
                    tempDF['id'] = meta_rDF.index.values
                    tempDF.set_index('id', inplace=True)
                    indContribDF = pd.merge(meta_rDF, tempDF, left_index=True, right_index=True, how='inner')
                    indContribDF.reset_index(drop=False, inplace=True)
                    indContribDF.rename(columns={'index': 'rank_id'}, inplace=True)
                else:
                    indContribDF = tempDF.copy()
                table = indContribDF.to_html(classes="table display")
                table = table.replace('border="1"', 'border="0"')
                finalDict['indContribDF'] = str(table)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                removeRIDPCA(RID)
                done[RID] = True
                return HttpResponse(res, content_type='application/json')

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "Error with pca!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = simplejson.dumps(myDict)
            removeRIDPCA(RID)
            return HttpResponse(res, content_type='application/json')


def getTaxaDF(selectAll, savedDF, metaDF, DepVar, RID, stops, PID):
    global base
    try:
        base[RID] = 'Step 2 of 4: Selecting your chosen taxa or KEGG level...'
        taxaDF = pd.DataFrame(columns=['sampleid', 'rank', 'rank_id', 'rank_name', 'rel_abund', 'abund_16S', 'rich', 'diversity'])

        if selectAll == 2:
            taxaDF = savedDF.loc[:, ['sampleid', 'phylaid', 'phylaName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'phylaid': 'rank_id', 'phylaName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Phyla'
        elif selectAll == 3:
            taxaDF = savedDF.loc[:, ['sampleid', 'classid', 'className', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'classid': 'rank_id', 'className': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Class'
        elif selectAll == 4:
            taxaDF = savedDF.loc[:, ['sampleid', 'orderid', 'orderName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'orderid': 'rank_id', 'orderName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Order'
        elif selectAll == 5:
            taxaDF = savedDF.loc[:, ['sampleid', 'familyid', 'familyName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'familyid': 'rank_id', 'familyName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Family'
        elif selectAll == 6:
            taxaDF = savedDF.loc[:, ['sampleid', 'genusid', 'genusName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'genusid': 'rank_id', 'genusName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Genus'
        elif selectAll == 7:
            taxaDF = savedDF.loc[:, ['sampleid', 'speciesid', 'speciesName', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'speciesid': 'rank_id', 'speciesName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Species'

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        taxaDF.drop('sampleid', axis=1, inplace=True)
        finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')

        wantedList = ['sampleid', 'rank', 'rank_name', 'rank_id']
        if DepVar == 1:
            finalDF = finalDF.groupby(wantedList)[['rel_abund']].sum()
        elif DepVar == 4:
            finalDF = finalDF.groupby(wantedList)[['abund_16S']].sum()
        elif DepVar == 2:
            finalDF = finalDF.groupby(wantedList)[['rich']].sum()
        elif DepVar == 3:
            finalDF = finalDF.groupby(wantedList)[['diversity']].sum()

        finalDF.reset_index(drop=False, inplace=True)
        return finalDF

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "Error with pca!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getKeggDF(keggAll, savedDF, tempDF, DepVar, RID, stops, PID):
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
        wanted = ['sampleid', 'speciesid', 'rel_abund', 'abund_16S', 'rich', 'diversity']
        profileDF = tempDF.loc[:, wanted]
        profileDF.set_index('speciesid', inplace=True)

        # get PICRUSt data for species
        speciesList = pd.unique(profileDF.index.ravel().tolist())
        qs = PICRUSt.objects.using('picrust').filter(speciesid__in=speciesList)
        picrustDF = read_frame(qs, fieldnames=['speciesid__speciesid', 'geneCount'])
        picrustDF.set_index('speciesid__speciesid', inplace=True)

        levelList = []
        for key in koDict:
            levelList.append(str(key))

        picrustDF = pd.concat([picrustDF, pd.DataFrame(columns=levelList)])
        picrustDF.fillna(0.0, inplace=True)
        picrustDF = sumKEGG(speciesList, picrustDF, koDict, RID, PID, stops)
        picrustDF.drop('geneCount', axis=1, inplace=True)
        picrustDF[picrustDF > 0.0] = 1.0

        # merge to get final gene counts for all selected samples
        taxaDF = pd.merge(profileDF, picrustDF, left_index=True, right_index=True, how='inner')

        for level in levelList:
            if DepVar == 1:
                taxaDF[level] = taxaDF['rel_abund'] * taxaDF[level]
            elif DepVar == 2:
                taxaDF[level] = np.where(taxaDF['rel_abund'] * taxaDF[level] > 0, 1, 0)
            elif DepVar == 3:
                taxaDF[level] = taxaDF['rel_abund'] * taxaDF[level]
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
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rel_abund')
        elif DepVar == 4:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund_16S')
        elif DepVar == 2:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rich')
        elif DepVar == 3:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='diversity')

        wanted = ['sampleid']
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
            myDict['error'] = "Error with pca!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getNZDF(nzAll, savedDF, tempDF, DepVar, RID, stops, PID):
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
        wanted = ['sampleid', 'speciesid', 'rel_abund', 'abund_16S', 'rich', 'diversity']
        profileDF = tempDF.loc[:, wanted]
        profileDF.set_index('speciesid', inplace=True)

        # get PICRUSt data for species
        speciesList = pd.unique(profileDF.index.ravel().tolist())
        qs = PICRUSt.objects.using('picrust').filter(speciesid__in=speciesList)
        picrustDF = read_frame(qs, fieldnames=['speciesid__speciesid', 'geneCount'])
        picrustDF.set_index('speciesid__speciesid', inplace=True)

        levelList = []
        for key in nzDict:
            levelList.append(str(key))

        picrustDF = pd.concat([picrustDF, pd.DataFrame(columns=levelList)])
        picrustDF.fillna(0.0, inplace=True)
        picrustDF = sumKEGG(speciesList, picrustDF, nzDict, RID, PID, stops)
        picrustDF.drop('geneCount', axis=1, inplace=True)
        picrustDF[picrustDF > 0.0] = 1.0

        # merge to get final gene counts for all selected samples
        taxaDF = pd.merge(profileDF, picrustDF, left_index=True, right_index=True, how='inner')
        for level in levelList:
            if DepVar == 1:
                taxaDF[level] = taxaDF['rel_abund'] * taxaDF[level]
            elif DepVar == 2:
                taxaDF[level] = np.where(taxaDF['rel_abund'] * taxaDF[level] > 0, 1, 0)
            elif DepVar == 3:
                taxaDF[level] = taxaDF['rel_abund'] * taxaDF[level]
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
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rel_abund')
        elif DepVar == 4:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund_16S')
        elif DepVar == 2:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rich')
        elif DepVar == 3:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='diversity')

        wanted = ['sampleid']
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
            myDict['error'] = "Error with pca!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def removePCAFiles(request):
    if request.is_ajax():
        RID = request.GET["all"]

        file = "myPhyloDB/media/temp/pca/Rplots/" + str(RID) + ".pca.pdf"
        if os.path.exists(file):
            os.remove(file)

        file = "myPhyloDB/media/temp/pca/" + str(RID) + ".pkl"
        if os.path.exists(file):
            os.remove(file)

        file = "myPhyloDB/media/temp/pca/" + str(RID) + ".csv"
        if os.path.exists(file):
            os.remove(file)

        return HttpResponse()


def getTabPCA(request):
    if request.is_ajax():
        RID = request.GET["all"]
        myDir = 'myPhyloDB/media/temp/pca/'
        fileName = str(myDir) + str(RID) + '.pkl'
        savedDF = pd.read_pickle(fileName)

        myDir = 'myPhyloDB/media/temp/pca/'
        fileName = str(myDir) + str(RID) + '.csv'
        savedDF.to_csv(fileName)

        myDict = {}
        myDir = '/myPhyloDB/media/temp/pca/'
        fileName = str(myDir) + str(RID) + '.csv'
        myDict['name'] = str(fileName)
        res = simplejson.dumps(myDict)

        return HttpResponse(res, content_type='application/json')


def getTabWGCNA(request):
    if request.is_ajax():
        RID = request.GET["all"]
        myDir = 'myPhyloDB/media/temp/wgcna/'
        fileName = str(myDir) + str(RID) + '.pkl'
        savedDF = pd.read_pickle(fileName)

        myDir = 'myPhyloDB/media/temp/wgcna/'
        fileName = str(myDir) + str(RID) + '.csv'
        savedDF.to_csv(fileName)

        myDict = {}
        myDir = '/myPhyloDB/media/temp/wgcna/'
        fileName = str(myDir) + str(RID) + '.csv'
        myDict['name'] = str(fileName)
        res = simplejson.dumps(myDict)

        return HttpResponse(res, content_type='application/json')


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
