import datetime
from django.http import HttpResponse
import logging
import pandas as pd
from pyper import *
import simplejson

from database.utils import multidict
from database.utils_kegg import getTaxaDF, getKeggDF, getNZDF
from database.utils_kegg import getFullTaxonomy, getFullKO, getFullNZ, insertTaxaInfo
import database.queue


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getPCA(request, stops, RID, PID):
    try:
        while True:
            if request.is_ajax():
                allJson = request.body.split('&')[0]
                all = simplejson.loads(allJson)
                database.queue.setBase(RID, 'Step 1 of 4 Selecting your chosen meta-variables...')
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
                        if idDictCat[key] not in catSampleIDs:
                            catSampleIDs.extend(idDictCat[key])

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

                finalSampleList = metaDF.sampleid.tolist()
                wantedList = allFields + ['sampleid']
                metaDF = metaDF[wantedList]
                metaDF.set_index('sampleid', drop=True, inplace=True)

                result += 'Categorical variables selected by user: ' + ", ".join(catFields) + '\n'
                result += 'Categorical variables removed from analysis (contains only 1 level): ' + ", ".join(removed) + '\n'
                result += '===============================================\n\n'

                button3 = int(all['button3'])
                DepVar = 1
                if button3 == 1:
                    DepVar = int(all["DepVar_taxa"])
                elif button3 == 2:
                    DepVar = int(all["DepVar_kegg"])
                elif button3 == 3:
                    DepVar = int(all["DepVar_nz"])

                if DepVar == 4:
                    rnaDF = savedDF.loc[savedDF['abund_16S'] != 0]
                    rows, cols = rnaDF.shape
                    if rows < 1:
                        myDict = {'error': "Error: no qPCR or 'rRNA gene copies' data were found for this dataset"}
                        res = simplejson.dumps(myDict)
                        return HttpResponse(res, content_type='application/json')

                    remSampleList = list(set(catSampleIDs) - set(finalSampleList))

                    result += str(len(remSampleList)) + " samples were removed from analysis (missing 'rRNA gene copies' data)\n"
                    result += '===============================================\n\n'

                database.queue.setBase(RID, 'Step 1 of 8: Selecting your chosen meta-variables...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 2 of 4: Selecting your chosen taxa or KEGG level...')

                finalDF = pd.DataFrame()
                if button3 == 1:
                    finalDF, missingList = getTaxaDF('rel_abund', selectAll, '', savedDF, metaDF, catFields_edit, DepVar, RID, stops, PID)
                    if selectAll == 8:
                        result += '\nThe following PGPRs were not detected: ' + ", ".join(missingList) + '\n'
                        result += '===============================================\n'

                if button3 == 2:
                    finalDF = getKeggDF('rel_abund', keggAll, '', savedDF, metaDF, catFields_edit, DepVar, RID, stops, PID)

                if button3 == 3:
                    finalDF = getNZDF('rel_abund', nzAll, '', savedDF, metaDF, catFields_edit, DepVar, RID, stops, PID)

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

                database.queue.setBase(RID, 'Step 2 of 4: Selecting your chosen taxa or KEGG level...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 3 of 4: Performing statistical test...')

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                r.assign("data", count_rDF)
                r.assign("cols", count_rDF.columns.values.tolist())
                r("colnames(data) <- cols")

                r.assign("meta", metaDF)
                r.assign("rows", metaDF.index.values.tolist())
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

                scale = all['scaled']
                if scale == 'yes':
                    r("res.pca <- PCA(data, ncp=(nrow(data)-1), scale.unit=TRUE, graph=FALSE)")
                else:
                    r("res.pca <- PCA(data, ncp=(nrow(data)-1), scale.unit=FALSE, graph=FALSE)")

                r("fviz_screeplot(res.pca)")

                r.assign("PC1", PC1)
                r.assign("PC2", PC2)

                addContrib = all['addContrib']
                if addContrib == 'yes':
                    contrib = int(all["contribVal"])
                    row, taxa = count_rDF.shape
                    samples, col = metaDF.shape
                    contrib = min(contrib, taxa, samples)
                    r.assign("contrib", contrib)
                if addContrib == 'no':
                    row, taxa = count_rDF.shape
                    samples, col = metaDF.shape
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

                rankDF = finalDF.drop_duplicates(subset='rank_id', take_last=True)
                rankDF.set_index('rank_id', inplace=True)
                r.assign('rankDF', rankDF['rank_name'])

                r("varDF <- data.frame(varnames=rownames(res.pca$var$coord), \
                    x=res.pca$var$coord[,PC1], \
                    y=res.pca$var$coord[,PC2]) \
                ")

                r('varDF <- merge(varDF, rankDF, by="row.names", all.x=FALSE)')

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
                        r("p <- p + guides(color=guide_legend('Colors'), shape=guide_legend('Symbols'))")
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
                    r("p <- p + stat_ellipse(aes(x=x, y=y, color=factor(Fill)), geom='polygon', level=0.95, alpha=0)")

                r("p <- p + geom_hline(aes(yintercept=0), linetype='dashed')")
                r("p <- p + geom_vline(aes(xintercept=0), linetype='dashed')")

                if addContrib2 == 'yes':
                    r("p <- p + geom_segment(data=varDF, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,'cm')), alpha=0.75, color='blue')")
                    r("p <- p + geom_text(data=varDF, aes(x=v1, y=v2, label=rank_name, vjust=ifelse(v2 >= 0, -1, 2)), size=2, position=position_jitter(width=0, height=0), color='blue')")

                r("p <- p + ggtitle('Biplot of variables and individuals')")
                r("p <- p + xlab(paste('Dim.', PC1, ' (', round(res.pca$eig[PC1,2], 1), '%)', sep=''))")
                r("p <- p + ylab(paste('Dim.', PC2, ' (', round(res.pca$eig[PC2,2], 1), '%)', sep=''))")
                r("print(p)")

                r("dev.off()")
                database.queue.setBase(RID, 'Step 3 of 4: Performing statistical test...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 4 of 4: Formatting graph data...')

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

                if button3 == 1:
                    zipped = getFullTaxonomy(list(varCoordDF['rank_id']))
                    insertTaxaInfo(button3, zipped, varCoordDF, pos=1)
                elif button3 == 2:
                    zipped = getFullKO(list(varCoordDF['rank_id']))
                    insertTaxaInfo(button3, zipped, varCoordDF, pos=1)
                elif button3 == 3:
                    zipped = getFullNZ(list(varCoordDF['rank_id']))
                    insertTaxaInfo(button3, zipped, varCoordDF, pos=1)

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

                if button3 == 1:
                    zipped = getFullTaxonomy(list(varContribDF['rank_id']))
                    insertTaxaInfo(button3, zipped, varContribDF, pos=1)
                elif button3 == 2:
                    zipped = getFullKO(list(varContribDF['rank_id']))
                    insertTaxaInfo(button3, zipped, varContribDF, pos=1)
                elif button3 == 3:
                    zipped = getFullNZ(list(varContribDF['rank_id']))
                    insertTaxaInfo(button3, zipped, varContribDF, pos=1)

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
                if not metaDF.empty:
                    tempDF['id'] = metaDF.index.values.tolist()
                    tempDF.set_index('id', inplace=True)
                    indCoordDF = pd.merge(metaDF, tempDF, left_index=True, right_index=True, how='inner')
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
                if not metaDF.empty:
                    tempDF['id'] = metaDF.index.values
                    tempDF.set_index('id', inplace=True)
                    indContribDF = pd.merge(metaDF, tempDF, left_index=True, right_index=True, how='inner')
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

