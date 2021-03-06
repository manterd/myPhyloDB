import datetime
from django.http import HttpResponse
import logging
import numpy as np
import pandas as pd
from pyper import *
import json

import functions


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getPCA(request, stops, RID, PID):
    try:
        while True:
            if request.is_ajax():
                allJson = request.body.split('&')[0]
                all = json.loads(allJson)
                functions.setBase(RID, 'Step 1 of 5: Reading normalized data file...')

                functions.setBase(RID, 'Step 2 of 5 Selecting your chosen meta-variables...')
                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])

                method = all["Method"]
                scale = all['scaled']
                constrain = all["constrain"]
                PC1 = int(all["PC1"])
                PC2 = int(all["PC2"])

                result = ''
                treeType = int(all['treeType'])
                if treeType == 1:
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
                    elif selectAll == 9:
                        result += 'Taxa level: OTU_99' + '\n'
                elif treeType == 2:
                    if keggAll == 1:
                        result += 'KEGG Pathway level: 1' + '\n'
                    elif keggAll == 2:
                        result += 'KEGG Pathway level: 2' + '\n'
                    elif keggAll == 3:
                        result += 'KEGG Pathway level: 3' + '\n'
                elif treeType == 3:
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
                metaValsQuant = all['metaValsQuant']
                metaIDsQuant = all['metaIDsQuant']

                treeType = int(all['treeType'])
                DepVar = int(all["DepVar"])

                # Create meta-variable DataFrame, final sample list, final category and quantitative field lists based on tree selections
                savedDF, metaDF, finalSampleIDs, catFields, remCatFields, quantFields, catValues, quantValues = functions.getMetaDF(request.user, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, DepVar)
                allFields = catFields + quantFields

                result = ''
                result += 'Categorical variables selected by user: ' + ", ".join(catFields + remCatFields) + '\n'
                result += 'Categorical variables not included in the statistical analysis (contains only 1 level): ' + ", ".join(remCatFields) + '\n'
                result += 'Quantitative variables selected by user: ' + ", ".join(quantFields) + '\n'
                result += '===============================================\n\n'

                functions.setBase(RID, 'Step 2 of 5: Selecting your chosen meta-variables...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 3 of 5: Selecting your chosen taxa or KEGG level...')

                # filter otus based on user settings
                remUnclass = all['remUnclass']
                remZeroes = all['remZeroes']
                perZeroes = int(all['perZeroes'])
                filterData = all['filterData']
                filterPer = int(all['filterPer'])
                filterMeth = int(all['filterMeth'])
                mapTaxa = 'no'

                finalDF = pd.DataFrame()
                if treeType == 1:
                    if selectAll != 8:
                        filteredDF = functions.filterDF(savedDF, DepVar, selectAll, remUnclass, remZeroes, perZeroes, filterData, filterPer, filterMeth)
                    else:
                        filteredDF = savedDF.copy()

                    finalDF, missingList = functions.getTaxaDF(selectAll, '', filteredDF, metaDF, allFields, DepVar, RID, stops, PID)

                    if selectAll == 8:
                        result += '\nThe following PGPRs were not detected: ' + ", ".join(missingList) + '\n'
                        result += '===============================================\n'

                if treeType == 2:
                    finalDF, allDF = functions.getKeggDF(keggAll, '', savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)

                if treeType == 3:
                    finalDF, allDF = functions.getNZDF(nzAll, '', savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)

                if finalDF.empty:
                    error = "Selected taxa were not found in your selected samples."
                    myDict = {'error': error}
                    res = json.dumps(myDict)
                    return HttpResponse(res, content_type='application/json')

                # make sure column types are correct
                finalDF[catFields] = finalDF[catFields].astype(str)

                # transform Y, if requested
                transform = int(all["transform"])
                finalDF = functions.transformDF(transform, DepVar, finalDF)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/pca/'
                if not os.path.exists(myDir):
                    os.makedirs(myDir)

                path = str(myDir) + str(RID) + '.biom'
                functions.imploding_panda(path, treeType, DepVar, finalSampleIDs, metaDF, finalDF)

                # STOP SHARED
                # START STATSGRAPH

                count_rDF = pd.DataFrame()
                if DepVar == 0:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='abund')
                elif DepVar == 1:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='rel_abund')
                elif DepVar == 2:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='rich')
                elif DepVar == 3:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='diversity')
                elif DepVar == 4:
                    count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='abund_16S')

                count_rDF.fillna(0, inplace=True)

                functions.setBase(RID, 'Step 3 of 5: Selecting your chosen taxa or KEGG level...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 4 of 5: Performing statistical test...')

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                functions.setBase(RID, 'Verifying R packages...missing packages are being installed')

                r("list.of.packages <- c('fpc', 'vegan', 'ggplot2')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                print r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

                functions.setBase(RID, 'Step 4 of 5: Performing statistical test...')

                r("options(width=5000)")
                print r('library(fpc)')
                print r('library(ggplot2)')
                print r('library(vegan)')
                print r('source("R/myFunctions/myFunctions.R")')

                count_rDF.sort_index(axis=0, inplace=True)
                r.assign("data", count_rDF)
                r.assign("cols", count_rDF.columns.values.tolist())
                r("colnames(data) <- unlist(cols)")

                metaDF.sort_values('sampleid', inplace=True)
                r.assign("meta", metaDF)
                r.assign("rows", metaDF.index.values.tolist())
                r("rownames(meta) <- unlist(rows)")

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                r.assign("PC1", PC1)
                r.assign("PC2", PC2)

                # Can only constrain if meta-variables have been selected
                constrain2 = 'no'
                if constrain == 'yes':
                    if catFields or quantFields:
                        constrain2 = 'yes'
                    else:
                        constrain2 = 'no'

                if not method == 'decorana':
                    if constrain2 == 'no':
                        if scale == 'yes':
                            pca_string = 'res.pca <- ' + method + '(data, scale=TRUE)'
                            r.assign("cmd", pca_string)
                            r("eval(parse(text=cmd))")
                        else:
                            pca_string = 'res.pca <- ' + method + '(data, scale=FALSE)'
                            r.assign("cmd", pca_string)
                            r("eval(parse(text=cmd))")

                    if constrain2 == 'yes':
                        if scale == 'yes':
                            pca_string = 'res.pca <- ' + method + '(data ~ ., data=meta, scale=TRUE)'
                            r.assign("cmd", pca_string)
                            r("eval(parse(text=cmd))")
                        else:
                            pca_string = 'res.pca <- ' + method + '(data ~ ., data=meta, scale=FALSE)'
                            r.assign("cmd", pca_string)
                            r("eval(parse(text=cmd))")

                if method == 'decorana':
                    pca_string = 'res.pca <- ' + method + '(data)'
                    r.assign("cmd", pca_string)
                    r("eval(parse(text=cmd))")

                res = r.get('res.pca')
                if res is None:
                    error = "Your analysis failed due to infinite or missing values.\nPlease try tranforming your data and/or selecting different samples."
                    myDict = {'error': error}
                    res = json.dumps(myDict)
                    return HttpResponse(res, content_type='application/json')

                result += str(r('print(res.pca)')) + '\n'
                result += '===============================================\n'

                addContrib1 = all['addContrib1']
                contribVal1 = float(all['contribVal1'])
                addContrib2 = all['addContrib2']
                contribVal2 = float(all['contribVal2'])

                # Use vegan to calculate regression between ord axes and quantFields
                if addContrib1 == 'yes':
                    r("ef1 <- envfit(res.pca, data, add=False)")

                if quantFields and addContrib2 == 'yes':
                    r.assign('quantFields', quantFields)
                    r('ef2 <- envfit(res.pca, meta[,paste(quantFields)], add=False)')

                # get scores from vegan
                r('sites <- scores(res.pca, display="sites", choices=c(PC1,PC2))')
                r('species <- scores(res.pca, display="species", choices=c(PC1,PC2))')

                ellipseVal = all['ellipseVal']
                if ellipseVal == 'None':
                    r("ellipseTrt <- c('All')")
                if ellipseVal == 'interaction':
                    r.assign("catFields", catFields)
                    r("ellipseTrt <- interaction(meta[,paste(catFields)])")
                if ellipseVal != 'None' and ellipseVal != 'k-means' and ellipseVal != 'interaction':
                    r.assign("ellipseVal", ellipseVal)
                    r("ellipseTrt <- as.factor(meta[,paste(ellipseVal)])")
                if ellipseVal != 'None' and ellipseVal == 'k-means':
                    r("pamk.best <- pamk(sites)")
                    r("km <- kmeans(sites, centers=pamk.best$nc)")
                    r("ellipseTrt <- as.factor(paste('k-cluster: ', km$cluster, sep=''))")
                r("if (!exists('ellipseTrt')) {ellipseTrt <- c('All')}")

                colorVal = all['colorVal']
                if colorVal == 'None':
                    r("colorTrt <- c('All')")
                if colorVal == 'interaction':
                    r.assign("catFields", catFields)
                    r("colorTrt <- interaction(meta[,paste(catFields)])")
                if colorVal != 'None' and colorVal != 'k-means' and colorVal != 'interaction':
                    r.assign("colorVal", colorVal)
                    r("colorTrt <- as.factor(meta[,paste(colorVal)])")
                if colorVal != 'None' and colorVal == 'k-means':
                    r("pamk.best <- pamk(sites)")
                    r("km <- kmeans(sites, centers=pamk.best$nc)")
                    r("colorTrt <- as.factor(paste('k-cluster: ', km$cluster, sep=''))")
                r("if (!exists('colorTrt')) {colorTrt <- c('All')}")

                shapeVal = all['shapeVal']
                if shapeVal == 'None':
                    r("shapeTrt <- 'All'")
                if shapeVal == 'interaction':
                    r.assign("catFields", catFields)
                    r("shapeTrt <- interaction(meta[,paste(catFields)])")
                if shapeVal != 'None' and shapeVal != 'k-means' and shapeVal != 'interaction':
                    r.assign("shapeVal", shapeVal)
                    r("shapeTrt <- as.factor(meta[,paste(shapeVal)])")
                if shapeVal != 'None' and shapeVal == 'k-means':
                    r("pamk.best <- pamk(sites)")
                    r("km <- kmeans(sites, centers=pamk.best$nc)")
                    r("shapeTrt <- as.factor(paste('k-cluster: ', km$cluster, sep=''))")
                r("if (!exists('shapeTrt')) {shapeTrt <- c('All')}")

                r("indDF <- data.frame( \
                    x=sites[,PC1], \
                    y=sites[,PC2], \
                    Color=colorTrt, \
                    Shape=shapeTrt, \
                    Fill=ellipseTrt) \
                ")

                gridVal_X = all['gridVal_X']
                if gridVal_X != 'None':
                    r.assign("gridVal_X", gridVal_X)
                    r("indDF$myGrid_X <- meta[,paste(gridVal_X)]")

                gridVal_Y = all['gridVal_Y']
                if gridVal_Y != 'None':
                    r.assign("gridVal_Y", gridVal_Y)
                    r("indDF$myGrid_Y <- meta[,paste(gridVal_Y)]")

                r("varDF <- data.frame( \
                    x=species[,PC1], \
                    y=species[,PC2]) \
                ")

                # get taxa rank names
                rankNameDF = finalDF.drop_duplicates(subset='rank_id', keep='last')
                rankNameDF.set_index('rank_id', inplace=True)
                rankNameDF['rank_name'] = rankNameDF['rank_name'].str.split('|').str[-1]
                r.assign('rankNameDF', rankNameDF['rank_name'])
                r('varDF <- merge(varDF, rankNameDF, by="row.names", all.x=TRUE)')

                # rescale
                r("mult <- min(max(indDF$x)-min(indDF$x)/(max(varDF$x)-min(varDF$x)), max(indDF$y)-min(indDF$y)/(max(varDF$y)-min(varDF$y)))")

                # Create biplot using ggplot
                r("p <- ggplot(indDF, aes(x,y))")

                if gridVal_X != 'None' and gridVal_Y == 'None':
                    r("p <- p + facet_grid(. ~ myGrid_X)")
                    r("p <- p + theme(strip.text.x=element_text(size=10, colour='blue', angle=0))")
                elif gridVal_X == 'None' and gridVal_Y != 'None':
                    r("p <- p + facet_grid(myGrid_Y ~ .)")
                    r("p <- p + theme(strip.text.y=element_text(size=10, colour='blue', angle=90))")
                elif gridVal_X != 'None' and gridVal_Y != 'None':
                    r("p <- p + facet_grid(myGrid_Y ~ myGrid_X)")
                    r("p <- p + theme(strip.text.x=element_text(size=10, colour='blue', angle=0))")
                    r("p <- p + theme(strip.text.y=element_text(size=10, colour='blue', angle=90))")

                myPalette = all['palette']
                r.assign("myPalette", myPalette)

                r('number <- nlevels(indDF$Shape)')
                r('shapes <- rep(c(21, 22, 23, 24, 25), length.out = number) ')

                if not colorVal == 'None':
                    if not shapeVal == 'None':
                        r("p <- p + geom_point(aes(fill=factor(Color), shape=factor(Shape)), size=4)")
                        r("p <- p + scale_fill_brewer(name='Symbol-colors', palette=myPalette, guide=guide_legend(override.aes=list(shape=21)))")
                        r("p <- p + scale_shape_manual(name='Symbol-shapes', values=shapes)")
                    else:
                        r("p <- p + geom_point(aes(fill=factor(Color)), shape=21, size=4)")
                        r("p <- p + scale_fill_brewer(name='Symbol-colors', palette=myPalette, guide=guide_legend(override.aes=list(shape=21)))")
                else:
                    if not shapeVal == 'None':
                        r("p <- p + geom_point(aes(shape=factor(Shape)), size=4)")
                        r("p <- p + scale_shape_manual(name='Symbol-shapes', values=shapes)")
                    else:
                        r("p <- p + geom_point(color='gray', size=4)")

                if not ellipseVal == 'None':
                    myCI = float(all["CI"])
                    r.assign("myCI", myCI)
                    r("p <- p + stat_ellipse(aes(color=factor(Fill)), geom='polygon', level=myCI, alpha=0)")
                    r("p <- p + scale_color_brewer(palette=myPalette)")
                    r("p <- p + guides(color=guide_legend('Ellipse-colors'))")

                r("p <- p + geom_hline(aes(yintercept=0), linetype='dashed')")
                r("p <- p + geom_vline(aes(xintercept=0), linetype='dashed')")

                if addContrib1 == 'yes':
                    r('efDF <- as.data.frame(ef1$vectors$arrows*ef1$vectors$r)')
                    r('efDF$p <- ef1$vectors$pvals')
                    r('pvals.adj <- round(p.adjust(efDF$p, method="BH"),3)')
                    r('efDF$p.adj <- pvals.adj')
                    r('efDF <- efDF[ order(row.names(efDF)), ]')
                    r('row.names(varDF) <- varDF$Row.names')
                    r('varDF <- varDF[ order(row.names(varDF)), ]')
                    r('efDF$label <- varDF$rank_name')

                    # scale and remove non-significant objects
                    r.assign("contribVal1", contribVal1)
                    r('efDF.adj <- efDF[efDF$p <= paste(contribVal1),]')
                    r('efDF.adj$v1 <- efDF.adj[,PC1] * mult * 0.7')
                    r('efDF.adj$v2 <- efDF.adj[,PC2] * mult * 0.7')
                    efDF_adj = r.get("efDF.adj")

                    # send data to result string
                    envfit = r("efDF")
                    result += 'Envfit results for species scores\n'
                    result += str(envfit) + '\n'
                    result += '===============================================\n'

                    if not efDF_adj.empty:
                        r("p <- p + geom_segment(data=efDF.adj, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,'cm')), alpha=0.75, color='blue')")
                        r("p <- p + geom_text(data=efDF.adj, aes(x=v1, y=v2, label=label, vjust=ifelse(v2 >= 0, -1, 2)), size=3, color='blue')")

                if quantFields and addContrib2 == 'yes':
                    r('efDF <- data.frame(ef2$vectors$arrows*sqrt(ef2$vectors$r))')
                    r('efDF$p <- ef2$vectors$pvals')
                    r('pvals.adj <- round(p.adjust(efDF$p, method="BH"),3)')
                    r('efDF$p.adj <- pvals.adj')
                    r('efDF$label <- unlist(quantFields)')

                    # scale and remove non-significant objects
                    r.assign("contribVal2", contribVal2)
                    r('efDF.adj <- efDF[efDF$p < paste(contribVal2),]')
                    r('efDF.adj$v1 <- efDF.adj[,PC1] * mult * 0.7')
                    r('efDF.adj$v2 <- efDF.adj[,PC2] * mult * 0.7')
                    efDF_adj = r.get('efDF.adj')

                    # send data to result string
                    envfit = r("efDF")
                    result += 'EnvFit for selected quantitative variables\n'
                    result += str(envfit) + '\n'
                    result += '===============================================\n'

                    if not efDF_adj.empty:
                        r("p <- p + geom_segment(data=efDF.adj, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,'cm')), alpha=0.75, color='red')")
                        r("p <- p + geom_text(data=efDF.adj, aes(x=v1, y=v2, label=label, vjust=ifelse(v2 >= 0, -1, 2)), size=3, color='red')")

                # add labels to plot
                r("p <- p + ggtitle('Biplot of variables and individuals')")
                if method != 'decorana':
                    r("eig <- eigenvals(res.pca)")
                else:
                    r("eig <- res.pca$evals")
                r("perExp <- eig / sum(eig) * 100")
                r("p <- p + xlab(paste('Axis', PC1, ' (', round(perExp[[PC1]], 1), '%)', sep=''))")
                r("p <- p + ylab(paste('Axis', PC2, ' (', round(perExp[[PC2]], 1), '%)', sep=''))")

                path = "myPhyloDB/media/temp/pca/Rplots"
                if not os.path.exists(path):
                    os.makedirs(path)

                r.assign("path", path)
                r.assign("RID", RID)
                r("file <- paste(path, '/', RID, '.pca.pdf', sep='')")
                r("p <- set_panel_size(p, height=unit(2.9, 'in'), width=unit(2.9, 'in'))")
                r("nlev <- nlevels(as.factor(indDF$myGrid_X))")
                r('if (nlev == 0) { \
                        myWidth <- 8 \
                    } else { \
                        myWidth <- 3*nlev+4 \
                }')
                r("nlev <- nlevels(as.factor(indDF$myGrid_Y))")
                r('if (nlev == 0) { \
                        myHeight <- 8 \
                    } else { \
                        myHeight <- 3*nlev+4 \
                }')
                r("ggsave(filename=file, plot=p, units='in', height=myHeight, width=myWidth, limitsize=F)")

                functions.setBase(RID, 'Step 4 of 5: Performing statistical test...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 5 of 5: Formatting graph data...')

                finalDict = {}
                r("options(width=5000)")
                finalDict['text'] = result

                ## variables
                nameDF = finalDF[['rank_id']].drop_duplicates(subset='rank_id', keep='last')
                nameDF.set_index('rank_id', inplace=True)

                r("df <- data.frame(species)")
                tempDF = r.get("df")
                IDs = r.get("row.names(df)")
                tempDF['id'] = IDs
                tempDF.set_index('id', inplace=True)
                varCoordDF = pd.merge(nameDF, tempDF, left_index=True, right_index=True, how='inner')
                varCoordDF.reset_index(drop=False, inplace=True)
                varCoordDF.rename(columns={'index': 'rank_id'}, inplace=True)

                if treeType == 1:
                    idList = functions.getFullTaxonomy(list(varCoordDF.rank_id.unique()))
                    varCoordDF['Taxonomy'] = varCoordDF['rank_id'].map(idList)
                elif treeType == 2:
                    idList = functions.getFullKO(list(varCoordDF.rank_id.unique()))
                    varCoordDF['Taxonomy'] = varCoordDF['rank_id'].map(idList)
                elif treeType == 3:
                    idList = functions.getFullNZ(list(varCoordDF.rank_id.unique()))
                    varCoordDF['Taxonomy'] = varCoordDF['rank_id'].map(idList)

                varCoordDF.replace(to_replace='N/A', value=np.nan, inplace=True)
                varCoordDF.dropna(axis=1, how='all', inplace=True)
                table = varCoordDF.to_html(classes="table display")
                table = table.replace('border="1"', 'border="0"')
                finalDict['varCoordDF'] = str(table)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                if ellipseVal == 'k-means' or colorVal == 'k-means' or shapeVal == 'k-means':
                    r("df <- data.frame(km$cluster, sites)")
                else:
                    r("df <- data.frame(sites)")

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

                finalDict['error'] = 'none'
                res = json.dumps(finalDict)
                return HttpResponse(res, content_type='application/json')

    except Exception as e:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "There was an error during your analysis:\nError: " + str(e.message) + "\nTimestamp: " + str(datetime.datetime.now())
            res = json.dumps(myDict)
            return HttpResponse(res, content_type='application/json')

