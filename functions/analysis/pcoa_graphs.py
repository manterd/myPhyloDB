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


def getPCoA(request, stops, RID, PID):
    try:
        while True:
            if request.is_ajax():
                allJson = request.body.split('&')[0]
                all = json.loads(allJson)
                functions.setBase(RID, 'Step 1 of 9: Reading normalized data file...')

                functions.setBase(RID, 'Step 2 of 9: Selecting your chosen meta-variables...')
                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])

                distance = int(all["distance"])
                PC1 = int(all["PC1"])
                PC2 = int(all["PC2"])
                test = int(all["test"])
                alpha = float(all["alpha"])
                perms = int(all["perms"])

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

                if not catFields:
                    error = "Selected categorical variable(s) contain only one level.\nPlease select different variable(s)."
                    myDict = {'error': error}
                    res = json.dumps(myDict)
                    return HttpResponse(res, content_type='application/json')

                if not finalSampleIDs:
                    error = "No valid samples were contained in your final dataset.\nPlease select different variable(s)."
                    myDict = {'error': error}
                    res = json.dumps(myDict)
                    return HttpResponse(res, content_type='application/json')

                result = ''
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

                if distance == 1:
                    result += 'Distance score: Manhattan' + '\n'
                elif distance == 2:
                    result += 'Distance score: Euclidean' + '\n'
                elif distance == 3:
                    result += 'Distance score: Canberra' + '\n'
                elif distance == 4:
                    result += 'Distance score: Bray-Curtis' + '\n'
                elif distance == 5:
                    result += 'Distance score: Kulczynski' + '\n'
                elif distance == 6:
                    result += 'Distance score: Jaccard' + '\n'
                elif distance == 7:
                    result += 'Distance score: Gower' + '\n'
                elif distance == 8:
                    result += 'Distance score: altGower' + '\n'
                elif distance == 9:
                    result += 'Distance score: Morisita' + '\n'
                elif distance == 10:
                    result += 'Distance score: Horn' + '\n'
                elif distance == 11:
                    result += 'Distance score: Mountford' + '\n'
                elif distance == 12:
                    result += 'Distance score: Binomial' + '\n'
                elif distance == 13:
                    result += 'Distance score: Chao' + '\n'
                elif distance == 14:
                    result += 'Distance score: Cao' + '\n'
                elif distance == 15:
                    result += 'Distance score: wOdum' + '\n'
                    result += 'alpha: ' + str(alpha) + '\n'

                result += 'Categorical variables selected by user: ' + ", ".join(catFields + remCatFields) + '\n'
                result += 'Categorical variables not included in the statistical analysis (contains only 1 level): ' + ", ".join(remCatFields) + '\n'
                result += 'Quantitative variables selected by user: ' + ", ".join(quantFields) + '\n'
                result += '===============================================\n\n'

                functions.setBase(RID, 'Step 2 of 9: Selecting your chosen meta-variables...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 3 of 9: Selecting your chosen taxa or KEGG level...')

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
                finalDF[quantFields] = finalDF[quantFields].astype(float)

                # transform Y, if requested
                transform = int(all["transform"])
                finalDF = functions.transformDF(transform, DepVar, finalDF)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/pcoa/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)

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

                functions.setBase(RID, 'Step 3 of 9: Selecting your chosen taxa...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 4 of 9: Calculating distance matrix...')

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                functions.setBase(RID, 'Verifying R packages...missing packages are being installed')

                r("list.of.packages <- c('vegan', 'ggplot2', 'data.table')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                print r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

                functions.setBase(RID, 'Step 4 of 9: Calculating distance matrix...')

                r("options(width=5000)")
                print r("library(vegan)")
                print r("library(ggplot2)")
                print r('library(data.table)')
                print r('source("R/myFunctions/myFunctions.R")')

                count_rDF.sort_index(axis=0, inplace=True)
                r.assign("data", count_rDF)
                r.assign("cols", count_rDF.columns.values.tolist())
                r("colnames(data) <- cols")

                if distance == 1:
                    r("dist <- vegdist(data, method='manhattan')")
                elif distance == 2:
                    r("dist <- vegdist(data, method='euclidean')")
                elif distance == 3:
                    r("dist <- vegdist(data, method='canberra')")
                elif distance == 4:
                    r("dist <- vegdist(data, method='bray')")
                elif distance == 5:
                    r("dist <- vegdist(data, method='kulczynski')")
                elif distance == 6:
                    r("dist <- vegdist(data, method='jaccard')")
                elif distance == 7:
                    r("dist <- vegdist(data, method='gower')")
                elif distance == 8:
                    r("dist <- vegdist(data, method='altGower')")
                elif distance == 9:
                    r("dist <- vegdist(data, method='morisita')")
                elif distance == 10:
                    r("dist <- vegdist(data, method='horn')")
                elif distance == 11:
                    r("dist <- vegdist(data, method='mountford')")
                elif distance == 12:
                    r("dist <- vegdist(data, method='binomial')")
                elif distance == 13:
                    r("dist <- vegdist(data, method='chao')")
                elif distance == 14:
                    r("dist <- vegdist(data, method='cao')")
                elif distance == 15:
                    datamtx = np.asarray(count_rDF)
                    dists = functions.wOdum(datamtx, alpha)
                    r.assign("dist", dists)
                    r("dist <- as.dist(dist)")

                r("mat <- as.matrix(dist, diag=TRUE, upper=TRUE)")
                mat = r.get("mat")

                metaDF.sort('sampleid', inplace=True)
                rowList = metaDF.sampleid.values.tolist()
                distDF = pd.DataFrame(mat, columns=[rowList], index=rowList)

                functions.setBase(RID, 'Step 4 of 9: Calculating distance matrix...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 5 of 9: Principal coordinates analysis...')

                trtLength = len(set(catValues))
                trtString = " * ".join(catFields)

                bigf = ''
                r.assign("PC1", PC1)
                r.assign("PC2", PC2)

                addContrib2 = all['addContrib2']
                contribVal2 = float(all['contribVal2'])

                r.assign("meta", metaDF)
                method = all['Method']
                if method == 'capscale':
                    if trtLength > 0:
                        pcoa_string = "ord <- capscale(dist ~ " + str(trtString) + ", meta)"
                        r.assign("cmd", pcoa_string)
                        r("eval(parse(text=cmd))")
                    else:
                        state = "Your selected variable(s) only have one treatment level, please select additional data!"
                        myDict = {}
                        myDict['error'] = state
                        res = json.dumps(myDict)
                        return HttpResponse(res, content_type='application/json')
                elif method == 'metaMDS':
                    r('ord <- metaMDS(dist, autotransform=FALSE, trace=FALSE)')
                elif method == 'wcmdscale':
                    print r('ord <- wcmdscale(dist, eig=TRUE)')

                result += str(r('print(ord)')) + '\n'
                result += '===============================================\n'

                r('sites <- scores(ord, display="sites")')
                r("pcoa <- data.frame(meta, sites)")
                pcoaDF = r.get("pcoa")

                pcoaDF.rename(columns={'sampleid': 'Sample ID'}, inplace=True)

                eigDF = pd.DataFrame()
                if method != 'metaMDS':
                    r("Stat <- c('Eigenvalue', 'Proportion Explained', 'Cumulative Proportion')")
                    r("res <- summary(ord)")
                    r("eig <- data.frame(Stat, res$cont$importance)")
                    eigDF = r.get("eig")

                if quantFields:
                    r.assign("quantFields", quantFields)
                    r("ef <- envfit(ord, meta[,paste(quantFields)], add=False)")

                    # create dataframe from envfit for export and adding to biplot
                    r('efDF <- as.data.frame(ef$vectors$arrows*ef$vectors$r)')
                    r('efDF$r2 <- ef$vectors$r')
                    r('efDF$p <- ef$vectors$pvals')
                    r('pvals.adj <- round(p.adjust(efDF$p, method="BH"),3)')
                    r('efDF$p.adj <- pvals.adj')

                    # send data to result string
                    envfit = r("efDF")
                    result += 'EnvFit for selected quantitative variables\n'
                    result += str(envfit) + '\n'
                    result += '===============================================\n'

                colorVal = all['colorVal']
                if colorVal == 'None':
                    r("colorTrt <- c('All')")
                if colorVal == 'interaction':
                    r.assign("catFields", catFields)
                    r("colorTrt <- interaction(meta[,paste(catFields)])")
                if colorVal != 'None' and colorVal != 'interaction':
                    r.assign("colorVal", colorVal)
                    r("colorTrt <- as.factor(meta[,paste(colorVal)])")
                r("if (!exists('colorTrt')) {colorTrt <- c('All')}")

                shapeVal = all['shapeVal']
                if shapeVal == 'None':
                    r("shapeTrt <- c('All')")
                if shapeVal == 'interaction':
                    r.assign("catFields", catFields)
                    r("shapeTrt <- interaction(meta[,paste(catFields)])")
                if shapeVal != 'None' and shapeVal != 'interaction':
                    r.assign("shapeVal", shapeVal)
                    r("shapeTrt <- as.factor(meta[,paste(shapeVal)])")
                r("if (!exists('shapeTrt')) {shapeTrt <- c('All')}")

                ellipseVal = all['ellipseVal']
                if ellipseVal == 'None':
                    r("ellipseTrt <- c('All')")
                if ellipseVal != 'None' and ellipseVal != 'interaction':
                    r.assign("ellipseVal", ellipseVal)
                    r("ellipseTrt <- as.factor(meta[,paste(ellipseVal)])")
                if ellipseVal == 'interaction':
                    r.assign("catFields", catFields)
                    r("ellipseTrt <- interaction(meta[,paste(catFields)])")
                r("if (!exists('ellipseTrt')) {ellipseTrt <- c('All')}")

                surfVal = all['surfVal']
                if surfVal != 'None':
                    r.assign("surfVal", surfVal)
                    r("quant <- meta[,paste(surfVal)]")
                    r("ordi <- ordisurf(ord ~ quant, add=FALSE)")
                    r("ordi.grid <- ordi$grid")
                    r("ordi.mat <- expand.grid(x=ordi.grid$x, y=ordi.grid$y)")
                    r("ordi.mat$z <- as.vector(ordi.grid$z)")
                    r("ordi.mat <- data.frame(na.omit(ordi.mat))")

                # extract data and create dataframe for plotting
                r("indDF <- data.frame( \
                    x=as.vector(scores(ord, choices=c(PC1), display=c('sites'))), \
                    y=as.vector(scores(ord, choices=c(PC2), display=c('sites'))), \
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

                # set up plot
                r("p <- ggplot(indDF, aes(x, y))")

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

                if not surfVal == 'None':
                    r("p <- p + stat_contour(data=ordi.mat, aes(x, y, z=z, label=..level..), color='red')")
                    # get the last element in p (i.e., the one with the contour lines)
                    r("p.data <- tail(ggplot_build(p)$data, n=1)")
                    r("DT <- as.data.table(p.data[[1]], n=1)")
                    r("tmp <- unique(DT, by='level', fromLast=TRUE)")
                    r("p <- p + geom_text(aes(label=level, z=NULL), data=tmp)")

                if quantFields and addContrib2 == 'yes':
                    # scale and remove non-significant objects from efDF
                    r('names(efDF) <- c("PC1", "PC2", "r2", "p", "p.adj")')
                    r('efDF$label <- row.names(efDF)')
                    r.assign("contribVal2", contribVal2)
                    r('efDF.adj <- efDF[efDF$p.adj <= contribVal2,]')
                    r("mult <- min( max(indDF$x)-min(indDF$x), max(indDF$y)-min(indDF$y) )")
                    r('efDF.adj$v1 <- efDF.adj[,PC1] * mult * 0.7')
                    r('efDF.adj$v2 <- efDF.adj[,PC2] * mult * 0.7')
                    sigVar = r.get("nrow(efDF.adj)")
                    if sigVar >= 1:
                        r("p <- p + geom_segment(data=efDF.adj, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,'cm')), alpha=0.75, color='red')")
                        r("p <- p + geom_text(data=efDF.adj, aes(x=v1, y=v2, label=label, vjust=ifelse(v2 >= 0, -1, 2)), size=3, color='red')")

                r("p <- p + geom_hline(aes(yintercept=0), linetype='dashed')")
                r("p <- p + geom_vline(aes(xintercept=0), linetype='dashed')")

                r("p <- p + ggtitle('Principal Coordinates Analysis')")

                if method != 'metaMDS':
                    r("eig <- eigenvals(ord)")
                    r("perExp <- eig / sum(eig) * 100")
                    r("p <- p + xlab(paste('Axis', PC1, ' (', round(perExp[[PC1]], 1), '%)', sep=''))")
                    r("p <- p + ylab(paste('Axis', PC2, ' (', round(perExp[[PC2]], 1), '%)', sep=''))")
                else:
                    r("p <- p + xlab(paste('Axis', PC1, sep=''))")
                    r("p <- p + ylab(paste('Axis', PC2, sep=''))")

                path = "myPhyloDB/media/temp/pcoa/Rplots"
                if not os.path.exists(path):
                    os.makedirs(path)

                r.assign("path", path)
                r.assign("RID", RID)
                r("file <- paste(path, '/', RID, '.pcoa.pdf', sep='')")
                r("p <- set_panel_size(p, height=unit(4, 'in'), width=unit(4, 'in'))")
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
                r("ggsave(filename=file, plot=p, units='in', height=myHeight, width=myWidth)")

                functions.setBase(RID, 'Step 5 of 9: Principal coordinates analysis...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 6 of 9: Performing perMANOVA...')

                if perms < 10:
                    bigf = 'A minimum of 10 permutations is required...'
                elif len(catFields) == 0:
                    bigf = 'No categorical variables are available for perMANOVA/betaDisper analysis'
                elif perms >= 10 and len(catFields) > 0:
                    if test == 1:
                        for i in catFields:
                            factor_string = str(i) + " <- factor(meta$" + str(i) + ")"
                            r.assign("cmd", factor_string)
                            r("eval(parse(text=cmd))")

                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                            if stops[PID] == RID:
                                res = ''
                                return HttpResponse(res, content_type='application/json')
                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                        r.assign("perms", perms)
                        trtString = " * ".join(catFields)
                        amova_string = "res <- adonis(dist ~ " + str(trtString) + ", perms=perms)"
                        r.assign("cmd", amova_string)
                        r("eval(parse(text=cmd))")

                        res_aov = r("res$aov.tab")

                        tempStuff = res_aov.split('\n')
                        for part in tempStuff:
                            if part != tempStuff[0]:
                                bigf += part + '\n'
                        functions.setBase(RID, 'Step 6 of 9: Performing perMANOVA...done!')

                    elif test == 2:
                        functions.setBase(RID, 'Step 5 of 9: Principal coordinates analysis...done!')
                        functions.setBase(RID, 'Step 6 of 9: Performing BetaDisper...')

                        for i in catFields:
                            factor_string = str(i) + " <- factor(meta$" + str(i) + ")"
                            r.assign("cmd", factor_string)
                            r("eval(parse(text=cmd))")

                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                            if stops[PID] == RID:
                                res = ''
                                return HttpResponse(res, content_type='application/json')
                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                        r.assign("perms", perms)
                        for i in catFields:
                            beta_string = "res <- betadisper(dist, " + str(i) + ")"
                            r.assign("cmd", beta_string)
                            r("eval(parse(text=cmd))")

                            r("something <- anova(res)")
                            beta = r("something")
                            tempStuff = beta.split('\n')
                            bigf += 'group: ' + str(i) + '\n'
                            for part in tempStuff:
                                if part != tempStuff[0]:
                                    bigf += part + '\n'

                            betaString = str(r('res'))
                            lines = betaString.split('\n')
                            for line in lines[1:]:
                                bigf += str(line) + '\n'

                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                            if stops[PID] == RID:
                                res = ''
                                return HttpResponse(res, content_type='application/json')
                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                            functions.setBase(RID, 'Step 6 of 9: Performing BetaDisper...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 7 of 9: Formatting graph data for display...')
                finalDict = {}
                seriesList = []
                xAxisDict = {}
                yAxisDict = {}

                CAP1 = PC1 + len(catFields) + len(quantFields) + 1
                CAP2 = PC2 + len(catFields) + len(quantFields) + 1

                if catFields:
                    grouped = pcoaDF.groupby(catFields)
                    for name, group in grouped:
                        if len(catFields) > 1:
                            trt = "; ".join(name)
                        else:
                            trt = name

                        dataList = []
                        for index, row in group.iterrows():
                            dataDict = {}
                            dataDict['name'] = row['Sample ID']
                            dataDict['x'] = float(row[CAP1])
                            dataDict['y'] = float(row[CAP2])
                            dataList.append(dataDict)

                        seriesDict = {}
                        seriesDict['name'] = str(trt)
                        seriesDict['data'] = dataList
                        seriesList.append(seriesDict)

                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                        if stops[PID] == RID:
                            res = ''
                            return HttpResponse(res, content_type='application/json')
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                xTitle = {}
                if method == 'capscale':
                    xTitle['text'] = 'Axis' + str(PC1) + " (" + str(round(eigDF.iloc[1][PC1] * 100, 1)) + "%)"
                else:
                    xTitle['text'] = "Axis" + str(PC1)
                xTitle['style'] = {'fontSize': '18px', 'fontWeight': 'bold'}
                xAxisDict['title'] = xTitle

                yTitle = {}
                if method == 'capscale':
                    yTitle['text'] = 'Axis' + str(PC2) + " (" + str(round(eigDF.iloc[1][PC2] * 100, 1)) + "%)"
                else:
                    yTitle['text'] = "Axis" + str(PC2)
                yTitle['style'] = {'fontSize': '18px', 'fontWeight': 'bold'}
                yAxisDict['title'] = yTitle

                styleDict = {'style': {'fontSize': '14px'}}
                xAxisDict['labels'] = styleDict
                yAxisDict['labels'] = styleDict

                finalDict['series'] = seriesList
                finalDict['xAxis'] = xAxisDict
                finalDict['yAxis'] = yAxisDict

                if test == 1:
                    result += 'perMANOVA results:' + '\n'
                if test == 2:
                    result += 'betaDisper results:' + '\n'

                if len(catFields) == 0:
                    result += 'test cannot be run...' + '\n'
                else:
                    bigf = bigf.decode('utf-8')
                    result += bigf + '\n'

                result += '===============================================\n\n\n'

                finalDict['text'] = result

                functions.setBase(RID, 'Step 7 of 9: Formatting graph data for display...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 8 of 9: Formatting PCoA table...')

                res_table = pcoaDF.to_html(classes="table display")
                res_table = res_table.replace('border="1"', 'border="0"')
                finalDict['res_table'] = str(res_table)

                functions.setBase(RID, 'Step 8 of 9: Formatting PCoA table...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 9 of 9: Formatting distance score table...')

                distDF.sort_index(axis=1, inplace=True)
                distDF.sort_index(axis=0, inplace=True)
                dist_table = distDF.to_html(classes="table display")
                dist_table = dist_table.replace('border="1"', 'border="0"')
                finalDict['dist_table'] = str(dist_table)

                functions.setBase(RID, 'Step 9 of 9: Formatting distance score table...done!')

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
