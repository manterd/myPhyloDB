import datetime
from django.http import HttpResponse
import logging
from natsort import natsorted
import numpy as np
import pandas as pd
from PyPDF2 import PdfFileReader, PdfFileMerger
from pyper import *
import simplejson

from database.utils import getMetaDF, transformDF
from database.utils_kegg import getTaxaDF, getKeggDF, getNZDF
from database.utils_kegg import getFullTaxonomy, getFullKO, getFullNZ, insertTaxaInfo, filterDF
import database.queue


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getRF(request, stops, RID, PID):
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
                savedDF, metaDF, finalSampleIDs, catFields, remCatFields, quantFields, catValues, quantValues = getMetaDF(savedDF, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, DepVar)
                allFields = catFields + quantFields

                result = ''
                result += 'Categorical variables selected by user: ' + ", ".join(catFields + remCatFields) + '\n'
                result += 'Categorical variables not included in the statistical analysis (contains only 1 level): ' + ", ".join(remCatFields) + '\n'
                result += 'Quantitative variables selected by user: ' + ", ".join(quantFields) + '\n'
                result += '===============================================\n\n'

                database.queue.setBase(RID, 'Step 1 of 8: Selecting your chosen meta-variables...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 2 of 4: Selecting your chosen taxa or KEGG level...')

                # filter phylotypes based on user settings
                remUnclass = all['remUnclass']
                remZeroes = all['remZeroes']
                perZeroes = int(all['perZeroes'])
                filterData = all['filterData']
                filterPer = int(all['filterPer'])
                filterMeth = int(all['filterMeth'])

                finalDF = pd.DataFrame()
                if treeType == 1:
                    if selectAll != 8:
                        filteredDF = filterDF(savedDF, DepVar, selectAll, remUnclass, remZeroes, perZeroes, filterData, filterPer, filterMeth)
                    else:
                        filteredDF = savedDF.copy()

                    finalDF, missingList = getTaxaDF(selectAll, '', filteredDF, metaDF, allFields, DepVar, RID, stops, PID)

                    if selectAll == 8:
                        result += '\nThe following PGPRs were not detected: ' + ", ".join(missingList) + '\n'
                        result += '===============================================\n'

                if treeType == 2:
                    finalDF, allDF = getKeggDF(keggAll, '', savedDF, metaDF, allFields, DepVar, RID, stops, PID)

                if treeType == 3:
                    finalDF, allDF = getNZDF(nzAll, '', savedDF, metaDF, allFields, DepVar, RID, stops, PID)

                # make sure column types are correct
                finalDF[catFields] = finalDF[catFields].astype(str)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/rf/'
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

                database.queue.setBase(RID, 'Verifying R packages...missing packages are being installed')

                # R packages from cran
                r("list.of.packages <- c('caret', 'randomForest', 'NeuralNetTools', 'e1071', 'stargazer')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                print r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

                database.queue.setBase(RID, 'Step 3 of 4: Performing statistical test...')

                print r('library(caret)')
                print r('library(reshape2)')
                print r('library(RColorBrewer)')
                print r('library(stargazer)')
                print r('source("database/myFunctions.R")')

                method = all['Method']
                if method == 'rf':
                    print r('library(randomForest)')
                elif method == 'nnet':
                    print r('library(NeuralNetTools)')
                elif method == 'svm':
                    print r('library(e1071)')

                # Wrangle data into R
                rankNameDF = finalDF.drop_duplicates(subset='rank_id', take_last=True)
                rankNameDF.sort(columns='rank_id', inplace=True)
                rankNameDF['name_id'] = rankNameDF[['rank_name', 'rank_id']].apply(lambda x: ' id: '.join(x), axis=1)
                r.assign('rankNames', rankNameDF.name_id.values)

                r.assign("data", count_rDF)
                r("names(data) <- rankNames")

                metaDF.drop('sample_name', axis=1, inplace=True)
                myList = list(metaDF.select_dtypes(include=['object']).columns)
                for i in myList:
                    metaDF[i] = metaDF[i].str.replace(' ', '_')

                r.assign("meta", metaDF)
                r.assign("rows", metaDF.index.values.tolist())
                r("rownames(meta) <- rows")

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                r("X = data")
                r("nzv_cols <- nearZeroVar(X)")
                r("if(length(nzv_cols > 0)) X <- X[,-nzv_cols]")
                r.assign("allFields", allFields)
                r("Y = meta[,allFields]")
                r("n <- names(data)")
                r("n.vars <- ncol(data)")
                r("myData <- data.frame(Y, X)")

                # Initialize R output to pdf
                path = 'myPhyloDB/media/temp/rf/Rplots/%s' % RID
                if not os.path.exists(path):
                    os.makedirs(path)

                r.assign("path", path)
                r.assign("RID", RID)
                r("pdf_counter <- 1")

                finalDict = {}

                # set up tuneGrid
                r("set.seed(1)")
                if method == 'rf':
                    r("method <- 'rf' ")
                    r("title <- 'Random Forest' ")
                    r("grid <- expand.grid(.mtry=c(1, 2, 5))")
                elif method == 'nnet':
                    r("method <- 'nnet' ")
                    r("grid <- expand.grid(.size=c(1, 5, 10), .decay=c(0, 0.05, 1, 2, 3, 4, 5))")
                    r("title <- 'Neural Network' ")
                elif method == 'svm':
                    r("method <- 'svmLinear2' ")
                    r("grid <- expand.grid(.cost=c(1, 5, 10))")
                    r("title <- 'Support Vector Machine' ")

                if catFields:
                    r("ctrl <- trainControl(method='cv', repeats=3, number=10, \
                        classProbs=T, savePredictions=T)")
                    r("fit <- train(Y ~ ., data=myData, method=method, linout=F, trace=F, trControl=ctrl, \
                        tuneGrid=grid, importance=T, preProcess=c('center', 'scale'))")
                if quantFields:
                    r("ctrl <- trainControl(method='cv', repeats=3, number=10, \
                        classProbs=F, savePredictions=T)")
                    r("fit <- train(Y ~ ., data=myData, method=method, linout=T, trace=F, trControl=ctrl, \
                        tuneGrid=grid, importance=T, preProcess=c('center', 'scale'))")

                r("predY <- predict(fit)")

                r("if (exists('predY')) {fitError <- FALSE} else {fitError <- TRUE}")
                fitError = r.get("fitError")
                if fitError:
                    myDict = {'error': "Model could not be fit:\nPlease try a different model"}
                    res = simplejson.dumps(myDict)
                    return HttpResponse(res, content_type='application/json')
                else:
                    result += str(r('print(fit)')) + '\n'
                    result += '===============================================\n'

                if catFields:
                    r("varDF <- filterVarImp(X, Y)")
                    r("rankDF <- apply(-abs(varDF), 2, rank, ties.method='random')")
                    r("rankDF <- (rankDF <= 6)")
                    r("rankDF <- rankDF * 1")
                    r("myFilter <- as.vector(rowSums(rankDF) > 0)")
                    r("fVarDF <- varDF[myFilter,]")
                    r("fVarDF['rank_id'] <- row.names(fVarDF)")
                    r("graphDF <- melt(fVarDF, id='rank_id')")
                    r("foo <- data.frame(do.call('rbind', strsplit(as.character(graphDF$rank_id), ' id: ', fixed=TRUE)))")
                    r("graphDF['taxa'] <- foo$X1")

                    r("pdf_counter <- pdf_counter + 1")
                    r("par(mar=c(2,2,1,1),family='serif')")
                    r("p <- ggplot(graphDF, aes(x=taxa, y=value, fill=variable))")
                    r("p <- p + geom_bar(stat='identity', alpha=0.9, colour='black', size=0.1)")
                    r("p <- p + facet_grid(variable ~ .)")
                    r("p <- p + theme(axis.ticks=element_line(size = 0.2))")
                    r("p <- p + theme(strip.text.y=element_text(size=5, colour='blue', angle=0))")
                    r("p <- p + theme(legend.position='none')")
                    r("p <- p + theme(axis.title.y=element_text(size=8))")
                    r("p <- p + theme(axis.text.x = element_text(size=5, angle = 90, hjust = 0))")
                    r("p <- p + theme(axis.text.y = element_text(size=4))")
                    r("p <- p + theme(plot.title = element_text(size=10))")
                    r("p <- p + theme(plot.subtitle = element_text(size=7))")
                    r("p <- p + labs(y='Importance', x='', \
                        title=title, \
                        subtitle='Relative importance (top 6 for each factor)')")

                    r("file <- paste(path, '/rf_temp', pdf_counter, '.pdf', sep='')")
                    r("panel.width <- nlevels(graphDF$taxa)*0.4")
                    r("p <- set_panel_size(p, height=unit(0.75, 'cm'), width=unit(panel.width, 'cm'))")
                    r("ggsave(filename=file, plot=p, units='cm', height=4+nlevels(graphDF$variable), width=5+panel.width)")

                    # graph probabilites
                    r("probY <- predict(fit, type='prob')")
                    r("newX <- X[,myFilter]")
                    r("tempDF <- cbind(newX, Y, probY)")
                    r("tempDF['sampleid'] <- row.names(tempDF)")
                    r("myFactors <- levels(Y)")
                    r("myTaxa <- names(newX)")

                    r("df1 <- melt(tempDF, id.vars=c('sampleid', 'Y', myTaxa), measure.vars=myFactors)")
                    r("names(df1) <- c('sampleid', 'Y', myTaxa, 'variable', 'prob')")
                    r("dfeq <- df1[df1$Y==df1$variable,]")

                    r("graphDF <- melt(dfeq, id.vars=c('variable', 'prob'), measure.vars=myTaxa)")
                    r("names(graphDF) <- c('trt', 'prob', 'rank_id', 'count')")
                    r("foo <- data.frame(do.call('rbind', strsplit(as.character(graphDF$rank_id), ' id: ', fixed=TRUE)))")
                    r("graphDF['taxa'] <- foo$X1")

                    r("pdf_counter <- pdf_counter + 1")
                    r("par(mar=c(2,2,1,1),family='serif')")
                    r("p <- ggplot(graphDF, aes(x=count, y=prob, colour=trt))")
                    r("p <- p + geom_point(size=0.5)")
                    r("p <- p + facet_grid(taxa ~ .)")
                    r("p <- p + theme(strip.text.y=element_text(size=5, colour='blue', angle=0))")
                    #r("p <- p + scale_x_log10()")
                    r("p <- p + theme(axis.title=element_text(size=8))")
                    r("p <- p + theme(axis.text.x = element_text(size=5, angle = 90, hjust = 0))")
                    r("p <- p + theme(axis.text.y = element_text(size=4))")
                    r("p <- p + theme(plot.title = element_text(size=10))")
                    r("p <- p + theme(plot.subtitle = element_text(size=7))")
                    r("p <- p + labs(y='Probability', x='Abundance', \
                        title=title, \
                        subtitle='Probability by factor')")

                    r("file <- paste(path, '/rf_temp', pdf_counter, '.pdf', sep='')")
                    r("p <- set_panel_size(p, height=unit(1.5, 'cm'), width=unit(4, 'cm'))")
                    r("ggsave(filename=file, plot=p, units='cm', height=1.5*length(myTaxa)+7, width=10)")

                    # confusion matrix
                    r("tab <- table(predY, Y)")
                    r("cm <- confusionMatrix(tab)")

                    result += str(r('print(cm)')) + '\n'
                    result += '===============================================\n'

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                if quantFields:
                    r("varDF <- filterVarImp(X, Y)")
                    r("rankDF <- apply(-abs(varDF), 2, rank, ties.method='random')")
                    r("rankDF <- (rankDF <= 10)")
                    r("rankDF <- rankDF * 1")
                    r("myFilter <- as.vector((rowSums(rankDF) > 0))")
                    r("Overall <- varDF[myFilter,]")
                    r("rank_id <- row.names(varDF)[myFilter]")
                    r("graphDF <- data.frame(rank_id, Overall)")
                    r("foo <- data.frame(do.call('rbind', strsplit(as.character(graphDF$rank_id), ' id: ', fixed=TRUE)))")
                    r("graphDF['taxa'] <- foo$X1")

                    r("pdf_counter <- pdf_counter + 1")
                    r("par(mar=c(2,2,1,1),family='serif')")
                    r("p <- ggplot(graphDF, aes(x=taxa, y=Overall))")
                    r("p <- p + geom_bar(stat='identity', alpha=0.9, , fill='blue', colour='black', size=0.1)")
                    r("p <- p + theme(axis.ticks=element_line(size = 0.2))")
                    r("p <- p + theme(strip.text.y=element_text(size=5, colour='blue', angle=0))")
                    r("p <- p + theme(legend.position='none')")
                    r("p <- p + theme(axis.title.y=element_text(size=8))")
                    r("p <- p + theme(axis.text.x = element_text(size=5, angle = 90, hjust = 0))")
                    r("p <- p + theme(axis.text.y = element_text(size=4))")
                    r("p <- p + theme(plot.title = element_text(size=10))")
                    r("p <- p + theme(plot.subtitle = element_text(size=7))")
                    r("p <- p + labs(y='Importance', x='', \
                        title=title, \
                        subtitle='Overall importance (top 10)')")

                    r("file <- paste(path, '/rf_temp', pdf_counter, '.pdf', sep='')")
                    r("panel.width <- nlevels(graphDF$taxa)*0.5")
                    r("p <- set_panel_size(p, height=unit(2, 'cm'), width=unit(panel.width, 'cm'))")
                    r("ggsave(filename=file, plot=p, units='cm', height=5, width=5+panel.width)")

                    r("pdf_counter <- pdf_counter + 1")
                    r("graphDF <- data.frame(x=Y, y=predY)")
                    r("p <- ggplot(graphDF, aes(x=x, y=y))")
                    r("p <- p + geom_point(shape=1)")
                    r("p <- p + geom_smooth(method='lm')")
                    r("p <- p + labs(y='Predicted', x='Observed', \
                        title=title, \
                        subtitle='Predicted vs Observed')")
                    r("file <- paste(path, '/rf_temp', pdf_counter, '.pdf', sep='')")
                    r("ggsave(filename=file, plot=p, units='in', height=4, width=4)")

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                r("myTable <- stargazer(varDF, type='text', summary=F, rownames=T)")
                result += 'Variable Importance\n'
                myString = r.get("myTable")
                result += str(myString)
                result += '\n===============================================\n'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                # Neural Network Plot
                #if method == 'nnet':
                #    r("pdf_counter <- pdf_counter + 1")
                #    r("pdf(paste(path, '/rf_temp', pdf_counter, '.pdf', sep=''), height= 5 + ncol(data)*0.1)")
                #    r("plotnet(fit, cex=0.6)")
                #    r("dev.off()")

                # Combining Pdf files
                finalFile = 'myPhyloDB/media/temp/rf/Rplots/' + str(RID) + '/rf_final.pdf'

                pdf_files = [f for f in os.listdir(path) if f.endswith("pdf")]
                pdf_files = natsorted(pdf_files, key=lambda y: y.lower())

                merger = PdfFileMerger()
                for filename in pdf_files:
                    merger.append(PdfFileReader(os.path.join(path, filename), 'rb'))

                merger.write(finalFile)

                database.queue.setBase(RID, 'Step 3 of 4: Performing statistical test...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 4 of 4: Formatting graph data...')
                r("options(width=5000)")
                finalDict['text'] = result

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

