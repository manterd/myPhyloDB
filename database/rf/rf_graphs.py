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

                # R packages from cran
                r("list.of.packages <- c('stargazer', 'e1071', 'randomForest', 'forestFloor', 'caret', 'NeuralNetTools', 'sparseLDA', 'pROC')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                print r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

                print r('library(caret)')
                print r('library(reshape2)')
                print r('library(stargazer)')
                print r('library(randomForest)')
                print r('library(NeuralNetTools)')
                print r('library(e1071)')

                #r('library(pROC)')
                #r('library(forestFloor)')

                # Wrangle data into R
                rankNameDF = finalDF.drop_duplicates(subset='rank_id', take_last=True)
                rankNameDF.set_index('rank_id', inplace=True)
                r.assign('rankNames', rankNameDF.rank_name.values)

                r.assign("data", count_rDF)
                r("names(data) <- rankNames")

                metaDF.drop('sample_name', axis=1, inplace=True)
                r.assign("meta", metaDF)
                r.assign("rows", metaDF.index.values.tolist())
                r("rownames(meta) <- rows")

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                r("X = data")
                r("nzv_cols <- X[,-nearZeroVar(X)]")
                r("if(length(nzv_col > 0)) X <- X[,-nzv_cols]")
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

                method = all['Method']
                finalDict = {}
                if method == 'rf':
                    r("set.seed(1)")
                    r("grid <- expand.grid(.mtry=c(1,2,5,10))")
                    if catFields:
                        r("ctrl <- trainControl(method='cv', repeats=3, number=10, \
                            classProbs=T, savePredictions=T)")
                        r("fit <- train(Y ~ ., data=myData, method='rf', linout=F, trace=F, trControl=ctrl, \
                            tuneGrid=grid, importance=T, preProcess=c('center', 'scale'))")
                    if quantFields:
                        r("ctrl <- trainControl(method='cv', repeats=3, number=10, \
                            classProbs=F, savePredictions=T)")
                        r("fit <- train(Y ~ ., data=myData, method='rf', linout=T, trace=F, trControl=ctrl, \
                            tuneGrid=grid, importance=T, preProcess=c('center', 'scale'))")

                    r("predY <- predict(fit)")

                    result += str(r('print(fit)')) + '\n'
                    result += '===============================================\n'

                    if catFields:
                        r("varDF <- filterVarImp(X, Y)")    # calculated from area under the ROC curve
                        r("rankDF <- apply(-abs(varDF), 2, rank, ties.method='min')")
                        r("rankDF <- (rankDF <= 6)")
                        r("rankDF <- rankDF * 1")
                        r("myFilter <- as.vector((rowSums(rankDF) > 0))")
                        r("fVarDF <- varDF[myFilter,]")

                        r("pdf_counter <- pdf_counter + 1")
                        r("pdf(paste(path, '/rf_temp', pdf_counter, '.pdf', sep=''), width=4+nrow(fVarDF)*0.2)")
                        r("par(mar=c(3,4,1,1),family='serif')")
                        r("fVarDF['taxa'] <- row.names(fVarDF)")
                        r("graphDF <- melt(fVarDF, id='taxa')")
                        r("p <- ggplot(graphDF, aes(x=taxa, y=value, fill=variable))")
                        r("p <- p + geom_bar(stat='identity')")
                        r("p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0))")
                        r("p <- p + labs(y='Relative Importance', x='', \
                            title='Relative Importance of Variables', \
                            subtitle='Random Forest')")
                        r("print(p)")
                        r("dev.off()")

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
                        r("rankDF <- apply(-abs(varDF), 2, rank, ties.method='min')")
                        r("rankDF <- (rankDF <= 6)")
                        r("rankDF <- rankDF * 1")
                        r("myFilter <- as.vector((rowSums(rankDF) > 0))")
                        r("Overall <- varDF[myFilter,]")
                        r("taxa <- row.names(varDF)[myFilter]")
                        r("graphDF <- data.frame(taxa, Overall)")

                        r("pdf_counter <- pdf_counter + 1")
                        r("pdf(paste(path, '/rf_temp', pdf_counter, '.pdf', sep=''), width=4+nrow(graphDF)*0.2)")
                        r("par(mar=c(3,4,1,1),family='serif')")
                        r("p <- ggplot(graphDF, aes(x=taxa, y=Overall))")
                        r("p <- p + geom_bar(stat='identity', fill='blue')")
                        r("p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0))")
                        r("p <- p + labs(y='Relative Importance', x='', \
                            title='Relative Importance of Variables', \
                            subtitle='Random Forest')")
                        r("print(p)")
                        r("dev.off()")

                        r("pdf_counter <- pdf_counter + 1")
                        r("pdf(paste(path, '/rf_temp', pdf_counter, '.pdf', sep=''))")
                        r("graphDF <- data.frame(x=Y, y=predY)")
                        r("p <- ggplot(graphDF, aes(x=x, y=y))")
                        r("p <- p + geom_point(shape=1)")
                        r("p <- p + geom_smooth(method='lm')")
                        r("p <- p + labs(y='Predicted', x='Observed', \
                            title='Predicted vs Observed', \
                            subtitle='Random Forest')")
                        r("print(p)")
                        r("dev.off()")

                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                        if stops[PID] == RID:
                            res = ''
                            return HttpResponse(res, content_type='application/json')
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                if method == 'nnet':
                    r("set.seed(1)")
                    r("grid <- expand.grid(.size=c(1, 5, 10), .decay=c(0, 0.05, 1, 2, 3, 4, 5))")
                    if catFields:
                        r("ctrl <- trainControl(method='cv', repeats=3, number=10, \
                            classProbs=T, savePredictions=T)")
                        r("fit <- train(Y ~ ., data=myData, method='nnet', linout=F, trace=F, trControl=ctrl, \
                            tuneGrid=grid, importance=T, preProcess=c('center', 'scale'))")
                    if quantFields:
                        r("ctrl <- trainControl(method='cv', repeats=3, number=10, \
                            classProbs=F, savePredictions=T)")
                        r("fit <- train(Y ~ ., data=myData, method='nnet', linout=T, trace=F, trControl=ctrl, \
                            tuneGrid=grid, importance=T, preProcess=c('center', 'scale'))")

                    r("predY <- predict(fit)")

                    result += str(r('print(fit)')) + '\n'
                    result += '===============================================\n'

                    if catFields:
                        r("varDF <- filterVarImp(X, Y)")    # calculated from area under the ROC curve
                        r("rankDF <- apply(-abs(varDF), 2, rank, ties.method='min')")
                        r("rankDF <- (rankDF <= 6)")
                        r("rankDF <- rankDF * 1")
                        r("myFilter <- as.vector((rowSums(rankDF) > 0))")
                        r("fVarDF <- varDF[myFilter,]")

                        # decision boundary
                        r("pdf_counter <- pdf_counter + 1")
                        r("pdf(paste(path, '/rf_temp', pdf_counter, '.pdf', sep=''), width=4+nrow(fVarDF)*0.2)")
                        r("hs <- 0.01")
                        r("x_min <- min(X)")
                        r("x_max <- max(X)")
                        r("y_min <- min(X)")
                        r("y_max <- max(X)")
                        r('grid <- as.matrix(expand.grid())')
                        r("dev.off()")

                        # relative importance
                        r("pdf_counter <- pdf_counter + 1")
                        r("pdf(paste(path, '/rf_temp', pdf_counter, '.pdf', sep=''), width=4+nrow(fVarDF)*0.2)")
                        r("par(mar=c(3,4,1,1),family='serif')")
                        r("fVarDF['taxa'] <- row.names(fVarDF)")
                        r("graphDF <- melt(fVarDF, id='taxa')")
                        r("p <- ggplot(graphDF, aes(x=taxa, y=value, fill=variable))")
                        r("p <- p + geom_bar(stat='identity')")
                        r("p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0))")
                        r("p <- p + labs(y='Relative Importance', x='', \
                            title='Relative Importance of Variables', \
                            subtitle='Neural Network')")
                        r("print(p)")
                        r("dev.off()")

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
                        r("rankDF <- apply(-abs(varDF), 2, rank, ties.method='min')")
                        r("rankDF <- (rankDF <= 6)")
                        r("rankDF <- rankDF * 1")
                        r("myFilter <- as.vector((rowSums(rankDF) > 0))")
                        r("Overall <- varDF[myFilter,]")
                        r("taxa <- row.names(varDF)[myFilter]")
                        r("graphDF <- data.frame(taxa, Overall)")

                        r("pdf_counter <- pdf_counter + 1")
                        r("pdf(paste(path, '/rf_temp', pdf_counter, '.pdf', sep=''), width=4+nrow(graphDF)*0.2)")
                        r("par(mar=c(3,4,1,1),family='serif')")
                        r("p <- ggplot(graphDF, aes(x=taxa, y=Overall))")
                        r("p <- p + geom_bar(stat='identity', fill='blue')")
                        r("p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0))")
                        r("p <- p + labs(y='Relative Importance', x='', \
                            title='Relative Importance of Variables', \
                            subtitle='Neural Network')")
                        r("print(p)")
                        r("dev.off()")

                        r("pdf_counter <- pdf_counter + 1")
                        r("pdf(paste(path, '/rf_temp', pdf_counter, '.pdf', sep=''))")
                        r("graphDF <- data.frame(x=Y, y=predY)")
                        r("p <- ggplot(graphDF, aes(x=x, y=y))")
                        r("p <- p + geom_point(shape=1)")
                        r("p <- p + geom_smooth(method='lm')")
                        r("p <- p + labs(y='Predicted', x='Observed', \
                            title='Predicted vs Observed', \
                            subtitle='Neural Network')")
                        r("print(p)")
                        r("dev.off()")

                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                        if stops[PID] == RID:
                            res = ''
                            return HttpResponse(res, content_type='application/json')
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    r("pdf_counter <- pdf_counter + 1")
                    r("pdf(paste(path, '/rf_temp', pdf_counter, '.pdf', sep=''), height= 5 + ncol(data)*0.1)")
                    r("plotnet(fit, cex=0.6)")
                    r("dev.off()")

                if method == 'svm':
                    r("set.seed(1)")
                    r("grid <- expand.grid(.cost=c(1))")
                    if catFields:
                        r("ctrl <- trainControl(method='cv', repeats=3, number=10, \
                            classProbs=T, savePredictions=T)")
                        r("fit <- train(Y ~ ., data=myData, method='svmLinear2', linout=F, trace=F, trControl=ctrl, \
                            tuneGrid=grid, importance=T, preProcess=c('center', 'scale'))")
                    if quantFields:
                        r("ctrl <- trainControl(method='cv', repeats=3, number=10, \
                            classProbs=F, savePredictions=T)")
                        r("fit <- train(Y ~ ., data=myData, method='svmLinear2', linout=T, trace=F, trControl=ctrl, \
                            tuneGrid=grid, importance=T, preProcess=c('center', 'scale'))")

                    r("predY <- predict(fit)")

                    result += str(r('print(fit)')) + '\n'
                    result += '===============================================\n'

                    if catFields:
                        r("varDF <- filterVarImp(X, Y)")    # calculated from area under the ROC curve
                        r("rankDF <- apply(-abs(varDF), 2, rank, ties.method='min')")
                        r("rankDF <- (rankDF <= 6)")
                        r("rankDF <- rankDF * 1")
                        r("myFilter <- as.vector((rowSums(rankDF) > 0))")
                        r("fVarDF <- varDF[myFilter,]")

                        r("pdf_counter <- pdf_counter + 1")
                        r("pdf(paste(path, '/rf_temp', pdf_counter, '.pdf', sep=''), width=4+nrow(fVarDF)*0.2)")
                        r("par(mar=c(3,4,1,1),family='serif')")
                        r("fVarDF['taxa'] <- row.names(fVarDF)")
                        r("graphDF <- melt(fVarDF, id='taxa')")
                        r("p <- ggplot(graphDF, aes(x=taxa, y=value, fill=variable))")
                        r("p <- p + geom_bar(stat='identity')")
                        r("p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0))")
                        r("p <- p + labs(y='Relative Importance', x='', \
                            title='Relative Importance of Variables', \
                            subtitle='Support Vector Machine')")
                        r("print(p)")
                        r("dev.off()")

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
                        r("rankDF <- apply(-abs(varDF), 2, rank, ties.method='min')")
                        r("rankDF <- (rankDF <= 6)")
                        r("rankDF <- rankDF * 1")
                        r("myFilter <- as.vector((rowSums(rankDF) > 0))")
                        r("Overall <- varDF[myFilter,]")
                        r("taxa <- row.names(varDF)[myFilter]")
                        r("graphDF <- data.frame(taxa, Overall)")

                        r("pdf_counter <- pdf_counter + 1")
                        r("pdf(paste(path, '/rf_temp', pdf_counter, '.pdf', sep=''), width=4+nrow(graphDF)*0.2)")
                        r("par(mar=c(3,4,1,1),family='serif')")
                        r("p <- ggplot(graphDF, aes(x=taxa, y=Overall))")
                        r("p <- p + geom_bar(stat='identity', fill='blue')")
                        r("p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0))")
                        r("p <- p + labs(y='Relative Importance', x='', \
                            title='Relative Importance of Variables', \
                            subtitle='Support Vector Machine')")
                        r("print(p)")
                        r("dev.off()")

                        r("pdf_counter <- pdf_counter + 1")
                        r("pdf(paste(path, '/rf_temp', pdf_counter, '.pdf', sep=''))")
                        r("graphDF <- data.frame(x=Y, y=predY)")
                        r("p <- ggplot(graphDF, aes(x=x, y=y))")
                        r("p <- p + geom_point(shape=1)")
                        r("p <- p + geom_smooth(method='lm')")
                        r("p <- p + labs(y='Predicted', x='Observed', \
                            title='Predicted vs Observed', \
                            subtitle='Support Vector Machine')")
                        r("print(p)")
                        r("dev.off()")

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

