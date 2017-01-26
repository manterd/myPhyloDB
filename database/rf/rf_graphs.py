import datetime
from django.http import HttpResponse
import logging
from natsort import natsorted
import pandas as pd
from PyPDF2 import PdfFileReader, PdfFileMerger
from pyper import *
import simplejson

from database.utils import getMetaDF, transformDF
from database.utils_kegg import getTaxaDF, getKeggDF, getNZDF
from database.utils_kegg import filterDF
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
                path = str(myDir) + 'usr_norm_data.pkl'
                savedDF = pd.read_pickle(path)

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

                # filter otus based on user settings
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

                # transform Y, if requested
                transform = int(all["transform"])
                finalDF = transformDF(transform, DepVar, finalDF)

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
                r("list.of.packages <- c('caret', 'randomForest', 'NeuralNetTools', 'e1071', 'stargazer', 'stringr')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                print r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

                database.queue.setBase(RID, 'Step 3 of 4: Performing statistical test...')

                print r('library(caret)')
                print r('library(reshape2)')
                print r('library(RColorBrewer)')
                print r("library(plyr)")
                print r('library(stargazer)')
                print r('library(stringr)')
                print r('source("R/myFunctions/myFunctions.R")')

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

                r.assign("treeType", treeType)
                r.assign("data", count_rDF)
                r("names(data) <- rankNames")

                myList = list(metaDF.select_dtypes(include=['object']).columns)
                for i in myList:
                    metaDF[i] = metaDF[i].str.replace(' ', '_')

                r.assign("meta_full", metaDF)
                r.assign("rows", metaDF.index.values.tolist())

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                # Predictors
                r("X_full = data")
                r("nzv_cols <- nearZeroVar(X_full)")
                r("if(length(nzv_cols > 0)) X_full <- X_full[,-nzv_cols]")
                r("n.vars <- ncol(data)")

                # Response
                r.assign("allFields", allFields)
                r("Y_full = meta_full")

                # Subset train data
                trainIDs = all['trainArray']
                r.assign("trainIDs", trainIDs)
                r("X <- X_full[row.names(X_full) %in% trainIDs,]")
                r("meta <- meta_full[row.names(meta_full) %in% trainIDs,]")
                r("Y <- Y_full[row.names(Y_full) %in% trainIDs,]")
                r("Y <- Y[,allFields]")
                r("myData <- data.frame(Y, X)")

                # Subset test data
                testIDs = list(all['testArray'])
                if testIDs:
                    r.assign("testIDs", testIDs)
                    r("X_test <- X_full[row.names(X_full) %in% testIDs,]")
                    r("meta_test <- meta_full[row.names(meta_full) %in% testIDs,]")
                    r("Y_test <- Y_full[row.names(Y_full) %in% testIDs,]")
                    r("Y_test <- Y_test[,allFields]")
                    r("myData_test <- data.frame(Y_test, X_test)")
                    r("nameVec <- c('Y_test', names(X_test))")
                    r("nameVec <- make.names(nameVec)")
                    r("names(myData_test) <- nameVec")

                # Initialize R output to pdf
                path = 'myPhyloDB/media/temp/rf/Rplots/%s' % RID
                if not os.path.exists(path):
                    os.makedirs(path)

                r.assign("path", path)
                r.assign("RID", RID)
                r("pdf_counter <- 1")

                finalDict = {}

                # set up tuneGrid
                if method == 'rf':
                    r("method <- 'rf' ")
                    r("title <- 'Random Forest' ")
                    r("grid <- expand.grid(.mtry=seq(1:10))")
                elif method == 'nnet':
                    r("method <- 'nnet' ")
                    r("grid <- expand.grid(.size=seq(1:5), .decay=seq(0, 2, 0.5))")
                    r("title <- 'Neural Network' ")
                elif method == 'svm':
                    r("method <- 'svmLinear2' ")
                    r("grid <- expand.grid(.cost=seq(1:10))")
                    r("title <- 'Support Vector Machine' ")

                trainMethod = all['trainMethod']
                r.assign("trainMethod", trainMethod)

                number1 = int(all['number1'])
                r.assign("number1", number1)
                number2 = int(all['number2'])
                r.assign("number2", number2)
                repeats = int(all['repeats'])
                r.assign("repeats", repeats)
                proportion = float(all['proportion'])
                r.assign("proportion", proportion)

                if trainMethod == 'boot':
                    if catFields:
                        r("ctrl <- trainControl(method='boot', number=number2, \
                            classProbs=T, savePredictions=T)")
                    if quantFields:
                        r("ctrl <- trainControl(method='boot', number=number2, \
                            classProbs=F, savePredictions=T)")
                elif trainMethod == 'cv':
                    if catFields:
                        r("ctrl <- trainControl(method='cv', number=number1, \
                            classProbs=T, savePredictions=T)")
                    if quantFields:
                        r("ctrl <- trainControl(method='cv', number=number1, \
                            classProbs=F, savePredictions=T)")
                elif trainMethod == 'repeatedcv':
                    if catFields:
                        r("ctrl <- trainControl(method='repeatedcv', number=number1, \
                            repeats=repeats, classProbs=T, savePredictions=T)")
                    if quantFields:
                        r("ctrl <- trainControl(method='repeatedcv', number=number1, \
                            repeats=repeats, classProbs=F, savePredictions=T)")
                elif trainMethod == 'LOOCV':
                    if catFields:
                        r("ctrl <- trainControl(method='LOOCV', number=number1, \
                            classProbs=T, savePredictions=T)")
                    if quantFields:
                        r("ctrl <- trainControl(method='LOOCV', number=number1, \
                            classProbs=F, savePredictions=T)")
                elif trainMethod == 'LGOCV':
                    if catFields:
                        r("ctrl <- trainControl(method='LGOCV', number=number1, \
                            p=proportion, classProbs=T, savePredictions=T)")
                    if quantFields:
                        r("ctrl <- trainControl(method='LGOCV', number=number1, \
                            p=proportion, classProbs=F, savePredictions=T)")

                if catFields:
                    r("fit <- train(Y ~ ., data=myData, method=method, linout=F, trace=F, trControl=ctrl, \
                        tuneGrid=grid, importance=T, preProcess=c('center', 'scale'))")
                if quantFields:
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
                    r("vi <- varImp(fit, scale=F)")
                    r("varDF = as.data.frame(vi[[1]])")
                    r("goodNameVec <- names(X)")
                    r("badNameVec <- names(myData)[2:length(names(myData))]")
                    r("row.names(varDF) <- mapvalues(row.names(varDF), from=badNameVec, to=goodNameVec)")

                    r("rankDF <- apply(-varDF, 2, rank, ties.method='random')")
                    r("rankDF <- (rankDF <= 6)")
                    r("rankDF <- rankDF * 1")
                    r("myFilter <- as.vector(rowSums(rankDF) > 0)")
                    r("fVarDF <- varDF[myFilter,]")
                    r("fVarDF['rank_id'] <- row.names(fVarDF)")
                    r("graphDF <- melt(fVarDF, id='rank_id')")
                    r("graphDF$rank_id <- gsub(' ', '.', graphDF$rank_id)")
                    r("graphDF$rank_id <- gsub(':', '.', graphDF$rank_id)")
                    r("myVec <- unlist(str_split_fixed(as.character(graphDF$rank_id), '\\.id\\.\\.', 2))")
                    r("graphDF$taxa <- myVec[,1]")
                    r("graphDF$id <- myVec[,2]")
                    r("graphDF <- graphDF[with(graphDF, order(taxa, id)),] ")

                    r("pdf_counter <- pdf_counter + 1")
                    r("p <- ggplot(graphDF, aes(x=variable, y=value, fill=rank_id))")
                    r("parse_labels <- function(value) { \
                        myVec <- unlist(str_split_fixed(value, '\\.id\\.\\.', 2)); \
                        myVec[,1]; \
                    }")
                    r("p <- p + facet_wrap(~rank_id, nc=4, labeller=as_labeller(parse_labels))")
                    r("p <- p + geom_bar(stat='identity', alpha=0.9, colour='black', size=0.1)")
                    r("p <- p + theme(axis.ticks=element_line(size = 0.2))")
                    r("p <- p + theme(strip.text.x=element_text(size=7, colour='blue', angle=0))")
                    r("p <- p + theme(legend.position='none')")
                    r("p <- p + theme(axis.title.y=element_text(size=10))")
                    r("p <- p + theme(axis.text.x=element_text(size=7, angle=90, hjust=1, vjust=0.5))")
                    r("p <- p + theme(axis.text.y=element_text(size=6))")
                    r("p <- p + theme(plot.title=element_text(size=12))")
                    r("p <- p + theme(plot.subtitle=element_text(size=9))")
                    r("p <- p + labs(y='Importance', x='', \
                        title=title, \
                        subtitle='Importance (top 6 for each factor)')")

                    r("file <- paste(path, '/rf_temp', pdf_counter, '.pdf', sep='')")
                    r("char.width <- max(nchar(as.character(graphDF$rank_id)))/35")
                    r("bar.width <- 0.2+(nlevels(as.factor(graphDF$rank_id))*0.05)")
                    r("panel.width <- max(char.width, bar.width)")
                    r("n_wrap <- ceiling(nlevels(as.factor(graphDF$rank_id))/4)")
                    r("p <- set_panel_size(p, height=unit(1.25, 'in'), width=unit(panel.width, 'in'))")
                    r("char.width <- max(nchar(as.character(graphDF$variable)))/35")
                    r("ggsave(filename=file, plot=p, units='in', height=2+(1.5*n_wrap)+char.width, width=2+(4*panel.width))")

                    # graph probabilites for each training sample
                    r("probY <- predict(fit, type='prob')")
                    r("myFactors <- levels(Y)")
                    r("nFactors <- nlevels(Y)")
                    r("tempDF <- cbind(meta, probY)")
                    r("tempDF['sampleid'] = row.names(meta)")
                    r("graphDF <- melt(tempDF, id.vars=c('sampleid', 'sample_name', allFields), measure.vars=myFactors)")
                    r("names(graphDF) <- c('sampleid', 'sample_name', 'obs', 'variable', 'value')")

                    r("pdf_counter <- pdf_counter + 1")
                    r("p <- ggplot(graphDF, aes(x=sampleid, y=value, fill=variable))")
                    r("p <- p + geom_bar(stat='identity', alpha=0.9, colour='black', size=0.1)")
                    r("p <- p + facet_wrap(~ obs, scales='free_x', nc=3)")
                    r("p <- p + scale_x_discrete(labels=element_blank())")
                    r("p <- p + theme(strip.text.x=element_text(size=7, colour='blue', angle=0))")
                    r("p <- p + theme(axis.ticks.x=element_blank())")
                    r("p <- p + theme(axis.ticks.y=element_line(size = 0.2))")
                    r("p <- p + theme(legend.title=element_blank())")
                    r("p <- p + theme(legend.text=element_text(size=6))")
                    r("p <- p + theme(axis.title.y=element_text(size=10))")
                    r("p <- p + theme(axis.text.y=element_text(size=6))")
                    r("p <- p + theme(plot.title=element_text(size=12))")
                    r("p <- p + theme(plot.subtitle=element_text(size=9))")
                    r("p <- p + labs(y='Probability', x='', \
                        title=title, \
                        subtitle='Training Dataset: probabilities')")

                    r("file <- paste(path, '/rf_temp', pdf_counter, '.pdf', sep='')")
                    r("p <- set_panel_size(p, height=unit(1, 'in'), width=unit(1.5, 'in'))")
                    r("n_wrap <- ceiling(nlevels(graphDF$obs)/3)")
                    r("ggsave(filename=file, plot=p, units='in', height=2+(1.25*n_wrap), width=2+6)")

                    # graph probabilites for each test sample
                    if testIDs:
                        r("probY_test <- predict(fit, myData_test, type='prob')")
                        r("myFactors <- levels(Y_test)")
                        r("nFactors <- nlevels(Y_test)")
                        r("tempDF <- cbind(meta_test, probY_test)")
                        r("tempDF['sampleid'] = row.names(meta_test)")
                        r("graphDF <- melt(tempDF, id.vars=c('sampleid', 'sample_name', allFields), measure.vars=myFactors)")
                        r("names(graphDF) <- c('sampleid', 'sample_name', 'obs', 'variable', 'value')")

                        r("pdf_counter <- pdf_counter + 1")
                        r("p <- ggplot(graphDF, aes(x=sampleid, y=value, fill=variable))")
                        r("p <- p + geom_bar(stat='identity', alpha=0.9, colour='black', size=0.1)")
                        r("p <- p + facet_wrap(~ obs, scales='free_x', nc=3)")
                        r("p <- p + scale_x_discrete(labels=element_blank())")
                        r("p <- p + theme(strip.text.x=element_text(size=7, colour='blue', angle=0))")
                        r("p <- p + theme(axis.ticks.x=element_blank())")
                        r("p <- p + theme(axis.ticks.y=element_line(size = 0.2))")
                        r("p <- p + theme(legend.title=element_blank())")
                        r("p <- p + theme(legend.text=element_text(size=6))")
                        r("p <- p + theme(axis.title.y=element_text(size=10))")
                        r("p <- p + theme(axis.text.y=element_text(size=6))")
                        r("p <- p + theme(plot.title=element_text(size=12))")
                        r("p <- p + theme(plot.subtitle=element_text(size=9))")
                        r("p <- p + labs(y='Probability', x='', \
                            title=title, \
                            subtitle='Test dataset: assignment probabilities')")

                        r("file <- paste(path, '/rf_temp', pdf_counter, '.pdf', sep='')")
                        r("p <- set_panel_size(p, height=unit(1, 'in'), width=unit(1.5, 'in'))")
                        r("n_wrap <- ceiling(nlevels(obs)/3)")
                        r("ggsave(filename=file, plot=p, units='in', height=2+(1.25*n_wrap), width=2+6)")

                    # graph probabilites by taxa
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
                    r("graphDF$rank_id <- gsub(' ', '.', graphDF$rank_id)")
                    r("graphDF$rank_id <- gsub(':', '.', graphDF$rank_id)")
                    r("myVec <- unlist(str_split_fixed(as.character(graphDF$rank_id), '\\.id\\.\\.', 2))")
                    r("graphDF$taxa <- myVec[,1]")
                    r("graphDF$id <- myVec[,2]")
                    r("graphDF <- graphDF[with(graphDF, order(taxa, id)),] ")

                    r("pdf_counter <- pdf_counter + 1")
                    r("par(mar=c(2,2,1,1),family='serif')")
                    r("p <- ggplot(graphDF, aes(x=count, y=prob, colour=trt))")
                    r("parse_labels <- function(value) { \
                        myVec <- unlist(str_split_fixed(value, '\\.id\\.\\.', 2)); \
                        myVec[,1]; \
                    }")
                    r("p <- p + facet_wrap(~rank_id, nc=4, labeller=as_labeller(parse_labels))")
                    r("p <- p + geom_point(size=0.5)")
                    r("p <- p + theme(strip.text.x=element_text(size=7, colour='blue', angle=0))")
                    r("p <- p + theme(legend.title=element_blank())")
                    r("p <- p + theme(legend.text=element_text(size=6))")
                    r("p <- p + theme(axis.title=element_text(size=10))")
                    r("p <- p + theme(axis.text.x=element_text(size=7, angle=0))")
                    r("p <- p + theme(axis.text.y=element_text(size=6))")
                    r("p <- p + theme(plot.title=element_text(size=12))")
                    r("p <- p + theme(plot.subtitle=element_text(size=9))")
                    r("p <- p + labs(y='Probability', x='Abundance', \
                        title=title, \
                        subtitle='Probability of correct sample assignment vs taxa abundance')")

                    r("file <- paste(path, '/rf_temp', pdf_counter, '.pdf', sep='')")
                    r("char.width <- max(nchar(as.character(graphDF$rank_id)))/35")
                    r("bar.width <- 1.5")
                    r("panel.width <- max(char.width, bar.width)")
                    r("n_wrap <- ceiling(nlevels(as.factor(graphDF$rank_id))/4)")
                    r("p <- set_panel_size(p, height=unit(1.25, 'in'), width=unit(panel.width, 'in'))")
                    r("char.width <- max(nchar(as.character(graphDF$count)))/35")
                    r("ggsave(filename=file, plot=p, units='in', height=2+(1.5*n_wrap)+char.width, width=2+(4*panel.width))")

                    # confusion matrix - train
                    r("tab <- table(Observed=Y, Predicted=predY)")
                    r("cm <- confusionMatrix(tab)")

                    result += '\nTraining Dataset\n'
                    result += str(r('print(cm)')) + '\n'
                    result += '===============================================\n'

                    # confusion matrix - test
                    if testIDs:
                        r("predY_test <- predict(fit, myData_test)")
                        r("tab <- table(Observed=Y_test, Predicted=predY_test)")
                        r("cm <- confusionMatrix(tab)")

                        result += '\nTest Dataset\n'
                        result += str(r('print(cm)')) + '\n'
                        result += '===============================================\n'

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                if quantFields:
                    r("vi <- varImp(fit, scale=F)")
                    r("varDF = as.data.frame(vi[[1]])")
                    r("goodNameVec <- names(X)")
                    r("badNameVec <- names(myData)[2:length(names(myData))]")
                    r("row.names(varDF) <- mapvalues(row.names(varDF), from=badNameVec, to=goodNameVec)")

                    r("rankDF <- apply(-varDF, 2, rank, ties.method='random')")
                    r("rankDF <- (rankDF <= 10)")
                    r("rankDF <- rankDF * 1")
                    r("myFilter <- as.vector((rowSums(rankDF) > 0))")
                    r("Overall <- varDF[myFilter,]")
                    r("rank_id <- row.names(varDF)[myFilter]")
                    r("graphDF <- data.frame(rank_id, Overall)")
                    r("myVec <- unlist(str_split_fixed(as.character(graphDF$rank_id), '\\.id\\.\\.', 2))")
                    r("graphDF$taxa <- myVec[,1]")
                    r("graphDF$id <- myVec[,2]")
                    r("unique <- make.unique(as.vector(graphDF$taxa))")
                    r("graphDF['taxa2'] <- unique")

                    r("pdf_counter <- pdf_counter + 1")
                    r("par(mar=c(2,2,1,1),family='serif')")
                    r("p <- ggplot(graphDF, aes(x=taxa2, y=Overall))")
                    r("p <- p + geom_bar(stat='identity', alpha=0.9, , fill='blue', colour='black', size=0.1)")
                    r("p <- p + theme(axis.ticks=element_line(size = 0.2))")
                    r("p <- p + theme(strip.text.y=element_text(size=7, colour='blue', angle=0))")
                    r("p <- p + theme(legend.position='none')")
                    r("p <- p + theme(axis.title.y=element_text(size=10))")
                    r("p <- p + theme(axis.text.x = element_text(size=7, angle = 90))")
                    r("p <- p + theme(axis.text.y = element_text(size=6))")
                    r("p <- p + theme(plot.title = element_text(size=12))")
                    r("p <- p + theme(plot.subtitle = element_text(size=9))")
                    r("p <- p + labs(y='Importance', x='', \
                        title=title, \
                        subtitle='Overall importance (top 10)')")

                    r("file <- paste(path, '/rf_temp', pdf_counter, '.pdf', sep='')")
                    r("panel.width <- nlevels(as.factor(graphDF$rank_id))*0.2")
                    r("p <- set_panel_size(p, height=unit(2, 'in'), width=unit(panel.width, 'in'))")
                    r("charWidth <- max(nchar(graphDF$taxa2))/12")
                    r("ggsave(filename=file, plot=p, units='in', height=1+(2)+charWidth, width=1+panel.width)")

                    r("pdf_counter <- pdf_counter + 1")
                    r("graphDF <- data.frame(x=Y, y=predY)")
                    r("p <- ggplot(graphDF, aes(x=x, y=y, color='blue'))")
                    r("p <- p + geom_abline(yintercept=0, slope=1, color='gray')")
                    r("p <- p + geom_point()")
                    r("p <- p + geom_smooth(method='lm', colour='black', se=F)")
                    r("p <- p + xlim(range(Y, predY))")
                    r("p <- p + ylim(range(Y, predY))")
                    r("p <- p + labs(y='Predicted', x='Observed', \
                        title=title, \
                        subtitle='Predicted vs Observed')")

                    # Add test data if available
                    if testIDs:
                        r("probY_test <- predict(fit, myData_test, type='prob')")
                        r("predY_test <- predict(fit, myData_test)")
                        r("graphDF2 <- data.frame(x=Y_test, y=predY_test)")
                        r("p <- p + geom_point(data=graphDF2, aes(x=x, y=y, color='red'))")
                        r("p <- p + geom_smooth(data=graphDF2, method='lm', colour='black', se=F)")

                    r("p <- p + theme(legend.position='right')")
                    r("p <- p + scale_color_manual(name='Dataset', values=c('blue', 'red'), labels=c('Train', 'Test'))")
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

                if catFields:
                    r("mergeDF <- cbind(meta, probY)")
                    r("myTable <- stargazer(mergeDF, type='text', summary=F, rownames=T)")
                    result += 'Train Dataset: Probabilities\n'
                    myString = r.get("myTable")
                    result += str(myString)
                    result += '\n===============================================\n'

                    r("mergeDF <- cbind(meta_test, probY_test)")
                    r("myTable <- stargazer(mergeDF, type='text', summary=F, rownames=T)")
                    result += 'Test Dataset: Probabilities\n'
                    myString = r.get("myTable")
                    result += str(myString)
                    result += '\n===============================================\n'

                if quantFields:
                    r("mergeDF <- cbind(meta, predY)")
                    r("myTable <- stargazer(mergeDF, type='text', summary=F, rownames=T)")
                    result += 'Train Dataset: Observed vs. Predicted\n'
                    myString = r.get("myTable")
                    result += str(myString)
                    result += '\n===============================================\n'

                    r("mergeDF <- cbind(meta_test, predY_test)")
                    r("myTable <- stargazer(mergeDF, type='text', summary=F, rownames=T)")
                    result += 'Test Dataset: Observed vs. Predicted\n'
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

