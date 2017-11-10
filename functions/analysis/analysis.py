# goals: create extendable class for all analyses,
# with functions covering duplicate sections (with arguments for variances)

# class objects created within queue call, destroyed after returning from main (more likely called RUN)
# almost all variables become local instance variables created in constructor
# this way all variables are visible to all functions (pseudo-global)
# PID and RID get passed in alongside Request as arguments for constructor
# each analysis's RUN function calls the sequence of helper methods involved

# make each analysis its own class inheriting from main analysis template (abstract) class
# ie myAnalysis = new Analysis(request, rid, pid)
# results = myAnalysis.runAnova(potential extra args)


# required sections:
#   request validation
#   variable checking (from request dict, which settings are active)
#   getTaxaDF (or equivalent) call
#
#   interface with R?
#   differentiate between quant and cat analyses

from abc import ABCMeta, abstractmethod
import json
import functions
from django.http import HttpResponse
from functions.utils.utils_df import recLabels
import zipfile
import pandas as pd
from pyper import *
import os




# for now, analysis on its own is basically anova. Meaning you should be able to create an Anova and use run()
# to complete the analysis entirely. At the same time, analysis should be split into functions to make it easy to
# overwrite the anova specific components where needed
class Analysis:  # abstract parent class, not to be run on its own. Instead, should be used as a template
    __metaclass__ = ABCMeta

    def __init__(self, iRequest, iRID, iStops, iPID):
        self.request = iRequest
        self.RID = iRID
        self.stopList = iStops
        self.PID = iPID

        # 'declare' all variables used in multiple steps here
        self.all, self.selectAll, self.treeType, self.savedDF, self.DepVar, self.metaDF, self.allFields = (None,)*7
        self.result, self.keggAll, self.nzAll, self.catFields, self.quantFields, self.finalSampleIDs = (None,)*6
        self.finalDict, self.seriesList, self.xAxisDict, self.yAxisDict, self.finalDF, self.pValDict = (None,)*6
        self.sig_only, self.transform, self.finalDict, self.zipFile, self.mapTaxa, self.allDF = (None,)*6

    @abstractmethod
    def run(self):
        print "Really shouldn't be able to see this"

    def validate(self):
        print "Validate!"
        # Get variables from web page
        allJson = self.request.body.split('&')[0]
        self.all = json.loads(allJson)
        functions.setBase(self.RID, 'Step 1 of 4: Selecting your chosen meta-variables...')

        self.selectAll = int(self.all["selectAll"])
        self.keggAll = int(self.all["keggAll"])
        self.nzAll = int(self.all["nzAll"])
        self.sig_only = int(self.all["sig_only"])

        metaValsCat = self.all['metaValsCat']
        metaIDsCat = self.all['metaIDsCat']
        metaValsQuant = self.all['metaValsQuant']
        metaIDsQuant = self.all['metaIDsQuant']

        self.treeType = int(self.all['treeType'])
        self.DepVar = int(self.all["DepVar"])

        # Create meta-variable DataFrame, final sample list, final category and quantitative field lists based on tree selections
        self.savedDF, self.metaDF, self.finalSampleIDs, self.catFields, remCatFields, self.quantFields, self.catValues, self.quantValues = functions.getMetaDF(self.request.user, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, self.DepVar, levelDep=True)
        self.allFields = self.catFields + self.quantFields

        if not self.catFields:
            error = "Selected categorical variable(s) contain only one level.\nPlease select different variable(s)."
            myDict = {'error': error}
            res = json.dumps(myDict)
            return HttpResponse(res, content_type='application/json')

        if not self.finalSampleIDs:
            error = "No valid samples were contained in your final dataset.\nPlease select different variable(s)."
            myDict = {'error': error}
            res = json.dumps(myDict)
            return HttpResponse(res, content_type='application/json')

        result = ''
        result += 'Categorical variables selected by user: ' + ", ".join(self.catFields + remCatFields) + '\n'
        result += 'Categorical variables not included in the statistical analysis (contains only 1 level): ' + ", ".join(remCatFields) + '\n'
        result += 'Quantitative variables selected by user: ' + ", ".join(self.quantFields) + '\n'
        result += '===============================================\n\n'

        functions.setBase(self.RID, 'Step 1 of 4: Selecting your chosen meta-variables...done')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        return 0

    def query(self):
        print "Query!"
        functions.setBase(self.RID, 'Step 2 of 4: Selecting your chosen taxa or KEGG level...')
        # filter otus based on user settings
        remUnclass = self.all['remUnclass']
        remZeroes = self.all['remZeroes']
        perZeroes = int(self.all['perZeroes'])
        filterData = self.all['filterData']
        filterPer = int(self.all['filterPer'])
        filterMeth = int(self.all['filterMeth'])
        self.mapTaxa = self.all['map_taxa']
        self.finalDF = pd.DataFrame()
        self.allDF = pd.DataFrame()
        if self.treeType == 1:
            if self.selectAll == 0 or self.selectAll == 8:
                taxaString = self.all["taxa"]
                taxaDict = json.JSONDecoder(object_pairs_hook=functions.multidict).decode(taxaString)
                filteredDF = self.savedDF.copy()
            else:
                taxaDict = ''
                filteredDF = functions.filterDF(self.savedDF, self.DepVar, self.selectAll, remUnclass, remZeroes, perZeroes, filterData, filterPer, filterMeth)
            self.finalDF, missingList = functions.getTaxaDF(self.selectAll, taxaDict, filteredDF, self.metaDF, self.allFields, self.DepVar, self.RID, self.stopList, self.PID)

            if self.selectAll == 8:
                self.result += '\nThe following PGPRs were not detected: ' + ", ".join(missingList) + '\n'
                self.result += '===============================================\n'
        if self.treeType == 2:
            keggDict = ''
            if self.keggAll == 0:
                keggString = self.all["kegg"]
                keggDict = json.JSONDecoder(object_pairs_hook=functions.multidict).decode(keggString)
            self.finalDF, self.allDF = functions.getKeggDF(self.keggAll, keggDict, self.savedDF, self.metaDF, self.DepVar, self.mapTaxa, self.RID, self.stopList, self.PID)
        if self.treeType == 3:
            keggDict = ''
            if self.nzAll == 0:
                keggString = self.all["nz"]
                keggDict = json.JSONDecoder(object_pairs_hook=functions.multidict).decode(keggString)
            self.finalDF, self.allDF = functions.getNZDF(self.nzAll, keggDict, self.savedDF, self.metaDF, self.DepVar, self.mapTaxa, self.RID, self.stopList, self.PID)
        if self.finalDF.empty:
            error = "Selected taxa were not found in your selected samples."
            myDict = {'error': error}
            res = json.dumps(myDict)
            return HttpResponse(res, content_type='application/json')

        # make sure column types are correct
        self.finalDF[self.catFields] = self.finalDF[self.catFields].astype(str)
        self.finalDF[self.quantFields] = self.finalDF[self.quantFields].astype(float)

        # transform Y, if requested
        transform = int(self.all["transform"])
        self.finalDF = functions.transformDF(transform, self.DepVar, self.finalDF)

        # save location info to session
        myDir = 'myPhyloDB/media/temp/anova/'
        if not os.path.exists(myDir):
            os.makedirs(myDir)

        path = str(myDir) + str(self.RID) + '.biom'
        functions.imploding_panda(path, self.treeType, self.DepVar, self.finalSampleIDs, self.metaDF, self.finalDF)
        functions.setBase(self.RID, 'Step 2 of 4: Selecting your chosen taxa or KEGG level...done')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        return 0

    def stats(self):
        print "Stats!"

        functions.setBase(self.RID, 'Step 3 of 4: Performing statistical test...')
        self.finalDict = {}
        self.seriesList = []
        self.xAxisDict = {}
        self.yAxisDict = {}

        if os.name == 'nt':
            r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
        else:
            r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

        functions.setBase(self.RID, 'Verifying R packages...missing packages are being installed')

        # R packages from cran
        r("list.of.packages <- c('lsmeans', 'ggplot2', 'RColorBrewer', 'ggthemes')")
        r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
        print r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

        functions.setBase(self.RID, 'Step 3 of 4: Performing statistical test...')

        print r("library(lsmeans)")
        print r("library(ggplot2)")
        print r("library(ggthemes)")
        print r("library(RColorBrewer)")
        print r('source("R/myFunctions/myFunctions.R")')
        # R graph
        r.assign('finalDF', self.finalDF)

        colorVal = self.all['colorVal']
        if colorVal != 'None':
            r.assign('colorVal', colorVal)
        else:
            r.assign('colorVal', 'rank_name')

        xVal = self.all['xVal']
        if xVal != 'None':
            r.assign('xVal', xVal)
        else:
            r.assign('xVal', self.catFields[0])
        print "Checkpoint"
        gridVal_X = self.all['gridVal_X']
        r.assign('gridVal_X', gridVal_X)

        gridVal_Y = self.all['gridVal_Y']
        r.assign('gridVal_Y', gridVal_Y)

        if self.DepVar == 0:
            r('DepVar <- "abund"')
        elif self.DepVar == 1:
            r('DepVar <- "rel_abund"')
        elif self.DepVar == 2:
            r('DepVar <- "rich"')
        elif self.DepVar == 3:
            r('DepVar <- "diversity"')
        elif self.DepVar == 4:
            r('DepVar <- "abund_16S"')
        if gridVal_X == 'None' and gridVal_Y == 'None':
            r("gDF <- data.frame(x=finalDF[,paste(xVal)], y=finalDF[,paste(DepVar)], \
                myFill=finalDF[,paste(colorVal)])")
            r("gDF <- aggregate(gDF[, 'y'], list(gDF$x, gDF$myFill), mean)")
            r("names(gDF) <- c('x', 'myFill', 'y')")
        elif gridVal_X != 'None' and gridVal_Y == 'None':
            r("gDF <- data.frame(x=finalDF[,paste(xVal)], y=finalDF[,paste(DepVar)], \
                gridVal_X=finalDF[,paste(gridVal_X)], \
                myFill=finalDF[,paste(colorVal)])")
            r("gDF <- aggregate(gDF[, 'y'], list(gDF$x, gDF$myFill, gDF$gridVal_X), mean)")
            r("names(gDF) <- c('x', 'myFill', 'gridVal_X', 'y')")
        elif gridVal_X == 'None' and gridVal_Y != 'None':
            r("gDF <- data.frame(x=finalDF[,paste(xVal)], y=finalDF[,paste(DepVar)], \
                gridVal_Y=finalDF[,paste(gridVal_Y)], \
                myFill=finalDF[,paste(colorVal)])")
            r("gDF <- aggregate(gDF[, 'y'], list(gDF$x, gDF$myFill, gDF$gridVal_Y), mean)")
            r("names(gDF) <- c('x', 'myFill', 'gridVal_Y', 'y')")
        elif gridVal_X != 'None' and gridVal_Y != 'None':
            r("gDF <- data.frame(x=finalDF[,paste(xVal)], y=finalDF[,paste(DepVar)], \
                gridVal_X=finalDF[,paste(gridVal_X)], gridVal_Y=finalDF[,paste(gridVal_Y)], \
                myFill=finalDF[,paste(colorVal)])")
            r("gDF <- aggregate(gDF[, 'y'], list(gDF$x, gDF$myFill, gDF$gridVal_X, gDF$gridVal_Y), mean)")
            r("names(gDF) <- c('x', 'myFill', 'gridVal_X', 'gridVal_Y', 'y')")

        r("p <- ggplot(gDF, aes(x=x, y=y, fill=myFill))")
        r("p <- p + geom_bar(stat='identity', position='stack')")

        if gridVal_X != 'None' and gridVal_Y == 'None':
            r("p <- p + facet_grid(. ~ gridVal_X)")
            r("p <- p + theme(strip.text.x=element_text(size=10, colour='blue', angle=0))")
        elif gridVal_X == 'None' and gridVal_Y != 'None':
            r("p <- p + facet_grid(gridVal_Y ~ .)")
            r("p <- p + theme(strip.text.y=element_text(size=10, colour='blue', angle=90))")
        elif gridVal_X != 'None' and gridVal_Y != 'None':
            r("p <- p + facet_grid(gridVal_Y ~ gridVal_X)")
            r("p <- p + theme(strip.text.x=element_text(size=10, colour='blue', angle=0))")
            r("p <- p + theme(strip.text.y=element_text(size=10, colour='blue', angle=90))")

        r("p <- p + theme(strip.text.x=element_text(size=10, colour='blue', angle=0))")
        r("p <- p + theme(strip.text.y=element_text(size=10, colour='blue', angle=90))")
        palette = self.all['palette']
        r.assign('palette', palette)
        if palette == 'gdocs':
            r('pal <- gdocs_pal()(20)')
        elif palette == 'hc':
            r('pal <- hc_pal()(10)')
        elif palette == 'Set1':
            r('pal <- brewer.pal(8, "Set1")')
        elif palette == 'Set2':
            r('pal <- brewer.pal(8, "Set2")')
        elif palette == 'Set3':
            r('pal <- brewer.pal(12, "Set3")')
        elif palette == 'Paired':
            r('pal <- brewer.pal(12, "Paired")')
        elif palette == 'Dark2':
            r('pal <- brewer.pal(12, "Dark2")')
        elif palette == 'Accent':
            r('pal <- brewer.pal(12, "Accent")')
        r('nColors <- length(pal)')

        r("p <- p + scale_fill_manual(values=rep(pal, ceiling(nlevels(gDF$myFill)/nColors)))")

        r("p <- p + theme(legend.text=element_text(size=7))")
        r("p <- p + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))")

        r("p <- p + theme(legend.title=element_blank())")
        r("p <- p + theme(legend.position='bottom')")
        r("p <- p + guides(fill=guide_legend(ncol=8))")

        if self.DepVar == 0:
            r("p <- p + ylab('Abundance') + xlab('')")
        elif self.DepVar == 1:
            r("p <- p + ylab('Relative Abundance') + xlab('')")
        elif self.DepVar == 2:
            r("p <- p + ylab('OTU Richness') + xlab('')")
        elif self.DepVar == 3:
            r("p <- p + ylab('OTU Diversity') + xlab('')")
        elif self.DepVar == 4:
            r("p <- p + ylab('Total Abundance') + xlab('')")

        path = "myPhyloDB/media/temp/anova/Rplots"
        if not os.path.exists(path):
            os.makedirs(path)

        r.assign("path", path)
        r.assign("RID", self.RID)
        r("file <- paste(path, '/', RID, '.anova.pdf', sep='')")
        r("p <- set_panel_size(p, height=unit(2.9, 'in'), width=unit(2.9, 'in'))")

        r("nlev <- nlevels(as.factor(gDF$gridVal_X))")

        r('if (nlev == 0) { \
                myWidth <- 11 \
            } else { \
                myWidth <- 3*nlev+4 \
        }')

        r("nlev <- nlevels(as.factor(gDF$gridVal_Y))")
        r("nRow <- ceiling(nlevels(gDF$myFill)/8)")
        r('if (nlev == 0) { \
                myHeight <- 8 \
            } else { \
                myHeight <- 2+(3*nlev)+(1*nRow) \
        }')

        r("ggsave(filename=file, plot=p, units='in', height=myHeight, width=myWidth, limitsize=F)")

        print "Checkpoint"

        # group DataFrame by each taxa level selected
        grouped1 = self.finalDF.groupby(['rank_name', 'rank_id'])
        pValDict = {}
        counter = 1
        for name1, group1 in grouped1:
            D = ''
            r.assign("df", group1)
            trtString = " * ".join(self.allFields)

            if self.DepVar == 0:
                anova_string = "fit <- aov(df$abund ~ " + str(trtString) + ", data=df)"
                r.assign("cmd", anova_string)
                r("eval(parse(text=cmd))")

            elif self.DepVar == 1:
                anova_string = "fit <- aov(df$rel_abund ~ " + str(trtString) + ", data=df)"
                r.assign("cmd", anova_string)
                r("eval(parse(text=cmd))")

            elif self.DepVar == 2:
                anova_string = "fit <- aov(df$rich ~ " + str(trtString) + ", data=df)"
                r.assign("cmd", anova_string)
                r("eval(parse(text=cmd))")

            elif self.DepVar == 3:
                anova_string = "fit <- aov(df$diversity ~ " + str(trtString) + ", data=df)"
                r.assign("cmd", anova_string)
                r("eval(parse(text=cmd))")

            elif self.DepVar == 4:
                anova_string = "fit <- aov(df$abund_16S ~ " + str(trtString) + ", data=df)"
                r.assign("cmd", anova_string)
                r("eval(parse(text=cmd))")

            aov = r("summary(fit)")
            pString = r("summary(fit)[[1]][['Pr(>F)']]")

            tempStuff = pString.split(' ')
            pList = []
            for part in tempStuff:
                try:
                    pList.append(float(part))
                except Exception:
                    pass

            if pList:
                p_val = min(pList)
                tempStuff = aov.split('\n')
                for part in tempStuff:
                    if part != tempStuff[0]:
                        D += part + '\n'

                fList = []
                for part in tempStuff:
                    if part != tempStuff[0] and part != tempStuff[1]:
                        part = part.replace('\r', '')
                        part1 = part.split(' ')
                        if part1[0] == 'Residuals':
                            break
                        if part1[0] not in fList:
                            fList.append(part1[0])

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if self.stopList[self.PID] == self.RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                D += "\nLSmeans & Tukey's HSD post-hoc test:\n\n"

                if len(self.quantFields) == 0:
                    for i in fList:
                        hsd_string = "lsm <- lsmeans(fit, list(pairwise ~ " + str(i) + "))"
                        r.assign("cmd", hsd_string)
                        r("eval(parse(text=cmd))")
                        r("options(width=5000)")
                        table = r("lsm")
                        tempStuff = table.split('\n')
                        for i in xrange(len(tempStuff)):
                            if i > 0:
                                D += tempStuff[i] + '\n'
                else:
                    for i in fList:
                        if i not in self.quantFields:
                            hsd_string = "lsm <- lsmeans(fit, list(pairwise ~ " + str(i) + "))"
                            r.assign("cmd", hsd_string)
                            r("eval(parse(text=cmd))")
                            r("options(width=5000)")
                            table = r("lsm")
                            tempStuff = table.split('\n')
                            for i in xrange(len(tempStuff)):
                                if i > 0:
                                    D += tempStuff[i] + '\n'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if self.stopList[self.PID] == self.RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            else:
                p_val = 1.0
                D = 'ANOVA cannot be performed, please check that you have more than one treatment level and appropriate replication.\n'

            pValDict[name1] = p_val

            self.result += 'Name: ' + str(name1[0]) + '\n'
            self.result += 'ID: ' + str(name1[1]) + '\n'
            if self.DepVar == 0:
                self.result += 'Dependent Variable: Abundance' + '\n'
            elif self.DepVar == 1:
                self.result += 'Dependent Variable: Relative Abundance' + '\n'
            elif self.DepVar == 2:
                self.result += 'Dependent Variable: OTU Richness' + '\n'
            elif self.DepVar == 3:
                self.result += 'Dependent Variable: OTU Diversity' + '\n'
            elif self.DepVar == 4:
                self.result += 'Dependent Variable: Total Abundance' + '\n'

            self.result += '\nANCOVA table:\n'
            D = D.decode('utf-8')
            self.result += D + '\n'
            self.result += '===============================================\n'
            self.result += '\n\n\n\n'

            taxa_no = len(grouped1)
            functions.setBase(self.RID, 'Step 3 of 4: Performing statistical test...taxa ' + str(counter) + ' of ' + str(taxa_no) + ' is complete!')
            counter += 1

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if self.stopList[self.PID] == self.RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        functions.setBase(self.RID, 'Step 3 of 4: Performing statistical test...done!')
        return 0

    def graph(self):
        print "Graph!"
        functions.setBase(self.RID, 'Step 4 of 4: Formatting graph data for display...')

        grouped1 = self.finalDF.groupby(['rank_name', 'rank_id'])
        for name1, group1 in grouped1:
            dataList = []
            errorList = []
            pValue = self.pValDict[name1]

            if self.sig_only == 0:
                if self.DepVar == 0:
                    mean = group1.groupby(self.catFields)['abund'].mean()
                    se = group1.groupby(self.catFields)['abund'].std()
                    se.fillna(0, inplace=True)
                    high = [x + y for x, y in zip(mean, se)]
                    low = [x - y for x, y in zip(mean, se)]
                    dataList = list(mean)
                    errorTuple = zip(low, high)
                    errorList = [list(elem) for elem in errorTuple]
                elif self.DepVar == 1:
                    mean = group1.groupby(self.catFields)['rel_abund'].mean()
                    se = group1.groupby(self.catFields)['rel_abund'].std()
                    se.fillna(0, inplace=True)
                    high = [x + y for x, y in zip(mean, se)]
                    low = [x - y for x, y in zip(mean, se)]
                    dataList = list(mean)
                    errorTuple = zip(low, high)
                    errorList = [list(elem) for elem in errorTuple]
                elif self.DepVar == 2:
                    mean = group1.groupby(self.catFields)['rich'].mean()
                    se = group1.groupby(self.catFields)['rich'].std()
                    se.fillna(0, inplace=True)
                    high = [x + y for x, y in zip(mean, se)]
                    low = [x - y for x, y in zip(mean, se)]
                    dataList = list(mean)
                    errorTuple = zip(low, high)
                    errorList = [list(elem) for elem in errorTuple]
                elif self.DepVar == 3:
                    mean = group1.groupby(self.catFields)['diversity'].mean()
                    se = group1.groupby(self.catFields)['diversity'].std()
                    se.fillna(0, inplace=True)
                    high = [x + y for x, y in zip(mean, se)]
                    low = [x - y for x, y in zip(mean, se)]
                    dataList = list(mean)
                    errorTuple = zip(low, high)
                    errorList = [list(elem) for elem in errorTuple]
                elif self.DepVar == 4:
                    mean = group1.groupby(self.catFields)['abund_16S'].mean()
                    se = group1.groupby(self.catFields)['abund_16S'].std()
                    se.fillna(0, inplace=True)
                    high = [x + y for x, y in zip(mean, se)]
                    low = [x - y for x, y in zip(mean, se)]
                    dataList = list(mean)
                    errorTuple = zip(low, high)
                    errorList = [list(elem) for elem in errorTuple]

                seriesDict = {}
                seriesDict['name'] = name1
                seriesDict['type'] = 'column'
                seriesDict['data'] = dataList
                self.seriesList.append(seriesDict)

                seriesDict = {}
                seriesDict['name'] = name1
                seriesDict['type'] = 'errorbar'
                seriesDict['visible'] = False
                seriesDict['data'] = errorList
                self.seriesList.append(seriesDict)

            elif self.sig_only == 1:
                if pValue < 0.05:
                    if self.DepVar == 0:
                        mean = group1.groupby(self.catFields)['abund'].mean()
                        se = group1.groupby(self.catFields)['abund'].std()
                        se.fillna(0, inplace=True)
                        high = [x + y for x, y in zip(mean, se)]
                        low = [x - y for x, y in zip(mean, se)]
                        dataList = list(mean)
                        errorTuple = zip(low, high)
                        errorList = [list(elem) for elem in errorTuple]
                    elif self.DepVar == 1:
                        mean = group1.groupby(self.catFields)['rel_abund'].mean()
                        se = group1.groupby(self.catFields)['rel_abund'].std()
                        se.fillna(0, inplace=True)
                        high = [x + y for x, y in zip(mean, se)]
                        low = [x - y for x, y in zip(mean, se)]
                        dataList = list(mean)
                        errorTuple = zip(low, high)
                        errorList = [list(elem) for elem in errorTuple]
                    elif self.DepVar == 2:
                        mean = group1.groupby(self.catFields)['rich'].mean()
                        se = group1.groupby(self.catFields)['rich'].std()
                        se.fillna(0, inplace=True)
                        high = [x + y for x, y in zip(mean, se)]
                        low = [x - y for x, y in zip(mean, se)]
                        dataList = list(mean)
                        errorTuple = zip(low, high)
                        errorList = [list(elem) for elem in errorTuple]
                    elif self.DepVar == 3:
                        mean = group1.groupby(self.catFields)['diversity'].mean()
                        se = group1.groupby(self.catFields)['diversity'].std()
                        se.fillna(0, inplace=True)
                        high = [x + y for x, y in zip(mean, se)]
                        low = [x - y for x, y in zip(mean, se)]
                        dataList = list(mean)
                        errorTuple = zip(low, high)
                        errorList = [list(elem) for elem in errorTuple]
                    elif self.DepVar == 4:
                        mean = group1.groupby(self.catFields)['abund_16S'].mean()
                        se = group1.groupby(self.catFields)['abund_16S'].std()
                        se.fillna(0, inplace=True)
                        high = [x + y for x, y in zip(mean, se)]
                        low = [x - y for x, y in zip(mean, se)]
                        dataList = list(mean)
                        errorTuple = zip(low, high)
                        errorList = [list(elem) for elem in errorTuple]

                    seriesDict = {}
                    seriesDict['name'] = name1
                    seriesDict['type'] = 'column'
                    seriesDict['data'] = dataList
                    self.seriesList.append(seriesDict)

                    seriesDict = {}
                    seriesDict['name'] = name1
                    seriesDict['type'] = 'errorbar'
                    seriesDict['visible'] = False
                    seriesDict['data'] = errorList
                    self.seriesList.append(seriesDict)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if self.stopList[self.PID] == self.RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            '''
            catFieldsList = []
            for i in catFields:
                catFieldsList.append(len(group1.groupby(i)))

            catFields = [x for (y, x) in sorted(zip(catFieldsList, catFields))]
            '''

            if self.DepVar == 0:
                grouped2 = group1.groupby(self.catFields)['abund'].mean()
            elif self.DepVar == 1:
                grouped2 = group1.groupby(self.catFields)['rel_abund'].mean()
            elif self.DepVar == 4:
                grouped2 = group1.groupby(self.catFields)['abund_16S'].mean()
            elif self.DepVar == 2:
                grouped2 = group1.groupby(self.catFields)['rich'].mean()
            elif self.DepVar == 3:
                grouped2 = group1.groupby(self.catFields)['diversity'].mean()
            else:
                raise Exception("Something went horribly wrong")

            if self.catFields.__len__() == 1:
                self.xAxisDict['categories'] = grouped2.index.values.tolist()
            else:
                g2indexvals = grouped2.index.values
                level = g2indexvals[0].__len__()
                labelTree = recLabels(g2indexvals, level)
                self.xAxisDict['categories'] = labelTree['categories']

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if self.stopList[self.PID] == self.RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        yTitle = {}
        if self.DepVar == 0:
            yTitle['text'] = 'Abundance'
        elif self.DepVar == 1:
            yTitle['text'] = 'Relative Abundance'
        elif self.DepVar == 2:
            yTitle['text'] = 'OTU Richness'
        elif self.DepVar == 3:
            yTitle['text'] = 'OTU Diversity'
        elif self.DepVar == 4:
            yTitle['text'] = 'Total Abundance'
        yTitle['style'] = {'fontSize': '18px', 'fontWeight': 'bold'}

        if self.transform != 0:
            tname = {
                '1': "Ln", '2': "Log10", '3': "Sqrt", '4': "Logit", '5': "Arcsin"
            }
            yTitle['text'] = tname[str(self.transform)] + "(" + str(yTitle['text']) + ")"

        self.yAxisDict['title'] = yTitle

        xStyleDict = {'style': {'fontSize': '14px'}, 'rotation': 0}
        self.xAxisDict['labels'] = xStyleDict
        yStyleDict = {'style': {'fontSize': '14px'}}
        self.yAxisDict['labels'] = yStyleDict

        self.finalDict['series'] = self.seriesList
        self.finalDict['xAxis'] = self.xAxisDict
        self.finalDict['yAxis'] = self.yAxisDict
        self.finalDict['text'] = self.result

        if not self.seriesList:
            self.finalDict['empty'] = 0
        else:
            self.finalDict['empty'] = 1

        functions.setBase(self.RID, 'Step 4 of 4: Formatting graph data for display...done!')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        # datatable of taxa mapped to selected kegg orthologies
        if not self.treeType == 1 and self.mapTaxa == 'yes':
            myDir = 'myPhyloDB/media/temp/anova/'
            fileName = str(myDir) + 'Mapped_Taxa.csv'
            self.allDF.to_csv(fileName)

            myDir = 'myPhyloDB/media/temp/anova/'
            fileName2 = str(myDir) + 'Mapped_Taxa.gz'
            zf = zipfile.ZipFile(fileName2, "w", zipfile.ZIP_DEFLATED, allowZip64=True)
            zf.write(fileName, 'Mapped_Taxa.csv')
            zf.close()

        self.finalDict['resType'] = 'res'
        self.finalDict['error'] = 'none'

        res = json.dumps(self.finalDict)
        return HttpResponse(res, content_type='application/json')


class Anova(Analysis):

    def run(self):
        print "I'm an anova. WIP"
        ret = self.validate()
        if ret == 0:
            ret = self.query()
            if ret == 0:
                ret = self.stats()
                if ret == 0:
                    return self.graph()
        print "Something went wrong with Anova"
        return ret
