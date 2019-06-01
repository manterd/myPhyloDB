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
import collections
import numpy as np
from scipy import stats
from database.models import Sample
from natsort import natsorted
from PyPDF2 import PdfFileReader, PdfFileMerger


# stop function for creating proper return message
def getStopDict():
    stopFinal = {}
    stopFinal['resType'] = 'res'
    stopFinal['error'] = 'Your request was stopped'
    return json.dumps(stopFinal)


# analysis class is a base for other analyses to overwrite (not intended for direct use)
# TODO create R based alternative functions, which accept data querysets and R scripts to produce PDF format graphs
class Analysis:  # abstract parent class, not to be run on its own. Instead, should be used as a template
    __metaclass__ = ABCMeta

    # Initializes the slightly obnoxious amount of instance variables
    # not all values are to be used at one time
    # if a value must cross multiple functions, declare and initialize it here (yes, 'self.' is required)
    def __init__(self, iRequest, iRID, iStops, iPID, debug=False):
        self.request = iRequest
        self.RID = iRID
        self.stopList = iStops
        self.PID = iPID
        self.debug = debug
        self.result = ""
        # 'declare' all variables used in multiple steps here
        self.all, self.selectAll, self.treeType, self.savedDF, self.DepVar, self.metaDF, self.allFields = (None,)*7
        self.keggAll, self.nzAll, self.catFields, self.quantFields, self.finalSampleIDs = (None,)*5
        self.finalDict, self.seriesList, self.xAxisDict, self.yAxisDict, self.finalDF, self.pValDict = (None,)*6
        self.sig_only, self.transform, self.finalDict, self.zipFile, self.mapTaxa, self.allDF = (None,)*6
        self.distance, self.alpha, = (None,)*2

    # construct 'run' in each analysis separately, as settings and sequencing vary
    @abstractmethod
    def run(self):
        pass

    # verify request was formatted correctly and all variables are accounted for, as well as configure future steps
    def validate(self, sig=True, metaCat=True, metaQuant=True, reqMultiLevel=True, dist=False, selAll=True, taxTree=True):  # supports flags for sig_only, metaValsCat
        if self.debug:
            print "Validate!"
        # Get variables from web page
        allJson = self.request.body.split('&')[0]
        self.all = json.loads(allJson)

        if self.debug:
            print "Analysis data: ", self.all

        functions.setBase(self.RID, 'Step 1 of 4: Selecting your chosen meta-variables...')

        if selAll:
            self.selectAll = int(self.all["selectAll"])
            self.keggAll = int(self.all["keggAll"])
            self.nzAll = int(self.all["nzAll"])

        if sig:  # check for sig_only support, use value if yes, treat as false if not
            self.sig_only = int(self.all["sig_only"])
        else:
            self.sig_only = 0
        if dist:
            self.distance = int(self.all["distance"])
            self.alpha = float(self.all["alpha"])
        else:
            self.distance = 0

        if metaCat:
            metaValsCat = self.all['metaValsCat']
            metaIDsCat = self.all['metaIDsCat']
        else:
            metaValsCat = []
            metaIDsCat = []

        if metaQuant:
            metaValsQuant = self.all['metaValsQuant']
            metaIDsQuant = self.all['metaIDsQuant']
        else:
            metaValsQuant = []
            metaIDsQuant = []

        if taxTree:
            self.treeType = int(self.all['treeType'])

        self.DepVar = int(self.all["DepVar"])

        if self.debug:
            print "First!"

        # Create meta-variable DataFrame, final sample list, final category and quantitative field lists based on tree selections
        self.savedDF, self.metaDF, self.finalSampleIDs, self.catFields, remCatFields, self.quantFields, self.catValues, self.quantValues = functions.getMetaDF(self.request.user, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, self.DepVar, levelDep=True)
        self.allFields = self.catFields + self.quantFields
        if reqMultiLevel:
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

        if self.debug:
            print "Second!"

        self.result = ''
        if self.treeType == 1:
            if self.selectAll == 1:
                self.result += 'Taxa level: Kingdom' + '\n'
            elif self.selectAll == 2:
                self.result += 'Taxa level: Phyla' + '\n'
            elif self.selectAll == 3:
                self.result += 'Taxa level: Class' + '\n'
            elif self.selectAll == 4:
                self.result += 'Taxa level: Order' + '\n'
            elif self.selectAll == 5:
                self.result += 'Taxa level: Family' + '\n'
            elif self.selectAll == 6:
                self.result += 'Taxa level: Genus' + '\n'
            elif self.selectAll == 7:
                self.result += 'Taxa level: Species' + '\n'
            elif self.selectAll == 9:
                self.result += 'Taxa level: OTU_99' + '\n'
        elif self.treeType == 2:
            if self.keggAll == 1:
                self.result += 'KEGG Pathway level: 1' + '\n'
            elif self.keggAll == 2:
                self.result += 'KEGG Pathway level: 2' + '\n'
            elif self.keggAll == 3:
                self.result += 'KEGG Pathway level: 3' + '\n'
        elif self.treeType == 3:
            if self.nzAll == 1:
                self.result += 'KEGG Enzyme level: 1' + '\n'
            elif self.nzAll == 2:
                self.result += 'KEGG Enzyme level: 2' + '\n'
            elif self.nzAll == 3:
                self.result += 'KEGG Enzyme level: 3' + '\n'
            elif self.nzAll == 4:
                self.result += 'KEGG Enzyme level: 4' + '\n'
            elif self.keggAll == 5:
                self.result += 'KEGG Enzyme level: GIBBs' + '\n'
            elif self.keggAll == 6:
                self.result += 'KEGG Enzyme level: Nitrogen cycle' + '\n'

        if self.distance == 1:
            self.result += 'Distance score: Manhattan' + '\n'
        elif self.distance == 2:
            self.result += 'Distance score: Euclidean' + '\n'
        elif self.distance == 3:
            self.result += 'Distance score: Canberra' + '\n'
        elif self.distance == 4:
            self.result += 'Distance score: Bray-Curtis' + '\n'
        elif self.distance == 5:
            self.result += 'Distance score: Kulczynski' + '\n'
        elif self.distance == 6:
            self.result += 'Distance score: Jaccard' + '\n'
        elif self.distance == 7:
            self.result += 'Distance score: Gower' + '\n'
        elif self.distance == 8:
            self.result += 'Distance score: altGower' + '\n'
        elif self.distance == 9:
            self.result += 'Distance score: Morisita' + '\n'
        elif self.distance == 10:
            self.result += 'Distance score: Horn' + '\n'
        elif self.distance == 11:
            self.result += 'Distance score: Mountford' + '\n'
        elif self.distance == 12:
            self.result += 'Distance score: Binomial' + '\n'
        elif self.distance == 13:
            self.result += 'Distance score: Chao' + '\n'
        elif self.distance == 14:
            self.result += 'Distance score: Cao' + '\n'
        elif self.distance == 15:
            self.result += 'Distance score: wOdum' + '\n'
            self.result += 'alpha: ' + str(self.alpha) + '\n'

        if self.debug:
            print "Third!"

        self.result += 'Categorical variables selected by user: ' + ", ".join(self.catFields + remCatFields) + '\n'
        self.result += 'Categorical variables not included in the statistical analysis (contains only 1 level): ' + ", ".join(remCatFields) + '\n'
        if metaQuant:
            self.result += 'Quantitative variables selected by user: ' + ", ".join(self.quantFields) + '\n'
        self.result += '===============================================\n\n'

        functions.setBase(self.RID, 'Step 1 of 4: Selecting your chosen meta-variables...done')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        if self.debug:
            print "Fourth!"

        return 0

    # get actual data from database and format it for actual analysis step based on settings
    def query(self, taxmap=True, usetransform=True, filterable=True):
        if self.debug:
            print "Query!"
        functions.setBase(self.RID, 'Step 2 of 4: Selecting your chosen taxa or KEGG level...')
        if filterable:
            # filter otus based on user settings
            remUnclass = self.all['remUnclass']
            remZeroes = self.all['remZeroes']
            perZeroes = int(self.all['perZeroes'])
            filterData = self.all['filterData']
            filterPer = int(self.all['filterPer'])
            filterMeth = int(self.all['filterMeth'])

        if self.debug:
            print "First!"

        if taxmap:
            self.mapTaxa = self.all['map_taxa']
        else:
            self.mapTaxa = "no"

        self.finalDF = pd.DataFrame()
        self.allDF = pd.DataFrame()
        if self.treeType == 1:
            if self.selectAll == 0 or self.selectAll == 8:
                try:
                    taxaString = self.all["taxa"]
                    taxaDict = json.JSONDecoder(object_pairs_hook=functions.multidict).decode(taxaString)
                except:
                    taxaDict = {}
                filteredDF = self.savedDF.copy()
            else:   # Filter lists being sent, need default values for unfiltered for gage? TODO More readable code
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

        if self.debug:
            print "Second!"

        # make sure column types are correct
        self.finalDF[self.catFields] = self.finalDF[self.catFields].astype(str)
        self.finalDF[self.quantFields] = self.finalDF[self.quantFields].astype(float)

        if usetransform:
            # transform Y, if requested
            self.transform = int(self.all["transform"])
            self.finalDF = functions.transformDF(self.transform, self.DepVar, self.finalDF)

        # save location info to session
        myDir = 'myPhyloDB/media/temp/analysis/'   # TODO this gets called by more than ANOVA, path should match caller
        if not os.path.exists(myDir):
            os.makedirs(myDir)

        if self.debug:
            print "Third!"

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        path = str(myDir) + str(self.RID) + '.biom'
        functions.imploding_panda(path, self.treeType, self.DepVar, self.finalSampleIDs, self.metaDF, self.finalDF)
        functions.setBase(self.RID, 'Step 2 of 4: Selecting your chosen taxa or KEGG level...done')

        if self.debug:
            print "Fourth!"

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        return 0

    # do the actual statistics side of the analysis
    def stats(self):
        if self.debug:
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

        if self.debug:
            print "First!"

        # R packages from cran
        r("list.of.packages <- c('lsmeans', 'ggplot2', 'RColorBrewer', 'ggthemes')")
        r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
        r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

        functions.setBase(self.RID, 'Step 3 of 4: Performing statistical test...')

        r("library(lsmeans)")
        r("library(ggplot2)")
        r("library(ggthemes)")
        r("library(RColorBrewer)")
        r('source("R/myFunctions/myFunctions.R")')
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

        # not sure this is worth splitting


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

        # group DataFrame by each taxa level selected
        grouped1 = self.finalDF.groupby(['rank_name', 'rank_id'])
        self.pValDict = {}
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
                    if self.debug:
                        print "Stopping!"
                    return HttpResponse(getStopDict(), content_type='application/json')
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
                    if self.debug:
                        print "Stopping!"
                    return HttpResponse(getStopDict(), content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            else:
                p_val = 1.0
                D = 'ANOVA cannot be performed, please check that you have more than one treatment level and appropriate replication.\n'

            self.pValDict[name1] = p_val
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
                if self.debug:
                    print "Stopping!"
                return HttpResponse(getStopDict(), content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        functions.setBase(self.RID, 'Step 3 of 4: Performing statistical test...done!')
        return 0

    # format results from stats into desired output, typically the last step before returning request
    def graph(self):
        if self.debug:
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
                if self.debug:
                    print "Stopping!"
                return HttpResponse(getStopDict(), content_type='application/json')
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
                if self.debug:
                    print "Stopping!"
                return HttpResponse(getStopDict(), content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        if self.debug:
            print "First!"

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
            science = tname[str(self.transform)]
            yTitle['text'] = science + "(" + str(yTitle['text']) + ")"
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

        if self.debug:
            print "Second!"

        functions.setBase(self.RID, 'Step 4 of 4: Formatting graph data for display...done!')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        # datatable of taxa mapped to selected kegg orthologies
        if not self.treeType == 1 and self.mapTaxa == 'yes':
            myDir = 'myPhyloDB/media/temp/anova/'  # anova specific in graph function because its used by anova only
            fileName = str(myDir) + 'Mapped_Taxa.csv'
            self.allDF.to_csv(fileName)

            myDir = 'myPhyloDB/media/temp/anova/'
            fileName2 = str(myDir) + 'Mapped_Taxa.gz'
            zf = zipfile.ZipFile(fileName2, "w", zipfile.ZIP_DEFLATED, allowZip64=True)
            zf.write(fileName, 'Mapped_Taxa.csv')
            zf.close()

        self.finalDict['resType'] = 'res'
        self.finalDict['error'] = 'none'

        if self.debug:
            print "Third!"

        res = json.dumps(self.finalDict)
        return HttpResponse(res, content_type='application/json')

    # smaller shared functions
    # prepare R for use
    def initializeR(self):
        if os.name == 'nt':
            r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
        else:
            r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

        # R packages from cran
        r("list.of.packages <- c('lsmeans', 'ggplot2', 'RColorBrewer', 'ggthemes', 'picante', 'FD')")
        r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
        r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

        r("library(lsmeans)")
        r("library(ggplot2)")
        r("library(ggthemes)")
        r("library(RColorBrewer)")
        r('source("R/myFunctions/myFunctions.R")')

        return r

    def fullRAnalysis(self, script):    # script is string for filename in R Scripts folder, R should do the most work
        if self.validate() == 0:
            if self.debug:
                print "Here we go"
            curUser = self.request.user
            # query data vs read from phyloseq.biom
            # need metadata to be passed in regardless, query would involve formatting for pandas only to send to R
            # should send metadata to R but drop Query for phyloseq reading
            pathToBiom = "myPhyloDB/media/usr_temp/"+str(curUser)
            r = self.initializeR()  # need to see an example of a full script to plug it in correctly
            if self.debug:
                print "R initialized"
            r.assign("pathToBiom", pathToBiom)
            # actually run R given script and data, not sure which is needed where
            r.assign("script", script)
            print r("getwd()")
            if self.debug:
                print "Assigned script"
            output = r("source(script)")
            if self.debug:
                print "Script run"
            print "Output: ", output
            print r("warnings()")
            # get output from R and send it back up
            return output
        else:
            print "Validate failed"

    def deprint(self, text):
        if self.debug:
            print text



# TODO trim these down and/or split them up further (more helper functions!)
class Anova(Analysis):  # base template, is essentially getCatUnivData

    # so quant Univ is huge... like >700 lines long.... implementing is one thing, but this could easily be split 6 ways
    def quantstats(self):       # for getQuantUnivData      WIP

        if self.debug:
            print "statsgraph! (AnovaQuant)"


        functions.setBase(self.RID, 'Step 3 of 4: Performing statistical test...!')

        finalDict = {}
        # group DataFrame by each taxa level selected
        shapes = ['circle', 'square', 'triangle', 'triangle-down', 'diamond']

        functions.setBase(self.RID, 'Verifying R packages...missing packages are being installed')

        r = self.initializeR()  # get R ready by checking installed packages and libraries (includes os check)

        functions.setBase(self.RID, 'Step 3 of 4: Performing statistical test...')

        # R graph
        r.assign('finalDF', self.finalDF)

        colorVal = self.all['colorVal']
        if colorVal == 'None':
            r("colorTrt <- c('All')")
        else:
            r.assign("colorVal", colorVal)
            r("colorTrt <- as.factor(finalDF[,paste(colorVal)])")

        r.assign('xVal', self.quantFields[0])

        gridVal_X = self.all['gridVal_X']
        if gridVal_X == 'None':
            r("gridTrt_X <- c('All')")
        else:
            r.assign("gridVal_X", gridVal_X)
            r("gridTrt_X <- as.factor(finalDF[,paste(gridVal_X)])")

        gridVal_Y = self.all['gridVal_Y']
        if gridVal_Y == 'None':
            r("gridTrt_Y <- c('All')")
        else:
            r.assign("gridVal_Y", gridVal_Y)
            r("gridTrt_Y <- as.factor(finalDF[,paste(gridVal_Y)])")

        shapeVal = self.all['shapeVal']
        if shapeVal == 'None':
            r("shapeTrt <- c('All')")
        else:
            r.assign("shapeVal", shapeVal)
            r("shapeTrt <- as.factor(finalDF[,paste(shapeVal)])")

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

        r("gDF <- data.frame(x=finalDF[,paste(xVal)], y=finalDF[,paste(DepVar)], \
            gridVal_X=gridTrt_X, gridVal_Y=gridTrt_Y, \
            myColor=colorTrt, myShape=shapeTrt)")

        r("p <- ggplot(gDF, aes(x=x, y=y, fill=factor(myColor), shape=factor(myShape)) )")
        r("p <- p + geom_point(size=4)")

        if gridVal_X != 'None' and gridVal_Y == 'None':
            r("p <- p + facet_grid(. ~ gridVal_X)")
        elif gridVal_X == 'None' and gridVal_Y != 'None':
            r("p <- p + facet_grid(gridVal_Y ~ .)")
        elif gridVal_X != 'None' and gridVal_Y != 'None':
            r("p <- p + facet_grid(gridVal_Y ~ gridVal_X)")

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

        r('number <- nlevels(gDF$myColor)')
        r('colors <- rep(pal, length.out=number) ')
        r("p <- p + scale_fill_manual(name='', values=colors, guide=guide_legend(override.aes=list(shape=21)))")

        r('number <- nlevels(gDF$myShape)')
        r('shapes <- rep(c(21, 22, 23, 24, 25), length.out=number) ')
        r("p <- p + scale_shape_manual(name='', values=shapes)")

        r("p <- p + theme(legend.text=element_text(size=7))")
        r("p <- p + theme(legend.position='bottom')")

        r("my.formula <- y ~ x")
        r("p <- p + geom_smooth(method='lm', se=T, color='black', formula=my.formula)")

        if self.DepVar == 0:
            r("p <- p + ylab('Abundance') + xlab(paste(xVal))")
        elif self.DepVar == 1:
            r("p <- p + ylab('Relative Abundance') + xlab(paste(xVal))")
        elif self.DepVar == 2:
            r("p <- p + ylab('OTU Richness') + xlab(paste(xVal))")
        elif self.DepVar == 3:
            r("p <- p + ylab('OTU Diversity') + xlab(paste(xVal))")
        elif self.DepVar == 4:
            r("p <- p + ylab('Total Abundance') + xlab(paste(xVal))")

        path = "myPhyloDB/media/temp/anova/Rplots"
        if not os.path.exists(path):
            os.makedirs(path)

        r.assign("path", path)
        r.assign("RID", self.RID)
        r("file <- paste(path, '/', RID, '.anova.pdf', sep='')")
        r("p <- set_panel_size(p, height=unit(2.9, 'in'), width=unit(2.9, 'in'))")

        r("nlev <- nlevels(as.factor(gDF$gridVal_X))")
        r('if (nlev == 0) { \
                myWidth <- 8 \
            } else { \
                myWidth <- min(3*nlev+4, 50) \
        }')

        r("nlev <- nlevels(as.factor(gDF$gridVal_Y))")
        r('if (nlev == 0) { \
                myHeight <- 8 \
            } else { \
                myHeight <- min(3*nlev+4, 50) \
        }')

        r("ggsave(filename=file, plot=p, units='in', height=myHeight, width=myWidth, limitsize=F)")

        pValDict = {}
        counter = 1
        catLevels = len(set(self.catValues))
        grouped1 = self.finalDF.groupby(['rank_name', 'rank_id'])
        for name1, group1 in grouped1:
            D = ''
            r.assign("df", group1)

            trtString = " * ".join(self.allFields)
            if self.DepVar == 0:
                anova_string = "fit <- lm(abund ~ " + str(trtString) + ", data=df)"
                r.assign("cmd", anova_string)
                r("eval(parse(text=cmd))")
            elif self.DepVar == 1:
                anova_string = "fit <- lm(rel_abund ~ " + str(trtString) + ", data=df)"
                r.assign("cmd", anova_string)
                r("eval(parse(text=cmd))")
            elif self.DepVar == 2:
                anova_string = "fit <- lm(rich ~ " + str(trtString) + ", data=df)"
                r.assign("cmd", anova_string)
                r("eval(parse(text=cmd))")
            elif self.DepVar == 3:
                anova_string = "fit <- lm(diversity ~ " + str(trtString) + ", data=df)"
                r.assign("cmd", anova_string)
                r("eval(parse(text=cmd))")
            elif self.DepVar == 4:
                anova_string = "fit <- lm(abund_16S ~ " + str(trtString) + ", data=df)"
                r.assign("cmd", anova_string)
                r("eval(parse(text=cmd))")

            # calculate predicted scores (full model)
            pred_string = "df$pred <- predict(fit, df)"
            r.assign("cmd", pred_string)
            r("eval(parse(text=cmd))")

            r("options(width=5000)")
            r("aov <- anova(fit)")
            aov = r("aov")

            tempStuff = aov.split('\n')
            for i in xrange(len(tempStuff)):
                if i >= 4:
                    D += tempStuff[i] + '\n'

            fit = r("summary(fit)")
            tempStuff = fit.split('\n')
            for i in xrange(len(tempStuff)):
                if i >= 8:
                    D += tempStuff[i] + '\n'
            r("p_vals <- summary(fit)$coefficients[,4]")
            p_vals = r.get("p_vals")

            if not np.isnan(p_vals).any():
                p_value = min(p_vals)
                pValDict[name1] = p_value
            else:
                pValDict[name1] = np.nan

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

            # stop check per loop iteration
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
            if self.stopList[self.PID] == self.RID:
                if self.debug:
                    print "Stopping!"
                return HttpResponse(getStopDict(), content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

        functions.setBase(self.RID, 'Step 3 of 4: Performing statistical test...done!')

        # stop check after step completion
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

        functions.setBase(self.RID, 'Step 4 of 4: Formatting graph data for display...')

        self.finalDF['sample_name'] = ''
        for index, row in self.finalDF.iterrows():
            val = Sample.objects.get(sampleid=row['sampleid']).sample_name
            self.finalDF.loc[index, 'sample_name'] = val

        shapes_idx = 0
        seriesList = []
        grouped1 = self.finalDF.groupby(['rank_name', 'rank_id'])
        for name1, group1 in grouped1:
            pValue = pValDict[name1]

            if self.sig_only == 0:
                if catLevels > 1:
                    grouped2 = group1.groupby(self.catFields)
                    for name2, group2 in grouped2:
                        dataList = []
                        x = []
                        y = []
                        if self.DepVar == 0:
                            x = group2[self.quantFields[0]].values.astype(float).tolist()
                            y = group2['abund'].values.astype(float).tolist()
                        elif self.DepVar == 1:
                            x = group2[self.quantFields[0]].values.astype(float).tolist()
                            y = group2['rel_abund'].values.astype(float).tolist()
                        elif self.DepVar == 2:
                            x = group2[self.quantFields[0]].values.astype(float).tolist()
                            y = group2['rich'].values.astype(float).tolist()
                        elif self.DepVar == 3:
                            x = group2[self.quantFields[0]].values.astype(float).tolist()
                            y = group2['diversity'].values.astype(float).tolist()
                        elif self.DepVar == 4:
                            x = group2[self.quantFields[0]].values.astype(float).tolist()
                            y = group2['abund_16S'].values.astype(float).tolist()

                        if self.DepVar == 0:
                            for index, row in group2.iterrows():
                                dataDict = {}
                                dataDict['name'] = row['sample_name']
                                dataDict['x'] = float(row[self.quantFields[0]])
                                dataDict['y'] = float(row['abund'])
                                dataList.append(dataDict)
                        elif self.DepVar == 1:
                            for index, row in group2.iterrows():
                                dataDict = {}
                                dataDict['name'] = row['sample_name']
                                dataDict['x'] = float(row[self.quantFields[0]])
                                dataDict['y'] = float(row['rel_abund'])
                                dataList.append(dataDict)
                        elif self.DepVar == 2:
                            for index, row in group2.iterrows():
                                dataDict = {}
                                dataDict['name'] = row['sample_name']
                                dataDict['x'] = float(row[self.quantFields[0]])
                                dataDict['y'] = float(row['rich'])
                                dataList.append(dataDict)
                        elif self.DepVar == 3:
                            for index, row in group2.iterrows():
                                dataDict = {}
                                dataDict['name'] = row['sample_name']
                                dataDict['x'] = float(row[self.quantFields[0]])
                                dataDict['y'] = float(row['diversity'])
                                dataList.append(dataDict)
                        elif self.DepVar == 4:
                            for index, row in group2.iterrows():
                                dataDict = {}
                                dataDict['name'] = row['sample_name']
                                dataDict['x'] = float(row[self.quantFields[0]])
                                dataDict['y'] = float(row['abund_16S'])
                                dataList.append(dataDict)

                        seriesDict = {}
                        seriesDict['turboThreshold'] = 0
                        seriesDict['type'] = 'scatter'
                        seriesDict['name'] = str(name1[0]) + ": " + str(name2)
                        seriesDict['data'] = dataList

                        markerDict = {}
                        markerDict['symbol'] = shapes[shapes_idx]
                        seriesDict['marker'] = markerDict
                        seriesDict['data'] = dataList

                        seriesList.append(seriesDict)

                        slp, inter, r_value, p, std_err = stats.linregress(x, y)
                        min_y = float(slp*min(x) + inter)
                        max_y = float(slp*max(x) + inter)
                        slope = "%0.3f" % slp
                        intercept = "%0.3f" % inter
                        r_sq = r_value * r_value
                        r_square = "%0.3f" % r_sq

                        regrList = []
                        regrList.append([float(min(x)), min_y])
                        regrList.append([float(max(x)), max_y])

                        regrDict = {}
                        regrDict['type'] = 'line'
                        sup2 = u"\u00B2"
                        regrDict['name'] = 'y = ' + str(slope) + 'x' + ' + ' + str(intercept) + '; R' + sup2 + ' = ' + str(r_square)
                        regrDict['data'] = regrList
                        regrDict['color'] = 'black'

                        markerDict = {}
                        markerDict['enabled'] = False
                        regrDict['marker'] = markerDict
                        seriesList.append(regrDict)

                        shapes_idx += 1
                        if shapes_idx >= len(shapes):
                            shapes_idx = 0

                else:   # if catLevel <=1
                    dataList = []
                    x = []
                    y = []
                    if self.DepVar == 0:
                        x = group1[self.quantFields[0]].values.astype(float).tolist()
                        y = group1['abund'].values.astype(float).tolist()
                    elif self.DepVar == 1:
                        x = group1[self.quantFields[0]].values.astype(float).tolist()
                        y = group1['rel_abund'].values.astype(float).tolist()
                    elif self.DepVar == 2:
                        x = group1[self.quantFields[0]].values.astype(float).tolist()
                        y = group1['rich'].values.astype(float).tolist()
                    elif self.DepVar == 3:
                        x = group1[self.quantFields[0]].values.astype(float).tolist()
                        y = group1['diversity'].values.astype(float).tolist()
                    elif self.DepVar == 4:
                        x = group1[self.quantFields[0]].values.astype(float).tolist()
                        y = group1['abund_16S'].values.astype(float).tolist()

                    if self.DepVar == 0:
                        for index, row in group1.iterrows():
                            dataDict = {}
                            dataDict['name'] = row['sample_name']
                            dataDict['x'] = float(row[self.quantFields[0]])
                            dataDict['y'] = float(row['abund'])
                            dataList.append(dataDict)
                    elif self.DepVar == 1:
                        for index, row in group1.iterrows():
                            dataDict = {}
                            dataDict['name'] = row['sample_name']
                            dataDict['x'] = float(row[self.quantFields[0]])
                            dataDict['y'] = float(row['rel_abund'])
                            dataList.append(dataDict)
                    elif self.DepVar == 2:
                        for index, row in group1.iterrows():
                            dataDict = {}
                            dataDict['name'] = row['sample_name']
                            dataDict['x'] = float(row[self.quantFields[0]])
                            dataDict['y'] = float(row['rich'])
                            dataList.append(dataDict)
                    elif self.DepVar == 3:
                        for index, row in group1.iterrows():
                            dataDict = {}
                            dataDict['name'] = row['sample_name']
                            dataDict['x'] = float(row[self.quantFields[0]])
                            dataDict['y'] = float(row['diversity'])
                            dataList.append(dataDict)
                    elif self.DepVar == 4:
                        for index, row in group1.iterrows():
                            dataDict = {}
                            dataDict['name'] = row['sample_name']
                            dataDict['x'] = float(row[self.quantFields[0]])
                            dataDict['y'] = float(row['abund_16S'])
                            dataList.append(dataDict)

                    seriesDict = {}
                    seriesDict['turboThreshold'] = 0
                    seriesDict['type'] = 'scatter'
                    seriesDict['name'] = str(name1[0])
                    seriesDict['data'] = dataList

                    markerDict = {}
                    markerDict['symbol'] = shapes[shapes_idx]
                    seriesDict['marker'] = markerDict
                    seriesDict['data'] = dataList

                    seriesList.append(seriesDict)

                    slp, inter, r_value, p, std_err = stats.linregress(x, y)
                    min_y = float(slp*min(x) + inter)
                    max_y = float(slp*max(x) + inter)
                    slope = "%0.3f" % slp
                    intercept = "%0.3f" % inter
                    r_sq = r_value * r_value
                    r_square = "%0.3f" % r_sq

                    regrList = []
                    regrList.append([float(min(x)), min_y])
                    regrList.append([float(max(x)), max_y])

                    regrDict = {}
                    regrDict['type'] = 'line'
                    sup2 = u"\u00B2"
                    regrDict['name'] = 'y = ' + str(slope) + 'x' + ' + ' + str(intercept) + '; R' + sup2 + ' = ' + str(r_square)
                    regrDict['data'] = regrList
                    regrDict['color'] = 'black'

                    markerDict = {}
                    markerDict['enabled'] = False
                    regrDict['marker'] = markerDict
                    seriesList.append(regrDict)

            elif self.sig_only == 1:
                if pValue < 0.05:
                    if catLevels > 1:
                        grouped2 = group1.groupby(self.catFields)
                        for name2, group2 in grouped2:
                            dataList = []
                            x = []
                            y = []
                            if self.DepVar == 0:
                                x = group2[self.quantFields[0]].values.astype(float).tolist()
                                y = group2['abund'].values.astype(float).tolist()
                            elif self.DepVar == 1:
                                x = group2[self.quantFields[0]].values.astype(float).tolist()
                                y = group2['rel_abund'].values.astype(float).tolist()
                            elif self.DepVar == 2:
                                x = group2[self.quantFields[0]].values.astype(float).tolist()
                                y = group2['rich'].values.astype(float).tolist()
                            elif self.DepVar == 3:
                                x = group2[self.quantFields[0]].values.astype(float).tolist()
                                y = group2['diversity'].values.astype(float).tolist()
                            elif self.DepVar == 4:
                                x = group2[self.quantFields[0]].values.astype(float).tolist()
                                y = group2['abund_16S'].values.astype(float).tolist()

                            if self.DepVar == 0:
                                for index, row in group2.iterrows():
                                    dataDict = {}
                                    dataDict['name'] = row['sample_name']
                                    dataDict['x'] = float(row[self.quantFields[0]])
                                    dataDict['y'] = float(row['abund'])
                                    dataList.append(dataDict)
                            elif self.DepVar == 1:
                                for index, row in group2.iterrows():
                                    dataDict = {}
                                    dataDict['name'] = row['sample_name']
                                    dataDict['x'] = float(row[self.quantFields[0]])
                                    dataDict['y'] = float(row['rel_abund'])
                                    dataList.append(dataDict)
                            elif self.DepVar == 2:
                                for index, row in group2.iterrows():
                                    dataDict = {}
                                    dataDict['name'] = row['sample_name']
                                    dataDict['x'] = float(row[self.quantFields[0]])
                                    dataDict['y'] = float(row['rich'])
                                    dataList.append(dataDict)
                            elif self.DepVar == 3:
                                for index, row in group2.iterrows():
                                    dataDict = {}
                                    dataDict['name'] = row['sample_name']
                                    dataDict['x'] = float(row[self.quantFields[0]])
                                    dataDict['y'] = float(row['diversity'])
                                    dataList.append(dataDict)
                            elif self.DepVar == 4:
                                for index, row in group2.iterrows():
                                    dataDict = {}
                                    dataDict['name'] = row['sample_name']
                                    dataDict['x'] = float(row[self.quantFields[0]])
                                    dataDict['y'] = float(row['abund_16S'])
                                    dataList.append(dataDict)

                            seriesDict = {}
                            seriesDict['turboThreshold'] = 0
                            seriesDict['type'] = 'scatter'
                            seriesDict['name'] = str(name1[0]) + ": " + str(name2)
                            seriesDict['data'] = dataList

                            markerDict = {}
                            markerDict['symbol'] = shapes[shapes_idx]
                            seriesDict['marker'] = markerDict
                            seriesDict['data'] = dataList

                            seriesList.append(seriesDict)

                            slp, inter, r_value, p, std_err = stats.linregress(x, y)
                            min_y = float(slp*min(x) + inter)
                            max_y = float(slp*max(x) + inter)
                            slope = "%0.3f" % slp
                            intercept = "%0.3f" % inter
                            r_sq = r_value * r_value
                            r_square = "%0.3f" % r_sq

                            regrList = []
                            regrList.append([float(min(x)), min_y])
                            regrList.append([float(max(x)), max_y])

                            regrDict = {}
                            regrDict['type'] = 'line'
                            regrDict['name'] = 'y = ' + str(slope) + 'x' + ' + ' + str(intercept) + '; R2 = ' + str(r_square)
                            regrDict['data'] = regrList
                            regrDict['color'] = 'black'

                            markerDict = {}
                            markerDict['enabled'] = False
                            regrDict['marker'] = markerDict
                            seriesList.append(regrDict)

                            shapes_idx += 1
                            if shapes_idx >= len(shapes):
                                shapes_idx = 0

                    else:   # if catLevel <=1
                        dataList = []
                        x = []
                        y = []
                        if self.DepVar == 0:
                            x = group1[self.quantFields[0]].values.astype(float).tolist()
                            y = group1['abund'].values.astype(float).tolist()
                        elif self.DepVar == 1:
                            x = group1[self.quantFields[0]].values.astype(float).tolist()
                            y = group1['rel_abund'].values.astype(float).tolist()
                        elif self.DepVar == 2:
                            x = group1[self.quantFields[0]].values.astype(float).tolist()
                            y = group1['rich'].values.astype(float).tolist()
                        elif self.DepVar == 3:
                            x = group1[self.quantFields[0]].values.astype(float).tolist()
                            y = group1['diversity'].values.astype(float).tolist()
                        elif self.DepVar == 4:
                            x = group1[self.quantFields[0]].values.astype(float).tolist()
                            y = group1['abund_16S'].values.astype(float).tolist()

                        if self.DepVar == 0:
                            for index, row in group1.iterrows():
                                dataDict = {}
                                dataDict['name'] = row['sample_name']
                                dataDict['x'] = float(row[self.quantFields[0]])
                                dataDict['y'] = float(row['abund'])
                                dataList.append(dataDict)
                        elif self.DepVar == 1:
                            for index, row in group1.iterrows():
                                dataDict = {}
                                dataDict['name'] = row['sample_name']
                                dataDict['x'] = float(row[self.quantFields[0]])
                                dataDict['y'] = float(row['rel_abund'])
                                dataList.append(dataDict)
                        elif self.DepVar == 2:
                            for index, row in group1.iterrows():
                                dataDict = {}
                                dataDict['name'] = row['sample_name']
                                dataDict['x'] = float(row[self.quantFields[0]])
                                dataDict['y'] = float(row['rich'])
                                dataList.append(dataDict)
                        elif self.DepVar == 3:
                            for index, row in group1.iterrows():
                                dataDict = {}
                                dataDict['name'] = row['sample_name']
                                dataDict['x'] = float(row[self.quantFields[0]])
                                dataDict['y'] = float(row['diversity'])
                                dataList.append(dataDict)
                        elif self.DepVar == 4:
                            for index, row in group1.iterrows():
                                dataDict = {}
                                dataDict['name'] = row['sample_name']
                                dataDict['x'] = float(row[self.quantFields[0]])
                                dataDict['y'] = float(row['abund_16S'])
                                dataList.append(dataDict)

                        seriesDict = {}
                        seriesDict['turboThreshold'] = 0
                        seriesDict['type'] = 'scatter'
                        seriesDict['name'] = str(name1[0])
                        seriesDict['data'] = dataList

                        markerDict = {}
                        markerDict['symbol'] = shapes[shapes_idx]
                        seriesDict['marker'] = markerDict
                        seriesDict['data'] = dataList

                        seriesList.append(seriesDict)

                        slp, inter, r_value, p, std_err = stats.linregress(x, y)
                        min_y = float(slp*min(x) + inter)
                        max_y = float(slp*max(x) + inter)
                        slope = "%0.3f" % slp
                        intercept = "%0.3f" % inter
                        r_sq = r_value * r_value
                        r_square = "%0.3f" % r_sq

                        regrList = []
                        regrList.append([float(min(x)), min_y])
                        regrList.append([float(max(x)), max_y])

                        regrDict = {}
                        regrDict['type'] = 'line'
                        regrDict['name'] = 'y = ' + str(slope) + 'x' + ' + ' + str(intercept) + '; R2 = ' + str(r_square)
                        regrDict['data'] = regrList
                        regrDict['color'] = 'black'

                        markerDict = {}
                        markerDict['enabled'] = False
                        regrDict['marker'] = markerDict
                        seriesList.append(regrDict)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #
            if self.stopList[self.PID] == self.RID:
                if self.debug:
                    print "Stopping!"
                return HttpResponse(getStopDict(), content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ #

        xAxisDict = {}
        xTitle = {}
        xTitle['text'] = self.quantFields[0]
        xTitle['style'] = {'fontSize': '18px', 'fontWeight': 'bold'}
        xAxisDict['title'] = xTitle

        yAxisDict = {}
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
        yAxisDict['title'] = yTitle

        if self.transform != 0:
            tname = {
                '1': "Ln", '2': "Log10", '3': "Sqrt", '4': "Logit", '5': "Arcsin"
            }
            yTitle['text'] = tname[str(self.transform)] + "(" + str(yTitle['text']) + ")"

        yTitle['style'] = {'fontSize': '18px', 'fontWeight': 'bold'}
        yAxisDict['title'] = yTitle

        styleDict = {'style': {'fontSize': '14px'}}
        xAxisDict['labels'] = styleDict
        yAxisDict['labels'] = styleDict

        finalDict['series'] = seriesList
        finalDict['xAxis'] = xAxisDict
        finalDict['yAxis'] = yAxisDict

        if not seriesList:
            finalDict['empty'] = 0
        else:
            finalDict['empty'] = 1

        functions.setBase(self.RID, 'Step 4 of 4: Formatting graph data for display...done!')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        # datatable of taxa mapped to selected kegg orthologies
        if not self.treeType == 1 and self.mapTaxa == 'yes':
            records = self.allDF.values.tolist()
            finalDict['taxData'] = json.dumps(records)
            columns = self.allDF.columns.values.tolist()
            finalDict['taxColumns'] = json.dumps(columns)

        finalDict['resType'] = 'res'
        finalDict['text'] = self.result

        finalDict['error'] = 'none'
        res = json.dumps(finalDict)
        return HttpResponse(res, content_type='application/json')   #
    # ^ WHAT IS THIS MONSTER OF A FUNCTION ^

    def run(self, quant=False):
        if self.debug:
            print "Running Anova. Quant = ", quant

        if quant:
            ret = self.validate(reqMultiLevel=False)
            if ret == 0:
                ret = self.query()
                if ret == 0:
                    return self.quantstats()
        else:
            ret = self.validate()
            if ret == 0:
                ret = self.query()
                if ret == 0:
                    ret = self.stats()
                    if ret == 0:
                        return self.graph()
            # experimental R substitution
            #return self.fullRAnalysis("functions/analysis/R_scripts/anova.r")

        if self.debug:
            if self.stopList[self.PID] == self.RID:
                print "Stopped Anova"
            else:
                print "Something went wrong with Anova"
        return ret


class Corr(Analysis):

    # overwrite stats and graph
    def statsGraph(self):
        print 'statsGraph'
        if self.debug:
            print ("statsgraph! (corr)")
        functions.setBase(self.RID, 'Step 3 of 4: Calculating Correlations Matrix...')
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
        self.result += '\n===============================================\n'

        count_rDF = pd.DataFrame()
        if self.DepVar == 0:
            count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='abund')
        elif self.DepVar == 1:
            count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='rel_abund')
        elif self.DepVar == 2:
            count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='rich')
        elif self.DepVar == 3:
            count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='diversity')
        elif self.DepVar == 4:
            count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='abund_16S')

        count_rDF.fillna(0, inplace=True)

        if os.name == 'nt':
            r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
        else:
            r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

        functions.setBase(self.RID, 'Verifying R packages...missing packages are being installed')

        # R packages from cran
        r("list.of.packages <- c('corrplot', 'RColorBrewer', 'WGCNA')")
        r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
        r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

        r("library(corrplot)")
        r("library(WGCNA)")
        r("library(RColorBrewer)")

        idList = count_rDF.columns.values.tolist()
        namesList = []
        if self.treeType == 1:
            namesDict = functions.getFullTaxonomy(idList)
            od = collections.OrderedDict(sorted(namesDict.items()))  # unexpected argument warning TODO resolve
            namesList = od.values()
            namesList = [i.split('|')[-1] for i in namesList]
        elif self.treeType == 2:
            namesDict = functions.getFullKO(idList)
            od = collections.OrderedDict(sorted(namesDict.items()))
            namesList = od.values()
            namesList = [i.split('|')[-1] for i in namesList]
        elif self.treeType == 3:
            if self.nzAll < 5:
                namesDict = functions.getFullNZ(idList)
                od = collections.OrderedDict(sorted(namesDict.items()))
                namesList = od.values()
                namesList = [i.split('|')[-1] for i in namesList]
            else:
                namesList = [i.split(':', 1)[0] for i in idList]

        count_rDF.sort_index(axis=0, inplace=True)
        self.metaDF.sort_values('sampleid', inplace=True)   # apparently inplace is considered bad practice

        r.assign("X", count_rDF)
        r('X <- X * 1.0')
        r.assign("names", namesList)
        r("colnames(X) <- names")

        if self.quantFields:
            r.assign("Y", self.metaDF[self.quantFields])
            r('M <- cor(as.matrix(Y), as.matrix(X))')
            r('M[is.na(M)] <- 0')
            r("M.p <- corPvalueFisher(M, nrow(X))")
        else:
            r('M <- cor(X)')
            r('M[is.na(M)] <- 0')
            r("M.p <- corPvalueFisher(M, nrow(X))")

        path = "myPhyloDB/media/temp/corr/Rplots/" + str(self.RID) + ".corr.pdf"
        if os.path.exists(path):
            os.remove(path)

        if not os.path.exists('myPhyloDB/media/temp/corr/Rplots'):
            os.makedirs('myPhyloDB/media/temp/corr/Rplots')

        row, col = count_rDF.shape
        width = 2 + col*0.15
        if self.quantFields:
            height = 4 + len(self.quantFields)*0.1
        else:
            height = width
        file = "pdf('myPhyloDB/media/temp/corr/Rplots/" + str(self.RID) + ".corr.pdf', height=" + str(height) + ", width=" + str(width) + ")"
        r.assign("cmd", file)
        r("eval(parse(text=cmd))")

        if self.quantFields:
            r('corrplot(M, p.mat=M.p, method="pie", sig.level=0.05, insig="blank", cl.length=5, tl.cex=0.7, cl.cex=0.5)')
        else:
            r('corrplot(M, p.mat=M.p, method="pie", sig.level=0.05, insig="blank", cl.length=5, tl.cex=0.7, cl.cex=0.5, order="hclust")')

        r("dev.off()")

        functions.setBase(self.RID, 'Step 3 of 4: Calculating correlation matrix...done!')


        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        functions.setBase(self.RID, 'Step 4 of 4: Formatting graph data for display...')

        finalDict = {}
        dat = r.get("M")
        if self.quantFields:
            coeffsDF = pd.DataFrame(dat, columns=namesList, index=self.quantFields)
        else:
            coeffsDF = pd.DataFrame(dat, columns=namesList, index=namesList)
        res_table = coeffsDF.to_html(classes="table display")
        res_table = res_table.replace('border="1"', 'border="0"')
        finalDict['coeff_table'] = str(res_table)

        dat = r.get("M.p")
        if self.quantFields:
            coeffsDF = pd.DataFrame(dat, columns=namesList, index=self.quantFields)
        else:
            coeffsDF = pd.DataFrame(dat, columns=namesList, index=namesList)
        res_table = coeffsDF.to_html(classes="table display")
        res_table = res_table.replace('border="1"', 'border="0"')
        finalDict['p_table'] = str(res_table)

        finalDict['text'] = self.result

        functions.setBase(self.RID, 'Step 4 of 4: Formatting graph data for display...done!')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        finalDict['error'] = 'none'
        res = json.dumps(finalDict)
        return HttpResponse(res, content_type='application/json')

    def run(self):
        if self.debug:
            print "Running Corr"

        ret = self.validate(sig=False, metaCat=False, metaQuant=True, reqMultiLevel=False)
        if ret == 0:
            ret = self.query(taxmap=False)
            if ret == 0:
                return self.statsGraph()
        if self.debug:
            if self.stopList[self.PID] == self.RID:
                print "Stopped Corr"
            else:
                print "Something went wrong with Corr"
        return ret


class diffAbund(Analysis):

    def statsGraph(self):
        if self.debug:
            print "statsGraph! (diffabund)"

        if len(self.catFields) > 1:
            for index, row in self.metaDF.iterrows():
                self.metaDF.loc[index, 'merge'] = ".".join(row[self.catFields])
        else:
            self.metaDF.loc[:, 'merge'] = self.metaDF.loc[:, self.catFields[0]]

        functions.setBase(self.RID, 'Step 3 of 4: Performing statistical test...')
        count_rDF = pd.DataFrame()
        if self.DepVar == 0:
            count_rDF = self.finalDF.pivot(index='rank_id', columns='sampleid', values='abund')
        elif self.DepVar == 4:
            count_rDF = self.finalDF.pivot(index='rank_id', columns='sampleid', values='abund_16S')

        if self.debug:
            print "Check!"

        count_rDF.fillna(0, inplace=True)

        self.finalDict = {}
        if os.name == 'nt':
            r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
        else:
            r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

        if self.debug:
            print "Check!"

        functions.setBase(self.RID, 'Verifying R packages...missing packages are being installed')

        r("list.of.packages <- c('edgeR')")
        r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
        r("if (length(new.packages)) source('http://bioconductor.org/biocLite.R')")
        r("if (length(new.packages)) biocLite(new.packages, type='source', suppressUpdate=T, dependencies=T)")

        functions.setBase(self.RID, 'Step 4 of 6: Performing statistical test...')

        if self.debug:
            print "Check!"

        r("library(edgeR)")

        if self.DepVar == 0:
            self.result += 'Dependent Variable: Abundance' + '\n'
        elif self.DepVar == 4:
            self.result += 'Dependent Variable: Total Abundance' + '\n'
        self.result += '\n===============================================\n\n\n'

        if self.debug:
            print "Check!"

        myList = list(self.metaDF.select_dtypes(include=['object']).columns)
        for i in myList:
            self.metaDF[i] = self.metaDF[i].str.replace(' ', '_')
            self.metaDF[i] = self.metaDF[i].str.replace('-', '.')
            self.metaDF[i] = self.metaDF[i].str.replace('(', '.')
            self.metaDF[i] = self.metaDF[i].str.replace(')', '.')

        self.metaDF.sort_values('sampleid', inplace=True)
        r.assign("metaDF", self.metaDF)
        r("trt <- factor(metaDF$merge)")

        r.assign("count", count_rDF)
        r('e <- DGEList(counts=count)')
        r('e <- calcNormFactors(e, method="none")')

        r('design <- model.matrix(~ 0 + trt)')    # here
        r('trtLevels <- levels(trt)')
        r('colnames(design) <- trtLevels')

        r('e <- estimateGLMCommonDisp(e, design)')
        r('e <- estimateGLMTrendedDisp(e, design)')
        r('e <- estimateGLMTagwiseDisp(e, design)')
        r('fit <- glmFit(e, design)')
        fit = r.get('fit')

        if self.debug:
            print "Check!"

        if not fit:
            error = "edgeR failed!\nUsually this is caused by one or more taxa having a negative disperion.\nTry filtering your data to remove problematic taxa (e.g. remove phylotypes with 50% or more zeros)."
            myDict = {'error': error}
            res = json.dumps(myDict)
            return HttpResponse(res, content_type='application/json')

        nTopTags = int(self.all['nTopTags'])
        r.assign('nTopTags', nTopTags)

        if self.debug:
            print "Check!"

        mergeList = self.metaDF['merge'].tolist()
        mergeSet = list(set(mergeList))
        nbinom_res = pd.DataFrame()
        for i, val in enumerate(mergeSet):
            start = i + 1
            stop = int(len(mergeSet))
            for j in range(start, stop):
                if i != j:
                    r.assign("trt1", mergeSet[i])
                    r.assign("trt2", mergeSet[j])

                    r('contVec <- sprintf("%s-%s", trt1, trt2)')
                    r('cont.matrix= makeContrasts(contVec, levels=design)')
                    r('lrt <- glmLRT(fit, contrast=cont.matrix)')
                    r("res <- as.data.frame(topTags(lrt, sort.by='PValue', n=min(nTopTags, nrow(fit$table))))")
                    r('res <- res[ order(row.names(res)), ]')
                    taxaIDs = r.get("row.names(res)")

                    baseMean = count_rDF.mean(axis=1)
                    baseMean = baseMean.loc[baseMean.index.isin(taxaIDs)]

                    listA = self.metaDF[self.metaDF['merge'] == mergeSet[i]].sampleid.tolist()
                    baseMeanA = count_rDF[listA].mean(axis=1)
                    baseMeanA = baseMeanA.loc[baseMeanA.index.isin(taxaIDs)]

                    listB = self.metaDF[self.metaDF['merge'] == mergeSet[j]].sampleid.tolist()
                    baseMeanB = count_rDF[listB].mean(axis=1)
                    baseMeanB = baseMeanB.loc[baseMeanB.index.isin(taxaIDs)]

                    r.assign("baseMean", baseMean)
                    r.assign("baseMeanA", baseMeanA)
                    r.assign("baseMeanB", baseMeanB)

                    r('baseMean <- baseMean[ order(as.numeric(row.names(baseMean))), ]')
                    r('baseMeanA <- baseMeanA[ order(as.numeric(row.names(baseMeanA))), ]')
                    r('baseMeanB <- baseMeanB[ order(as.numeric(row.names(baseMeanB))), ]')

                    r("df <- data.frame(rank_id=rownames(res), baseMean=baseMean, baseMeanA=baseMeanA, \
                         baseMeanB=baseMeanB, logFC=-res$logFC, logCPM=res$logCPM, \
                         LR=res$LR, pval=res$PValue, FDR=res$FDR)")
                    df = r.get("df")

                    if df is None:
                        myDict = {'error': "edgeR failed!\nPlease try a different data combination."}
                        res = json.dumps(myDict)
                        return HttpResponse(res, content_type='application/json')

                    # remove taxa that failed (i.e., both trts are zero or log2FoldChange is NaN)
                    df = df.loc[pd.notnull(df[' logFC '])]

                    if self.treeType == 1:
                        idList = functions.getFullTaxonomy(list(df.rank_id.unique()))
                        df['Taxonomy'] = df['rank_id'].map(idList)
                    elif self.treeType == 2:
                        idList = functions.getFullKO(list(df.rank_id.unique()))
                        df['Taxonomy'] = df['rank_id'].map(idList)
                    elif self.treeType == 3:
                        if self.nzAll < 5:
                            idList = functions.getFullNZ(list(df.rank_id.unique()))
                            df['Taxonomy'] = df['rank_id'].map(idList)

                    df.rename(columns={'rank_id': 'Rank ID'}, inplace=True)

                    iterationName = str(mergeSet[i]) + ' vs ' + str(mergeSet[j])
                    df.insert(1, 'Comparison', iterationName)
                    df.rename(columns={' baseMean ': 'baseMean'}, inplace=True)
                    df.rename(columns={' baseMeanA ': 'baseMeanA'}, inplace=True)
                    df.rename(columns={' baseMeanB ': 'baseMeanB'}, inplace=True)
                    df.rename(columns={' logFC ': 'logFC'}, inplace=True)
                    df['logFC'] = df['logFC'].round(3).astype(float)
                    df.rename(columns={' logCPM ': 'logCPM'}, inplace=True)
                    df.rename(columns={' LR ': 'Likelihood Ratio'}, inplace=True)
                    df['Likelihood Ratio'] = df['Likelihood Ratio'].round(3).astype(float)
                    df.rename(columns={' pval ': 'p-value'}, inplace=True)
                    df.rename(columns={' FDR ': 'False Discovery Rate'}, inplace=True)

                    functions.setBase(self.RID, 'Step 3 of 4: Performing statistical test...' + str(iterationName) + ' is done!')

                    if nbinom_res.empty:
                        nbinom_res = df.copy()
                    else:
                        nbinom_res = nbinom_res.append(df)

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if self.stopList[self.PID] == self.RID:
                        if self.debug:
                            print "Stopping!"
                        return HttpResponse(getStopDict(), content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if self.stopList[self.PID] == self.RID:
                if self.debug:
                    print "Stopping!"
                return HttpResponse(getStopDict(), content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        if self.debug:
            print "Check!"

        functions.setBase(self.RID, 'Step 3 of 4: Performing statistical test...done!')
        functions.setBase(self.RID, 'Step 4 of 4: Formatting graph data for display...')

        seriesList = []
        xAxisDict = {}
        yAxisDict = {}

        nbinom_res.fillna(value=1.0, inplace=True)
        nbinom_res.reset_index(drop=True, inplace=True)
        grouped = nbinom_res.groupby('Comparison')

        listOfShapes = ['circle', 'square', 'triangle', 'triangle-down', 'diamond']
        shapeIterator = 0

        if self.debug:
            print "Check!"

        FdrVal = float(self.all['FdrVal'])
        for name, group in grouped:
            nosigDF = group[group["False Discovery Rate"] > FdrVal]
            nosigData = []
            for index, row in nosigDF.iterrows():
                dataDict = {}
                dataDict['name'] = 'ID: ' + row['Rank ID']
                dataDict['x'] = float(row['logCPM'])
                dataDict['y'] = float(row['logFC'])
                nosigData.append(dataDict)

            seriesDict = {}
            seriesDict['name'] = "NotSig: " + str(name)
            seriesDict['data'] = nosigData

            markerDict = {}
            markerDict['symbol'] = listOfShapes[shapeIterator]
            seriesDict['marker'] = markerDict
            seriesList.append(seriesDict)

            sigDF = group[group["False Discovery Rate"] <= FdrVal]
            sigData = []
            for index, row in sigDF.iterrows():
                dataDict = {}
                dataDict['name'] = row['Rank ID']
                dataDict['x'] = float(row['logCPM'])
                dataDict['y'] = float(row['logFC'])
                sigData.append(dataDict)

            seriesDict = {}
            seriesDict['name'] = "Sig: " + str(name)
            seriesDict['data'] = sigData

            markerDict = {}
            markerDict['symbol'] = listOfShapes[shapeIterator]
            seriesDict['marker'] = markerDict
            seriesList.append(seriesDict)

            shapeIterator += 1
            if shapeIterator >= len(listOfShapes):
                shapeIterator = 0

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if self.stopList[self.PID] == self.RID:
                if self.debug:
                    print "Stopping!"
                return HttpResponse(getStopDict(), content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        if self.debug:
            print "Check!"

        xTitle = {}
        xTitle['text'] = "logCPM"
        xTitle['style'] = {'fontSize': '18px', 'fontWeight': 'bold'}
        xAxisDict['title'] = xTitle
        xAxisDict['type'] = 'linear'

        yTitle = {}
        yTitle['text'] = "logFC"
        yTitle['style'] = {'fontSize': '18px', 'fontWeight': 'bold'}
        yAxisDict['title'] = yTitle
        yAxisDict['type'] = 'linear'

        styleDict = {'style': {'fontSize': '14px'}}
        xAxisDict['labels'] = styleDict
        yAxisDict['labels'] = styleDict

        if self.debug:
            print "Check!"

        self.finalDict['series'] = seriesList
        self.finalDict['xAxis'] = xAxisDict
        self.finalDict['yAxis'] = yAxisDict

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        nbinom_res.replace(to_replace='N/A', value=np.nan, inplace=True)
        nbinom_res.dropna(axis=1, how='all', inplace=True)
        res_table = nbinom_res.to_html(classes="table display")
        res_table = res_table.replace('border="1"', 'border="0"')
        self.finalDict['res_table'] = str(res_table)

        self.finalDict['text'] = self.result

        functions.setBase(self.RID, 'Step 4 of 4: Formatting results for display...done!')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        self.finalDict['error'] = 'none'
        res = json.dumps(self.finalDict)

        return HttpResponse(res, content_type='application/json')

    def run(self):
        if self.debug:
            print "Running diffAbund"
        ret = self.validate(sig=False, metaQuant=False, reqMultiLevel=False)
        if ret == 0:
            ret = self.query(taxmap=False, usetransform=False)
            if ret == 0:
                return self.statsGraph()

        if self.debug:
            print "Something went wrong with diffAbund"
        return ret


class PCA(Analysis):

    def statsGraph(self):

        if self.debug:
            print "statsgraph! (PCA)"

        method = self.all["Method"]
        scale = self.all['scaled']
        constrain = self.all["constrain"]
        PC1 = int(self.all["PC1"])
        PC2 = int(self.all["PC2"])

        count_rDF = pd.DataFrame()
        if self.DepVar == 0:
            count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='abund')
        elif self.DepVar == 1:
            count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='rel_abund')
        elif self.DepVar == 2:
            count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='rich')
        elif self.DepVar == 3:
            count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='diversity')
        elif self.DepVar == 4:
            count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='abund_16S')

        count_rDF.fillna(0, inplace=True)

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        functions.setBase(self.RID, 'Step 3 of 4: Performing statistical test...')

        if os.name == 'nt':
            r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
        else:
            r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

        functions.setBase(self.RID, 'Verifying R packages...missing packages are being installed')

        r("list.of.packages <- c('fpc', 'vegan', 'ggplot2')")
        r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
        r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")


        r("options(width=5000)")
        r('library(fpc)')
        r('library(ggplot2)')
        r('library(vegan)')
        r('source("R/myFunctions/myFunctions.R")')

        count_rDF.sort_index(axis=0, inplace=True)
        r.assign("data", count_rDF)
        r.assign("cols", count_rDF.columns.values.tolist())
        r("colnames(data) <- unlist(cols)")

        self.metaDF.sort_values('sampleid', inplace=True)
        r.assign("meta", self.metaDF)
        r.assign("rows", self.metaDF.index.values.tolist())
        r("rownames(meta) <- unlist(rows)")

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        r.assign("PC1", PC1)
        r.assign("PC2", PC2)

        # Can only constrain if meta-variables have been selected
        constrain2 = 'no'
        if constrain == 'yes':
            if self.catFields or self.quantFields:
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

        self.result += str(r('print(res.pca)')) + '\n'
        self.result += '===============================================\n'

        addContrib1 = self.all['addContrib1']
        contribVal1 = float(self.all['contribVal1'])
        addContrib2 = self.all['addContrib2']
        contribVal2 = float(self.all['contribVal2'])

        # Use vegan to calculate regression between ord axes and quantFields
        if addContrib1 == 'yes':
            r("ef1 <- envfit(res.pca, data, add=False)")

        if self.quantFields and addContrib2 == 'yes':
            r.assign('quantFields', self.quantFields)
            r('ef2 <- envfit(res.pca, meta[,paste(quantFields)], add=False)')

        # get scores from vegan
        r('sites <- scores(res.pca, display="sites", choices=c(PC1,PC2))')
        r('species <- scores(res.pca, display="species", choices=c(PC1,PC2))')

        ellipseVal = self.all['ellipseVal']
        if ellipseVal == 'None':
            r("ellipseTrt <- c('All')")
        if ellipseVal == 'interaction':
            r.assign("catFields", self.catFields)
            r("ellipseTrt <- interaction(meta[,paste(catFields)])")
        if ellipseVal != 'None' and ellipseVal != 'k-means' and ellipseVal != 'interaction':
            r.assign("ellipseVal", ellipseVal)
            r("ellipseTrt <- as.factor(meta[,paste(ellipseVal)])")
        if ellipseVal != 'None' and ellipseVal == 'k-means':
            r("pamk.best <- pamk(sites)")
            r("km <- kmeans(sites, centers=pamk.best$nc)")
            r("ellipseTrt <- as.factor(paste('k-cluster: ', km$cluster, sep=''))")
        r("if (!exists('ellipseTrt')) {ellipseTrt <- c('All')}")

        colorVal = self.all['colorVal']
        if colorVal == 'None':
            r("colorTrt <- c('All')")
        if colorVal == 'interaction':
            r.assign("catFields", self.catFields)
            r("colorTrt <- interaction(meta[,paste(catFields)])")
        if colorVal != 'None' and colorVal != 'k-means' and colorVal != 'interaction':
            r.assign("colorVal", colorVal)
            r("colorTrt <- as.factor(meta[,paste(colorVal)])")
        if colorVal != 'None' and colorVal == 'k-means':
            r("pamk.best <- pamk(sites)")
            r("km <- kmeans(sites, centers=pamk.best$nc)")
            r("colorTrt <- as.factor(paste('k-cluster: ', km$cluster, sep=''))")
        r("if (!exists('colorTrt')) {colorTrt <- c('All')}")

        shapeVal = self.all['shapeVal']
        if shapeVal == 'None':
            r("shapeTrt <- 'All'")
        if shapeVal == 'interaction':
            r.assign("catFields", self.catFields)
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

        gridVal_X = self.all['gridVal_X']
        if gridVal_X != 'None':
            r.assign("gridVal_X", gridVal_X)
            r("indDF$myGrid_X <- meta[,paste(gridVal_X)]")

        gridVal_Y = self.all['gridVal_Y']
        if gridVal_Y != 'None':
            r.assign("gridVal_Y", gridVal_Y)
            r("indDF$myGrid_Y <- meta[,paste(gridVal_Y)]")

        r("varDF <- data.frame( \
            x=species[,PC1], \
            y=species[,PC2]) \
        ")

        # get taxa rank names
        rankNameDF = self.finalDF.drop_duplicates(subset='rank_id', keep='last')
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

        myPalette = self.all['palette']
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
            myCI = float(self.all["CI"])
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
            self.result += 'Envfit results for species scores\n'
            self.result += str(envfit) + '\n'
            self.result += '===============================================\n'

            if not efDF_adj.empty:
                r("p <- p + geom_segment(data=efDF.adj, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,'cm')), alpha=0.75, color='blue')")
                r("p <- p + geom_text(data=efDF.adj, aes(x=v1, y=v2, label=label, vjust=ifelse(v2 >= 0, -1, 2)), size=3, color='blue')")

        if self.quantFields and addContrib2 == 'yes':
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
            self.result += 'EnvFit for selected quantitative variables\n'
            self.result += str(envfit) + '\n'
            self.result += '===============================================\n'

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
        r.assign("RID", self.RID)
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

        functions.setBase(self.RID, 'Step 3 of 4: Performing statistical test...done')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        functions.setBase(self.RID, 'Step 4 of 4: Formatting graph data...')

        finalDict = {}
        r("options(width=5000)")
        finalDict['text'] = self.result

        ## variables
        nameDF = self.finalDF[['rank_id']].drop_duplicates(subset='rank_id', keep='last')
        nameDF.set_index('rank_id', inplace=True)

        r("df <- data.frame(species)")
        tempDF = r.get("df")
        IDs = r.get("row.names(df)")
        tempDF['id'] = IDs
        tempDF.set_index('id', inplace=True)
        varCoordDF = pd.merge(nameDF, tempDF, left_index=True, right_index=True, how='inner')
        varCoordDF.reset_index(drop=False, inplace=True)
        varCoordDF.rename(columns={'index': 'rank_id'}, inplace=True)

        if self.treeType == 1:
            idList = functions.getFullTaxonomy(list(varCoordDF.rank_id.unique()))
            varCoordDF['Taxonomy'] = varCoordDF['rank_id'].map(idList)
        elif self.treeType == 2:
            idList = functions.getFullKO(list(varCoordDF.rank_id.unique()))
            varCoordDF['Taxonomy'] = varCoordDF['rank_id'].map(idList)
        elif self.treeType == 3:
            idList = functions.getFullNZ(list(varCoordDF.rank_id.unique()))
            varCoordDF['Taxonomy'] = varCoordDF['rank_id'].map(idList)

        varCoordDF.replace(to_replace='N/A', value=np.nan, inplace=True)
        varCoordDF.dropna(axis=1, how='all', inplace=True)
        table = varCoordDF.to_html(classes="table display")
        table = table.replace('border="1"', 'border="0"')
        finalDict['varCoordDF'] = str(table)

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        if ellipseVal == 'k-means' or colorVal == 'k-means' or shapeVal == 'k-means':
            r("df <- data.frame(km$cluster, sites)")
        else:
            r("df <- data.frame(sites)")

        tempDF = r.get("df")
        if not self.metaDF.empty:
            tempDF['id'] = self.metaDF.index.values.tolist()
            tempDF.set_index('id', inplace=True)
            indCoordDF = pd.merge(self.metaDF, tempDF, left_index=True, right_index=True, how='inner')
            indCoordDF.reset_index(drop=False, inplace=True)
            indCoordDF.rename(columns={'index': 'rank_id', ' km.cluster ': 'k-means cluster'}, inplace=True)
        else:
            indCoordDF = tempDF.copy()
            indCoordDF.rename(columns={' km.cluster ': 'k-means cluster'}, inplace=True)
        table = indCoordDF.to_html(classes="table display")
        table = table.replace('border="1"', 'border="0"')
        finalDict['indCoordDF'] = str(table)

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        finalDict['error'] = 'none'
        res = json.dumps(finalDict)
        return HttpResponse(res, content_type='application/json')

    def run(self):
        if self.debug:
            print "Running PCA+"
        ret = self.validate(sig=False)
        if ret == 0:
            ret = self.query(taxmap=False)
            if ret == 0:
                return self.statsGraph()

        if self.debug:
            print "Something went wrong with PCA+"
        return ret


class PCoA(Analysis):

    def statsGraph(self):

        if self.debug:
            print "statsgraph! (PCoA)"

        PC1 = int(self.all["PC1"])
        PC2 = int(self.all["PC2"])
        test = int(self.all["test"])
        perms = int(self.all["perms"])

        count_rDF = pd.DataFrame()
        if self.DepVar == 0:
            count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='abund')
        elif self.DepVar == 1:
            count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='rel_abund')
        elif self.DepVar == 2:
            count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='rich')
        elif self.DepVar == 3:
            count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='diversity')
        elif self.DepVar == 4:
            count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='abund_16S')

        count_rDF.fillna(0, inplace=True)

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        functions.setBase(self.RID, 'Step 3 of 9: Calculating distance matrix...')

        if os.name == 'nt':
            r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
        else:
            r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

        functions.setBase(self.RID, 'Verifying R packages...missing packages are being installed')

        r("list.of.packages <- c('vegan', 'ggplot2', 'data.table')")
        r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
        r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

        r("options(width=5000)")
        r("library(vegan)")
        r("library(ggplot2)")
        r('library(data.table)')
        r('source("R/myFunctions/myFunctions.R")')

        count_rDF.sort_index(axis=0, inplace=True)
        r.assign("data", count_rDF)
        r.assign("cols", count_rDF.columns.values.tolist())
        r("colnames(data) <- cols")

        if self.distance == 1:
            r("dist <- vegdist(data, method='manhattan')")
        elif self.distance == 2:
            r("dist <- vegdist(data, method='euclidean')")
        elif self.distance == 3:
            r("dist <- vegdist(data, method='canberra')")
        elif self.distance == 4:
            r("dist <- vegdist(data, method='bray')")
        elif self.distance == 5:
            r("dist <- vegdist(data, method='kulczynski')")
        elif self.distance == 6:
            r("dist <- vegdist(data, method='jaccard')")
        elif self.distance == 7:
            r("dist <- vegdist(data, method='gower')")
        elif self.distance == 8:
            r("dist <- vegdist(data, method='altGower')")
        elif self.distance == 9:
            r("dist <- vegdist(data, method='morisita')")
        elif self.distance == 10:
            r("dist <- vegdist(data, method='horn')")
        elif self.distance == 11:
            r("dist <- vegdist(data, method='mountford')")
        elif self.distance == 12:
            r("dist <- vegdist(data, method='binomial')")
        elif self.distance == 13:
            r("dist <- vegdist(data, method='chao')")
        elif self.distance == 14:
            r("dist <- vegdist(data, method='cao')")
        elif self.distance == 15:
            datamtx = np.asarray(count_rDF)
            dists = functions.wOdum(datamtx, self.alpha)
            r.assign("dist", dists)
            r("dist <- as.dist(dist)")

        r("mat <- as.matrix(dist, diag=TRUE, upper=TRUE)")
        mat = r.get("mat")

        if len(self.catFields) > 1:
            for index, row in self.metaDF.iterrows():
                self.metaDF.loc[index, 'merge'] = ".".join(row[self.catFields])
        else:
            self.metaDF.loc[:, 'merge'] = self.metaDF.loc[:, self.catFields[0]]

        self.metaDF.sort_values('sampleid', inplace=True)
        rowList = self.metaDF.sampleid.values.tolist()
        distDF = pd.DataFrame(mat, columns=[rowList], index=rowList)

        functions.setBase(self.RID, 'Step 3 of 9: Calculating distance matrix...done!')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        functions.setBase(self.RID, 'Step 4 of 9: Principal coordinates analysis...')

        trtLength = len(set(self.catValues))
        trtString = " * ".join(self.catFields)

        bigf = ''
        r.assign("PC1", PC1)
        r.assign("PC2", PC2)

        addContrib2 = self.all['addContrib2']
        contribVal2 = float(self.all['contribVal2'])
        r.assign("meta", self.metaDF)
        method = self.all['Method']
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
            r('ord <- wcmdscale(dist, eig=TRUE)')

        self.result += str(r('print(ord)')) + '\n'
        self.result += '===============================================\n'

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

        if self.quantFields:
            r.assign("quantFields", self.quantFields)
            r("ef <- envfit(ord, meta[,paste(quantFields)], add=False)")

            # create dataframe from envfit for export and adding to biplot
            r('efDF <- as.data.frame(ef$vectors$arrows*ef$vectors$r)')
            r('efDF$r2 <- ef$vectors$r')
            r('efDF$p <- ef$vectors$pvals')
            r('pvals.adj <- round(p.adjust(efDF$p, method="BH"),3)')
            r('efDF$p.adj <- pvals.adj')

            # send data to result string
            envfit = r("efDF")
            self.result += 'EnvFit for selected quantitative variables\n'
            self.result += str(envfit) + '\n'
            self.result += '===============================================\n'

        colorVal = self.all['colorVal']
        if colorVal == 'None':
            r("colorTrt <- c('All')")
        if colorVal == 'interaction':
            r.assign("catFields", self.catFields)
            r("colorTrt <- interaction(meta[,paste(catFields)])")
        if colorVal != 'None' and colorVal != 'interaction':
            r.assign("colorVal", colorVal)
            r("colorTrt <- as.factor(meta[,paste(colorVal)])")
        r("if (!exists('colorTrt')) {colorTrt <- c('All')}")

        shapeVal = self.all['shapeVal']
        if shapeVal == 'None':
            r("shapeTrt <- c('All')")
        if shapeVal == 'interaction':
            r.assign("catFields", self.catFields)
            r("shapeTrt <- interaction(meta[,paste(catFields)])")
        if shapeVal != 'None' and shapeVal != 'interaction':
            r.assign("shapeVal", shapeVal)
            r("shapeTrt <- as.factor(meta[,paste(shapeVal)])")
        r("if (!exists('shapeTrt')) {shapeTrt <- c('All')}")

        ellipseVal = self.all['ellipseVal']
        if ellipseVal == 'None':
            r("ellipseTrt <- c('All')")
        if ellipseVal != 'None' and ellipseVal != 'interaction':
            r.assign("ellipseVal", ellipseVal)
            r("ellipseTrt <- as.factor(meta[,paste(ellipseVal)])")
        if ellipseVal == 'interaction':
            r.assign("catFields", self.catFields)
            r("ellipseTrt <- interaction(meta[,paste(catFields)])")
        r("if (!exists('ellipseTrt')) {ellipseTrt <- c('All')}")

        surfVal = self.all['surfVal']
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

        gridVal_X = self.all['gridVal_X']
        if gridVal_X != 'None':
            r.assign("gridVal_X", gridVal_X)
            r("indDF$myGrid_X <- meta[,paste(gridVal_X)]")

        gridVal_Y = self.all['gridVal_Y']
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

        myPalette = self.all['palette']
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
            myCI = float(self.all["CI"])
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

        if self.quantFields and addContrib2 == 'yes':
            # scale and remove non-significant objects from efDF
            r('names(efDF) <- c("PC1", "PC2", "r2", "p", "p.adj")')
            r('efDF$label <- unlist(quantFields)')
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
        r.assign("RID", self.RID)
        r("file <- paste(path, '/', RID, '.pcoa.pdf', sep='')")
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

        functions.setBase(self.RID, 'Step 4 of 9: Principal coordinates analysis...done!')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        functions.setBase(self.RID, 'Step 5 of 9: Performing perMANOVA...')

        if perms < 10:
            bigf = 'A minimum of 10 permutations is required...'
        elif len(self.catFields) == 0 or trtLength <= 1:
            bigf = 'No categorical variables are available for perMANOVA/betaDisper analysis'
        elif perms >= 10 and len(self.catFields) > 0:
            if test == 1:
                if len(self.catFields) == 1:
                    trtString = str(self.catFields[0])
                    factor_string = trtString + " <- factor(meta$" + trtString + ")"
                    r.assign("cmd", factor_string)
                    r("eval(parse(text=cmd))")
                    amova_string = "res <- adonis(dist ~ " + str(trtString) + ", perms=perms)"
                else:
                    for i in self.catFields:
                        factor_string = str(i) + " <- factor(meta$" + str(i) + ")"
                        r.assign("cmd", factor_string)
                        r("eval(parse(text=cmd))")
                    trtString = " * ".join(self.catFields)
                    amova_string = "res <- adonis(dist ~ " + str(trtString) + ", perms=perms, data)"

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if self.stopList[self.PID] == self.RID:
                        if self.debug:
                            print "Stopping!"
                        return HttpResponse(getStopDict(), content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                r.assign("perms", perms)
                r.assign("cmd", amova_string)
                r("eval(parse(text=cmd))")

                res_aov = r("res$aov.tab")

                tempStuff = res_aov.split('\n')
                for part in tempStuff:
                    if part != tempStuff[0]:
                        bigf += part + '\n'
                functions.setBase(self.RID, 'Step 5 of 9: Performing perMANOVA...done!')

            elif test == 2:
                functions.setBase(self.RID, 'Step 5 of 9: Principal coordinates analysis...done!')
                functions.setBase(self.RID, 'Step 6 of 9: Performing BetaDisper...')

                r.assign("perms", perms)

                r("res <- betadisper(dist, meta$merge)")

                r("something <- anova(res)")
                beta = r("something")
                tempStuff = beta.split('\n')
                for part in tempStuff:
                    if part != tempStuff[0]:
                        bigf += part + '\n'

                betaString = str(r('res'))
                lines = betaString.split('\n')
                for line in lines[1:]:
                    bigf += str(line) + '\n'

                # get centroids
                bigf += 'Treatment centroids for first 6 PC axes:\n'
                r('centroids <- data.frame(res$centroids[,1:6])')
                centroids = str(r('centroids'))
                lines = centroids.split('\n')
                for line in lines[1:]:
                    bigf += str(line) + '\n'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if self.stopList[self.PID] == self.RID:
                    if self.debug:
                        print "Stopping!"
                    return HttpResponse(getStopDict(), content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(self.RID, 'Step 6 of 9: Performing BetaDisper...done!')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        functions.setBase(self.RID, 'Step 7 of 9: Formatting graph data for display...')
        finalDict = {}
        seriesList = []
        xAxisDict = {}
        yAxisDict = {}

        CAP1 = PC1 + len(self.catFields) + len(self.quantFields) + 2
        CAP2 = PC2 + len(self.catFields) + len(self.quantFields) + 2

        if self.catFields:
            grouped = pcoaDF.groupby(self.catFields)
            for name, group in grouped:
                if len(self.catFields) > 1:
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
                if self.stopList[self.PID] == self.RID:
                    if self.debug:
                        print "Stopping!"
                    return HttpResponse(getStopDict(), content_type='application/json')
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
            self.result += 'perMANOVA results:' + '\n'
        if test == 2:
            self.result += 'betaDisper results:' + '\n'

        if len(self.catFields) == 0:
            self.result += 'test cannot be run...' + '\n'
        else:
            bigf = bigf.decode('utf-8')
            self.result += bigf + '\n'

        self.result += '===============================================\n\n\n'

        finalDict['text'] = self.result

        functions.setBase(self.RID, 'Step 7 of 9: Formatting graph data for display...done!')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        functions.setBase(self.RID, 'Step 8 of 9: Formatting PCoA table...')

        res_table = pcoaDF.to_html(classes="table display")
        res_table = res_table.replace('border="1"', 'border="0"')
        finalDict['res_table'] = str(res_table)

        functions.setBase(self.RID, 'Step 8 of 9: Formatting PCoA table...done!')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        functions.setBase(self.RID, 'Step 9 of 9: Formatting distance score table...')

        distDF.sort_index(axis=1, inplace=True)
        distDF.sort_index(axis=0, inplace=True)
        dist_table = distDF.to_html(classes="table display")
        dist_table = dist_table.replace('border="1"', 'border="0"')
        finalDict['dist_table'] = str(dist_table)

        functions.setBase(self.RID, 'Step 9 of 9: Formatting distance score table...done!')

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if self.stopList[self.PID] == self.RID:
            if self.debug:
                print "Stopping!"
            return HttpResponse(getStopDict(), content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        finalDict['error'] = 'none'
        res = json.dumps(finalDict)
        return HttpResponse(res, content_type='application/json')

    def run(self):
        if self.debug:
            print "Running PCoA"
        ret = self.validate(sig=False, dist=True)
        if ret == 0:
            ret = self.query(taxmap=False)
            if ret == 0:
                return self.statsGraph()

        if self.debug:
            print "Something went wrong with PCoA"
        return ret


class Caret(Analysis):

    def statsGraph(self):
        try:
            count_rDF = pd.DataFrame()
            if self.DepVar == 0:
                count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='abund')
            elif self.DepVar == 1:
                count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='rel_abund')
            elif self.DepVar == 2:
                count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='rich')
            elif self.DepVar == 3:
                count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='diversity')
            elif self.DepVar == 4:
                count_rDF = self.finalDF.pivot(index='sampleid', columns='rank_id', values='abund_16S')

            count_rDF.fillna(0, inplace=True)
            r = self.initializeR()
            functions.setBase(self.RID, 'Verifying R packages...missing packages are being installed')

            # R packages from cran
            r("list.of.packages <- c('caret', 'randomForest', 'NeuralNetTools', 'e1071', 'stargazer', 'stringr', 'ROCR')")
            r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
            r(
                "if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

            functions.setBase(self.RID, 'Step 4 of 5: Performing statistical test...')

            r('library(caret)')
            r('library(reshape2)')
            r('library(RColorBrewer)')
            r('library(ROCR)')
            r("library(plyr)")
            r('library(stargazer)')
            r('library(stringr)')
            r('source("R/myFunctions/myFunctions.R")')

            method = self.all['Method']
            if method == 'rf':
                r('library(randomForest)')
            elif method == 'nnet':
                r('library(NeuralNetTools)')
            elif method == 'svm':
                r('library(e1071)')

            # Wrangle data into R
            rankNameDF = self.finalDF.drop_duplicates(subset='rank_id', keep='last')
            rankNameDF = rankNameDF.sort_values('rank_id')  # SettingWithCopyWarning when using inplace, fixed by manually setting
            if self.treeType == 3 and self.nzAll >= 5:
                rankNameDF.loc[:, 'name_id'] = rankNameDF['rank_name'].str.split(': ').str[0]
            else:
                rankNameDF.loc[:, 'name_id'] = rankNameDF[['rank_name', 'rank_id']].apply(lambda x: ' id: '.join(x), axis=1)
            r.assign('rankNames', rankNameDF.name_id.values)

            count_rDF.sort_index(axis=0, inplace=True)
            r.assign("treeType", self.treeType)
            r.assign("data", count_rDF)
            r("names(data) <- rankNames")

            myList = list(self.metaDF.select_dtypes(include=['object']).columns)
            for i in myList:
                self.metaDF[i] = self.metaDF[i].str.replace(' ', '_')
                self.metaDF[i] = self.metaDF[i].str.replace('-', '.')
                self.metaDF[i] = self.metaDF[i].str.replace('(', '.')
                self.metaDF[i] = self.metaDF[i].str.replace(')', '.')

            self.metaDF.sort_values('sampleid', inplace=True)
            self.metaDF.set_index('sampleid', inplace=True)
            r.assign("meta_full", self.metaDF)
            r.assign("rows", self.metaDF.index.values.tolist())

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if self.stopList[self.PID] == self.RID:
                if self.debug:
                    print "Stopping!"
                return HttpResponse(getStopDict(), content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            # Predictors
            r("X_full = data")
            r("nzv_cols <- nearZeroVar(X_full)")
            r("if(length(nzv_cols > 0)) X_full <- X_full[,-nzv_cols]")
            r("n.vars <- ncol(data)")

            # Response
            r.assign("allFields", self.allFields)
            r("Y_full = meta_full")

            # Subset train data
            trainIDs = self.all['trainArray']
            r.assign("trainIDs", trainIDs)

            r("X <- X_full[row.names(X_full) %in% trainIDs,]")
            r("meta <- meta_full[row.names(meta_full) %in% trainIDs,]")
            r("Y <- Y_full[row.names(Y_full) %in% trainIDs,]")
            r("Y <- Y[,paste(allFields)]")
            r("myData <- data.frame(Y, X)")

            # Subset test data
            testIDs = list(self.all['testArray'])
            if testIDs:
                r.assign("testIDs", testIDs)
                r("X_test <- X_full[row.names(X_full) %in% testIDs,]")
                r("meta_test <- meta_full[row.names(meta_full) %in% testIDs,]")
                r("Y_test <- Y_full[row.names(Y_full) %in% testIDs,]")
                r("Y_test <- Y_test[,paste(allFields)]")
                r("myData_test <- data.frame(Y_test, X_test)")
                r("nameVec <- c('Y_test', names(X_test))")
                r("nameVec <- make.names(nameVec)")
                r("names(myData_test) <- nameVec")

            # Initialize R output to pdf
            path = 'myPhyloDB/media/temp/rf/Rplots/%s' % self.RID
            if not os.path.exists(path):
                os.makedirs(path)

            r.assign("path", path)
            r.assign("RID", self.RID)
            r("pdf_counter <- 1")

            finalDict = {}

            # set up tuneGrid
            if method == 'rf':
                r("method <- 'rf' ")
                r("title <- 'Random Forest' ")
                r("grid <- expand.grid(.mtry=seq(1, nrow(myData), by=ceiling(nrow(myData)/50) ))")
            elif method == 'nnet':
                r("method <- 'nnet' ")
                r("grid <- expand.grid(.size=seq(1:5), .decay=seq(0, 2, 0.5))")
                r("title <- 'Neural Network' ")
            elif method == 'svm':
                r("method <- 'svmLinear2' ")
                r("grid <- expand.grid(.cost=seq(1:10))")
                r("title <- 'Support Vector Machine' ")

            trainMethod = self.all['trainMethod']
            r.assign("trainMethod", trainMethod)

            number1 = int(self.all['number1'])
            r.assign("number1", number1)
            number2 = int(self.all['number2'])
            r.assign("number2", number2)
            repeats = int(self.all['repeats'])
            r.assign("repeats", repeats)
            proportion = float(self.all['proportion'])
            r.assign("proportion", proportion)

            if trainMethod == 'boot':
                if self.catFields:
                    r("ctrl <- trainControl(method='boot', number=number2, \
                                    classProbs=T, savePredictions=T)")
                if self.quantFields:
                    r("ctrl <- trainControl(method='boot', number=number2, \
                                    classProbs=F, savePredictions=T)")
            elif trainMethod == 'cv':
                if self.catFields:
                    r("ctrl <- trainControl(method='cv', number=number1, \
                                    classProbs=T, savePredictions=T)")
                if self.quantFields:
                    r("ctrl <- trainControl(method='cv', number=number1, \
                                    classProbs=F, savePredictions=T)")
            elif trainMethod == 'repeatedcv':
                if self.catFields:
                    r("ctrl <- trainControl(method='repeatedcv', number=number1, \
                                    repeats=repeats, classProbs=T, savePredictions=T)")
                if self.quantFields:
                    r("ctrl <- trainControl(method='repeatedcv', number=number1, \
                                    repeats=repeats, classProbs=F, savePredictions=T)")
            elif trainMethod == 'LOOCV':
                if self.catFields:
                    r("ctrl <- trainControl(method='LOOCV', number=number1, \
                                    classProbs=T, savePredictions=T)")
                if self.quantFields:
                    r("ctrl <- trainControl(method='LOOCV', number=number1, \
                                    classProbs=F, savePredictions=T)")
            elif trainMethod == 'LGOCV':
                if self.catFields:
                    r("ctrl <- trainControl(method='LGOCV', number=number1, \
                                    p=proportion, classProbs=T, savePredictions=T)")
                if self.quantFields:
                    r("ctrl <- trainControl(method='LGOCV', number=number1, \
                                    p=proportion, classProbs=F, savePredictions=T)")

            if self.catFields:
                r("fit <- train(Y ~ ., data=myData, method=method, linout=F, trace=F, trControl=ctrl, \
                                tuneGrid=grid, importance=T, preProcess=c('center', 'scale'))")
            if self.quantFields:
                r("fit <- train(Y ~ ., data=myData, method=method, linout=T, trace=F, trControl=ctrl, \
                                tuneGrid=grid, importance=T, preProcess=c('center', 'scale'))")

            r("predY <- predict(fit)")

            r("if (exists('predY')) {fitError <- FALSE} else {fitError <- TRUE}")
            fitError = r.get("fitError")
            if fitError:
                myDict = {'error': "Model could not be fit:\nPlease try a different model"}
                res = json.dumps(myDict)
                return HttpResponse(res, content_type='application/json')
            else:
                self.result += str(r('print(fit)')) + '\n'
                self.result += '===============================================\n'

            if self.catFields:
                r("vi <- varImp(fit, scale=F)")
                r("varDF = as.data.frame(vi[[1]])")
                r("goodNameVec <- names(X)")
                r("badNameVec <- names(myData)[2:length(names(myData))]")
                r("row.names(varDF) <- mapvalues(row.names(varDF), from=badNameVec, to=goodNameVec)")

                r("rankDF <- apply(-abs(varDF), 2, rank, ties.method='random')")
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

                # graph probabilities for each test sample
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

                # graph probabilities by taxa
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
                r("p <- p + scale_x_log10()")
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

                self.result += '\nTraining Dataset\n'
                self.result += str(r('print(cm)')) + '\n'
                self.result += '===============================================\n'

                # confusion matrix - test
                if testIDs:
                    r("predY_test <- predict(fit, myData_test)")
                    r("tab <- table(Observed=Y_test, Predicted=predY_test)")
                    r("cm <- confusionMatrix(tab)")

                    self.result += '\nTest Dataset\n'
                    self.result += str(r('print(cm)')) + '\n'
                    self.result += '===============================================\n'

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if self.stopList[self.PID] == self.RID:
                    if self.debug:
                        print "Stopping!"
                    return HttpResponse(getStopDict(), content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            if self.quantFields:
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
                if self.stopList[self.PID] == self.RID:
                    if self.debug:
                        print "Stopping!"
                    return HttpResponse(getStopDict(), content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #


            r("myTable <- stargazer(varDF, type='text', summary=F, rownames=T)")
            self.result += 'Variable Importance\n'
            myString = r.get("myTable")
            self.result += str(myString)
            self.result += '\n===============================================\n'

            if self.catFields:
                r("mergeDF <- cbind(meta, probY)")
                r("myTable <- stargazer(mergeDF, type='text', summary=F, rownames=T)")
                self.result += 'Train Dataset: Probabilities\n'
                myString = r.get("myTable")
                self.result += str(myString)
                self.result += '\n===============================================\n'

                r("mergeDF <- cbind(meta_test, probY_test)")
                r("myTable <- stargazer(mergeDF, type='text', summary=F, rownames=T)")
                self.result += 'Test Dataset: Probabilities\n'
                myString = r.get("myTable")
                self.result += str(myString)
                self.result += '\n===============================================\n'

            if self.quantFields:
                r("mergeDF <- cbind(meta, predY)")
                r("myTable <- stargazer(mergeDF, type='text', summary=F, rownames=T)")
                self.result += 'Train Dataset: Observed vs. Predicted\n'
                myString = r.get("myTable")
                self.result += str(myString)
                self.result += '\n===============================================\n'

                r("mergeDF <- cbind(meta_test, predY_test)")
                r("myTable <- stargazer(mergeDF, type='text', summary=F, rownames=T)")
                self.result += 'Test Dataset: Observed vs. Predicted\n'
                myString = r.get("myTable")
                self.result += str(myString)
                self.result += '\n===============================================\n'

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if self.stopList[self.PID] == self.RID:
                if self.debug:
                    print "Stopping!"
                return HttpResponse(getStopDict(), content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            # Combining Pdf files
            finalFile = 'myPhyloDB/media/temp/rf/Rplots/' + str(self.RID) + '/rf_final.pdf'

            pdf_files = [f for f in os.listdir(path) if f.endswith("pdf")]
            pdf_files = natsorted(pdf_files, key=lambda y: y.lower())

            merger = PdfFileMerger()
            for filename in pdf_files:
                merger.append(PdfFileReader(file(os.path.join(path, filename), 'rb')))

            merger.write(finalFile)

            functions.setBase(self.RID, 'Step 4 of 5: Performing statistical test...done')

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if self.stopList[self.PID] == self.RID:
                if self.debug:
                    print "Stopping!"
                return HttpResponse(getStopDict(), content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            functions.setBase(self.RID, 'Step 5 of 5: Formatting graph data...')
            r("options(width=5000)")
            finalDict['text'] = self.result

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if self.stopList[self.PID] == self.RID:
                if self.debug:
                    print "Stopping!"
                return HttpResponse(getStopDict(), content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            self.deprint("Returning Caret results")

            finalDict['error'] = 'none'
            res = json.dumps(finalDict)
            return HttpResponse(res, content_type='application/json')

        except Exception as e:
            print "Error during Caret:", e.message
            print e
            if not self.stopList[self.PID] == self.RID:
                myDict = {}
                myDict['error'] = "There was an error during your analysis:\nError: " + str(
                    e.message) + "\n"
                res = json.dumps(myDict)
                return HttpResponse(res, content_type='application/json')

        return 0

    def run(self):
        if self.debug:
            print "Running Caret"
        ret = self.validate(sig=False, reqMultiLevel=False)
        if ret == 0:
            ret = self.query(taxmap=False)
            if ret == 0:
                return self.statsGraph()

        if self.debug:
            print "Something went wrong with Caret"
        return ret


class Gage(Analysis):

    def statsGraph(self):
        try:
            RID = self.RID
            all = self.all
            r = self.initializeR()
            functions.setBase(RID, 'Verifying R packages...missing packages are being installed')

            # R packages from biocLite
            r("list.of.packages <- c('gage', 'edgeR', 'pathview')")
            r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
            r("if (length(new.packages)) source('http://bioconductor.org/biocLite.R')")
            r("if (length(new.packages)) biocLite(new.packages, type='source', suppressUpdate=T, dependencies=T)")

            # R packages from cran
            r("list.of.packages <- c('png', 'grid', 'plyr')")
            r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
            r(
                "if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

            functions.setBase(RID, 'Step 3 of 6: Mapping phylotypes to KEGG pathways...')

            r("library(gage)")
            r("library(edgeR)")
            r("library(pathview)")
            r("library(png)")
            r("library(grid)")
            r("library(plyr)")

            keggString = all["kegg"]
            keggDict = json.JSONDecoder(object_pairs_hook=functions.multidict).decode(keggString)
            nameList = []
            for value in keggDict.itervalues():
                if isinstance(value, list):
                    nameList.extend(value)
                else:
                    nameList.append(value)

            # Enable this only if you want to update gage data and pathways
            '''
            r("list.of.packages <- c('gageData')")
            r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
            r("if (length(new.packages)) source('http://bioconductor.org/biocLite.R')")
            print r("if (length(new.packages)) biocLite(new.packages)")

            r('library(gageData)')
            r("data(kegg.sets.ko)")
            r("save(kegg.sets.ko, file='myPhyloDB/media/kegg/kegg.sets.ko.RData')")
            r("for (name in names(kegg.sets.ko)) { \
                id = substr(name, 3, 7); \
                download.kegg(pathway.id=id, species='ko', kegg.dir='myPhyloDB/media/kegg/pathways', file.type=c('xml', 'png')) \
                } \
            ")
            '''

            r("load('myPhyloDB/media/kegg/kegg.sets.ko.RData')")

            keggDict = {}
            r("selPaths <- vector()")
            for i in nameList:
                pathStr = i.split('[PATH:')[1].split(']')[0]
                r.assign("pathStr", pathStr)
                r("selPath <- kegg.sets.ko[grepl(paste(pathStr), names(kegg.sets.ko))]")
                key = r.get("names(selPath)")
                value = r.get("selPath$ko")
                keggDict[key] = value.tolist()
                r("selPaths <- append(selPaths, names(selPath))")

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if self.stopList[self.PID] == RID:
                    if self.debug:
                        print "Stopping!"
                    return HttpResponse(getStopDict(), content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            keggAll = 4
            mapTaxa = 'no'
            finalDF, junk = functions.getKeggDF(keggAll, '', self.savedDF, self.metaDF, self.DepVar, mapTaxa, RID, self.stopList, self.PID)

            # make sure column types are correct
            finalDF[self.catFields] = finalDF[self.catFields].astype(str)

            functions.setBase(RID, 'Step 4 of 6: Performing GAGE analysis...')

            # save location info to session
            myDir = 'myPhyloDB/media/temp/gage/'
            if not os.path.exists(myDir):
                os.makedirs(myDir)

            path = str(myDir) + str(RID) + '.biom'
            functions.imploding_panda(path, 2, self.DepVar, self.finalSampleIDs, self.metaDF, finalDF)

            count_rDF = pd.DataFrame()
            if self.DepVar == 0:
                count_rDF = finalDF.pivot(index='rank_id', columns='sampleid', values='abund')
            elif self.DepVar == 4:
                count_rDF = finalDF.pivot(index='rank_id', columns='sampleid', values='abund_16S')

            # need to change rank_id to kegg orthologies for gage analysis
            count_rDF.reset_index(drop=False, inplace=True)
            count_rDF.rename(columns={'index': 'rank_id'}, inplace=True)
            idList = count_rDF.rank_id.tolist()
            idDict = {}
            for id in idList:
                entry = self.ko_entry.objects.using('picrust').get(ko_lvl4_id=id).ko_orthology
                idDict[id] = entry
            count_rDF['ko'] = count_rDF['rank_id'].map(idDict)
            count_rDF.drop('rank_id', axis=1, inplace=True)
            count_rDF.drop_duplicates(keep='last', inplace=True)  # remove dups - KOs mapped to multiple pathways
            count_rDF.set_index('ko', drop=True, inplace=True)

            # make metaDF R compatible, remove offending characters in categorical variables
            for cat in self.catFields:
                self.metaDF[cat] = self.metaDF[cat].str.replace('-', '_')
                self.metaDF[cat] = self.metaDF[cat].str.replace(' ', '_')
                self.metaDF[cat] = self.metaDF[cat].str.replace('(', '_')
                self.metaDF[cat] = self.metaDF[cat].str.replace(')', '_')

            # Create combined metadata column
            if len(self.catFields) > 1:
                for index, row in self.metaDF.iterrows():
                    self.metaDF.loc[index, 'merge'] = ".".join(row[self.catFields])
            else:
                self.metaDF.loc[:, 'merge'] = self.metaDF.loc[:, self.catFields[0]]

            wantedList = ['merge', 'sampleid', 'sample_name']
            metaDF = self.metaDF.loc[:, wantedList]

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if self.stopList[self.PID] == RID:
                if self.debug:
                    print "Stopping!"
                return HttpResponse(getStopDict(), content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            finalDict = {}
            metaDF.sort_values('sampleid', inplace=True)
            r.assign("metaDF", metaDF)
            r("trt <- factor(metaDF$merge)")

            r.assign("count", count_rDF)
            r.assign("sampleIDs", count_rDF.columns.values.tolist())
            r("names(count) <- sampleIDs")

            r('e <- DGEList(counts=count)')
            r('e <- calcNormFactors(e, method="none")')

            r('design <- model.matrix(~ 0 + trt)')
            r('trtLevels <- levels(trt)')
            r('colnames(design) <- trtLevels')

            r('e <- estimateGLMCommonDisp(e, design)')
            r('e <- estimateGLMTrendedDisp(e, design)')
            r('e <- estimateGLMTagwiseDisp(e, design)')
            r('fit <- glmFit(e, design)')
            fit = r.get('fit')

            if not fit:
                error = "edgeR failed!\nUsually this is caused by one or more taxa having a negative disperion.\nTry filtering your data to remove problematic taxa (e.g. remove phylotypes with 50% or more zeros)."
                myDict = {'error': error}
                res = json.dumps(myDict)
                return HttpResponse(res, content_type='application/json')

            if self.DepVar == 0:
                self.result += 'Dependent Variable: Abundance' + '\n'
            elif self.DepVar == 4:
                self.result += 'Dependent Variable: Total Abundance' + '\n'
            self.result += '\n===============================================\n\n\n'

            levels = list(set(metaDF['merge'].tolist()))
            levels = natsorted(levels, key=lambda y: y.lower())

            r("pdf_counter <- 1")

            path = os.path.join('myPhyloDB', 'media', 'temp', 'gage', 'Rplots', RID)
            if not os.path.exists(path):
                os.makedirs(path)

            r.assign("path", path)
            r("setwd(path)")
            r("options(width=5000)")
            r.assign("RID", RID)

            gageDF = pd.DataFrame(
                columns=['comparison', 'pathway', ' p.geomean ', ' stat.mean ', ' p.val ', ' q.val ', ' set.size '])
            diffDF = pd.DataFrame(
                columns=['comparison', 'kegg', ' baseMean ', ' baseMeanA ', ' baseMeanB ', ' logFC ', ' logCPM ', ' LR ',
                         ' pval ', ' FDR '])

            mergeList = metaDF['merge'].tolist()
            mergeSet = list(set(mergeList))
            for i in xrange(len(levels) - 1):
                for j in xrange(i + 1, len(levels)):
                    trt1 = levels[i]
                    trt2 = levels[j]
                    r.assign("trt1", trt1)
                    r.assign("trt2", trt2)

                    '''
                    Error in makeContrasts(contVec, levels = design) :
                    The levels must by syntactically valid names in R, see help(make.names).  Non-valid names: A-pinene,B-caryophyllene
                    # potential fix on line 177
                    '''

                    r('contVec <- sprintf("%s-%s", trt1, trt2)')
                    r('cont.matrix= makeContrasts(contVec, levels=design)')
                    r('lrt <- glmLRT(fit, contrast=cont.matrix)')
                    r("res <- as.data.frame(topTags(lrt, n=nrow(lrt$table)))")
                    r('res <- res[ order(row.names(res)), ]')
                    r('res')
                    taxaIDs = r.get("row.names(res)")

                    r("change <- -res$logFC")
                    r("names(change) <- row.names(res)")

                    baseMean = count_rDF.mean(axis=1)
                    baseMean = baseMean.loc[baseMean.index.isin(taxaIDs)]

                    listA = metaDF[metaDF['merge'] == mergeSet[i]].sampleid.tolist()
                    baseMeanA = count_rDF[listA].mean(axis=1)
                    baseMeanA = baseMeanA.loc[baseMeanA.index.isin(taxaIDs)]

                    listB = metaDF[metaDF['merge'] == mergeSet[j]].sampleid.tolist()
                    baseMeanB = count_rDF[listB].mean(axis=1)
                    baseMeanB = baseMeanB.loc[baseMeanB.index.isin(taxaIDs)]

                    r.assign("baseMean", baseMean)
                    r.assign("baseMeanA", baseMeanA)
                    r.assign("baseMeanB", baseMeanB)

                    r('baseMean <- baseMean[ order(as.numeric(row.names(baseMean))), ]')
                    r('baseMeanA <- baseMeanA[ order(as.numeric(row.names(baseMeanA))), ]')
                    r('baseMeanB <- baseMeanB[ order(as.numeric(row.names(baseMeanB))), ]')

                    # output DiffAbund to DataTable
                    r("df <- data.frame(kegg=row.names(res), baseMean=baseMean, baseMeanA=baseMeanA, \
                                   baseMeanB=baseMeanB, logFC=-res$logFC, logCPM=res$logCPM, \
                                   LR=res$LR, pval=res$PValue, FDR=res$FDR) \
                               ")

                    nbinom_res = r.get("df")
                    nbinom_res.fillna(value=1.0, inplace=True)

                    if nbinom_res is None:
                        myDict = {'error': "edgeR failed!\nPlease try a different data combination."}
                        res = json.dumps(myDict)
                        return HttpResponse(res, content_type='application/json')

                    comparison = str(trt1) + ' vs. ' + str(trt2)
                    nbinom_res.insert(0, 'comparison', comparison)
                    diffDF = diffDF.append(nbinom_res, ignore_index=True)

                    ### GAGE analysis on all pathways...
                    r("gage.res <- gage(change, gsets=kegg.sets.ko, species='ko', same.dir=FALSE)")
                    r("df2 <- data.frame(pathway=row.names(gage.res$greater), p.geomean=gage.res$greater[, 1], stat.mean=gage.res$greater[, 2], \
                                   p.val=gage.res$greater[, 3], q.val=gage.res$greater[, 4], \
                                   set.size=gage.res$greater[, 5])")

                    compDF = r.get("df2")
                    compDF.insert(0, 'comparison', comparison)
                    gageDF = gageDF.append(compDF, ignore_index=True)

                    ### Get data way for pathview
                    # merge sign and sig to get vector (1=sig. positive, 0=not sig., -1=sig. negative)
                    r("binary <- change / abs(change)")
                    r("sig <- as.vector((res$PValue <= 0.05))")
                    r("sig <- sig * 1")

                    r("sig <- sig * binary")
                    r("names(sig) <- row.names(res)")

                    for key in keggDict.iterkeys():
                        r.assign("pathway", key)
                        r("pid <- substr(pathway, start=1, stop=7)")
                        r("pv <- pathview(gene.data=sig, pathway.id=pid, species='ko', kegg.dir='../../../../kegg/pathways', \
                                       kegg.native=T,  multi.state=F, same.layer=T, low='red', mid='gray', high='green')")

                        # convert to pdf
                        r("pdf(paste('gage_temp', pdf_counter, '.pdf', sep=''))")
                        r("plot.new()")
                        r("pngRaster <- readPNG(paste(pid, 'pathview.png', sep='.'))")
                        r("grid.raster(pngRaster)")
                        r("mtext(paste(trt1, ' vs ', trt2, sep=''), side=3, line=3, col='blue')")
                        r("dev.off()")
                        r("pdf_counter <- pdf_counter + 1")

                    functions.setBase(RID,
                                      'Step 4 of 6: Performing GAGE Analysis...\nComparison: ' + str(trt1) + ' vs ' + str(
                                          trt2) + ' is done!')

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if self.stopList[self.PID] == RID:
                        if self.debug:
                            print "Stopping!"
                        return HttpResponse(getStopDict(), content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if self.stopList[self.PID] == RID:
                    if self.debug:
                        print "Stopping!"
                    return HttpResponse(getStopDict(), content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            functions.setBase(RID, 'Step 5 of 6: Pooling pdf files for display...')

            # Combining Pdf files
            finalFile = 'myPhyloDB/media/temp/gage/Rplots/' + str(RID) + '/gage_final.pdf'
            pdf_files = [f for f in os.listdir(path) if f.endswith("pdf")]
            if pdf_files:
                pdf_files = natsorted(pdf_files, key=lambda y: y.lower())

                merger = PdfFileMerger()
                for filename in pdf_files:
                    merger.append(PdfFileReader(os.path.join(path, filename), 'rb'))

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if self.stopList[self.PID] == RID:
                        if self.debug:
                            print "Stopping!"
                        return HttpResponse(getStopDict(), content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                merger.write(finalFile)

            functions.setBase(RID, 'Step 6 of 6: Formatting result tables...this may take several minutes')
            # Export tables to html
            gage_table = gageDF.to_html(classes="table display")
            gage_table = gage_table.replace('border="1"', 'border="0"')
            finalDict['gage_table'] = str(gage_table)

            diff_table = diffDF.to_html(classes="table display")
            diff_table = diff_table.replace('border="1"', 'border="0"')
            finalDict['diff_table'] = str(diff_table)

            finalDict['text'] = self.result
            finalDict['error'] = 'none'
            res = json.dumps(finalDict)
            return HttpResponse(res, content_type='application/json')

        except Exception as e:
            print "Error during Gage:", e.message
            print e
            if not self.stopList[self.PID] == self.RID:
                myDict = {}
                myDict['error'] = "There was an error during your analysis:\nError: " + str(
                    e.message) + "\n"
                res = json.dumps(myDict)
                return HttpResponse(res, content_type='application/json')

    def run(self):
        if self.debug:
            print "Running Gage"
        ret = self.validate(sig=False, reqMultiLevel=False, selAll=False, metaQuant=False, taxTree=False)
        if ret == 0:
            ret = self.query(taxmap=False, filterable=False)
            if ret == 0:
                return self.statsGraph()

        if self.debug:
            print "Something went wrong with Gage"
        return ret