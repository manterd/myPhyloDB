import datetime
from django.http import HttpResponse
import logging
import pandas as pd
from pyper import *
import json

import functions


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getSpAC(request, stops, RID, PID):
    try:
        while True:
            if request.is_ajax():
                allJson = request.body.split('&')[0]
                all = json.loads(allJson)
                functions.setBase(RID, 'Step 1 of 5: Reading normalized data file...')

                functions.setBase(RID, 'Step 2 of 5: Selecting your chosen meta-variables...')
                selectAll = int(all["selectAll"])

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

                # Select samples and meta-variables from savedDF
                metaValsCat = all['metaValsCat']
                metaIDsCat = all['metaIDsCat']
                metaValsQuant = ''
                metaIDsQuant = ''

                treeType = 1
                DepVar = 0

                # Create meta-variable DataFrame, final sample list, final category and quantitative field lists based on tree selections
                savedDF, metaDF, finalSampleIDs, catFields, remCatFields, quantFields, catValues, quantValues = functions.getMetaDF(request.user, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, DepVar)
                allFields = catFields + quantFields

                result += 'Categorical variables selected by user: ' + ", ".join(catFields + remCatFields) + '\n'
                result += 'Categorical variables not included in the statistical analysis (contains only 1 level): ' + ", ".join(remCatFields) + '\n'
                result += 'Quantitative variables selected by user: ' + ", ".join(quantFields) + '\n'

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
                DepVar = 0

                finalDF = pd.DataFrame()
                if treeType == 1:
                    if selectAll != 8:
                        filteredDF = functions.filterDF(savedDF, DepVar, selectAll, remUnclass, remZeroes, perZeroes, filterData, filterPer, filterMeth)
                    else:
                        filteredDF = savedDF.copy()

                    finalDF, missingList = functions.getTaxaDF(selectAll, '', filteredDF, metaDF, allFields, DepVar, RID, stops, PID)

                # make sure column types are correct
                finalDF[catFields] = finalDF[catFields].astype(str)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/spac/'
                if not os.path.exists(myDir):
                    os.makedirs(myDir)

                path = str(myDir) + str(RID) + '.biom'
                functions.imploding_panda(path, treeType, DepVar, finalSampleIDs, metaDF, finalDF)

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

                functions.setBase(RID, 'Step 3 of 5: Selecting your chosen taxa...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 4 of 5: Calculating OTU Accumulation Curves...')

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                functions.setBase(RID, 'Verifying R packages...missing packages are being installed')

                r("list.of.packages <- c('vegan', 'ggplot2', 'ggthemes', 'RColorBrewer', 'reshape2')")
                r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
                print r("if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)")

                functions.setBase(RID, 'Step 4 of 5: Calculating OTU Accumulation Curves...')

                r("options(width=5000)")
                print r("library(vegan)")
                print r("library(reshape2)")
                print r("library(ggplot2)")
                print r("library(ggthemes)")
                print r("library(RColorBrewer)")
                print r('source("R/myFunctions/myFunctions.R")')

                r.assign('data', count_rDF)
                metaDF.sort_values('sampleid', axis=0, inplace=True)

                if len(catFields) > 1:
                    for index, row in metaDF.iterrows():
                       metaDF.loc[index, 'merge'] = ".".join(row[catFields])
                else:
                    metaDF.loc[:, 'merge'] = metaDF.loc[:, catFields[0]]

                if allFields:
                    grouped = metaDF.groupby('merge')
                    counter = 0
                    total = len(grouped)

                    skipped = []

                    for name, group in grouped:
                        if len(group) >= 3:
                            r.assign('name', name)
                            IDs = group['sampleid'].tolist()
                            data = count_rDF.ix[IDs]
                            r.assign("data", data)

                            r("x <- poolaccum(data, minsize=2)")

                            r("boot <- summary(x, display='boot')")
                            r("boot <- data.frame(boot$boot)")
                            r("chao <- summary(x, display='chao')")
                            r("chao <- data.frame(chao$chao)")
                            r("jack1 <- summary(x, display='jack1')")
                            r("jack1 <- data.frame(jack1$jack1)")
                            r("jack2 <- summary(x, display='jack2')")
                            r("jack2 <- data.frame(jack2$jack2)")
                            r("S <- summary(x, display='S')")
                            r("S <- data.frame(S$S)")

                            r("df <- data.frame(N.Samples=boot$N, \
                                 Observed=S$S, Observed.SD=S$Std.Dev, \
                                 Chao=chao$Chao, Chao.SD=chao$Std.Dev, \
                                 Jackknife.1=jack1$Jackknife.1, Jackknife.1.SD=jack1$Std.Dev, \
                                 Jackknife.2=jack2$Jackknife.2, Jackknife.2.SD=jack2$Std.Dev, \
                                 Bootstrap=boot$Bootstrap, Bootstrap.SD=boot$Std.Dev \
                            )")
                            r('df$Treatment <- name')

                            if counter == 0:
                                r('gDF <- cbind(df)')
                            else:
                                r('gDF <- rbind(gDF, df)')

                            counter += 1
                            myStr = 'Sample ' + str(counter) + ' out of ' + str(total) + ' is done.'
                            functions.setBase(RID, myStr)

                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                            if stops[PID] == RID:
                                res = ''
                                return HttpResponse(res, content_type='application/json')
                            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                        else:
                            skipped.append(name)
                            finalDict = {"error": "Error: Not enough replication"}
                            res = json.dumps(finalDict)
                            return HttpResponse(res, content_type='application/json')

                    if len(skipped) > 0:
                        result += "\n"+str(len(skipped)) + " groups ignored: \n\t- "
                        for thing in skipped:
                            result += thing + " - "
                    result += "\n\n"

                result += 'Quantitative variables selected by user: ' + ", ".join(quantFields) + '\n'
                result += '===============================================\n\n'

                finalDict = {}
                gDF = r.get('gDF')
                table = gDF.to_html(classes="table display")
                table = table.replace('border="1"', 'border="0"')
                finalDict['gDF'] = str(table)

                # create graph here
                r('gDF.val <- melt(gDF, id.vars=c("N.Samples", "Treatment"), measure.vars=c("Observed", "Chao", \
                    "Jackknife.1", "Jackknife.2", "Bootstrap"))')
                r('gDF.sd <- melt(gDF, id.vars=c("N.Samples", "Treatment"), measure.vars=c("Observed.SD", "Chao.SD", \
                    "Jackknife.1.SD", "Jackknife.2.SD", "Bootstrap.SD"))')
                r('gDF.final <- data.frame(N.Samples=gDF.val$N.Samples, Treatment=gDF.val$Treatment, \
                    variable=gDF.val$variable, value=gDF.val$value, sd=gDF.sd$value)')

                r('p <- ggplot(gDF.final, aes(x=N.Samples, y=value, color=factor(Treatment)))')
                r("p <- p + geom_point() + geom_line()")

                std_error = int(all["std_error"])
                if std_error == 1:
                    r("p <- p + geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0.1)")

                r("p <- p + facet_wrap(~ variable, ncol=3)")
                r("p <- p + theme(strip.text.x=element_text(size=10, colour='blue', angle=0))")

                palette = all['palette']
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

                r('number <- nlevels(as.factor(gDF$Treatment))')
                r('colors <- rep(pal, length.out=number) ')
                r("p <- p + scale_color_manual(values=rep(pal, ceiling(nlevels(as.factor(gDF$Treatment))/nColors)))")

                r("p <- p + theme(legend.title=element_blank())")
                r("p <- p + theme(legend.position='bottom')")
                r("p <- p + ylab('Richness Estimator') + xlab('Number of Samples')")

                path = "myPhyloDB/media/temp/spac/Rplots"
                if not os.path.exists(path):
                    os.makedirs(path)

                r.assign("path", path)
                r.assign("RID", RID)
                r("file <- paste(path, '/', RID, '.spac.pdf', sep='')")
                r("p <- set_panel_size(p, height=unit(2.9, 'in'), width=unit(2.9, 'in'))")
                r("ggsave(filename=file, plot=p, units='in', height=10, width=10)")

                functions.setBase(RID, 'Step 4 of 5: Calculating OTU Accumulation Curves...done!')
                functions.setBase(RID, 'Step 5 of 5: Sending files for display...')

                finalDict['text'] = result
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
