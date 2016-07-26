import datetime
from django.http import HttpResponse
import logging
import numpy as np
import pandas as pd
from pyper import *
import simplejson

from database.utils import multidict, wOdum
from database.utils_kegg import getTaxaDF, getKeggDF, getNZDF
import database.queue


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getPCoA(request, stops, RID, PID):
    try:
        while True:
            if request.is_ajax():
                allJson = request.body.split('&')[0]
                all = simplejson.loads(allJson)
                database.queue.setBase(RID, 'Step 1 of 8: Selecting your chosen meta-variables...')
                myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
                path = str(myDir) + 'usr_norm_data.csv'

                with open(path, 'rb') as f:
                    savedDF = pd.read_csv(f, index_col=0, sep=',')

                selectAll = int(all["selectAll"])
                keggAll = int(all["keggAll"])
                nzAll = int(all["nzAll"])

                distance = int(all["distance"])
                PC1 = int(all["PC1"])
                PC2 = int(all["PC2"])
                test = int(all["test"])
                alpha = float(all["alpha"])
                perms = int(all["perms"])

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

                if not catFields_edit:
                    myDict = {}
                    myDict['error'] = "Selected meta data only has one level.\nPlease select a different variable."
                    res = simplejson.dumps(myDict)
                    return HttpResponse(res, content_type='application/json')

                catSampleIDs = []
                if metaIDsCat:
                    idDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaIDsCat)
                    for key in sorted(idDictCat):
                        catSampleIDs.extend(idDictCat[key])

                metaValsQuant = all['metaValsQuant']
                metaIDsQuant = all['metaIDsQuant']

                metaDictQuant = {}
                quantFields = []
                quantValues = []
                if metaValsQuant:
                    metaDictQuant = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaValsQuant)
                    for key in sorted(metaDictQuant):
                        quantFields.append(key)
                        quantValues.extend(metaDictQuant[key])

                quantSampleIDs = []
                if metaIDsQuant:
                    idDictQuant = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaIDsQuant)
                    for key in sorted(idDictQuant):
                        quantSampleIDs.extend(idDictQuant[key])

                allSampleIDs = catSampleIDs + quantSampleIDs
                allFields = catFields_edit + quantFields

                # Removes samples (rows) that are not in our samplelist
                tempDF = savedDF.loc[savedDF['sampleid'].isin(allSampleIDs)]

                if metaDictCat:
                    for key in metaDictCat:
                        tempDF = tempDF.loc[tempDF[key].isin(metaDictCat[key])]

                if metaDictQuant:
                    for key in metaDictQuant:
                        valueList = [float(x) for x in metaDictQuant[key]]
                        tempDF = tempDF.loc[tempDF[key].isin(valueList)]

                wantedList = allFields + ['sampleid']

                metaDF = tempDF[allFields]

                result += 'Categorical variables selected by user: ' + ", ".join(catFields) + '\n'
                result += 'Categorical variables removed from analysis (contains only 1 level): ' + ", ".join(removed) + '\n'
                result += 'Quantitative variables selected by user: ' + ", ".join(quantFields) + '\n'
                result += '===============================================\n'

                database.queue.setBase(RID, 'Step 1 of 8: Selecting your chosen meta-variables...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 2 of 8: Selecting your chosen taxa or KEGG level...')

                DepVar = 1
                finalDF = pd.DataFrame()
                if button3 == 1:
                    DepVar = int(all["DepVar_taxa"])
                    finalDF, missingList = getTaxaDF('rel_abund', selectAll, '', savedDF, metaDF, allFields, DepVar, RID, stops, PID)
                    result += 'The following PGPRs were not detected: ' + ", ".join(missingList) + '\n'
                    result += '===============================================\n'

                if button3 == 2:
                    DepVar = int(all["DepVar_kegg"])
                    finalDF = getKeggDF('rel_abund', keggAll, '', savedDF, tempDF, allFields, DepVar, RID, stops, PID)

                if button3 == 3:
                    DepVar = int(all["DepVar_nz"])
                    finalDF = getNZDF('rel_abund', nzAll, '', savedDF, tempDF, allFields, DepVar, RID, stops, PID)

                # save location info to session
                myDir = 'myPhyloDB/media/temp/pcoa/'
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

                if metaDictQuant:
                    for key in metaDictQuant:
                        valueList = [float(x) for x in metaDictQuant[key]]
                        meta_rDF = meta_rDF.loc[meta_rDF[key].isin(valueList)]

                wantedList = allFields + ['sampleid', 'sample_name']
                meta_rDF = meta_rDF[wantedList]
                meta_rDF.set_index('sampleid', drop=True, inplace=True)

                database.queue.setBase(RID, 'Step 2 of 8: Selecting your chosen taxa...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 3 of 8: Calculating distance matrix...')

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                r.assign("data", count_rDF)
                r("library(vegan)")

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
                    dists = wOdum(datamtx, alpha)
                    r.assign("dist", dists)
                    r("dist <- as.dist(dist)")

                r("mat <- as.matrix(dist, diag=TRUE, upper=TRUE)")
                mat = r.get("mat")

                rowList = meta_rDF.sample_name.values.tolist()
                distDF = pd.DataFrame(mat, columns=[rowList], index=rowList)

                database.queue.setBase(RID, 'Step 3 of 8: Calculating distance matrix...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 4 of 8: Principal coordinates analysis...')

                anovaFields = []
                for i in metaDictCat:
                    levels = len(set(metaDictCat[i]))
                    if levels > 1:
                        anovaFields.append(i)

                trtLength = len(set(catValues))
                trtString = " + ".join(anovaFields)

                pcoaDF = pd.DataFrame()
                eigDF = pd.DataFrame()
                bigf = ''
                envFit = ''
                if trtLength > 0:
                    r.assign("meta", meta_rDF)
                    pcoa_string = "ord <- capscale(dist ~ " + str(trtString) + ", meta)"
                    r.assign("cmd", pcoa_string)
                    r("eval(parse(text=cmd))")

                    path = "myPhyloDB/media/temp/pcoa/Rplots/" + str(RID) + ".pcoa.pdf"
                    if os.path.exists(path):
                        os.remove(path)

                    if not os.path.exists('myPhyloDB/media/temp/pcoa/Rplots'):
                        os.makedirs('myPhyloDB/media/temp/pcoa/Rplots')

                    file = "pdf('myPhyloDB/media/temp/pcoa/Rplots/" + str(RID) + ".pcoa.pdf', height=4, width=6)"
                    r.assign("cmd", file)
                    r("eval(parse(text=cmd))")

                    r("par(mfrow=c(1,2))")
                    r("layout(matrix(c(1,2), 1), widths=c(4,2), heights=c(4,4))")
                    r("ordiplot(ord, type='n')")

                    ellipseVal = all['ellipseVal']
                    if ellipseVal != 'None':
                        myStr = "cat <- factor(meta$" + str(ellipseVal) + ")"
                        r.assign("cmd", myStr)
                        r("eval(parse(text=cmd))")
                        r("pl <- ordiellipse(ord, cat, kind='sd', conf=0.95, draw='polygon', border='black')")


                    surfVal = all['surfVal']
                    if surfVal != 'None':
                        myStr = "quant <- meta$" + str(surfVal)
                        r.assign("cmd", myStr)
                        r("eval(parse(text=cmd))")
                        r("ordisurf(ord, quant, cex=1, labcex=0.6, add=TRUE)")

                    colorVal = all['colorVal']
                    if colorVal != 'None':
                        myStr = "cat <- factor(meta$" + str(colorVal) + ")"
                        r.assign("cmd", myStr)
                        r("eval(parse(text=cmd))")
                        r("points(ord, display='sites', pch=15, col=cat)")
                        r("par(mar=c(0,0,4,0))")
                        r("plot(0, type='n', xlab='', ylab='', axes=FALSE)")
                        r("legend('topleft', legend=levels(cat), bty='n', pch=15, col=1:length(cat))")
                    else:
                        r("points(ord, display='sites', pch=15)")
                        r("plot(0, type='n', xlab='', ylab='', axes=FALSE)")

                    r("dev.off()")

                    if len(quantFields) > 0:
                        trtString = " + ".join(quantFields)
                        envfit_str = "fit <- envfit(ord ~ " + str(trtString) + ", meta, choices=c(" + str(PC1) + "," + str(PC2) + "))"
                        r.assign("cmd", envfit_str)
                        r("eval(parse(text=cmd))")

                        fit_out = r("fit")
                        fit_out = fit_out.replace('try({fit})', '')
                        fit_out = fit_out.replace('***VECTORS', '')
                        tempStuff = fit_out.split('\n')
                        for part in tempStuff:
                            if part > tempStuff[1]:
                                envFit += part + '\n'

                    r("res <- summary(ord)")
                    r("id <- rownames(meta)")
                    r("pcoa <- data.frame(id, meta, res$sites)")
                    pcoaDF = r.get("pcoa")

                    pcoaDF.rename(columns={'id': 'Sample ID'}, inplace=True)
                    pcoaDF.rename(columns={'sample_name': 'Sample Name'}, inplace=True)

                    r("Stat <- c('Eigenvalue', 'Proportion Explained', 'Cumulative Proportion')")
                    r("eig <- data.frame(Stat, res$cont$importance)")
                    eigDF = r.get("eig")

                    database.queue.setBase(RID, 'Step 4 of 8: Principal coordinates analysis...done!')

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                    database.queue.setBase(RID, 'Step 5 of 8: Performing perMANOVA...')

                    if perms < 10:
                        bigf = 'A minimum of 10 permutations is required...'
                    elif len(catFields_edit) == 0:
                        bigf = 'No categorical variables are available for perMANOVA/betaDisper analysis'
                    elif perms >= 10 and len(catFields_edit) > 0:
                        if test == 1:
                            for i in catFields_edit:
                                factor_string = str(i) + " <- factor(meta$" + str(i) + ")"
                                r.assign("cmd", factor_string)
                                r("eval(parse(text=cmd))")

                                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                                if stops[PID] == RID:
                                    res = ''
                                    return HttpResponse(res, content_type='application/json')
                                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                            r.assign("perms", perms)
                            trtString = " * ".join(catFields_edit)
                            amova_string = "res <- adonis(dist ~ " + str(trtString) + ", perms=perms)"
                            r.assign("cmd", amova_string)
                            r("eval(parse(text=cmd))")

                            res_aov = r("res$aov.tab")

                            tempStuff = res_aov.split('\n')
                            for part in tempStuff:
                                if part != tempStuff[0]:
                                    bigf += part + '\n'
                            database.queue.setBase(RID, 'Step 5 of 8: Performing perMANOVA...done!')

                        elif test == 2:
                            database.queue.setBase(RID, 'Step 4 of 8: Principal coordinates analysis...done!')
                            database.queue.setBase(RID, 'Step 5 of 8: Performing BetaDisper...')

                            for i in catFields_edit:
                                factor_string = str(i) + " <- factor(meta$" + str(i) + ")"
                                r.assign("cmd", factor_string)
                                r("eval(parse(text=cmd))")

                                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                                if stops[PID] == RID:
                                    res = ''
                                    return HttpResponse(res, content_type='application/json')
                                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                            r.assign("perms", perms)
                            for i in catFields_edit:
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
                                database.queue.setBase(RID, 'Step 5 of 8: Performing BetaDisper...done!')

                                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                                if stops[PID] == RID:
                                    res = ''
                                    return HttpResponse(res, content_type='application/json')
                                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                else:
                    state = "Your selected variable(s) only have one treatment level, please select additional data!"
                    myDict = {}
                    myDict['error'] = state
                    res = simplejson.dumps(myDict)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 6 of 8: Formatting graph data for display...')
                finalDict = {}
                seriesList = []
                xAxisDict = {}
                yAxisDict = {}

                CAP1 = PC1 + len(catFields_edit) + len(quantFields) + 1
                CAP2 = PC2 + len(catFields_edit) + len(quantFields) + 1

                colors = [
                    "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                    "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                    "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                    "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                    "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                    "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                    "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                    "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                    "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                    "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
                    "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
                    "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
                    "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58",
                    "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
                    "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
                    "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5",
                    "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
                    "#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01",
                    "#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966",
                    "#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
                    "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#C6DC99", "#203B3C",
                    "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
                    "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
                    "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
                    "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F",
                    "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
                    "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
                    "#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00",
                    "#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
                    "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25",
                    "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"
                ]

                colors_idx = 0
                if catFields_edit:
                    grouped = pcoaDF.groupby(catFields_edit)
                    for name, group in grouped:
                        if len(catFields_edit) > 1:
                            trt = "; ".join(name)
                        else:
                            trt = name

                        dataList = []
                        for index, row in group.iterrows():
                            dataDict = {}
                            dataDict['name'] = row['Sample Name']
                            dataDict['x'] = float(row[CAP1])
                            dataDict['y'] = float(row[CAP2])
                            dataList.append(dataDict)

                        seriesDict = {}
                        seriesDict['name'] = str(trt)
                        seriesDict['data'] = dataList
                        seriesDict['color'] = colors[colors_idx]
                        seriesList.append(seriesDict)
                        colors_idx += 1

                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                        if stops[PID] == RID:
                            res = ''
                            return HttpResponse(res, content_type='application/json')
                        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                xTitle = {}
                xTitle['text'] = str(eigDF.columns.values.tolist()[PC1]) + " (" + str(eigDF.iloc[1][PC1] * 100) + "%)"
                xTitle['style'] = {'color': 'black', 'fontSize': '18px', 'fontWeight': 'bold'}
                xAxisDict['title'] = xTitle

                yTitle = {}
                yTitle['text'] = str(eigDF.columns.values.tolist()[PC2]) + " (" + str(eigDF.iloc[1][PC2] * 100) + "%)"
                yTitle['style'] = {'color': 'black', 'fontSize': '18px', 'fontWeight': 'bold'}
                yAxisDict['title'] = yTitle

                styleDict = {'style': {'color': 'black', 'fontSize': '14px'}}
                xAxisDict['labels'] = styleDict
                yAxisDict['labels'] = styleDict

                finalDict['series'] = seriesList
                finalDict['xAxis'] = xAxisDict
                finalDict['yAxis'] = yAxisDict

                if test == 1:
                    result += 'perMANOVA results:' + '\n'
                if test == 2:
                    result += 'betaDisper results:' + '\n'

                if len(catFields_edit) == 0:
                    result += 'test cannot be run...' + '\n'
                else:
                    bigf = bigf.decode('utf-8')
                    result += bigf + '\n'
                result += '===============================================\n'

                if len(catFields_edit) > 0:
                    result += '\nenvfit results:\n'
                    envFit = envFit.decode('utf-8')
                    result += envFit
                result += '===============================================\n'

                result += '\nEigenvalues\n'
                eigStr = eigDF.to_string()
                result += str(eigStr) + '\n'
                result += '===============================================\n\n\n\n'

                finalDict['text'] = result

                database.queue.setBase(RID, 'Step 6 of 8: Formatting graph data for display...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 7 of 8: Formatting PCoA table...')

                res_table = pcoaDF.to_html(classes="table display")
                res_table = res_table.replace('border="1"', 'border="0"')
                finalDict['res_table'] = str(res_table)

                database.queue.setBase(RID, 'Step 7 of 8: Formatting PCoA table...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                database.queue.setBase(RID, 'Step 8 of 8: Formatting distance score table...')

                distDF.sort_index(axis=1, inplace=True)
                distDF.sort_index(axis=0, inplace=True)
                dist_table = distDF.to_html(classes="table display")
                dist_table = dist_table.replace('border="1"', 'border="0"')
                finalDict['dist_table'] = str(dist_table)

                database.queue.setBase(RID, 'Step 8 of 8: Formatting distance score table...done!')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                finalDict['error'] = 'none'
                res = simplejson.dumps(finalDict)
                return HttpResponse(res, content_type='application/json')

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "Error with PCoA!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')
