from django.http import HttpResponse
from django.db.models import Sum
import numpy as np
import pandas as pd
import pickle
from pyper import *
from scipy.spatial.distance import *
import simplejson

from database.pcoa.pcoa_DF import PCoAMetaDF, normalizePCoA
from database.models import Sample, Profile
from stats.distance import wOdum
from database.utils import multidict, taxaProfileDF


base = {}
stage = {}
time1 = {}
time2 = {}
TimeDiff = {}


def statusPCoA(request):
    global base, stage, time1, time2, TimeDiff
    if request.is_ajax():
        RID = request.GET["all"]
        time2[RID] = time.time()
        try:
            TimeDiff[RID] = time2[RID] - time1[RID]
        except:
            TimeDiff[RID] = 0
        myDict = {}
        try:
            stage[RID] = str(base[RID]) + '<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
        except:
            stage[RID] = '<br>Analysis has been running for %.1f seconds' % TimeDiff[RID]
        myDict['stage'] = stage[RID]
        json_data = simplejson.dumps(myDict, encoding="Latin-1")
        return HttpResponse(json_data, content_type='application/json')


def removeRIDPCOA(request):
    global base, stage, time1, time2, TimeDiff
    try:
        if request.is_ajax():
            RID = request.GET["all"]
            base.pop(RID, None)
            stage.pop(RID, None)
            time1.pop(RID, None)
            time2.pop(RID, None)
            TimeDiff.pop(RID, None)
            return True
        else:
            return False
    except:
        return False


def getPCoAData(request):
    try:
        global base, time1, TimeDiff
        samples = Sample.objects.all()
        samples.query = pickle.loads(request.session['selected_samples'])
        selected = samples.values_list('sampleid')
        qs1 = Sample.objects.all().filter(sampleid__in=selected)

        if request.is_ajax():
            allJson = request.GET["all"]
            all = simplejson.loads(allJson)

            RID = str(all["RID"])
            time1[RID] = time.time()
            base[RID] = 'Step 1 of 8: Querying database...'

            DepVar = int(all["DepVar"])
            taxaLevel = int(all["taxa"])
            distance = int(all["distance"])
            PC1 = int(all["PC1"])
            PC2 = int(all["PC2"])
            test = int(all["test"])
            alpha = float(all["alpha"])
            perms = int(all["perms"])
            NormMeth = int(all["NormMeth"])
            Iters = int(all["Iters"])
            NormVal = all["NormVal"]
            size = int(all["MinSize"])

            countList = []
            for sample in qs1:
                total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                if total['count__sum'] is not None:
                    countList.append(total['count__sum'])

            minSize = int(min(countList))
            medianSize = int(np.median(np.array(countList)))
            maxSize = int(max(countList))

            if NormVal == "min":
                NormReads = minSize
            elif NormVal == "median":
                NormReads = medianSize
            elif NormVal == "max":
                NormReads = maxSize
            elif NormVal == "none":
                NormReads = -1
            else:
                NormReads = int(all["NormVal"])

            newList = []
            result = ''
            if taxaLevel == 1:
                result = result + 'Taxa level: Kingdom' + '\n'
            elif taxaLevel == 2:
                result = result + 'Taxa level: Phyla' + '\n'
            elif taxaLevel == 3:
                result = result + 'Taxa level: Class' + '\n'
            elif taxaLevel == 4:
                result = result + 'Taxa level: Order' + '\n'
            elif taxaLevel == 5:
                result = result + 'Taxa level: Family' + '\n'
            elif taxaLevel == 6:
                result = result + 'Taxa level: Genus' + '\n'
            elif taxaLevel == 7:
                result = result + 'Taxa level: Species' + '\n'

            if distance == 1:
                result = result + 'Distance score: Manhattan' + '\n'
            elif distance == 2:
                result = result + 'Distance score: Euclidean' + '\n'
            elif distance == 3:
                result = result + 'Distance score: Canberra' + '\n'
            elif distance == 4:
                result = result + 'Distance score: Bray-Curtis' + '\n'
            elif distance == 5:
                result = result + 'Distance score: Kulczynski' + '\n'
            elif distance == 6:
                result = result + 'Distance score: Jaccard' + '\n'
            elif distance == 7:
                result = result + 'Distance score: Gower' + '\n'
            elif distance == 8:
                result = result + 'Distance score: altGower' + '\n'
            elif distance == 9:
                result = result + 'Distance score: Morisita' + '\n'
            elif distance == 10:
                result = result + 'Distance score: Horn' + '\n'
            elif distance == 11:
                result = result + 'Distance score: Mountford' + '\n'
            elif distance == 12:
                result = result + 'Distance score: Binomial' + '\n'
            elif distance == 13:
                result = result + 'Distance score: Chao' + '\n'
            elif distance == 14:
                result = result + 'Distance score: Cao' + '\n'
            elif distance == 15:
                result = result + 'Distance score: wOdum' + '\n'

            metaStrCat = all["metaValsCat"]
            fieldListCat = []
            valueListCat = []
            try:
                metaDictCat = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStrCat)
                for key in sorted(metaDictCat):
                    fieldListCat.append(key)
                    valueListCat.append(metaDictCat[key])
            except:
                placeholder = ''

            metaStrQuant = all["metaValsQuant"]
            fieldListQuant = []
            valueListQuant = []
            try:
                metaDictQuant = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStrQuant)
                for key in sorted(metaDictQuant):
                    fieldListQuant.append(key)
                    valueListQuant.extend(metaDictQuant[key])

            except:
                placeholder = ''

            metaStr = all["metaVals"]
            fieldList = []
            valueList = []
            idDict = {}
            try:
                metaDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(metaStr)
                for key in sorted(metaDict):
                    fieldList.append(key)
                    valueList.append(metaDict[key])

                idStr = all["metaIDs"]
                idDict = simplejson.JSONDecoder(object_pairs_hook=multidict).decode(idStr)

            except:
                placeholder = ''

            if DepVar == 4:
                idList = []
                for i in xrange(len(selected)):
                    idList.append(selected[i][0])
                idDict['rRNA_copies'] = idList

            result += 'Categorical variables selected: ' + ", ".join(fieldListCat) + '\n'
            result += 'Quantitative variables selected: ' + ", ".join(fieldListQuant) + '\n'
            result += '===============================================\n'
            result += '\nData Normalization:\n'

            # Limit reads to max value
            if NormMeth == 1:
                for sample in qs1:
                    total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                    if total['count__sum'] is not None:
                        id = sample.sampleid
                        newList.append(id)

            elif NormMeth == 2 or NormMeth == 3:
                if NormReads > maxSize:
                    NormReads = medianSize
                    result += 'The subsample size was too high and automatically reset to the median value...\n'

                for sample in qs1:
                    total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                    if NormMeth == 2:
                        if total['count__sum'] is not None and int(total['count__sum']) >= NormReads:
                            id = sample.sampleid
                            newList.append(id)
                    else:
                        if total['count__sum'] is not None:
                            id = sample.sampleid
                            newList.append(id)

                # If user set reads too high sample list will be blank
                if not newList:
                    NormReads = medianSize
                    for sample in qs1:
                        total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                        if total['count__sum'] is not None and int(total['count__sum']) >= NormReads:
                            id = sample.sampleid
                            newList.append(id)

            elif NormMeth == 4 or NormMeth == 5 or NormMeth == 6:
                if size > maxSize:
                    size = medianSize
                    result += 'The minimum sample size was too high and automatically reset to the median value...\n'
                for sample in qs1:
                    total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                    if total['count__sum'] is not None and int(total['count__sum']) >= size:
                        id = sample.sampleid
                        newList.append(id)

                # If user set reads too high sample list will be blank
                if not newList:
                    size = medianSize
                    for sample in qs1:
                        total = Profile.objects.filter(sampleid=sample.sampleid).aggregate(Sum('count'))
                        if total['count__sum'] is not None and int(total['count__sum']) >= size:
                            id = sample.sampleid
                            newList.append(id)

            metaDF = PCoAMetaDF(idDict)

            lenA, col = metaDF.shape

            metaDF = metaDF.ix[newList]
            metaDF.dropna(inplace=True)
            lenB, col = metaDF.shape

            selectRem = len(selected) - lenA
            normRem = lenA - lenB

            result += str(lenB) + ' selected samples were included in the final analysis.\n'
            if normRem > 0:
                result += str(normRem) + ' samples did not met the desired normalization criteria.\n'
            if selectRem:
                result += str(selectRem) + ' samples were deselected by the user.\n'

            # Create unique list of samples in meta dataframe (may be different than selected samples)
            myList = metaDF.index.values.tolist()

            taxaDF = taxaProfileDF(myList)

            base[RID] = 'Step 1 of 8: Querying database...done!'

            try:
                base[RID] = 'Step 2 of 8: Normalizing data...'

                # Select only the taxa of interest if user used the selectAll button
                taxaDict = {}
                if taxaLevel == 1:
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('kingdomid', flat='True').distinct()
                    taxaDict['Kingdom'] = qs3
                elif taxaLevel == 2:
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('phylaid', flat='True').distinct()
                    taxaDict['Phyla'] = qs3
                elif taxaLevel == 3:
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('classid', flat='True').distinct()
                    taxaDict['Class'] = qs3
                elif taxaLevel == 4:
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('orderid', flat='True').distinct()
                    taxaDict['Order'] = qs3
                elif taxaLevel == 5:
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('familyid', flat='True').distinct()
                    taxaDict['Family'] = qs3
                elif taxaLevel == 6:
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('genusid', flat='True').distinct()
                    taxaDict['Genus'] = qs3
                elif taxaLevel == 7:
                    qs3 = Profile.objects.all().filter(sampleid__in=myList).values_list('speciesid', flat='True').distinct()
                    taxaDict['Species'] = qs3

                normDF, DESeq_error = normalizePCoA(taxaDF, taxaDict, myList, NormMeth, NormReads, metaDF, Iters)

                finalDict = {}
                if NormMeth == 1:
                    result += 'No normalization was performed...\n'
                elif NormMeth == 2 or NormMeth == 3:
                    result = result + 'Data were rarefied to ' + str(NormReads) + ' sequence reads...\n'
                elif NormMeth == 4:
                    result += 'Data were normalized by the total number of sequence reads...\n'
                elif NormMeth == 5 and DESeq_error == 'no':
                    result += 'Data were normalized by DESeq2...\n'
                elif NormMeth == 5 and DESeq_error == 'yes':
                    result += 'DESeq2 cannot run estimateSizeFactors...\n'
                    result += 'Analysis was run without normalization...\n'
                    result += 'To try again, please select fewer samples or another normalization method...\n'
                elif NormMeth == 6 and DESeq_error == 'no':
                    result += 'Data were normalized by DESeq2 with variance stabilization...\n'
                elif NormMeth == 6 and DESeq_error == 'yes':
                    result += 'DESeq2 cannot run estimateSizeFactors...\n'
                    result += 'Analysis was run without normalization...\n'
                    result += 'To try again, please select fewer samples or another normalization method...\n'
                result += '===============================================\n\n'

                normDF.set_index('sampleid', inplace=True)

                finalDF = pd.merge(metaDF, normDF, left_index=True, right_index=True)
                if DepVar == 4:
                    finalDF['copies'] = finalDF.abund / NormReads * finalDF.rRNA_copies
                    finalDF[['abund', 'copies', 'rich', 'diversity']] = finalDF[['abund', 'copies', 'rich', 'diversity']].astype(float)
                else:
                    finalDF[['abund', 'rich', 'diversity']] = finalDF[['abund', 'rich', 'diversity']].astype(float)

                finalDF.reset_index(inplace=True)

                if DepVar == 1:
                    normDF = finalDF.pivot(index='sampleid', columns='taxa_id', values='abund')
                if DepVar == 2:
                    normDF = finalDF.pivot(index='sampleid', columns='taxa_id', values='rich')
                if DepVar == 3:
                    normDF = finalDF.pivot(index='sampleid', columns='taxa_id', values='diversity')
                if DepVar == 4:
                    normDF = finalDF.pivot(index='sampleid', columns='taxa_id', values='copies')

                base[RID] = 'Step 2 of 8: Normalizing data...done!'

            except:
                myDict = {}
                myDict['error'] = 'Your selections resulted in no valid observations'
                res = simplejson.dumps(myDict)
                return HttpResponse(res, content_type='application/json')

            base[RID] = 'Step 3 of 8: Calculating distance matrix...'

            if os.name == 'nt':
                r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
            else:
                r = R(RCMD="R/R-Linux/bin/R")

            r.assign("data", normDF)
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
                datamtx = np.asarray(normDF)
                dists = wOdum(datamtx, alpha)
                r.assign("dist", dists)

            r("mat <- as.matrix(dist, diag=TRUE, upper=TRUE)")
            mat = r.get("mat")

            rowList = metaDF['sample_name'].tolist()
            distDF = pd.DataFrame(mat, columns=[rowList], index=rowList)

            base[RID] = 'Step 3 of 8: Calculating distance matrix...done!'
            base[RID] = 'Step 4 of 8: Principal coordinates analysis...'

            trtLength = 0
            for i in valueList:
                if len(i) > trtLength:
                    trtLength = len(i)

            if trtLength > 1:
                r.assign("meta", metaDF)
                trtString = " + ".join(fieldListCat)
                pcoa_string = "ord <- capscale(mat ~ " + str(trtString) + ", meta)"
                r.assign("cmd", pcoa_string)
                r("eval(parse(text=cmd))")

                #r("usr_cat1 <- factor(meta$usr_cat1)")
                #r("usr_quant1 <- meta$usr_quant1")

                #r("my_col=c(1:2)[as.factor(usr_cat1)]")

                #user = request.user
                #myStr = "jpeg('media/Rplots/" + str(user) + ".pcoa.jpg')"
                #r.assign("cmd", myStr)
                #r("eval(parse(text=cmd))")

                #r("plot(ord, type='n')")
                #r("points(ord, display='sites', pch=15, col=my_col)")
                #r("pl <- ordiellipse(ord, usr_cat1, kind='sd', conf=0.95, draw='polygon', border='black')")
                #r("ordisurf(ord, usr_quant1, add=TRUE)")
                #r("dev.off()")

                envFit = ''
                if len(fieldListQuant) > 0:
                    trtString = " + ".join(fieldListQuant)
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

                bigf = ''
                if perms <= 10:
                    bigf = 'Please increase the number of permutations...'
                elif len(fieldListCat) == 0:
                    bigf = 'No categorical variables are available for perMANOVA/betaDisper analysis'
                elif perms > 10 and len(fieldListCat) > 0:
                    if test == 1:

                        base[RID] = 'Step 4 of 8: Principal coordinates analysis...done!'
                        base[RID] = 'Step 5 of 8: Performing perMANOVA...'

                        for i in fieldListCat:
                            factor_string = str(i) + " <- factor(meta$" + str(i) + ")"
                            r.assign("cmd", factor_string)
                            r("eval(parse(text=cmd))")

                        r.assign("perms", perms)
                        trtString = " + ".join(fieldListCat)
                        amova_string = "res <- adonis(dist ~ " + str(trtString) + ", perms=perms)"
                        r.assign("cmd", amova_string)
                        r("eval(parse(text=cmd))")

                        res_aov = r("res$aov.tab")

                        tempStuff = res_aov.split('\n')
                        for part in tempStuff:
                            if part != tempStuff[0]:
                                bigf += part + '\n'
                        base[RID] = 'Step 5 of 8: Performing perMANOVA...done!'

                    elif test == 2:
                        base[RID] = 'Step 4 of 8: Principal coordinates analysis...done!'
                        base[RID] = 'Step 5 of 8: Performing BetaDisper...'

                        for i in fieldList:
                            factor_string = str(i) + " <- factor(meta$" + str(i) + ")"
                            r.assign("cmd", factor_string)
                            r("eval(parse(text=cmd))")

                        r.assign("perms", perms)
                        for i in fieldListCat:
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
                            base[RID] = 'Step 5 of 8: Performing BetaDisper...done!'
            else:
                state = "Your selected variable(s) only have one treatment level, please select additional data!"
                myDict = {}
                myDict['error'] = state
                res = simplejson.dumps(myDict)
                return HttpResponse(res, content_type='application/json')

            base[RID] = 'Step 6 of 8: Formatting graph data for display...'

            seriesList = []
            xAxisDict = {}
            yAxisDict = {}

            CAP1 = PC1 + int(len(fieldList)) + 1
            CAP2 = PC2 + int(len(fieldList)) + 1

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
                "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", "#C6DC99", "#203B3C",
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
            if fieldListCat:
                grouped = pcoaDF.groupby(fieldListCat)
                for name, group in grouped:
                    dataList = group.icol([CAP1, CAP2]).values.astype(float).tolist()
                    if len(fieldListCat) > 1:
                        trt = "; ".join(name)
                    else:
                        trt = name
                    seriesDict = {}
                    seriesDict['name'] = str(trt)
                    seriesDict['data'] = dataList
                    seriesDict['color'] = colors[colors_idx]
                    seriesList.append(seriesDict)
                    colors_idx += 1

            xTitle = {}
            xTitle['text'] = str(eigDF.columns.values.tolist()[PC1]) + " (" + str(eigDF.iloc[1][PC1] * 100) + "%)"
            xAxisDict['title'] = xTitle

            yTitle = {}
            yTitle['text'] = str(eigDF.columns.values.tolist()[PC2]) + " (" + str(eigDF.iloc[1][PC2] * 100) + "%)"
            yAxisDict['title'] = yTitle

            finalDict['series'] = seriesList
            finalDict['xAxis'] = xAxisDict
            finalDict['yAxis'] = yAxisDict

            if test == 1:
                result = result + 'perMANOVA results:' + '\n'
            if test == 2:
                result = result + 'betaDisper results:' + '\n'

            if len(fieldListCat) == 0:
                result = result + 'test cannot be run...' + '\n'
            else:
                result = result + str(bigf) + '\n'
            result += '===============================================\n'

            if len(fieldListCat) > 0:
                result = result + '\nenvfit results:\n'
                result = result + str(envFit)
            result += '===============================================\n'

            result = result + '\nEigenvalues\n'
            eigStr = eigDF.to_string()
            result = result + str(eigStr) + '\n'
            result += '===============================================\n\n\n\n'

            finalDict['text'] = result

            base[RID] = 'Step 6 of 8: Formatting graph data for display...done!'
            base[RID] = 'Step 7 of 8: Formatting PCoA table...'

            pcoaDF.reset_index(drop=True, inplace=True)
            res_table = pcoaDF.to_html(classes="table display")
            res_table = res_table.replace('border="1"', 'border="0"')
            finalDict['res_table'] = str(res_table)

            base[RID] = 'Step 7 of 8: Formatting PCoA table...done!'
            base[RID] = 'Step 8 of 8: Formatting distance score table...'

            dist_table = distDF.to_html(classes="table display")
            dist_table = dist_table.replace('border="1"', 'border="0"')
            finalDict['dist_table'] = str(dist_table)

            base[RID] = 'Step 8 of 8: Formatting distance score table...done!'
            finalDict['error'] = 'none'

            res = simplejson.dumps(finalDict)
            return HttpResponse(res, content_type='application/json')

    except Exception as e:
        print "Error with PCoA CAT: ", e
        state = "Error with PCoA CAT: " + str(e)

        myDict = {}
        myDict['error'] = state
        res = simplejson.dumps(myDict)
        return HttpResponse(res, content_type='application/json')
