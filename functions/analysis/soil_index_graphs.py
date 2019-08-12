import datetime
from django.http import HttpResponse
from django_pandas.io import read_frame
import logging
from natsort import natsorted
import numpy as np
import pandas as pd
from pyper import *
from PyPDF2 import PdfFileReader, PdfFileMerger
from scipy.stats import lognorm, norm
import json

from database.models import PICRUSt

import functions
from functions.utils.debug import debug


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getsoil_index(request, stops, RID, PID):
    try:
        while True:
            if request.is_ajax():
                allJson = request.body.split('&')[0]
                all = json.loads(allJson)
                functions.setBase(RID, 'Step 1 of 5: Reading normalized data file...')

                result = ''
                # Select samples and meta-variables from savedDF
                functions.setBase(RID, 'Step 2 of 5: Selecting your chosen meta-variables...')
                metaValsCat = all['metaValsCat']
                metaIDsCat = all['metaIDsCat']
                metaValsQuant = ['soil_texture_sand', 'soil_texture_silt', 'soil_texture_clay',
                        'soil_organic_matter_rating', 'soil_active_C_rating',
                        'soil_water_cap_rating', 'soil_agg_stability_rating',
                        'soil_ACE_protein_index_rating', 'soil_respiration_four_day_rating',
                        'soil_pH_rating', 'soil_p_rating',
                        'soil_k_rating', 'soil_minor_elements_rating',
                        'CASH_SHI_rating']
                metaIDsQuant = []

                # Get meanVals for plotting
                meanVals = all['meanVals']
                sdVals = all['sdVals']
                locMax = all['locMax']

                # Create meta-variable DataFrame, final sample list, final category and quantitative field lists based on tree selections
                savedDF, metaDF, finalSampleIDs, catFields, remCatFields, quantFields, catValues, quantValues = functions.getMetaDF(request.user, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, 10, True)
                avail = metaDF.columns.values
                allFields = list(set(avail) & set(catFields)) + list(set(avail) & set(quantFields))
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

                functions.setBase(RID, 'Step 3 of 5: Selecting your chosen taxa or KEGG level...\n')

                selectAll = 9
                DepVar = 1
                finalDF, missingList = functions.getTaxaDF(selectAll, '', savedDF, metaDF, allFields, DepVar, RID, stops, PID, soilHealth=True)

                nzAll = 10
                DepVar = 1
                mapTaxa = 'no'
                metaDF.set_index('sampleid', drop=True, inplace=True)  # getTaxaDF resets index of metaDF
                debug("SoilHealth: NZDF")
                keggDF, mtDF = functions.getNZDF(nzAll, '', savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)
                # TODO 1.3 check for empties and return http error msg or somesuch
                # TODO 1.3 move soil_index to analysis

                # make sure column types are correct
                #finalDF[catFields] = finalDF[catFields].astype(str)

                #count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='abund_16S')
                #count_rDF.fillna(0, inplace=True)
                #count_rDF.sort_index(axis=0, inplace=True)


                #want = allFields + ['sampleid']
                #metaDF.reset_index(drop=False, inplace=True)
                #bioDF = metaDF[want].copy()
                #print 'bioDF:', bioDF.columns.values

                debug("SoilHealth: CASH")
                # get CASH data
                cashDF = metaDF.sort_values('sampleid')
                for name in metaValsQuant:
                    if name not in cashDF.columns.values:
                        cashDF[name] = 0.0
                for col in cashDF:
                    if cashDF[col].dtype.name == "category":
                        if '0.0' not in cashDF[col].cat.categories:
                            cashDF[col] = cashDF[col].cat.add_categories('0.0')
                cashDF.fillna('0.0', inplace=True)
                cashDF.set_index('sampleid', drop=True, inplace=True)

                debug("SoilHealth: RNA")
                # get rRNA copy numbers data
                otuList = finalDF['rank_id'].tolist()   # TODO 1.3 final error here, after taxaDF fails groupby sum()
                otuList = list(set(otuList))
                qs = PICRUSt.objects.using('picrust').filter(otuid__in=otuList)
                rnaDF = read_frame(qs, fieldnames=['otuid__otuid', 'rRNACount'])
                rnaDF.rename(columns={'otuid__otuid': 'rank_id'}, inplace=True)

                rnaDF = pd.merge(finalDF, rnaDF, left_on='rank_id', right_on='rank_id', how='inner')

                rnaDF = rnaDF.replace(r'\s+', np.nan, regex=True)
                rnaDF.fillna(0.0, axis=1, inplace=True)
                rnaDF['bins'] = 'bin1'
                rnaDF.ix[(rnaDF.rRNACount == 1), 'bins'] = 'bin1'
                rnaDF.ix[(rnaDF.rRNACount == 2), 'bins'] = 'bin2'
                rnaDF.ix[(rnaDF.rRNACount == 3), 'bins'] = 'bin3'
                rnaDF.ix[(rnaDF.rRNACount == 4), 'bins'] = 'bin4'
                rnaDF.ix[(rnaDF.rRNACount == 5), 'bins'] = 'bin5'
                rnaDF.ix[(rnaDF.rRNACount == 6), 'bins'] = 'bin6'
                rnaDF.ix[rnaDF.rRNACount > 6, 'bins'] = 'bin7'
                myBins = list(set(rnaDF['bins'].tolist()))

                binDF = rnaDF.groupby(['sampleid', 'bins'])['rel_abund'].sum()
                binDF = binDF.unstack(level=-1)
                binDF.fillna(0.0, axis=1, inplace=True)
                binDF.reset_index(drop=False, inplace=True)

                binDF['stability'] = (binDF['bin1'] * 7) + (binDF['bin2'] * 6) + (binDF['bin3'] * 5) + \
                    (binDF['bin4'] * 4) + (binDF['bin5'] * 3) + (binDF['bin6'] * 2) + (binDF['bin7'] * 1)
                binDF.set_index('sampleid', drop=True, inplace=True)
                cashDF = pd.merge(cashDF, binDF, left_index=True, right_index=True, how='inner')

                myList = [
                    # GIBBs
                    # C decomposition
                    u'bglX: beta-glucosidase',
                    u'bglB: beta-glucosidase',
                    u'E3.2.1.21: beta-glucosidase',
                    #
                    # N decomposition
                    u'amiE: amidase',
                    u'ureC: urease alpha subunit',
                    #
                    # P decomposition
                    u'phoA: alkaline phosphatase',
                    u'phoD: alkaline phosphatase',
                    u'E3.1.3.2: acid phosphatase',
                    u'appA: acid phosphatase',
                    u'phoN: acid phosphatase',
                    #
                    # S decomposition
                    u'aslA: arylsulfatase',
                    #
                    # N-fixation
                    u'nifH: nitrogenase Fe protein',
                    u'nifDK: nitrogenase Mo-Fe protein',
                    #
                    # P solubility
                    u'pqqC: pyrroloquinoline-quinone synthase',
                    #
                    # Biocontrol
                    u'hcnA: glycine dehydrogenase (cyanide-forming)',
                    u'budA: acetolactate decarboxylase',
                    u'budC: (S,S)-butanediol dehydrogenase',
                    u'E3.2.1.6: endo-1,3(4)-beta-gulcanase',
                    u'E3.2.1.14: chitinase',
                    u'prnD: aminopyrrolnitrin oxygenase',
                    u'phlD: phloroglucinol synthase',
                    u'ituA: iturin family lipopeptide synthetase A',
                    u'fenA: fengycin family lipopeptide synthetase D',
                    u'srfAA: surfactin family lipopeptide synthetase A',
                    u'rifM: AHBA synthesis associated protein',
                    u'phzE: 2-amino-4-deoxychorismate synthase',
                    #
                    # Root growth
                    u'ipdC: indolepyruvate decarboxylase',
                    u'acdS: 1-aminocyclopropane-1-carboxylate deaminase',
                    #
                    # Siderophore
                    u'mbtI: salicylate synthetase',
                    u'entA: 2,3-dihydro-2,3-dihydroxybenzoate dehydrogenase',
                    u'pchB: isochorismate pyruvate lysase',
                    #
                    # N-fixation
                    # nitrogenase -- already in list
                    #
                    # Ammonification
                    # amidase -- already in list
                    # urease -- already in list
                    #
                    # Nitrificiation
                    u'pmoA-amoA: methane/ammonia monooxygenase',
                    u'hao: hydroxylamine dehydrogenase',
                    u'narGH: nitrate reductase',
                    #
                    # DNRA
                    u'nirBD: nitrite reductase (NADH)',
                    u'nrfA: nitrite reductase (cytochrome c-552)',
                    #
                    # ANRA
                    u'nirA: ferrodoxin-nitrite reductase',
                    u'NIT-6: nitrite reductase (NAD(P)H)',
                    #
                    # Denitrification
                    u'nirK: nitrite reductase (NO forming)',
                    u'nirS:  nitrite reductase (NO forming)',
                    u'norBC: nitric oxide reductase',
                    u'nosZ: nitrous-oxide reductase'
                ]

                keggDF['rank_name'] = keggDF['rank_name'].str.split("|").str[-1]
                keggDF = keggDF.pivot(index='sampleid', columns='rank_name', values='rel_abund')
                for i in myList:
                    if i not in keggDF:
                        keggDF[i] = 0.0

                keggDF.rename(columns={
                    # GIBBs
                    # C decomposition
                    u'bglX: beta-glucosidase': 'bglX',
                    u'bglB: beta-glucosidase': 'bglB',
                    u'E3.2.1.21: beta-glucosidase': 'E3.2.1.21',
                    #
                    # N decomposition
                    u'amiE: amidase': 'amiE',
                    u'ureC: urease alpha subunit': 'ureC',
                    #
                    # P decomposition
                    u'phoA: alkaline phosphatase': 'phoA',
                    u'phoD: alkaline phosphatase': 'phoD',
                    u'E3.1.3.2: acid phosphatase': 'E3.1.3.2',
                    u'appA: acid phosphatase': 'appA',
                    u'phoN: acid phosphatase': 'phoN',
                    #
                    # S decomposition
                    u'aslA: arylsulfatase': 'aslA',
                    #
                    # N-fixation
                    u'nifH: nitrogenase Fe protein': 'nifH',
                    u'nifDK: nitrogenase Mo-Fe protein': 'nifDK',
                    #
                    # P solubility
                    u'pqqC: pyrroloquinoline-quinone synthase': 'pqqC',
                    #
                    # Biocontrol
                    u'hcnA: glycine dehydrogenase (cyanide-forming)': 'hcnA',
                    u'budA: acetolactate decarboxylase': 'budA',
                    u'budC: (S,S)-butanediol dehydrogenase': 'budC',
                    u'E3.2.1.6: endo-1,3(4)-beta-gulcanase': 'E3.2.1.6',
                    u'E3.2.1.14: chitinase': 'E3.2.1.14',
                    u'prnD: aminopyrrolnitrin oxygenase': 'prnD',
                    u'phlD: phloroglucinol synthase': 'phlD',
                    u'ituA: iturin family lipopeptide synthetase A': 'ituA',
                    u'fenA: fengycin family lipopeptide synthetase D': 'fenA',
                    u'srfAA: surfactin family lipopeptide synthetase A': 'srfAA',
                    u'rifM: AHBA synthesis associated protein': 'rifM',
                    u'phzE: 2-amino-4-deoxychorismate synthase': 'phzE',
                    #
                    # Root growth
                    u'ipdC: indolepyruvate decarboxylase': 'ipdC',
                    u'acdS: 1-aminocyclopropane-1-carboxylate deaminase': 'acdS',
                    #
                    # Siderophore
                    u'mbtI: salicylate synthetase': 'mbtI',
                    u'entA: 2,3-dihydro-2,3-dihydroxybenzoate dehydrogenase': 'entA',
                    u'pchB: isochorismate pyruvate lysase': 'pchB',
                    #
                    # N-fixation
                    # nitrogenase -- already in list
                    #
                    # Ammonification
                    # amidase -- already in list
                    # urease -- already in list
                    #
                    # Nitrificiation
                    u'pmoA-amoA: methane/ammonia monooxygenase': 'pmoA-amoA',
                    u'hao: hydroxylamine dehydrogenase': 'hao',
                    u'narGH: nitrate reductase': 'narGH',
                    #
                    # DNRA
                    u'nirBD: nitrite reductase (NADH)': 'nirBD',
                    u'nrfA: nitrite reductase (cytochrome c-552)': 'nrfA',
                    #
                    # ANRA
                    u'nirA: ferrodoxin-nitrite reductase': 'nirA',
                    u'NIT-6: nitrite reductase (NAD(P)H)': 'NIT-6',
                    #
                    # Denitrification
                    u'nirK: nitrite reductase (NO forming)': 'nirK',
                    u'nirS:  nitrite reductase (NO forming)': 'nirS',
                    u'norBC: nitric oxide reductase': 'norBC',
                    u'nosZ: nitrous-oxide reductase': 'nosZ'
                }, inplace=True)

                debug("SoilHealth: Merge")
                cashDF = pd.merge(cashDF, keggDF, left_index=True, right_index=True, how='inner')

                geneList = [
                    u'bglX',
                    u'bglB',
                    u'E3.2.1.21',
                    #
                    # N decomposition
                    u'amiE',
                    u'ureC',
                    #
                    # P decomposition
                    u'phoA',
                    u'phoD',
                    u'E3.1.3.2',
                    u'appA',
                    u'phoN',
                    #
                    # S decomposition
                    u'aslA',
                    #
                    # N-fixation
                    u'nifH',
                    u'nifDK',
                    #
                    # P solubility
                    u'pqqC',
                    #
                    # Biocontrol
                    u'hcnA',
                    u'budA',
                    u'budC',
                    u'E3.2.1.6',
                    u'E3.2.1.14',
                    u'prnD',
                    u'phlD',
                    u'ituA',
                    u'fenA',
                    u'srfAA',
                    u'rifM',
                    u'phzE',
                    #
                    # Root growth
                    u'ipdC',
                    u'acdS',
                    #
                    # Siderophore
                    u'mbtI',
                    u'entA',
                    u'pchB',
                    #
                    # N cycle
                    # N-fixation
                    # nitrogenase -- already in list
                    #
                    # Ammonification
                    # amidase -- already in list
                    # urease -- already in list
                    #
                    # Nitrificiation
                    u'pmoA-amoA',
                    u'hao',
                    u'narGH',
                    #
                    # DNRA
                    u'nirBD',
                    u'nrfA',
                    #
                    # ANRA
                    u'nirA',
                    u'NIT-6',
                    #
                    # Denitrification
                    u'nirK',
                    u'nirS',
                    u'norBC',
                    u'nosZ'
                ]
                cashList = ['soil_active_C_rating', 'soil_organic_matter_rating', 'soil_texture_sand',
                            'soil_k_rating', 'CASH_SHI_rating',
                            'soil_pH_rating', 'soil_texture_clay', 'soil_agg_stability_rating',
                            'soil_ACE_protein_index_rating', 'soil_texture_silt', 'soil_p_rating',
                            'soil_minor_elements_rating', 'soil_water_cap_rating', 'soil_respiration_four_day_rating']
                cashDF[cashList] = cashDF[cashList].astype(np.float64)

                cashDF_scale = cashDF.copy()
                for key in geneList:
                    if locMax:
                        meanVal = cashDF_scale.mean()[key]
                        sdVal = cashDF_scale.std()[key]
                    else:
                        meanVal = float(meanVals[key])
                        sdVal = float(sdVals[key])

                    ''' get value from cumulative distribution function'''
                    if key == 'Biomass':
                        cashDF_scale[key] = lognorm.logcdf(cashDF_scale[key], loc=meanVal, scale=sdVal)
                    else:
                        cashDF_scale[key] = norm.cdf(cashDF_scale[key], loc=meanVal, scale=sdVal)

                debug("SoilHealth: CASH Groupby")
                cashDF_means = cashDF.groupby(catFields).mean()
                cashDF_means.reset_index(drop=False, inplace=True)

                if len(catFields) > 1:
                    for index, row in cashDF_means.iterrows():
                        cashDF_means.loc[index, 'Treatment'] = " & ".join(row[catFields])
                else:
                    cashDF_means.loc[:, 'Treatment'] = cashDF_means.loc[:, catFields[0]]
                    cashDF_means.set_index('Treatment', inplace=True)
                cashDF_scale_means = cashDF_scale.groupby(catFields).mean()
                cashDF_scale_means.reset_index(drop=False, inplace=True)
                if len(catFields) > 1:
                    for index, row in cashDF_scale_means.iterrows():
                        cashDF_scale_means.loc[index, 'Treatment'] = " & ".join(row[catFields])
                else:
                    cashDF_scale_means.loc[:, 'Treatment'] = cashDF_scale_means.loc[:, catFields[0]]
                cashDF_scale_means.set_index('Treatment', inplace=True)

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                path = os.path.join('myPhyloDB', 'media', 'temp', 'soil_index', 'Rplots', RID)
                if not os.path.exists(path):
                    os.makedirs(path)

                r.assign("path", path)
                r("setwd(path)")
                r("options('width'=5000)")
                r.assign("RID", RID)

                r.assign('data', cashDF_scale_means)
                rowcount = len(cashDF_scale_means)
                r.assign('rowcount', rowcount)
                r.assign('myBins', myBins)  # load bins

                r("pdf_counter <- 1")
                for row in range(1, rowcount+1):
                    r("pdf(paste('soil_index_temp', pdf_counter, '.pdf', sep=''), height=7.5, width=8)")
                    r('par(omi=c(0.25,0.25,0.75,0.25))')
                    r.assign('off', row)
                    r('curName <- row.names(data)[off]')
                    r('layout(matrix(c(1,2,3,3,4,4), 3, 2, byrow=T), widths=c(3,3), heights=c(3,1,3))')

                    # GIBBs
                    r("GIBBsList <- c( \
                        'bglX', \
                        'bglB', \
                        'E3.2.1.21', \
                        'phoA', \
                        'phoD', \
                        'E3.1.3.2', \
                        'appA', \
                        'phoN', \
                        'aslA', \
                        'pqqC', \
                        'hcnA', \
                        'budA', \
                        'budC', \
                        'prnD', \
                        'phlD', \
                        'ituA', \
                        'fenA', \
                        'srfAA', \
                        'rifM', \
                        'phzE', \
                        'E3.2.1.6', \
                        'E3.2.1.14', \
                        'ipdC', \
                        'acdS', \
                        'mbtI', \
                        'entA', \
                        'pchB' \
                     )")

                    r("legList <- c( \
                        'bglX', \
                        'bglB', \
                        'E3.2.1.21', \
                        'phoA', \
                        'phoD', \
                        'E3.1.3.2', \
                        'appA', \
                        'phoN', \
                        'aslA', \
                        'pqqC', \
                        'hcnA', \
                        'budA', \
                        'budC', \
                        'prnD', \
                        'phlD', \
                        'ituA', \
                        'fenA', \
                        'srfAA', \
                        'rifM', \
                        'phzE', \
                        'E3.2.1.6', \
                        'E3.2.1.14', \
                        'ipdC', \
                        'acdS', \
                        'mbtI', \
                        'entA', \
                        'pchB' \
                    )")

                    r("gibbs <- data[,GIBBsList]")

                    r('stars(matrix(1, ncol=27, nrow=1), \
                        draw.segments=T, col.segments=rep(adjustcolor("white", alpha=1), 27), \
                        scale=F, full=T, labels=NULL, len=1, axes=F, \
                        cex=0.7, mar=c(1,1,2,1), add=F, lty=2, \
                        key.xpd=T, key.loc=c(2.1, 2.1), key.labels=legList)')

                    r('mtext("GIBBs", side=3, line=0, at=0.25, cex=1, outer=T, font=2)')
                    r('mtext("N cycle", side=3, line=0, at=0.75, cex=1, outer=T, font=2)')
                    r('mtext(curName, side=3, line=3.5, cex=1.2, outer=T, font=2)')

                    r("colGIBBs <- c('#3366CC', '#3366CC', '#3366CC', \
                        '#DC3912', '#DC3912', '#DC3912', '#DC3912', '#DC3912', \
                        '#FF9900', \
                        '#109618', \
                        '#990099', '#990099', '#990099', '#990099', '#990099', '#990099', \
                        '#990099', '#990099', '#990099', '#990099', '#990099', '#990099', \
                        '#0099C6', '#0099C6', \
                        '#DD4477', '#DD4477', '#DD4477' \
                    )")

                    r('stars(gibbs[off,], \
                        draw.segments=T, col.segments=colGIBBs, \
                        ncol=1, scale=F, full=T, labels=NULL, \
                        cex=0.5, add=T, lty=1, key.xpd=F)')

                    r('stars(matrix(0.75, ncol=27, nrow=1), \
                        draw.segments=T, col.segments=rep(adjustcolor("white", alpha=0), 27), \
                        scale=F, full=T, labels=NULL, len=1, axes=F, \
                        cex=0.5, add=T, lty=2, key.xpd=F)')

                    r('stars(matrix(0.5, ncol=27, nrow=1), \
                        draw.segments=T, col.segments=rep(adjustcolor("white", alpha=0), 27), \
                        scale=F, full=T, labels=NULL, len=1, axes=F, \
                        cex=0.5, add=T, lty=2, key.xpd=F)')

                    r('stars(matrix(0.25, ncol=27, nrow=1), \
                        draw.segments=T, col.segments=rep(adjustcolor("white", alpha=0), 27), \
                        scale=F, full=T, labels=NULL, len=1, axes=F, \
                        cex=0.5, add=T, lty=2, key.xpd=F)')

                    # N cycle
                    r("NList <- c( \
                        'nifH', \
                        'nifDK', \
                        'amiE', \
                        'ureC', \
                        'pmoA.amoA', \
                        'hao', \
                        'narGH', \
                        'nirBD', \
                        'nrfA', \
                        'nirA', \
                        'NIT.6', \
                        'nirK', \
                        'nirS', \
                        'norBC', \
                        'nosZ' \
                     )")

                    r("legList <- c( \
                        'nifH', \
                        'nifDK', \
                        'amiE', \
                        'ureC', \
                        'pmoA.amoA', \
                        'hao', \
                        'narGH', \
                        'nirBD', \
                        'nrfA', \
                        'nirA', \
                        'NIT.6', \
                        'nirK', \
                        'nirS', \
                        'norBC', \
                        'nosZ' \
                    )")

                    r("ncycle <- data[,NList]")

                    r('stars(matrix(1, ncol=15, nrow=1), \
                        draw.segments=T, col.segments=rep(adjustcolor("white", alpha=0), 15), \
                        scale=F, full=T, labels=NULL, len=1, axes=F, \
                        cex=0.7, mar=c(1,1,2,1), add=F, lty=2, \
                        key.xpd=T, key.loc=c(2.1, 2.1), key.labels=legList)')

                    r("colN <- c('#3366CC', '#3366CC', \
                        '#DC3912', '#DC3912', \
                        '#FF9900', '#FF9900', '#FF9900', \
                        '#109618', '#109618', \
                        '#990099', '#990099', \
                        '#0099C6', '#0099C6', '#0099C6', '#0099C6' \
                    )")

                    r('stars(ncycle[off,], \
                        draw.segments=T, col.segments=colN, \
                        ncol=1, scale=F, full=T, labels=NULL, \
                        cex=0.5, add=T, lty=1, key.xpd=F)')

                    r('stars(matrix(0.75, ncol=15, nrow=1), \
                        draw.segments=T, col.segments=rep(adjustcolor("white", alpha=0), 15), \
                        scale=F, full=T, labels=NULL, len=1, axes=F, \
                        cex=0.5, add=T, lty=2, key.xpd=F)')

                    r('stars(matrix(0.5, ncol=15, nrow=1), \
                        draw.segments=T, col.segments=rep(adjustcolor("white", alpha=0), 15), \
                        scale=F, full=T, labels=NULL, len=1, axes=F, \
                        cex=0.5, add=T, lty=2, key.xpd=F)')

                    r('stars(matrix(0.25, ncol=15, nrow=1), \
                        draw.segments=T, col.segments=rep(adjustcolor("white", alpha=0), 15), \
                        scale=F, full=T, labels=NULL, len=1, axes=F, \
                        cex=0.5, add=T, lty=2, key.xpd=F)')

                    # RNA bins
                    r('vals <- as.matrix(t(data[off, rev(myBins)]))')
                    r('par(mai=c(0.5,3,0,3))')
                    r('bar <- barplot(vals, xlim=c(0,1), ylim=c(0,1), \
                        width=0.2, horiz=T, space=0, xpd=T, \
                        cex.names=0.8, las=2, axes=F, col=rainbow(10), \
                        beside=F, yaxs="i", xaxs="i")')
                    r('axis(1, at=c(0, 0.25, 0.5, 0.75, 1), lwd.ticks=0.2, cex.axis=0.8)')
                    r('title("rRNA Copy Number", line=-1, xpd=T)')

                    # CASH
                    r.assign('vars', metaValsQuant)
                    r('vals <- as.numeric(data[off, vars])')
                    r('par(pin=c(3,0.15*length(vals)))')
                    r('bar <- barplot(vals, space=0, horiz=T, cex.names=0.8, \
                        names.arg=c("Sand", "Silt", "Clay", "SOM", "activeC", "water_cap", "agg_stab", "ACE", "resp", \
                                    "pH", "P", "K", "minor_elements", "SHI"), \
                        xpd=T, xlim=c(0,100), las=2, axes=T, col=seq(1:14), beside=T, yaxs="i", xaxs="i")')
                    r('grid(col="black", lwd=1)')
                    r('title("CASH: Comprehensive Assessment of Soil Health", line=1)')

                    r('dev.off()')
                    r("pdf_counter <- pdf_counter + 1")

                functions.setBase(RID, 'Step 3 of 5: Selecting your chosen taxa...done')

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                functions.setBase(RID, 'Step 4 of 5: Pooling pdf files for display...')

                # Combining Pdf files
                finalFile = 'myPhyloDB/media/temp/soil_index/Rplots/' + str(RID) + '/soil_index_final.pdf'

                pdf_files = [f for f in os.listdir(path) if f.endswith("pdf")]
                pdf_files = natsorted(pdf_files, key=lambda y: y.lower())

                merger = PdfFileMerger()
                for filename in pdf_files:
                    merger.append(PdfFileReader(os.path.join(path, filename), 'rb'))

                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                merger.write(finalFile)
                functions.setBase(RID, 'Step 4 of 5: Creating tables for display...')

                finalDict = {}

                cashDF.reset_index(inplace=True)
                res_table = cashDF.to_html(classes="table display")
                res_table = res_table.replace('border="1"', 'border="0"')
                finalDict['res_table'] = str(res_table)

                cashDF_scale.reset_index(inplace=True)
                scale_table = cashDF_scale.to_html(classes="table display")
                scale_table = scale_table.replace('border="1"', 'border="0"')
                finalDict['scale_table'] = str(scale_table)

                cashDF_means.reset_index(inplace=True)
                mean_table = cashDF_means.to_html(classes="table display")
                mean_table = mean_table.replace('border="1"', 'border="0"')
                finalDict['mean_table'] = str(mean_table)

                cashDF_scale_means.reset_index(inplace=True)
                scale_mean_table = cashDF_scale_means.to_html(classes="table display")
                scale_mean_table = scale_mean_table.replace('border="1"', 'border="0"')
                finalDict['scale_mean_table'] = str(scale_mean_table)

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
