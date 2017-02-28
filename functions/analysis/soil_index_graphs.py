import datetime
from django.http import HttpResponse
from django_pandas.io import read_frame
import logging
from natsort import natsorted
import numpy as np
import pandas as pd
from pyper import *
from PyPDF2 import PdfFileReader, PdfFileMerger
import json

from database.models import PICRUSt

import functions


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
                metaValsQuant = []
                metaIDsQuant = []

                # Get maxvals for plotting
                maxVals = all['maxVals']
                locMax = all['locMax']

                # Create meta-variable DataFrame, final sample list, final category and quantitative field lists based on tree selections
                savedDF, metaDF, finalSampleIDs, catFields, remCatFields, quantFields, catValues, quantValues = functions.getMetaDF(request.user, metaValsCat, metaIDsCat, metaValsQuant, metaIDsQuant, '')
                allFields = catFields + quantFields

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
                DepVar = 10
                finalDF, missingList = functions.getTaxaDF(selectAll, '', savedDF, metaDF, allFields, DepVar, RID, stops, PID)

                nzAll = 10
                DepVar = 4
                mapTaxa = 'no'
                metaDF.set_index('sampleid', drop=True, inplace=True)  # getTaxaDF resets index of metaDF
                keggDF, mtDF = functions.getNZDF(nzAll, '', savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID)

                # make sure column types are correct
                finalDF[catFields] = finalDF[catFields].astype(str)

                '''
                # save location info to session
                myDir = 'myPhyloDB/media/temp/soil_index/'
                path = str(myDir) + str(RID) + '.pkl'

                # now save file to computer
                if not os.path.exists(myDir):
                    os.makedirs(myDir)
                finalDF.to_pickle(path)
                '''

                # Calculate Good's Coverage Estimate
                singletons = finalDF[finalDF.abund == 1].groupby('sampleid')['abund'].count()
                sDF = pd.DataFrame(singletons)
                total = finalDF.groupby('sampleid')['abund'].sum()
                tDF = pd.DataFrame(total)
                coverage = 1 - sDF / tDF
                cDF = pd.DataFrame(coverage)
                cDF.rename(columns={'abund': 'coverage'}, inplace=True)
                cDF.fillna(1, inplace=True)

                count_rDF = finalDF.pivot(index='sampleid', columns='rank_id', values='abund_16S')
                count_rDF.fillna(0, inplace=True)

                # Add coverage to datagrames -- samples may be missing in cDF (use outer merge)
                count_rDF = pd.merge(count_rDF, cDF, left_index=True, right_index=True, how='outer')

                count_rDF.sort_index(axis=0, inplace=True)
                metaDF.sort('sampleid', inplace=True)

                want = catFields + ['sampleid']
                metaDF.reset_index(drop=False, inplace=True)
                bioDF = metaDF[want].copy()
                chemDF = metaDF[want].copy()
                physDF = metaDF[want].copy()

                # Biological
                if 'microbial_respiration' in metaDF.columns:
                    bioDF['microbial_respiration'] = metaDF['microbial_respiration']
                else:
                    bioDF['microbial_respiration'] = 0.0

                if 'soil_active_c' in metaDF.columns:
                    bioDF['soil_active_c'] = metaDF['soil_active_c']
                else:
                    bioDF['soil_active_c'] = 0.0

                if 'soil_ACE_protein' in metaDF.columns:
                    bioDF['soil_ACE_protein'] = metaDF['soil_ACE_protein']
                else:
                    bioDF['soil_ACE_protein'] = 0.0

                # Chemical
                if 'soil_pH' in metaDF.columns:
                    chemDF['soil_pH'] = metaDF['soil_pH']
                else:
                    chemDF['soil_pH'] = 0.0
                if 'soil_C' in metaDF.columns:
                    chemDF['soil_C'] = metaDF['soil_C']
                else:
                    chemDF['soil_C'] = 0.0
                if 'soil_OM' in metaDF.columns:
                    chemDF['soil_OM'] = metaDF['soil_OM']
                else:
                    chemDF['soil_OM'] = 0.0
                if 'soil_N' in metaDF.columns:
                    chemDF['soil_N'] = metaDF['soil_N']
                else:
                    chemDF['soil_N'] = 0.0
                if 'soil_P' in metaDF.columns:
                    chemDF['soil_P'] = metaDF['soil_P']
                else:
                    chemDF['soil_P'] = 0.0
                if 'soil_K' in metaDF.columns:
                    chemDF['soil_K'] = metaDF['soil_K']
                else:
                    chemDF['soil_K'] = 0.0

                # Physical
                if 'water_content_soil' in metaDF.columns:
                    physDF['water_content_soil'] = metaDF['water_content_soil']
                else:
                    physDF['water_content_soil'] = 0.0
                if 'bulk_density' in metaDF.columns:
                    physDF['bulk_density'] = metaDF['bulk_density']
                else:
                    physDF['bulk_density'] = 0.0
                if 'porosity' in metaDF.columns:
                    physDF['porosity'] = metaDF['porosity']
                else:
                    physDF['porosity'] = 0.0
                if 'soil_EC' in metaDF.columns:
                    physDF['soil_EC'] = metaDF['soil_EC']
                else:
                    physDF['soil_EC'] = 0.0
                if 'soil_surf_hard' in metaDF.columns:
                    physDF['soil_surf_hard'] = metaDF['soil_surf_hard']
                else:
                    physDF['soil_surf_hard'] = 0.0
                if 'soil_subsurf_hard' in metaDF.columns:
                    physDF['soil_subsurf_hard'] = metaDF['soil_subsurf_hard']
                else:
                    physDF['soil_subsurf_hard'] = 0.0
                if 'microbial_biomass_C' in metaDF.columns:
                    physDF['microbial_biomass_C'] = metaDF['microbial_biomass_C']
                else:
                    physDF['microbial_biomass_C'] = 0.0
                if 'microbial_biomass_N' in metaDF.columns:
                    physDF['microbial_biomass_N'] = metaDF['microbial_biomass_N']
                else:
                    physDF['microbial_biomass_N'] = 0.0

                if 'soil_agg_stability' in metaDF.columns:
                    physDF['soil_agg_stability'] = metaDF['soil_agg_stability']
                else:
                    physDF['soil_agg_stability'] = 0.0

                if 'soil_water_cap' in metaDF.columns:
                    physDF['soil_water_cap'] = metaDF['soil_water_cap']
                else:
                    physDF['soil_water_cap'] = 0.0

                # replace NaN's with 0's
                bioDF = bioDF.fillna(0)
                chemDF = chemDF.fillna(0)
                physDF = physDF.fillna(0)

                dataDF = pd.merge(chemDF, bioDF, on=want, how='inner')
                dataDF = pd.merge(dataDF, physDF, on=want, how='inner')
                dataDF.set_index('sampleid', drop=True, inplace=True)

                # get means
                bioDF = bioDF.groupby(catFields).mean()
                chemDF = chemDF.groupby(catFields).mean()
                physDF = physDF.groupby(catFields).mean()

                otuList = finalDF['rank_id'].tolist()
                otuList = list(set(otuList))
                qs = PICRUSt.objects.using('picrust').filter(otuid__in=otuList)
                df = read_frame(qs, fieldnames=['otuid__otuid', 'rRNACount'])
                df.rename(columns={'otuid__otuid': 'rank_id'}, inplace=True)
                finalDF = pd.merge(finalDF, df, left_on='rank_id', right_on='rank_id', how='inner')

                rnaDF = finalDF.replace(r'\s+', np.nan, regex=True)
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

                sumDF1 = binDF.groupby('sampleid')[myBins].sum()
                sumDF1['stability'] = (sumDF1['bin1'] * 1) + (sumDF1['bin2'] * 2) + (sumDF1['bin3'] * 3) + \
                    (sumDF1['bin4'] * 4) + (sumDF1['bin5'] * 5) + (sumDF1['bin6'] * 6) + (sumDF1['bin7'] * 7)

                sumDF2 = finalDF.groupby('sampleid')['abund_16S', 'rich', 'diversity'].sum()

                dataDF = pd.merge(dataDF, sumDF1, left_index=True, right_index=True, how='inner')
                dataDF = pd.merge(dataDF, sumDF2, left_index=True, right_index=True, how='inner')

                metaDF.set_index('sampleid', drop=True, inplace=True)
                wantedList = allFields + ['sample_name']
                metaDF = metaDF[wantedList]

                allDF = pd.merge(sumDF2, sumDF1, left_index=True, right_index=True, how='inner')
                allDF = pd.merge(allDF, metaDF, left_index=True, right_index=True, how='inner')
                allDF = pd.merge(allDF, count_rDF, left_index=True, right_index=True, how='inner')

                wantList = ['abund_16S', 'rich', 'diversity', 'coverage'] + myBins + ['stability']
                bytrt1 = allDF.groupby(catFields)[wantList]
                df1 = bytrt1.mean()

                # merge with other df's before sending to R
                df1 = pd.merge(df1, bioDF, left_index=True, right_index=True, how='inner')
                df1 = pd.merge(df1, chemDF, left_index=True, right_index=True, how='inner')
                df1 = pd.merge(df1, physDF, left_index=True, right_index=True, how='inner')

                df1.reset_index(drop=False, inplace=True)
                if len(catFields) > 1:
                    for index, row in df1.iterrows():
                        df1.loc[index, 'Treatment'] = " & ".join(row[catFields])
                else:
                    df1.loc[:, 'Treatment'] = df1.loc[:, catFields[0]]
                df1.set_index('Treatment', inplace=True)

                myList = [
                    # GIBBs
                    # C decomposition
                    u'3.2.1.4  cellulase',
                    u'3.2.1.91  cellulose 1,4-beta-cellobiosidase (non-reducing end)',
                    u'3.2.1.21  beta-glucosidase',
                    u'3.2.1.8  endo-1,4-beta-xylanase',
                    u'3.2.1.37  xylan 1,4-beta-xylosidase',
                    #
                    # N decomposition
                    u'3.5.1.4  amidase',
                    u'3.5.1.5  urease',
                    #
                    # P decomposition
                    u'3.1.3.1  alkaline phosphatase',
                    u'3.1.3.2  acid phosphatase',
                    #
                    # S decomposition
                    u'3.1.6.1  arylsulfatase',
                    #
                    # N-fixation
                    u'1.18.6.1  nitrogenase',
                    #
                    # P solubility
                    u'1.3.3.11  pyrroloquinoline-quinone synthase',
                    #
                    # Biocontrol
                    u'1.4.99.5  glycine dehydrogenase (cyanide-forming)',
                    u'4.1.1.5  acetolactate decarboxylase',
                    u'1.1.1.76  (S,S)-butanediol dehydrogenase',
                    u'3.2.1.6  endo-1,3(4)-beta-glucanase',
                    u'3.2.1.14  chitinase',
                    #
                    # Root growth
                    u'4.1.1.74  indolepyruvate decarboxylase',
                    u'3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase',
                    #
                    # Siderophore
                    u'6.3.2.39  aerobactin synthase',
                    u'6.3.2.14  enterobactin synthase',
                    u'mycobactin salicyl-AMP ligase [EC:6.3.2.-]',
                    u'yersiniabactin salicyl-AMP ligase [EC:6.3.2.-]',
                    #
                    # N-fixation
                    # nitrogenase -- already in list
                    #
                    # Ammonification
                    # amidase -- already in list
                    # urease -- already in list
                    #
                    # Nitrificiation
                    u'1.14.99.39  ammonia monooxygenase',
                    u'1.7.2.6  hydroxylamine dehydrogenase',
                    u'1.7.99.4  nitrate reductase',
                    #
                    # DNRA
                    # nitrate reductase -- already in list
                    u'1.7.5.1  nitrate reductase (quinone)',
                    u'1.7.1.15  nitrite reductase (NADH)',
                    u'1.7.2.2  nitrite reductase (cytochrome; ammonia-forming)',
                    #
                    # ANRA
                    u'1.7.7.1  ferredoxin---nitrite reductase',
                    u'1.7.7.2  ferredoxin---nitrite reductase',
                    # nitrate reductase -- already in list
                    u'1.7.1.4  nitrite reductase [NAD(P)H]',
                    u'1.7.1.1  nitrate reductase (NADH)',
                    #
                    # Denitrification
                    # nitrate reductase -- already in list
                    u'1.7.2.1  nitrite reductase (NO-forming)',
                    u'1.7.2.5  nitric oxide reductase (cytochrome c)',
                    u'1.7.2.4  nitrous-oxide reductase',
                    #
                    # Anammox
                    u'1.7.2.7  hydrazine synthase subunit',
                    u'1.7.2.8  hydrazine dehydrogenase'
                ]

                # shorten rank names
                keggDF['rank_name'] = keggDF['rank_name'].str.split("|").str[-1]
                starDF = keggDF.pivot(index='sampleid', columns='rank_name', values='abund_16S')

                dataDF = pd.merge(dataDF, starDF, left_index=True, right_index=True, how='inner')
                dataDF.reset_index(drop=False, inplace=True)

                # Add potentially missing columns
                for i in myList:
                    if i not in dataDF:
                        dataDF[i] = 0.0

                bytrt2 = dataDF.groupby(catFields)
                df2 = bytrt2[myList].mean()
                df2.rename(columns={
                    # GIBBs
                    # C decomposition
                    u'3.2.1.4  cellulase': 'cellulase',
                    u'3.2.1.91  cellulose 1,4-beta-cellobiosidase (non-reducing end)': 'cellulose 1,4-beta-cellobiosidase (non-reducing end)',
                    u'3.2.1.21  beta-glucosidase': 'beta-glucosidase',
                    u'3.2.1.8  endo-1,4-beta-xylanase': 'endo-1,4-beta-xylanase',
                    u'3.2.1.37  xylan 1,4-beta-xylosidase': 'xylan 1,4-beta-xylosidase',
                    #
                    # N decomposition
                    u'3.5.1.4  amidase': 'amidase',
                    u'3.5.1.5  urease': 'urease',
                    #
                    # P decomposition
                    u'3.1.3.1  alkaline phosphatase': 'alkaline phosphatase',
                    u'3.1.3.2  acid phosphatase': 'acid phosphatase',
                    #
                    # S decomposition
                    u'3.1.6.1  arylsulfatase': 'arylsulfatase',
                    #
                    # N-fixation
                    u'1.18.6.1  nitrogenase': 'nitrogenase',
                    #
                    # P solubility
                    u'1.3.3.11  pyrroloquinoline-quinone synthase': 'pyrroloquinoline-quinone synthase',
                    #
                    # Biocontrol
                    u'1.4.99.5  glycine dehydrogenase (cyanide-forming)': 'glycine dehydrogenase (cyanide-forming)',
                    u'4.1.1.5  acetolactate decarboxylase': 'acetolactate decarboxylase',
                    u'1.1.1.76  (S,S)-butanediol dehydrogenase': '(S,S)-butanediol dehydrogenase',
                    u'3.2.1.6  endo-1,3(4)-beta-glucanase': 'endo-1,3(4)-beta-glucanase',
                    u'3.2.1.14  chitinase': 'chitinase',
                    #
                    # Root growth
                    u'4.1.1.74  indolepyruvate decarboxylase': 'indolepyruvate decarboxylase',
                    u'3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase': '1-aminocyclopropane-1-carboxylate deaminase',
                    #
                    # Siderophore
                    u'6.3.2.39  aerobactin synthase': 'aerobactin synthase',
                    u'6.3.2.14  enterobactin synthase': 'enterobactin synthase',
                    u'mycobactin salicyl-AMP ligase [EC:6.3.2.-]': 'mycobactin salicyl-AMP ligase',
                    u'yersiniabactin salicyl-AMP ligase [EC:6.3.2.-]': 'yersiniabactin salicyl-AMP ligase',
                    #
                    # N-fixation
                    # nitrogenase -- already in list
                    #
                    # Ammonification
                    # amidase -- already in list
                    # urease -- already in list
                    #
                    # Nitrificiation
                    u'1.14.99.39  ammonia monooxygenase': 'ammonia monooxygenase',
                    u'1.7.2.6  hydroxylamine dehydrogenase': 'hydroxylamine dehydrogenase',
                    u'1.7.99.4  nitrate reductase': 'nitrate reductase',
                    #
                    # DNRA
                    # nitrate reductase -- already in list
                    u'1.7.5.1  nitrate reductase (quinone)': 'nitrate reductase (quinone)',
                    u'1.7.1.15  nitrite reductase (NADH)': 'nitrite reductase (NADH)',
                    u'1.7.2.2  nitrite reductase (cytochrome; ammonia-forming)': 'nitrite reductase (cytochrome; ammonia-forming)',
                    #
                    # ANRA
                    u'1.7.7.1  ferredoxin---nitrite reductase': 'ferredoxin---nitrite reductase 1',
                    u'1.7.7.2  ferredoxin---nitrite reductase': 'ferredoxin---nitrite reductase 2',
                    # nitrate reductase -- already in list
                    u'1.7.1.4  nitrite reductase [NAD(P)H]': 'nitrite reductase [NAD(P)H]',
                    u'1.7.1.1  nitrate reductase (NADH)': 'nitrate reductase (NADH)',
                    #
                    # Denitrification
                    # nitrate reductase -- already in list
                    u'1.7.2.1  nitrite reductase (NO-forming)': 'nitrite reductase (NO-forming)',
                    u'1.7.2.5  nitric oxide reductase (cytochrome c)': 'nitric oxide reductase (cytochrome c)',
                    u'1.7.2.4  nitrous-oxide reductase': 'nitrous-oxide reductase',
                    #
                    # Anammox
                    u'1.7.2.7  hydrazine synthase subunit': 'hydrazine synthase subunit',
                    u'1.7.2.8  hydrazine dehydrogenase': 'hydrazine dehydrogenase'
                }, inplace=True)

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

                df2.reset_index(drop=False, inplace=True)
                if len(catFields) > 1:
                    for index, row in df2.iterrows():
                        df2.loc[index, 'Treatment'] = " & ".join(row[catFields])
                else:
                    df2.loc[:, 'Treatment'] = df2.loc[:, catFields[0]]

                newList = [
                    # GIBBs
                    # C decomposition
                    u'cellulase',
                    u'cellulose 1,4-beta-cellobiosidase (non-reducing end)',
                    u'beta-glucosidase',
                    u'endo-1,4-beta-xylanase',
                    u'xylan 1,4-beta-xylosidase',
                    #
                    # N decomposition
                    u'amidase',
                    u'urease',
                    #
                    # P decomposition
                    u'alkaline phosphatase',
                    u'acid phosphatase',
                    #
                    # S decomposition
                    u'arylsulfatase',
                    #
                    # N-fixation
                    u'nitrogenase',
                    #
                    # P solubility
                    u'pyrroloquinoline-quinone synthase',
                    #
                    # Biocontrol
                    u'glycine dehydrogenase (cyanide-forming)',
                    u'acetolactate decarboxylase',
                    u'(S,S)-butanediol dehydrogenase',
                    u'endo-1,3(4)-beta-glucanase',
                    u'chitinase',
                    #
                    # Root growth
                    u'indolepyruvate decarboxylase',
                    u'1-aminocyclopropane-1-carboxylate deaminase',
                    #
                    # Siderophore
                    u'aerobactin synthase',
                    u'enterobactin synthase',
                    u'mycobactin salicyl-AMP ligase',
                    u'yersiniabactin salicyl-AMP ligase',
                    #
                    # N-fixation
                    # nitrogenase -- already in list
                    #
                    # Ammonification
                    # amidase -- already in list
                    # urease -- already in list
                    #
                    # Nitrificiation
                    u'ammonia monooxygenase',
                    u'hydroxylamine dehydrogenase',
                    u'nitrate reductase',
                    #
                    # DNRA
                    # nitrate reductase -- already in list
                    u'nitrate reductase (quinone)',
                    u'nitrite reductase (NADH)',
                    u'nitrite reductase (cytochrome; ammonia-forming)',
                    #
                    # ANRA
                    u'ferredoxin---nitrite reductase 1',
                    u'ferredoxin---nitrite reductase 2',
                    # nitrate reductase -- already in list
                    u'nitrite reductase [NAD(P)H]',
                    u'nitrate reductase (NADH)',
                    #
                    # Denitrification
                    # nitrate reductase -- already in list
                    u'nitrite reductase (NO-forming)',
                    u'nitric oxide reductase (cytochrome c)',
                    u'nitrous-oxide reductase',
                    #
                    # Anammox
                    u'hydrazine synthase subunit',
                    u'hydrazine dehydrogenase'
                ]

                df2.set_index('Treatment', inplace=True)
                df2 = df2[newList]

                maxDict = {
                    #GIBBs
                    # C decomposition
                    u'cellulase': 1e9,
                    u'cellulose 1,4-beta-cellobiosidase (non-reducing end)': 1e9,
                    u'beta-glucosidase': 1e9,
                    u'endo-1,4-beta-xylanase': 1e9,
                    u'xylan 1,4-beta-xylosidase': 1e9,
                    #
                    # N decomposition
                    u'amidase': 1e9,
                    u'urease': 1e9,
                    #
                    # P decomposition
                    u'alkaline phosphatase': 1e9,
                    u'acid phosphatase': 1e9,
                    #
                    # S decomposition
                    u'arylsulfatase': 1e9,
                    #
                    # N-fixation
                    u'nitrogenase': 1e9,
                    #
                    # P solubility
                    u'pyrroloquinoline-quinone synthase': 1e9,
                    #
                    # Biocontrol
                    u'glycine dehydrogenase (cyanide-forming)': 1e9,
                    u'acetolactate decarboxylase': 1e9,
                    u'(S,S)-butanediol dehydrogenase': 1e9,
                    u'endo-1,3(4)-beta-glucanase': 1e9,
                    u'chitinase': 1e9,
                    #
                    # Root growth
                    u'indolepyruvate decarboxylase': 1e9,
                    u'1-aminocyclopropane-1-carboxylate deaminase': 1e9,
                    #
                    # Siderophore
                    u'aerobactin synthase': 1e9,
                    u'enterobactin synthase': 1e9,
                    u'mycobactin salicyl-AMP ligase': 1e9,
                    u'yersiniabactin salicyl-AMP ligase': 1e9,
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
                    u'ammonia monooxygenase': 1e9,
                    u'hydroxylamine dehydrogenase': 1e9,
                    u'nitrate reductase': 1e9,
                    #
                    # DNRA
                    # nitrate reductase -- already in list
                    u'nitrate reductase (quinone)': 1e9,
                    u'nitrite reductase (NADH)': 1e9,
                    u'nitrite reductase (cytochrome; ammonia-forming)': 1e9,
                    #
                    # ANRA
                    u'ferredoxin---nitrite reductase 1': 1e9,
                    u'ferredoxin---nitrite reductase 2': 1e9,
                    # nitrate reductase -- already in list
                    u'nitrite reductase [NAD(P)H]': 1e9,
                    u'nitrate reductase (NADH)': 1e9,
                    #
                    # Denitrification
                    # nitrate reductase -- already in list
                    u'nitrite reductase (NO-forming)': 1e9,
                    u'nitric oxide reductase (cytochrome c)': 1e9,
                    u'nitrous-oxide reductase': 1e9,
                    #
                    # Anammox
                    u'hydrazine synthase subunit': 1e9,
                    u'hydrazine dehydrogenase': 1e9
                }

                scaleDF = df2.copy()
                for key in maxDict:
                    if locMax:
                        maxVal = max(scaleDF[key])
                    else:
                        maxVal = float(maxVals[key])
                    if maxVal == 0:
                        maxVal = 1
                    scaleDF[key] /= maxVal

                r.assign('data', scaleDF)
                r('trt <- row.names(data)')
                r.assign('enzymes', newList)

                rowcount = len(scaleDF)
                r.assign('rowcount', rowcount)

                r.assign('dat', df1)  # use dat for displaying rounded data
                r('odat <- dat')  # use odat for calculating graphs with more precision

                r('dat$abund_16S <- round(dat$abund_16S, 0)')
                r('dat$rich <- signif(dat$rich, 3)')
                r('dat$diversity <- signif(dat$diversity, 3)')
                r('dat$coverage <- signif(dat$coverage, 3)')

                # load max vals from maxVals dictionary into R
                if locMax:
                    # bio
                    r('maxAbund <- max(odat$abund_16S)')
                    r('if(maxAbund==0) maxAbund=1')
                    r('maxRich <- max(odat$rich)')
                    r('if(maxRich==0) maxRich=1')
                    r('maxDiversity <- max(odat$diversity)')
                    r('if(maxDiversity==0) maxDiversity=1')
                    r('maxACE <- max(odat$soil_ACE_protein)')
                    r('if(maxACE==0) maxACE=1')
                    r('maxActC <- max(odat$soil_active_c)')
                    r('if(maxActC==0) maxActC=1')
                    r('maxResp <- max(odat$microbial_respiration)')
                    r('if(maxResp==0) maxResp=1')

                    # chem
                    r('maxPH <- 14')
                    r('maxC <- max(odat$soil_C)')
                    r('if(maxC==0) maxC=1')
                    r('maxOM <- max(odat$soil_OM)')
                    r('if(maxOM==0) maxOM=1')
                    r('maxP <- max(odat$soil_P)')
                    r('if(maxP==0) maxP=1')
                    r('maxK <- max(odat$soil_K)')
                    r('if(maxK==0) maxK=1')
                    r('maxN <- max(odat$soil_N)')
                    r('if(maxN==0) maxN=1')

                    # phys
                    r('maxWCS <- max(odat$water_content_soil)')
                    r('if(maxWCS==0) maxWCS=1')
                    r('maxBD <- max(odat$bulk_density)')
                    r('if(maxBD==0) maxBD=1')
                    r('maxPORO <- max(odat$porosity)')
                    r('if(maxPORO==0) maxPORO=1')
                    r('maxSH <- max(odat$soil_surf_hard)')
                    r('if(maxSH==0) maxSH=1')
                    r('maxSSH <- max(odat$soil_subsurf_hard)')
                    r('if(maxSSH==0) maxSSH=1')
                    r('maxEC <- max(odat$soil_EC)')
                    r('if(maxEC==0) maxEC=1')
                    r('maxAGGSTAB <- max(odat$soil_agg_stability)')
                    r('if(maxAGGSTAB==0) maxAGGSTAB=1')
                    r('maxWC <- max(odat$soil_water_cap)')
                    r('if(maxWC==0) maxWC=1')
                    r('maxBIOC <- max(odat$microbial_biomass_C)')
                    r('if(maxBIOC==0) maxBIOC=1')
                    r('maxBION <- max(odat$microbial_biomass_N)')
                    r('if(maxBION==0) maxBION=1')

                else:
                    # bio
                    r.assign('maxAbund', float(maxVals['Abundance']))
                    r.assign('maxRich', float(maxVals['Richness']))
                    r.assign('maxDiversity', float(maxVals['Diversity']))
                    r.assign('maxACE', float(maxVals['ACE Soil Protein Index']))
                    r.assign('maxActC', float(maxVals['Active Carbon']))
                    r.assign('maxResp', float(maxVals['Respiration']))
                    # chem
                    r.assign('maxC', float(maxVals['Carbon']))
                    r.assign('maxOM', float(maxVals['Organic Matter']))
                    r.assign('maxPH', float(maxVals['pH']))
                    r.assign('maxP', float(maxVals['Phosphorus']))
                    r.assign('maxK', float(maxVals['Potassium']))
                    r.assign('maxN', float(maxVals['Nitrogen']))
                    # physical
                    r.assign('maxWCS', float(maxVals['Water Content']))
                    r.assign('maxBD', float(maxVals['Bulk Density']))
                    r.assign('maxPORO', float(maxVals['Porosity']))
                    r.assign('maxEC', float(maxVals['Electric Conductivity']))
                    r.assign('maxBIOC', float(maxVals['Microbial Biomass C']))
                    r.assign('maxBION', float(maxVals['Microbial Biomass N']))
                    r.assign('maxWC', float(maxVals['Available Water Capacity']))
                    r.assign('maxSH', float(maxVals['Surface Hardness']))
                    r.assign('maxSSH', float(maxVals['Subsurface Hardness']))
                    r.assign('maxAGGSTAB', float(maxVals['Aggregate Stability']))

                r.assign('myBins', myBins)  # load bins

                # bio
                r('odat$abund_16S <- odat$abund_16S / maxAbund')
                r('odat$rich <- odat$rich / maxRich')
                r('odat$diversity <- odat$diversity / maxDiversity')
                r('odat$soil_ACE_protein <- odat$soil_ACE_protein / maxACE')
                r('odat$soil_active_c <- odat$soil_active_c / maxActC')
                r('odat$microbial_respiration <- odat$microbial_respiration / maxResp')

                # chem
                r('odat$soil_pH <- odat$soil_pH / maxPH')
                r('odat$soil_C <- odat$soil_C / maxC')
                r('odat$soil_OM <- odat$soil_OM / maxOM')
                r('odat$soil_P <- odat$soil_P / maxP')
                r('odat$soil_K <- odat$soil_K / maxK')
                r('odat$soil_N <- odat$soil_N / maxN')

                # phys
                r('odat$water_content_soil <- odat$water_content_soil / maxWCS')
                r('odat$bulk_density <- odat$bulk_density / maxBD')
                r('odat$porosity <- odat$porosity / maxPORO')
                r('odat$soil_surf_hard <- odat$odat$soil_surf_hard / maxSH')
                r('odat$soil_subsurf_hard <- odat$soil_subsurf_hard / maxSSH')
                r('odat$soil_EC <- odat$soil_EC / maxEC')
                r('odat$soil_agg_stability <- odat$soil_agg_stability / maxAGGSTAB')
                r('odat$soil_water_cap <- odat$soil_water_cap / maxWC')
                r('odat$microbial_biomass_C <- odat$microbial_biomass_C / maxBIOC')
                r('odat$microbial_biomass_N <- odat$microbial_biomass_N / maxBION')

                r("pdf_counter <- 1")
                for row in range(1, rowcount+1):
                    r("pdf(paste('soil_index_temp', pdf_counter, '.pdf', sep=''), height=7, width=8)")
                    r('par(omi=c(0.25,0.25,0.25,0.25))')
                    r.assign('off', row)
                    r('curName <- row.names(data)[off]')
                    r('layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5,6,6,6,6,6,6), 3, 6, byrow=T), widths=c(1,1,1,1,1,1), heights=c(3,3,1))')

                    # GIBBs
                    r("GIBBsList <- c( \
                        'cellulase', \
                        'cellulose.1.4.beta.cellobiosidase..non.reducing.end.', \
                        'beta.glucosidase', \
                        'endo.1.4.beta.xylanase', \
                        'xylan.1.4.beta.xylosidase', \
                        'amidase', \
                        'urease', \
                        'alkaline.phosphatase', \
                        'acid.phosphatase', \
                        'arylsulfatase', \
                        'nitrogenase', \
                        'pyrroloquinoline.quinone.synthase', \
                        'glycine.dehydrogenase..cyanide.forming.', \
                        'acetolactate.decarboxylase', \
                        'X.S.S..butanediol.dehydrogenase', \
                        'endo.1.3.4..beta.glucanase', \
                        'chitinase', \
                        'indolepyruvate.decarboxylase', \
                        'X1.aminocyclopropane.1.carboxylate.deaminase', \
                        'aerobactin.synthase', \
                        'enterobactin.synthase', \
                        'mycobactin.salicyl.AMP.ligase', \
                        'yersiniabactin.salicyl.AMP.ligase' \
                    )")

                    r("legList <- c( \
                        'EC3.2.1.4', \
                        'EC3.2.1.91', \
                        'EC3.2.1.21', \
                        'EC3.2.1.8', \
                        'EC3.2.1.37', \
                        'EC3.5.1.4', \
                        'EC3.5.1.5', \
                        'EC3.1.3.1', \
                        'EC3.1.3.2', \
                        'EC1.18.6.1', \
                        'EC1.3.3.11', \
                        'EC1.4.99.5', \
                        'EC4.1.1.5', \
                        'EC1.1.1.76', \
                        'EC3.2.1.6', \
                        'EC3.2.1.14', \
                        'EC4.1.1.74', \
                        'EC3.5.99.7', \
                        'EC6.3.2.39', \
                        'EC6.3.2.14', \
                        'K04787', \
                        'K04783' \
                    )")

                    r("gibbs <- data[,GIBBsList]")

                    r('stars(matrix(1, ncol=23, nrow=1), \
                        draw.segments=T, col.segments=rep(adjustcolor("white", alpha=1), 23), \
                        scale=F, full=T, labels=NULL, len=1, axes=F, \
                        cex=0.7, mar=c(1,1,2,1), add=F, lty=2, \
                        key.vpd=T, key.loc=c(2.1, 2.1), key.labels=legList)')

                    r('mtext("GIBBs", side=3, line=0, at=0.25, cex=1, outer=T, font=2)')
                    r('mtext("N cycle", side=3, line=0, at=0.75, cex=1, outer=T, font=2)')
                    r('mtext(curName, side=2, cex=1.2, outer=T, font=2)')

                    r("colGIBBs <- c('blue', 'blue', 'blue', 'blue', 'blue', \
                        'green', 'green', 'red', 'red', 'turquoise', \
                        'magenta', 'yellow', 'salmon', 'salmon', 'salmon', \
                        'salmon', 'salmon', 'darkgray', 'darkgray', \
                        'brown', 'brown', 'brown', 'brown')")

                    r('stars(gibbs[off,], \
                        draw.segments=T, col.segments=colGIBBs, \
                        ncol=1, scale=F, full=T, labels=NULL, \
                        cex=0.5, add=T, lty=1, key.xpd=F)')

                    r('stars(matrix(0.75, ncol=23, nrow=1), \
                        draw.segments=T, col.segments=rep(adjustcolor("white", alpha=0), 23), \
                        scale=F, full=T, labels=NULL, len=1, axes=F, \
                        cex=0.5, add=T, lty=2, key.vpd=F)')

                    r('stars(matrix(0.5, ncol=23, nrow=1), \
                        draw.segments=T, col.segments=rep(adjustcolor("white", alpha=0), 23), \
                        scale=F, full=T, labels=NULL, len=1, axes=F, \
                        cex=0.5, add=T, lty=2, key.vpd=F)')

                    r('stars(matrix(0.25, ncol=23, nrow=1), \
                        draw.segments=T, col.segments=rep(adjustcolor("white", alpha=0), 23), \
                        scale=F, full=T, labels=NULL, len=1, axes=F, \
                        cex=0.5, add=T, lty=2, key.vpd=F)')

                    # N cycle
                    r("NList <- c( \
                        'nitrogenase', \
                        'amidase', \
                        'urease', \
                        'ammonia.monooxygenase', \
                        'hydroxylamine.dehydrogenase', \
                        'nitrate.reductase', \
                        'nitrate.reductase', \
                        'nitrate.reductase..quinone.', \
                        'nitrite.reductase..NADH.', \
                        'nitrite.reductase..cytochrome..ammonia.forming.', \
                        'ferredoxin...nitrite.reductase.1', \
                        'ferredoxin...nitrite.reductase.2', \
                        'nitrate.reductase', \
                        'nitrite.reductase..NAD.P.H.', \
                        'nitrate.reductase..NADH.', \
                        'nitrate.reductase', \
                        'nitrite.reductase..NO.forming.', \
                        'nitric.oxide.reductase..cytochrome.c.', \
                        'nitrous.oxide.reductase', \
                        'hydrazine.synthase.subunit', \
                        'hydrazine.dehydrogenase' \
                     )")

                    r("legList <- c( \
                        'EC1.18.6.1', \
                        'EC3.5.1.4', \
                        'EC3.5.1.5', \
                        'EC1.14.99.39', \
                        'EC1.7.2.6', \
                        'EC1.7.99.4', \
                        'EC1.7.99.4', \
                        'EC1.7.5.1', \
                        'EC1.7.1.15', \
                        'EC1.7.2.2', \
                        'EC1.7.7.1', \
                        'EC1.7.7.2', \
                        'EC1.7.99.4', \
                        'EC1.7.1.4', \
                        'EC1.7.1.1', \
                        'EC1.7.99.4', \
                        'EC1.7.2.1', \
                        'EC1.7.2.5', \
                        'EC1.7.2.4', \
                        'EC1.7.2.7', \
                        'EC1.7.2.8' \
                    )")

                    r("ncycle <- data[,NList]")

                    r('stars(matrix(1, ncol=21, nrow=1), \
                        draw.segments=T, col.segments=rep(adjustcolor("white", alpha=0), 21), \
                        scale=F, full=T, labels=NULL, len=1, axes=F, \
                        cex=0.7, mar=c(1,1,2,1), add=F, lty=2, \
                        key.vpd=T, key.loc=c(2.1, 2.1), key.labels=legList)')

                    r("colN <- c('blue', 'green', 'green', \
                        'red', 'red', 'red', 'turquoise', 'turquoise', \
                        'turquoise', 'turquoise', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', \
                        'yellow', 'yellow', 'yellow', 'yellow', 'salmon', 'salmon')")

                    r('stars(ncycle[off,], \
                        draw.segments=T, col.segments=colN, \
                        ncol=1, scale=F, full=T, labels=NULL, \
                        cex=0.5, add=T, lty=1, key.xpd=F)')

                    r('stars(matrix(0.75, ncol=21, nrow=1), \
                        draw.segments=T, col.segments=rep(adjustcolor("white", alpha=0), 21), \
                        scale=F, full=T, labels=NULL, len=1, axes=F, \
                        cex=0.5, add=T, lty=2, key.vpd=F)')

                    r('stars(matrix(0.5, ncol=21, nrow=1), \
                        draw.segments=T, col.segments=rep(adjustcolor("white", alpha=0), 21), \
                        scale=F, full=T, labels=NULL, len=1, axes=F, \
                        cex=0.5, add=T, lty=2, key.vpd=F)')

                    r('stars(matrix(0.25, ncol=21, nrow=1), \
                        draw.segments=T, col.segments=rep(adjustcolor("white", alpha=0), 21), \
                        scale=F, full=T, labels=NULL, len=1, axes=F, \
                        cex=0.5, add=T, lty=2, key.vpd=F)')

                    # Biological
                    r('vals <- as.numeric(odat[off, c("microbial_respiration", "soil_active_c", "soil_ACE_protein", "coverage", "diversity", "rich", "abund_16S")])')
                    r('par(pin=c(1,0.15*length(vals)))')
                    r('bar <- barplot(vals, space=0, horiz=T, cex.names=0.8, \
                        names.arg=c("Respiration", "Active Carbon", "ACE Soil Protein Index", "Coverage", "Diversity", "Richness", "Abundance"),\
                        xpd=T, xlim=c(0,1), las=2, axes=T, col=seq(1:10), beside=T, yaxs="i", xaxs="i")')
                    r('grid(col="black", lwd=1)')
                    r('values <- c(dat[off, c("microbial_respiration", "soil_active_c", "soil_ACE_protein", "coverage", "diversity", "rich", "abund_16S")])')
                    r('text(x=1, y=bar, labels=values, pos=4, cex=0.8, xpd=TRUE)')
                    r('title("Biological", line=1)')

                    # Chemical
                    r('vals <- as.numeric(odat[off, c("soil_N", "soil_K", "soil_P", "soil_pH", "soil_C", "soil_OM")])')
                    r('par(pin=c(1,0.15*length(vals)))')
                    r('bar <- barplot(vals, space=0, horiz=T, cex.names=0.8, \
                        names.arg=c("Nitrogen", "Potassium", "Phosphorus", "pH", "Carbon", "Organic Matter"),\
                        xpd=T, xlim=c(0,1), las=2, axes=T, col=seq(1:10), beside=T, yaxs="i", xaxs="i")')
                    r('grid(col="black", lwd=1)')
                    r('values <- dat[off, c("soil_N", "soil_K", "soil_P", "soil_pH", "soil_C", "soil_OM")]')
                    r('text(x=1, y=bar, labels=signif(values,3), pos=4, cex=0.8, xpd=TRUE)')
                    r('title("Chemical", line=1)')

                    # Physical
                    r('vals <- as.numeric(odat[off, c("soil_agg_stability", "soil_water_cap", "microbial_biomass_N", "microbial_biomass_C", "soil_EC", "soil_subsurf_hard", "soil_surf_hard", "water_content_soil", "bulk_density", "porosity")])')
                    r('par(pin=c(1,0.15*length(vals)))')
                    r('bar <- barplot(vals, space=0, horiz=T, cex.names=0.8, \
                        names.arg=c("Aggregate Stability", "Available Water Capacity", "Microbial Biomass N", "Microbial Biomass C", "Electric Conductivity", "Subsurface Hardness", "Surface Hardness", "Water Content", "Bulk Density", "Porosity"),\
                        xpd=T, xlim=c(0,1), las=2, axes=T, col=seq(1:10), beside=T, yaxs="i", xaxs="i")')
                    r('grid(col="black", lwd=1)')
                    r('values <- dat[off, c("soil_agg_stability", "soil_water_cap", "microbial_biomass_N", "microbial_biomass_C", "soil_EC", "soil_subsurf_hard", "soil_surf_hard", "water_content_soil", "bulk_density", "porosity")]')
                    r('text(x=1, y=bar, labels=values, pos=4, cex=0.8, xpd=TRUE)')
                    r('title("Physical", line=1)')

                    # RNA bins
                    r('vals <- as.matrix(t(odat[off, rev(myBins)]))')
                    r('par(mai=c(0.5,3,0,3))')
                    r('bar <- barplot(vals, xlim=c(0,1), ylim=c(0,1), \
                        width=0.2, horiz=T, space=0, xpd=T, \
                        cex.names=0.8, las=2, axes=F, col=rainbow(10), \
                        beside=F, yaxs="i", xaxs="i")')
                    r('axis(1, at=c(0, 0.25, 0.5, 0.75, 1), lwd.ticks=0.2, cex.axis=0.8)')
                    #r('text(x=c(0.2, 0.5, 0.8), y=bar+0.3, labels=round(vals, 3), cex=0.8, xpd=T)')
                    r('title("rRNA Copy Number", line=-0.8, xpd=T)')
                    #r('legend(1.1, 0.5, c("x <= 2", "2 < x < 4", "x >= 4"), cex=0.8, bty="n", fill=c("red", "darkgray", "green"), xpd=T)')

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

                res_table = dataDF.to_html(classes="table display")
                res_table = res_table.replace('border="1"', 'border="0"')
                finalDict['res_table'] = str(res_table)

                mean1 = r.get('odat')
                mean2 = r.get('data')

                mean2.rename(columns={
                    ' cellulase ': '3.2.1.4  cellulase',
                    ' cellulose.1.4.beta.cellobiosidase..non.reducing.end. ': '3.2.1.91  cellulose 1,4-beta-cellobiosidase (non-reducing end)',
                    ' beta.glucosidase ': '3.2.1.21  beta-glucosidase',
                    ' endo.1.4.beta.xylanase ': '3.2.1.8  endo-1,4-beta-xylanase',
                    ' xylan.1.4.beta.xylosidase ': '3.2.1.37  xylan 1,4-beta-xylosidase',
                    ' amidase ': '3.5.1.4  amidase',
                    ' urease ': '3.5.1.5  urease',
                    ' alkaline.phosphatase ': '3.1.3.1  alkaline phosphatase',
                    ' acid.phosphatase ': '3.1.3.2  acid phosphatase',
                    ' arylsulfatase ': '3.1.6.1  arylsulfatase',
                    ' nitrogenase ': '1.18.6.1  nitrogenase',
                    ' pyrroloquinoline.quinone.synthase ': '1.3.3.11  pyrroloquinoline-quinone synthase',
                    ' glycine.dehydrogenase..cyanide.forming. ': '1.4.99.5  glycine dehydrogenase (cyanide-forming)',
                    ' acetolactate.decarboxylase ': '4.1.1.5  acetolactate decarboxylase',
                    ' X.S.S..butanediol.dehydrogenase ': '1.1.1.76  (S,S)-butanediol dehydrogenase',
                    ' endo.1.3.4..beta.glucanase': '3.2.1.6  endo-1,3(4)-beta-glucanase',
                    ' chitinase ': '3.2.1.14  chitinase',
                    ' indolepyruvate.decarboxylase ': '4.1.1.74  indolepyruvate decarboxylase',
                    ' X1.aminocyclopropane.1.carboxylate.deaminase ': '3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase',
                    ' aerobactin.synthase ': '6.3.2.39  aerobactin synthase',
                    ' enterobactin.synthase ': '6.3.2.14  enterobactin synthase',
                    ' mycobactin.salicyl.AMP.ligase ': 'mycobactin salicyl-AMP ligase [EC:6.3.2.-]',
                    ' yersiniabactin.salicyl.AMP.ligase ': 'yersiniabactin salicyl-AMP ligase [EC:6.3.2.-]',
                    ' ammonia.monooxygenase ': '1.14.99.39  ammonia monooxygenase',
                    ' hydroxylamine.dehydrogenase ': '1.7.2.6  hydroxylamine dehydrogenase',
                    ' nitrate.reductase ': '1.7.99.4  nitrate reductase',
                    ' nitrate.reductase..quinone. ': '1.7.5.1  nitrate reductase (quinone)',
                    ' nitrite.reductase..NADH. ': '1.7.1.15  nitrite reductase (NADH)',
                    ' nitrite.reductase..cytochrome..ammonia.forming. ': '1.7.2.2  nitrite reductase (cytochrome; ammonia-forming)',
                    ' ferredoxin...nitrite.reductase.1 ': '1.7.7.1  ferredoxin---nitrite reductase',
                    ' ferredoxin...nitrite.reductase.2 ': '1.7.7.2  ferredoxin---nitrite reductase',
                    ' nitrite.reductase..NAD.P.H. ': '1.7.1.4  nitrite reductase [NAD(P)H]',
                    ' nitrate.reductase..NADH. ': '1.7.1.1  nitrate reductase (NADH)',
                    ' nitrite.reductase..NO.forming. ': '1.7.2.1  nitrite reductase (NO-forming)',
                    ' nitric.oxide.reductase..cytochrome.c. ': '1.7.2.5  nitric oxide reductase (cytochrome c)',
                    ' nitrous.oxide.reductase ': '1.7.2.4  nitrous-oxide reductase',
                    ' hydrazine.synthase.subunit ': '1.7.2.7  hydrazine synthase subunit',
                    ' hydrazine.dehydrogenase ': '1.7.2.8  hydrazine dehydrogenase'
                    },
                inplace=True)

                scaleDF = pd.merge(mean1, mean2, left_index=True, right_index=True)
                scale_table = scaleDF.to_html(classes="table display")
                scale_table = scale_table.replace('border="1"', 'border="0"')
                finalDict['scale_table'] = str(scale_table)

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
