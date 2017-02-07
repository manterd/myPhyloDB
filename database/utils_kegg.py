import ast
import datetime
from django import db
from django.http import HttpResponse
from django_pandas.io import read_frame
import logging
import numpy as np
import pandas as pd
import ujson

from database.models import Kingdom, Phyla, Class, Order, Family, Genus, Species, OTU_99
from database.models import PICRUSt
from database.models import ko_lvl1, ko_lvl2, ko_lvl3, ko_entry
from database.models import nz_lvl1, nz_lvl2, nz_lvl3, nz_lvl4, nz_entry


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getTaxaDF(selectAll, taxaDict, savedDF, metaDF, allFields, DepVar, RID, stops, PID):
    try:
        missingList = []
        taxaDF = pd.DataFrame(columns=['sampleid', 'rank', 'rank_id', 'rank_name', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity'])
        if selectAll == 0:
            for key in taxaDict:
                taxaList = taxaDict[key]
                if isinstance(taxaList, str):
                    if key == 'Kingdom':
                        tempDF = savedDF.loc[savedDF['kingdomid'] == taxaList]
                        tempDF = tempDF[['sampleid', 'kingdomid', 'kingdomName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'kingdomid': 'rank_id', 'kingdomName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Kingdom'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Phyla':
                        tempDF = savedDF.loc[savedDF['phylaid'] == taxaList]
                        tempDF = tempDF[['sampleid', 'phylaid', 'phylaName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'phylaid': 'rank_id', 'phylaName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Phyla'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Class':
                        tempDF = savedDF.loc[savedDF['classid'] == taxaList]
                        tempDF = tempDF[['sampleid', 'classid', 'className', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'classid': 'rank_id', 'className': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Class'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Order':
                        tempDF = savedDF.loc[savedDF['orderid'] == taxaList]
                        tempDF = tempDF[['sampleid', 'orderid', 'orderName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'orderid': 'rank_id', 'orderName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Order'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Family':
                        tempDF = savedDF.loc[savedDF['familyid'] == taxaList]
                        tempDF = tempDF[['sampleid', 'familyid', 'familyName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'familyid': 'rank_id', 'familyName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Family'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Genus':
                        tempDF = savedDF.loc[savedDF['genusid'] == taxaList]
                        tempDF = tempDF[['sampleid', 'genusid', 'genusName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'genusid': 'rank_id', 'genusName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Genus'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Species':
                        tempDF = savedDF.loc[savedDF['speciesid'] == taxaList]
                        tempDF = tempDF[['sampleid', 'speciesid', 'speciesName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'speciesid': 'rank_id', 'speciesName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Species'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'OTU_99':
                        tempDF = savedDF.loc[savedDF['otuid'] == taxaList]
                        tempDF = tempDF[['sampleid', 'otuid', 'otuName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'otuid': 'rank_id', 'otuName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'OTU_99'
                        taxaDF = taxaDF.append(tempDF)
                else:
                    if key == 'Kingdom':
                        tempDF = savedDF.loc[savedDF['kingdomid'].isin(taxaList)]
                        tempDF = tempDF[['sampleid', 'kingdomid', 'kingdomName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'kingdomid': 'rank_id', 'kingdomName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Kingdom'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Phyla':
                        tempDF = savedDF.loc[savedDF['phylaid'].isin(taxaList)]
                        tempDF = tempDF[['sampleid', 'phylaid', 'phylaName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'phylaid': 'rank_id', 'phylaName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Phyla'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Class':
                        tempDF = savedDF.loc[savedDF['classid'].isin(taxaList)]
                        tempDF = tempDF[['sampleid', 'classid', 'className', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'classid': 'rank_id', 'className': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Class'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Order':
                        tempDF = savedDF.loc[savedDF['orderid'].isin(taxaList)]
                        tempDF = tempDF[['sampleid', 'orderid', 'orderName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'orderid': 'rank_id', 'orderName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Order'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Family':
                        tempDF = savedDF.loc[savedDF['familyid'].isin(taxaList)]
                        tempDF = tempDF[['sampleid', 'familyid', 'familyName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'familyid': 'rank_id', 'familyName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Family'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Genus':
                        tempDF = savedDF.loc[savedDF['genusid'].isin(taxaList)]
                        tempDF = tempDF[['sampleid', 'genusid', 'genusName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'genusid': 'rank_id', 'genusName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Genus'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Species':
                        tempDF = savedDF.loc[savedDF['speciesid'].isin(taxaList)]
                        tempDF = tempDF[['sampleid', 'speciesid', 'speciesName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'speciesid': 'rank_id', 'speciesName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'Species'
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'OTU_99':
                        tempDF = savedDF.loc[savedDF['otuid'].isin(taxaList)]
                        tempDF = tempDF[['sampleid', 'otuid', 'otuName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
                        tempDF.rename(columns={'otuid': 'rank_id', 'otuName': 'rank_name'}, inplace=True)
                        tempDF.loc[:, 'rank'] = 'OTU_99'
                        taxaDF = taxaDF.append(tempDF)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif selectAll == 1:
            taxaDF = savedDF.loc[:, ['sampleid', 'kingdomid', 'kingdomName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'kingdomid': 'rank_id', 'kingdomName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Kingdom'
        elif selectAll == 2:
            taxaDF = savedDF.loc[:, ['sampleid', 'phylaid', 'phylaName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'phylaid': 'rank_id', 'phylaName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Phyla'
        elif selectAll == 3:
            taxaDF = savedDF.loc[:, ['sampleid', 'classid', 'className', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'classid': 'rank_id', 'className': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Class'
        elif selectAll == 4:
            taxaDF = savedDF.loc[:, ['sampleid', 'orderid', 'orderName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'orderid': 'rank_id', 'orderName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Order'
        elif selectAll == 5:
            taxaDF = savedDF.loc[:, ['sampleid', 'familyid', 'familyName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'familyid': 'rank_id', 'familyName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Family'
        elif selectAll == 6:
            taxaDF = savedDF.loc[:, ['sampleid', 'genusid', 'genusName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'genusid': 'rank_id', 'genusName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Genus'
        elif selectAll == 7:
            taxaDF = savedDF.loc[:, ['sampleid', 'speciesid', 'speciesName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'speciesid': 'rank_id', 'speciesName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Species'
        elif selectAll == 9:
            taxaDF = savedDF.loc[:, ['sampleid', 'otuid', 'otuName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'otuid': 'rank_id', 'otuName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'OTU_99'
        elif selectAll == 8:
            taxaDF = savedDF.loc[:, ['sampleid', 'genusid', 'genusName', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'genusid': 'rank_id', 'genusName': 'rank_name'}, inplace=True)
            pgprList = ['Actinomyces', 'Arthrobacter', 'Azospirillum', 'Bacillus', 'Burkholderia', 'Enterobacter', 'Mitsuaria', 'Pasteuria', 'Pseudomonas', 'Rhizobium']
            taxaDF = taxaDF.loc[taxaDF['rank_name'].isin(pgprList)]
            taxaDF.loc[:, 'rank'] = 'Genus'
            foundList = list(taxaDF['rank_name'])
            missingList = list(set(pgprList) - set(foundList))

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        metaDF.reset_index(drop=False, inplace=True)
        finalDF = pd.merge(metaDF, taxaDF, left_on='sampleid', right_on='sampleid', how='inner')
        finalDF.reset_index(drop=False, inplace=True)

        wantedList = allFields + ['sampleid', 'rank', 'rank_name', 'rank_id']
        if DepVar == 0:
            finalDF = finalDF.groupby(wantedList)[['abund']].sum()
        elif DepVar == 1:
            finalDF = finalDF.groupby(wantedList)[['rel_abund']].sum()
        elif DepVar == 4:
            finalDF = finalDF.groupby(wantedList)[['abund_16S']].sum()
            finalDF['abund_16S'] = finalDF['abund_16S'].astype(int)
        elif DepVar == 2:
            finalDF = finalDF.groupby(wantedList)[['rich']].sum()
        elif DepVar == 3:
            finalDF = finalDF.groupby(wantedList)[['diversity']].sum()
        elif DepVar == 10:
            finalDF = finalDF.groupby(wantedList)[['abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']].sum()
        finalDF.reset_index(drop=False, inplace=True)
        return finalDF, missingList

    except Exception:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "There was an error with your analysis!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = ujson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getKeggDF(keggAll, keggDict, savedDF, metaDF, allFields, DepVar, RID, stops, PID):
    try:
        koDict = {}
        if keggAll == 0:
            for key in keggDict:
                keggList = keggDict[key]
                if isinstance(keggList, str):
                    if key == 'Level1':
                        koList = ko_entry.objects.using('picrust').filter(ko_lvl1_id_id=keggList).values_list('ko_orthology', flat=True)
                        if koList:
                            koDict[keggList] = koList
                    elif key == 'Level2':
                        koList = ko_entry.objects.using('picrust').filter(ko_lvl2_id_id=keggList).values_list('ko_orthology', flat=True)
                        if koList:
                            koDict[keggList] = koList
                    elif key == 'Level3':
                        koList = ko_entry.objects.using('picrust').filter(ko_lvl3_id_id=keggList).values_list('ko_orthology', flat=True)
                        if koList:
                            koDict[keggList] = koList
                    elif key == 'Level4':
                        koList = ko_entry.objects.using('picrust').filter(ko_lvl4_id=keggList).values_list('ko_orthology', flat=True)
                        if koList:
                            koDict[keggList] = koList
                else:
                    for i in keggList:
                        if key == 'Level1':
                            koList = ko_entry.objects.using('picrust').filter(ko_lvl1_id_id=i).values_list('ko_orthology', flat=True)
                            if koList:
                                koDict[i] = koList
                        elif key == 'Level2':
                            koList = ko_entry.objects.using('picrust').filter(ko_lvl2_id_id=i).values_list('ko_orthology', flat=True)
                            if koList:
                                koDict[i] = koList
                        elif key == 'Level3':
                            koList = ko_entry.objects.using('picrust').filter(ko_lvl3_id_id=i).values_list('ko_orthology', flat=True)
                            if koList:
                                koDict[i] = koList
                        elif key == 'Level4':
                            koList = ko_entry.objects.using('picrust').filter(ko_lvl4_id=i).values_list('ko_orthology', flat=True)
                            if koList:
                                koDict[i] = koList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif keggAll == 1:
            keys = ko_lvl1.objects.using('picrust').values_list('ko_lvl1_id', flat=True)
            for key in keys:
                koList = ko_entry.objects.using('picrust').filter(ko_lvl1_id_id=key).values_list('ko_orthology', flat=True)
                if koList:
                    koDict[key] = koList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif keggAll == 2:
            keys = ko_lvl2.objects.using('picrust').values_list('ko_lvl2_id', flat=True)
            for key in keys:
                koList = ko_entry.objects.using('picrust').filter(ko_lvl2_id_id=key).values_list('ko_orthology', flat=True)
                if koList:
                    koDict[key] = koList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif keggAll == 3:
            keys = ko_lvl3.objects.using('picrust').values_list('ko_lvl3_id', flat=True)
            for key in keys:
                koList = ko_entry.objects.using('picrust').filter(ko_lvl3_id_id=key).values_list('ko_orthology', flat=True)
                if koList:
                    koDict[key] = koList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif keggAll == 4:
            keys = ko_entry.objects.using('picrust').values_list('ko_lvl4_id', flat=True)
            for key in keys:
                koList = ko_entry.objects.using('picrust').filter(ko_lvl4_id=key).values_list('ko_orthology', flat=True)
                if koList:
                    koDict[key] = koList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        # create sample and otu lists based on meta data selection
        wanted = ['sampleid', 'otuid', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']
        profileDF = savedDF.loc[:, wanted]
        profileDF.set_index('otuid', inplace=True)

        # get PICRUSt data for otu
        otuList = pd.unique(profileDF.index.ravel().tolist())
        qs = PICRUSt.objects.using('picrust').filter(otuid__in=otuList)
        picrustDF = read_frame(qs, fieldnames=['otuid__otuid', 'geneCount'])
        picrustDF.set_index('otuid__otuid', inplace=True)

        levelList = []
        for key in koDict:
            levelList.append(str(key))
            picrustDF[key] = 0.0

        sumKEGG(otuList, picrustDF, koDict, RID, PID, stops)
        picrustDF.drop('geneCount', axis=1, inplace=True)
        picrustDF[picrustDF > 0.0] = 1.0

        # merge to get final gene counts for all selected samples
        taxaDF = pd.merge(profileDF, picrustDF, left_index=True, right_index=True, how='inner')
        taxaDF.dropna(axis=0, how='any', inplace=True)

        for level in levelList:
            if DepVar == 0:
                taxaDF[level] = taxaDF['abund'] * taxaDF[level]
            elif DepVar == 1:
                taxaDF[level] = taxaDF['rel_abund'] * taxaDF[level]
            elif DepVar == 2:
                taxaDF[level] = taxaDF['rich'] * taxaDF[level]
            elif DepVar == 3:
                taxaDF[level] = taxaDF['diversity'] * taxaDF[level]
            elif DepVar == 4:
                taxaDF[level] = taxaDF['abund_16S'] * taxaDF[level]

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        # create dataframe with data by otu
        wanted = levelList[:]
        wanted.insert(0, 'sampleid')
        allDF = taxaDF[wanted]

        for level in levelList:
            name = ''
            if ko_lvl1.objects.using('picrust').filter(ko_lvl1_id=level).exists():
                name = ko_lvl1.objects.using('picrust').get(ko_lvl1_id=level).ko_lvl1_name
            elif ko_lvl2.objects.using('picrust').filter(ko_lvl2_id=level).exists():
                name = ko_lvl2.objects.using('picrust').get(ko_lvl2_id=level).ko_lvl2_name
            elif ko_lvl3.objects.using('picrust').filter(ko_lvl3_id=level).exists():
                name = ko_lvl3.objects.using('picrust').get(ko_lvl3_id=level).ko_lvl3_name
            elif ko_entry.objects.using('picrust').filter(ko_lvl4_id=level).exists():
                name = ko_entry.objects.using('picrust').get(ko_lvl4_id=level).ko_desc
            allDF.rename(columns={level: name}, inplace=True)

        # get taxonomy names
        idList = list(set(allDF.index.tolist()))
        nameDict = getFullTaxonomy(idList)
        allDF.index.name = 'otuid'
        allDF.reset_index(drop=False, inplace=True)
        allDF['Taxonomy'] = allDF['otuid'].map(nameDict)
        allDF.set_index('otuid', inplace=True)

        allDF.set_index('sampleid', drop=True, inplace=True)
        allDF = pd.merge(metaDF, allDF, left_index=True, right_index=True, how='inner')
        allDF.reset_index(drop=False, inplace=True)
        allDF = allDF.loc[(allDF.sum(axis=1) != 0)]

        # sum all otu
        taxaDF = taxaDF.groupby('sampleid')[levelList].agg('sum')
        taxaDF.reset_index(drop=False, inplace=True)

        if DepVar == 0:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund')
        elif DepVar == 1:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rel_abund')
        elif DepVar == 2:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rich')
        elif DepVar == 3:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='diversity')
        elif DepVar == 4:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund_16S')

        taxaDF.set_index('sampleid', drop=True, inplace=True)
        finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')
        finalDF.reset_index(drop=False, inplace=True)

        finalDF['rank'] = ''
        finalDF['rank_name'] = ''
        for index, row in finalDF.iterrows():
            if ko_lvl1.objects.using('picrust').filter(ko_lvl1_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl1'
                finalDF.loc[index, 'rank_name'] = ko_lvl1.objects.using('picrust').get(ko_lvl1_id=row['rank_id']).ko_lvl1_name
            elif ko_lvl2.objects.using('picrust').filter(ko_lvl2_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl2'
                finalDF.loc[index, 'rank_name'] = ko_lvl2.objects.using('picrust').get(ko_lvl2_id=row['rank_id']).ko_lvl2_name
            elif ko_lvl3.objects.using('picrust').filter(ko_lvl3_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl3'
                finalDF.loc[index, 'rank_name'] = ko_lvl3.objects.using('picrust').get(ko_lvl3_id=row['rank_id']).ko_lvl3_name
            elif ko_entry.objects.using('picrust').filter(ko_lvl4_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl4'
                finalDF.loc[index, 'rank_name'] = ko_entry.objects.using('picrust').get(ko_lvl4_id=row['rank_id']).ko_name

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        return finalDF, allDF

    except Exception:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "There was an error with your analysis!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = ujson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getNZDF(nzAll, myDict, savedDF, metaDF, allFields, DepVar, RID, stops, PID):
    try:
        nzDict = {}
        if nzAll == 0:
            for key in myDict:
                keggList = myDict[key]
                if isinstance(keggList, str):
                    if key == 'Level1':
                        nzList = nz_entry.objects.using('picrust').filter(nz_lvl1_id_id=keggList).values_list('nz_orthology', flat=True)
                        if nzList:
                            nzDict[keggList] = nzList
                    elif key == 'Level2':
                        nzList = nz_entry.objects.using('picrust').filter(nz_lvl2_id_id=keggList).values_list('nz_orthology', flat=True)
                        if nzList:
                            nzDict[keggList] = nzList
                    elif key == 'Level3':
                        nzList = nz_entry.objects.using('picrust').filter(nz_lvl3_id_id=keggList).values_list('nz_orthology', flat=True)
                        if nzList:
                            nzDict[keggList] = nzList
                    elif key == 'Level4':
                        nzList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=keggList).values_list('nz_orthology', flat=True)
                        if nzList:
                            nzDict[keggList] = nzList
                    elif key == 'Level5':
                        nzList = nz_entry.objects.using('picrust').filter(nz_lvl5_id=keggList).values_list('nz_orthology', flat=True)
                        if nzList:
                            nzDict[keggList] = nzList
                else:
                    for i in keggList:
                        if key == 'Level1':
                            nzList = nz_entry.objects.using('picrust').filter(nz_lvl1_id_id=i).values_list('nz_orthology', flat=True)
                            if nzList:
                                nzDict[i] = nzList
                        elif key == 'Level2':
                            nzList = nz_entry.objects.using('picrust').filter(nz_lvl2_id_id=i).values_list('nz_orthology', flat=True)
                            if nzList:
                                nzDict[i] = nzList
                        elif key == 'Level3':
                            nzList = nz_entry.objects.using('picrust').filter(nz_lvl3_id_id=i).values_list('nz_orthology', flat=True)
                            if nzList:
                                nzDict[i] = nzList
                        elif key == 'Level4':
                            nzList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=i).values_list('nz_orthology', flat=True)
                            if nzList:
                                nzDict[i] = nzList
                        elif key == 'Level5':
                            nzList = nz_entry.objects.using('picrust').filter(nz_lvl5_id=i).values_list('nz_orthology', flat=True)
                            if nzList:
                                nzDict[i] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 1:
            keys = nz_lvl1.objects.using('picrust').values_list('nz_lvl1_id', flat=True)
            for key in keys:
                nzList = nz_entry.objects.using('picrust').filter(nz_lvl1_id_id=key).values_list('nz_orthology', flat=True)
                if nzList:
                    nzDict[key] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 2:
            keys = nz_lvl2.objects.using('picrust').values_list('nz_lvl2_id', flat=True)
            for key in keys:
                nzList = nz_entry.objects.using('picrust').filter(nz_lvl2_id_id=key).values_list('nz_orthology', flat=True)
                if nzList:
                    nzDict[key] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 3:
            keys = nz_lvl3.objects.using('picrust').values_list('nz_lvl3_id', flat=True)
            for key in keys:
                nzList = nz_entry.objects.using('picrust').filter(nz_lvl3_id_id=key).values_list('nz_orthology', flat=True)
                if nzList:
                    nzDict[key] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 4:
            keys = nz_lvl4.objects.using('picrust').values_list('nz_lvl4_id', flat=True)
            for key in keys:
                nzList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=key).values_list('nz_orthology', flat=True)
                if nzList:
                    nzDict[key] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 5:
            ### N-fixation
            # 1.18.6.1  nitrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.18.6.1  nitrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### PO4 Solubility
            # 1.3.3.11  pyrroloquinoline-quinone synthase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.3.3.11  pyrroloquinoline-quinone synthase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### Biocontrol
            # 1.4.99.5  glycine dehydrogenase (cyanide-forming)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.4.99.5  glycine dehydrogenase (cyanide-forming)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 4.1.1.5  acetolactate decarboxylase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='4.1.1.5  acetolactate decarboxylase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.1.1.76  (S,S)-butanediol dehydrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.1.1.76  (S,S)-butanediol dehydrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.6  endo-1,3(4)-beta-glucanase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='3.2.1.6  endo-1,3(4)-beta-glucanase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.14  chitinase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='3.2.1.14  chitinase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### Root growth (drought/salt stress)
            # 4.1.1.74  indolepyruvate decarboxylase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='4.1.1.74  indolepyruvate decarboxylase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### Siderophores
            # 6.3.2.39  aerobactin synthase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='6.3.2.39  aerobactin synthase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 6.3.2.14  enterobactin synthase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='6.3.2.14  enterobactin synthase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # K04787 mycobactin salicyl-AMP ligase [EC:6.3.2.-]
            id = nz_entry.objects.using('picrust').get(nz_orthology='K04787').nz_lvl5_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl5_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # K04783 yersiniabactin salicyl-AMP ligase [EC:6.3.2.-]
            id = nz_entry.objects.using('picrust').get(nz_orthology='K04783').nz_lvl5_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl5_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### C decomposition
            # 3.2.1.4  cellulase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.4  cellulase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.91  cellulose 1,4-beta-cellobiosidase (non-reducing end)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.91  cellulose 1,4-beta-cellobiosidase (non-reducing end)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.21  beta-glucosidase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.21  beta-glucosidase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.8  endo-1,4-beta-xylanase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.8  endo-1,4-beta-xylanase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.37  xylan 1,4-beta-xylosidase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.37  xylan 1,4-beta-xylosidase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### N decomposition
            # 3.5.1.4  amidase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.5.1.4  amidase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.5.1.5  urease
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.5.1.5  urease').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### P decomposition
            # 3.1.3.1  alkaline phosphatase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.1.3.1  alkaline phosphatase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.1.3.2  acid phosphatase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.1.3.2  acid phosphatase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### S decomposition
            # 3.1.6.1  arylsulfatase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.1.6.1  arylsulfatase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

        elif nzAll == 6:
            ### N-fixation
            # 1.18.6.1  nitrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.18.6.1  nitrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### Ammonification
            # 3.5.1.4  amidase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.5.1.4  amidase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.5.1.5  urease
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='3.5.1.5  urease').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### Nitrification
            # 1.14.99.39  ammonia monooxygenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.14.99.39  ammonia monooxygenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.6  hydroxylamine dehydrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.6  hydroxylamine dehydrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.99.4  nitrate reductase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.99.4  nitrate reductase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### Dissimilatory nitrate reduction
            # 1.7.5.1  nitrate reductase (quinone)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.5.1  nitrate reductase (quinone)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.99.4  nitrate reductase

            # 1.7.1.15  nitrite reductase (NADH)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.1.15  nitrite reductase (NADH)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.2  nitrite reductase (cytochrome; ammonia-forming)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.2  nitrite reductase (cytochrome; ammonia-forming)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### Assimilatory nitrate reduction
            # 1.7.7.1  ferredoxin---nitrite reductase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.7.1  ferredoxin---nitrite reductase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.7.2  ferredoxin---nitrate reductase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.7.2  ferredoxin---nitrate reductase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.99.4  nitrate reductase

            # 1.7.1.4  nitrite reductase [NAD(P)H]
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.1.4  nitrite reductase [NAD(P)H]').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.1.1  nitrate reductase (NADH)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.1.1  nitrate reductase (NADH)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### Denitrification
            # 1.7.99.4  nitrate reductase

            # 1.7.2.1  nitrite reductase (NO-forming)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.1  nitrite reductase (NO-forming)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.5  nitric oxide reductase (cytochrome c)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.5  nitric oxide reductase (cytochrome c)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.4  nitrous-oxide reductase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.4  nitrous-oxide reductase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### Anammox
            # 1.7.2.7  hydrazine synthase subunit
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.7  hydrazine synthase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.8  hydrazine dehydrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.8  hydrazine dehydrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

        elif nzAll == 10:
            ### N-fixation
            # 1.18.6.1  nitrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.18.6.1  nitrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### PO4 Solubility
            # 1.3.3.11  pyrroloquinoline-quinone synthase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.3.3.11  pyrroloquinoline-quinone synthase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### Biocontrol
            # 1.4.99.5  glycine dehydrogenase (cyanide-forming)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.4.99.5  glycine dehydrogenase (cyanide-forming)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 4.1.1.5  acetolactate decarboxylase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='4.1.1.5  acetolactate decarboxylase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.1.1.76  (S,S)-butanediol dehydrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.1.1.76  (S,S)-butanediol dehydrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.6  endo-1,3(4)-beta-glucanase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='3.2.1.6  endo-1,3(4)-beta-glucanase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.14  chitinase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='3.2.1.14  chitinase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### Root growth (drought/salt stress)
            # 4.1.1.74  indolepyruvate decarboxylase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='4.1.1.74  indolepyruvate decarboxylase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### Siderophores
            # 6.3.2.39  aerobactin synthase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='6.3.2.39  aerobactin synthase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 6.3.2.14  enterobactin synthase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='6.3.2.14  enterobactin synthase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # K04787 mycobactin salicyl-AMP ligase [EC:6.3.2.-]
            id = nz_entry.objects.using('picrust').get(nz_orthology='K04787').nz_lvl5_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl5_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # K04783 yersiniabactin salicyl-AMP ligase [EC:6.3.2.-]
            id = nz_entry.objects.using('picrust').get(nz_orthology='K04783').nz_lvl5_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl5_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### C decomposition
            # 3.2.1.4  cellulase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.4  cellulase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.91  cellulose 1,4-beta-cellobiosidase (non-reducing end)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.91  cellulose 1,4-beta-cellobiosidase (non-reducing end)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.21  beta-glucosidase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.21  beta-glucosidase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.8  endo-1,4-beta-xylanase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.8  endo-1,4-beta-xylanase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.37  xylan 1,4-beta-xylosidase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.2.1.37  xylan 1,4-beta-xylosidase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### N decomposition
            # 3.5.1.4  amidase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.5.1.4  amidase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.5.1.5  urease
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.5.1.5  urease').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### P decomposition
            # 3.1.3.1  alkaline phosphatase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.1.3.1  alkaline phosphatase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.1.3.2  acid phosphatase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.1.3.2  acid phosphatase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### S decomposition
            # 3.1.6.1  arylsulfatase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.1.6.1  arylsulfatase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### N-fixation
            # 1.18.6.1  nitrogenase

            ### Ammonification
            # 3.5.1.4  amidase

            # 3.5.1.5  urease

            ### Nitrification
            # 1.14.99.39  ammonia monooxygenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.14.99.39  ammonia monooxygenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.6  hydroxylamine dehydrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.6  hydroxylamine dehydrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.99.4  nitrate reductase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.99.4  nitrate reductase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### Dissimilatory nitrate reduction
            # 1.7.5.1  nitrate reductase (quinone)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.5.1  nitrate reductase (quinone)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.99.4  nitrate reductase

            # 1.7.1.15  nitrite reductase (NADH)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.1.15  nitrite reductase (NADH)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.2  nitrite reductase (cytochrome; ammonia-forming)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.2  nitrite reductase (cytochrome; ammonia-forming)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### Assimilatory nitrate reduction
            # 1.7.7.1  ferredoxin---nitrite reductase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.7.1  ferredoxin---nitrite reductase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.7.2  ferredoxin---nitrate reductase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.7.2  ferredoxin---nitrate reductase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.99.4  nitrate reductase

            # 1.7.1.4  nitrite reductase [NAD(P)H]
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.1.4  nitrite reductase [NAD(P)H]').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.1.1  nitrate reductase (NADH)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.1.1  nitrate reductase (NADH)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### Denitrification
            # 1.7.99.4  nitrate reductase

            # 1.7.2.1  nitrite reductase (NO-forming)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.1  nitrite reductase (NO-forming)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.5  nitric oxide reductase (cytochrome c)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.5  nitric oxide reductase (cytochrome c)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.4  nitrous-oxide reductase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.4  nitrous-oxide reductase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            ### Anammox
            # 1.7.2.7  hydrazine synthase subunit
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.7  hydrazine synthase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.7.2.8  hydrazine dehydrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.7.2.8  hydrazine dehydrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

        # create sample and otu lists based on meta data selection
        wanted = ['sampleid', 'otuid', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']
        profileDF = savedDF.loc[:, wanted]
        profileDF.set_index('otuid', inplace=True)

        # get PICRUSt data for otu
        otuList = pd.unique(profileDF.index.ravel().tolist())
        qs = PICRUSt.objects.using('picrust').filter(otuid__in=otuList)
        picrustDF = read_frame(qs, fieldnames=['otuid__otuid', 'geneCount'], verbose=False)
        picrustDF.set_index('otuid__otuid', inplace=True)

        levelList = []
        for key in nzDict:
            levelList.append(str(key))
            picrustDF[key] = 0.0

        sumKEGG(otuList, picrustDF, nzDict, RID, PID, stops)
        picrustDF.drop('geneCount', axis=1, inplace=True)
        picrustDF[picrustDF > 0.0] = 1.0

        # merge to get final gene counts for all selected samples
        taxaDF = pd.merge(profileDF, picrustDF, left_index=True, right_index=True, how='inner')
        taxaDF.dropna(axis=0, how='any', inplace=True)

        for level in levelList:
            if DepVar == 0:
                taxaDF[level] = taxaDF['abund'] * taxaDF[level]
            elif DepVar == 1:
                taxaDF[level] = taxaDF['rel_abund'] * taxaDF[level]
            elif DepVar == 2:
                taxaDF[level] = taxaDF['rich'] * taxaDF[level]
            elif DepVar == 3:
                taxaDF[level] = taxaDF['diversity'] * taxaDF[level]
            elif DepVar == 4:
                taxaDF[level] = taxaDF['abund_16S'] * taxaDF[level]

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        # create dataframe with data by species
        wanted = list(levelList)
        wanted.insert(0, 'sampleid')
        allDF = taxaDF[wanted]

        for level in levelList:
            name = ''
            if nz_lvl1.objects.using('picrust').filter(nz_lvl1_id=level).exists():
                name = nz_lvl1.objects.using('picrust').get(nz_lvl1_id=level).nz_lvl1_name
            elif nz_lvl2.objects.using('picrust').filter(nz_lvl2_id=level).exists():
                name = nz_lvl2.objects.using('picrust').get(nz_lvl2_id=level).nz_lvl2_name
            elif nz_lvl3.objects.using('picrust').filter(nz_lvl3_id=level).exists():
                name = nz_lvl3.objects.using('picrust').get(nz_lvl3_id=level).nz_lvl3_name
            elif nz_lvl4.objects.using('picrust').filter(nz_lvl4_id=level).exists():
                name = nz_lvl4.objects.using('picrust').get(nz_lvl4_id=level).nz_lvl4_name
            elif nz_entry.objects.using('picrust').filter(nz_lvl5_id=level).exists():
                name = nz_entry.objects.using('picrust').get(nz_lvl5_id=level).nz_desc
            allDF.rename(columns={level: name}, inplace=True)

        # get taxonomy names
        idList = list(set(allDF.index.tolist()))
        nameDict = getFullTaxonomy(idList)
        allDF.index.name = 'otuid'
        allDF.reset_index(drop=False, inplace=True)
        allDF['Taxonomy'] = allDF['otuid'].map(nameDict)
        allDF.set_index('otuid', inplace=True)

        allDF.set_index('sampleid', drop=True, inplace=True)
        allDF = pd.merge(metaDF, allDF, left_index=True, right_index=True, how='inner')
        allDF.reset_index(drop=False, inplace=True)
        allDF = allDF.loc[(allDF.sum(axis=1) != 0)]

        # sum all otu
        taxaDF = taxaDF.groupby('sampleid')[levelList].agg('sum')
        taxaDF.reset_index(drop=False, inplace=True)

        if DepVar == 0:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund')
        elif DepVar == 1:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rel_abund')
        elif DepVar == 2:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rich')
        elif DepVar == 3:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='diversity')
        elif DepVar == 4:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund_16S')

        taxaDF.set_index('sampleid', drop=True, inplace=True)
        finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')
        finalDF.reset_index(drop=False, inplace=True)

        finalDF['rank'] = ''
        finalDF['rank_name'] = ''

        for index, row in finalDF.iterrows():
            if nz_lvl1.objects.using('picrust').filter(nz_lvl1_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl1'
                finalDF.loc[index, 'rank_name'] = nz_lvl1.objects.using('picrust').get(nz_lvl1_id=row['rank_id']).nz_lvl1_name
            elif nz_lvl2.objects.using('picrust').filter(nz_lvl2_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl2'
                finalDF.loc[index, 'rank_name'] = nz_lvl2.objects.using('picrust').get(nz_lvl2_id=row['rank_id']).nz_lvl2_name
            elif nz_lvl3.objects.using('picrust').filter(nz_lvl3_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl3'
                finalDF.loc[index, 'rank_name'] = nz_lvl3.objects.using('picrust').get(nz_lvl3_id=row['rank_id']).nz_lvl3_name
            elif nz_lvl4.objects.using('picrust').filter(nz_lvl4_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl4'
                finalDF.loc[index, 'rank_name'] = nz_lvl4.objects.using('picrust').get(nz_lvl4_id=row['rank_id']).nz_lvl4_name
            elif nz_entry.objects.using('picrust').filter(nz_lvl5_id=row['rank_id']).exists():
                finalDF.loc[index, 'rank'] = 'Lvl5'
                finalDF.loc[index, 'rank_name'] = nz_entry.objects.using('picrust').get(nz_lvl5_id=row['rank_id']).nz_desc

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        return finalDF, allDF

    except Exception:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "There was an error with your analysis!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = ujson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def sumKEGG(otuList, picrustDF, keggDict, RID, PID, stops):
    db.close_old_connections()
    for otu in otuList:
        try:
            cell = picrustDF.at[otu, 'geneCount']
            d = ast.literal_eval(cell)

            for key in keggDict:
                sum = 0.0
                myList = keggDict[key]

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                for k in myList:
                    if k in d:
                        sum += d[k]
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                    if stops[PID] == RID:
                        res = ''
                        return HttpResponse(res, content_type='application/json')
                    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

                picrustDF.at[otu, key] = sum
        except Exception:
            pass


def getFullTaxonomy(idList):
    recordDict = {}
    for id in idList:
        if Kingdom.objects.all().filter(kingdomid=id).exists():
            qs = Kingdom.objects.all().filter(kingdomid=id).values_list('kingdomName')
        else:
            if Phyla.objects.all().filter(phylaid=id).exists():
                qs = Phyla.objects.all().filter(phylaid=id).values_list('kingdomid_id__kingdomName', 'phylaName')
            else:
                if Class.objects.all().filter(classid=id).exists():
                    qs = Class.objects.all().filter(classid=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'className')
                else:
                    if Order.objects.all().filter(orderid=id).exists():
                        qs = Order.objects.all().filter(orderid=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderName')
                    else:
                        if Family.objects.all().filter(familyid=id).exists():
                            qs = Family.objects.all().filter(familyid=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyName')
                        else:
                            if Genus.objects.all().filter(genusid=id).exists():
                                qs = Genus.objects.all().filter(genusid=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyid_id__familyName', 'genusName')
                            else:
                                if Species.objects.all().filter(speciesid=id).exists():
                                    qs = Species.objects.all().filter(speciesid=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyid_id__familyName', 'genusid_id__genusName', 'speciesName')
                                else:
                                    if OTU_99.objects.all().filter(otuid=id).exists():
                                        qs = OTU_99.objects.all().filter(otuid=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyid_id__familyName', 'genusid_id__genusName', 'speciesid_id__speciesName', 'otuName')
                                    else:
                                        qs = ('Not found', )
        record = ';'.join(qs[0])
        recordDict[id] = record
    return recordDict


def getFullKO(idList):
    recordDict = {}
    for id in idList:
        if ko_lvl1.objects.using('picrust').all().filter(ko_lvl1_id=id).exists():
            qs = ko_lvl1.objects.using('picrust').all().filter(ko_lvl1_id=id).values_list('ko_lvl1_name')
        else:
            if ko_lvl2.objects.using('picrust').all().filter(ko_lvl2_id=id).exists():
                qs = ko_lvl2.objects.using('picrust').all().filter(ko_lvl2_id=id).values_list('ko_lvl1_id_id__ko_lvl1_name', 'ko_lvl2_name')
            else:
                if ko_lvl3.objects.using('picrust').all().filter(ko_lvl3_id=id).exists():
                    qs = ko_lvl3.objects.using('picrust').all().filter(ko_lvl3_id=id).values_list('ko_lvl1_id_id__ko_lvl1_name', 'ko_lvl2_id_id__ko_lvl2_name', 'ko_lvl3_name')
                else:
                    qs = ('Not found', )
        record = ';'.join(qs[0])
        recordDict[id] = record
    return recordDict


def getFullNZ(idList):
    recordDict = {}
    for id in idList:
        if nz_lvl1.objects.using('picrust').all().filter(nz_lvl1_id=id).exists():
            qs = nz_lvl1.objects.using('picrust').all().filter(nz_lvl1_id=id).values_list('nz_lvl1_name')
        else:
            if nz_lvl2.objects.using('picrust').all().filter(nz_lvl2_id=id).exists():
                qs = nz_lvl2.objects.using('picrust').all().filter(nz_lvl2_id=id).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_name')
            else:
                if nz_lvl3.objects.using('picrust').all().filter(nz_lvl3_id=id).exists():
                    qs = nz_lvl3.objects.using('picrust').all().filter(nz_lvl3_id=id).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_name')
                else:
                    if nz_lvl4.objects.using('picrust').all().filter(nz_lvl4_id=id).exists():
                        qs = nz_lvl4.objects.using('picrust').all().filter(nz_lvl4_id=id).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_id_id__nz_lvl3_name', 'nz_lvl4_name')
                    else:
                        qs = ('Not found', )
        record = ';'.join(qs[0])
        recordDict[id] = record
    return recordDict


def insertTaxaInfo(treeType, zipped, DF, pos=1):
    if treeType == 1:
        k, p, c, o, f, g, s, otu = map(None, *zipped)
        DF.insert(pos, 'Kingdom', k)
        DF.insert(pos+1, 'Phyla', p)
        DF.insert(pos+2, 'Class', c)
        DF.insert(pos+3, 'Order', o)
        DF.insert(pos+4, 'Family', f)
        DF.insert(pos+5, 'Genus', g)
        DF.insert(pos+6, 'Species', s)
        DF.insert(pos+7, 'OTU_99', otu)
    elif treeType == 2:
        L1, L2, L3 = map(None, *zipped)
        DF.insert(pos, 'Level_1', L1)
        DF.insert(pos+1, 'Level_2', L2)
        DF.insert(pos+2, 'Level_3', L3)
    elif treeType == 3:
        L1, L2, L3, L4 = map(None, *zipped)
        DF.insert(pos, 'Level_1', L1)
        DF.insert(pos+1, 'Level_2', L2)
        DF.insert(pos+2, 'Level_3', L3)
        DF.insert(pos+3, 'Level_4', L4)

    DF.fillna(value='N/A', inplace=True)


def filterDF(savedDF, DepVar, level, remUnclass, remZeroes, perZeroes, filterData, filterPer, filterMeth):

    numSamples = len(savedDF['sampleid'].unique())
    myVar = ''
    if DepVar == 0:
        myVar = 'abund'
    elif DepVar == 1:
        myVar = 'rel_abund'
    elif DepVar == 2:
        myVar = 'rich'
    elif DepVar == 3:
        myVar = 'diversity'
    elif DepVar == 4:
        myVar = 'abund_16S'

    myLevel = 'otuName'
    myID = ''
    numTaxa = 0
    if level == 1:
        myLevel = 'kingdomName'
        myID = 'kingdomid'
        numTaxa = len(savedDF['kingdomid'].unique())
    elif level == 2:
        myLevel = 'phylaName'
        myID = 'phylaid'
        numTaxa = len(savedDF['phylaid'].unique())
    elif level == 3:
        myLevel = 'className'
        myID = 'classid'
        numTaxa = len(savedDF['classid'].unique())
    elif level == 4:
        myLevel = 'orderName'
        myID = 'orderid'
        numTaxa = len(savedDF['orderid'].unique())
    elif level == 5:
        myLevel = 'familyName'
        myID = 'familyid'
        numTaxa = len(savedDF['familyid'].unique())
    elif level == 6:
        myLevel = 'genusName'
        myID = 'genusid'
        numTaxa = len(savedDF['genusid'].unique())
    elif level == 7:
        myLevel = 'speciesName'
        myID = 'speciesid'
        numTaxa = len(savedDF['speciesid'].unique())
    elif level == 9:
        myLevel = 'otuName'
        myID = 'otuid'
        numTaxa = len(savedDF['otuid'].unique())

    if remUnclass == 'yes':
        # check if selecting based on level first, create else statement for tree usage
        savedDF = savedDF[~savedDF[myLevel].str.contains('unclassified')]

    if remZeroes == 'yes' and perZeroes > 0:
        threshold = int(perZeroes / 100.0 * numSamples)
        bytag = savedDF.groupby(myID).aggregate(np.count_nonzero)
        tags = bytag[bytag[myVar] >= threshold].index.tolist()
        savedDF = savedDF[savedDF[myID].isin(tags)]

    if filterData == 'yes' and filterPer < 100 and not level == 8:
        threshold = int(filterPer / 100.0 * numTaxa)
        if filterMeth == 1:
            pass
        elif filterMeth == 2:
            bytag_q3 = savedDF.groupby(myID)[myVar].quantile(0.75)
            bytag_q1 = savedDF.groupby(myID)[myVar].quantile(0.25)
            bytag = bytag_q3 - bytag_q1
            bytag.sort(axis=0, ascending=False, inplace=True)
            tags = bytag[:threshold].index.tolist()
            savedDF = savedDF[savedDF[myID].isin(tags)]
        elif filterMeth == 3:
            bytag_sd = savedDF.groupby(myID)[myVar].std()
            bytag_mean = savedDF.groupby(myID)[myVar].mean()
            bytag = bytag_sd / bytag_mean
            bytag.sort(axis=0, ascending=False, inplace=True)
            tags = bytag[:threshold].index.tolist()
            savedDF = savedDF[savedDF[myID].isin(tags)]
        elif filterMeth == 4:
            bytag = savedDF.groupby(myID)[myVar].std()
            bytag.sort(axis=0, ascending=False, inplace=True)
            tags = bytag[:threshold].index.tolist()
            savedDF = savedDF[savedDF[myID].isin(tags)]
        elif filterMeth == 5:
            bytag = savedDF.groupby(myID)[myVar].mean()
            bytag.sort(axis=0, ascending=False, inplace=True)
            tags = bytag[:threshold].index.tolist()
            savedDF = savedDF[savedDF[myID].isin(tags)]
        elif filterMeth == 6:
            bytag = savedDF.groupby(myID)[myVar].median()
            bytag.sort(axis=0, ascending=False, inplace=True)
            tags = bytag[:threshold].index.tolist()
            savedDF = savedDF[savedDF[myID].isin(tags)]

    return savedDF


