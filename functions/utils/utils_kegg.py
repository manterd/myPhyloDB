import datetime
from django.http import HttpResponse
import logging
import numpy as np
import pandas as pd
import json

from database.models import Kingdom, Phyla, Class, Order, Family, Genus, Species, OTU_99, \
    PICRUSt, \
    ko_lvl1, ko_lvl2, ko_lvl3, ko_entry, \
    nz_lvl1, nz_lvl2, nz_lvl3, nz_lvl4, nz_entry

import functions


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getTaxaDF(selectAll, taxaDict, savedDF, metaDF, allFields, DepVar, RID, stops, PID):
    try:
        missingList = []
        taxaDF = pd.DataFrame(columns=['sampleid', 'rank', 'rank_id', 'rank_name', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity'])
        if selectAll == 0:
            for key in taxaDict:
                taxaList = taxaDict[key]
                if isinstance(taxaList, unicode):
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
            if taxaDF.empty:
                return pd.DataFrame(), pgprList

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

        if 'sample_name' not in allFields:
            wantedList = allFields + ['sample_name', 'sampleid', 'rank', 'rank_name', 'rank_id']
        else:
            wantedList = allFields + ['sampleid', 'rank', 'rank_name', 'rank_id']

        if DepVar == 0:
            finalDF = finalDF.groupby(wantedList)[['abund']].sum()
        elif DepVar == 1:
            finalDF = finalDF.groupby(wantedList)[['rel_abund']].sum()
        elif DepVar == 4:
            finalDF = finalDF.groupby(wantedList)[['abund_16S']].sum()
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
            res = json.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getKeggDF(keggAll, keggDict, savedDF, metaDF, DepVar, mapTaxa, RID, stops, PID):
    try:
        koDict = {}
        if keggAll == 0:
            for key in keggDict:
                keggList = keggDict[key]
                if isinstance(keggList, unicode):
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
        profileDF.set_index('otuid', drop=True, inplace=True)

        # get PICRUSt data for otu
        otuList = pd.unique(profileDF.index.ravel().tolist())
        qs = PICRUSt.objects.using('picrust').filter(otuid__in=otuList)
        piList = []
        for item in qs:
            piList.append([item.otuid_id, item.geneList])
        picrustDF = pd.DataFrame(piList, columns=['otuid', 'geneList'])
        picrustDF.set_index('otuid', drop=True, inplace=True)

        curStep = functions.getBase(RID)
        sumKEGG(picrustDF, koDict, 0, RID, PID, stops)
        functions.setBase(RID, curStep)

        picrustDF.drop('geneList', axis=1, inplace=True)
        picrustDF.fillna(value=0, inplace=True)
        levelList = picrustDF.columns.values.tolist()

        # convert profile to index (sampleid) and columns (keggid) and values (depvar)
        profileDF.reset_index(drop=False, inplace=True)
        if DepVar == 0:
            profileDF = profileDF.pivot(index='otuid', columns='sampleid', values='abund')
        elif DepVar == 1:
            profileDF = profileDF.pivot(index='otuid', columns='sampleid', values='rel_abund')
        elif DepVar == 2:
            profileDF = profileDF.pivot(index='otuid', columns='sampleid', values='rich')
        elif DepVar == 3:
            profileDF = profileDF.pivot(index='otuid', columns='sampleid', values='diversity')
        elif DepVar == 4:
            profileDF = profileDF.pivot(index='otuid', columns='sampleid', values='abund_16S')

        sampleList = profileDF.columns.values.tolist()
        taxaDF = pd.DataFrame()
        total = len(sampleList)
        counter = 1
        for i in sampleList:
            tempDF = pd.DataFrame(index=profileDF.index)
            for j in levelList:
                tempDF[j] = profileDF[i] * picrustDF[j]

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            tempDF = tempDF[levelList].sum()
            tempDF['sampleid'] = i
            taxaDF = taxaDF.append(tempDF, ignore_index=True)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            functions.setBase(RID, 'Calculating KEGG pathway abundances...sample ' + str(counter) + ' out of ' + str(total) + ' is finished!')
            counter += 1

        functions.setBase(RID, curStep)

        # df of mapped taxa to selected kegg orthologies
        if mapTaxa == 'yes':
            namesDict = {}
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
                else:
                    print str(level) + ' not found in database'
                namesDict[level] = name

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            allDF = picrustDF.reset_index(drop=False, inplace=False)
            allDF.rename(columns={'index': 'otuid'}, inplace=True)
            allDF.dropna(axis=1, how='any', inplace=False)
            allDF = allDF[(allDF.sum(axis=1) != 0)]

            for i in levelList:
                allDF[i] = allDF[i] / allDF[i]
            allDF.fillna(value=0, inplace=True)
            recordDict = {}
            for id in allDF.otuid.tolist():
                try:
                    qs = OTU_99.objects.all().filter(otuid=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyid_id__familyName', 'genusid_id__genusName', 'speciesid_id__speciesName', 'otuName')
                    record = '|'.join(qs[0])
                    recordDict[id] = record
                except:
                    recordDict[id] = 'No data'
            allDF['Taxonomy'] = allDF['otuid'].map(recordDict)
            order = ['otuid', 'Taxonomy'] + levelList
            allDF = allDF[order]
            allDF.rename(columns=namesDict, inplace=True)
        else:
            allDF = ''

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        #taxaDF.reset_index(drop=False, inplace=True)
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

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        idList = functions.getFullKO(list(finalDF.rank_id.unique()))
        finalDF['rank_name'] = finalDF['rank_id'].map(idList)

        # required to match getTaxaDF output
        metaDF.reset_index(drop=False, inplace=True)

        return finalDF, allDF

    except Exception:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "There was an error with your analysis!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = json.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getNZDF(nzAll, myDict, savedDF, metaDF,  DepVar, mapTaxa, RID, stops, PID):
    try:
        nzDict = {}
        if nzAll == 0:
            for key in myDict:
                keggList = myDict[key]
                if isinstance(keggList, unicode):
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
            id = 'nifH: nitrogenase Fe protein'
            idList = ['K02588']
            nzDict[id] = idList

            id = 'nifDK: nitrogenase Mo-Fe protein'
            idList = ['K02591', 'K02586']
            nzDict[id] = idList

            ### PO4 Solubility
            id = 'pqqC: pyrroloquinoline-quinone synthase'
            idList = ['K06137']
            nzDict[id] = idList

            ### Biocontrol
            id = 'hcnABC: glycine dehydrogenase (cyanide-forming)'
            idList = ['K10814', 'K10815', 'K10816']
            nzDict[id] = idList

            id = 'budA: acetolactate decarboxylase'
            idList = ['K01575']
            nzDict[id] = idList

            id = 'budC: (S,S)-butanediol dehydrogenase'
            idList = ['K18009']
            nzDict[id] = idList

            id = 'E3.2.1.6: endo-1,3(4)-beta-gulcanase'
            idList = ['K01180']
            nzDict[id] = idList

            id = 'E3.2.1.14: chitinase'
            idList = ['K01183']
            nzDict[id] = idList

            ### Root growth (drought/salt stress)
            id = 'ipdC: indolepyruvate decarboxylase'
            idList = ['K04103']
            nzDict[id] = idList

            id = 'acdS: 1-aminocyclopropane-1-carboxylate deaminase'
            idList = ['K01505']
            nzDict[id] = idList

            ### Siderophores
            id = 'iucC: aerobactin synthase'
            idList = ['K03895']
            nzDict[id] = idList

            id = 'entDF: enterobactin synthase'
            idList = ['K02362', 'K02364']
            nzDict[id] = idList

            id = 'mbtA: mycobactin salicyl-AMP ligase'
            idList = ['K04787']
            nzDict[id] = idList

            id = 'ybtE: yersiniabactin salicyl-AMP ligase'
            idList = ['K04783']
            nzDict[id] = idList

            ### C decomposition
            id = 'bglX: beta-glucosidase'
            idList = ['K05349']
            nzDict[id] = idList

            id = 'bglB: beta-glucosidase'
            idList = ['K05350']
            nzDict[id] = idList

            id = 'E3.2.1.21: beta-glucosidase'
            idList = ['K01188']
            nzDict[id] = idList

            ### N decomposition
            id = 'amiE: amidase'
            idList = ['K01426']
            nzDict[id] = idList

            id = 'ureC: urease alpha subunit'
            idList = ['K01428']
            nzDict[id] = idList

            ### P decomposition
            id = 'phoA: alkaline phosphatase'
            idList = ['K01077']
            nzDict[id] = idList

            id = 'phoD: alkaline phosphatase'
            idList = ['K01113']
            nzDict[id] = idList

            id = 'E3.1.3.2: acid phosphatase'
            idList = ['K01078']
            nzDict[id] = idList

            id = 'appA: acid phosphatase'
            idList = ['K01093']
            nzDict[id] = idList

            id = 'phoN: acid phosphatase'
            idList = ['K09474']
            nzDict[id] = idList

            ### S decomposition
            id = 'aslA: arylsulfatase'
            idList = ['K01130']
            nzDict[id] = idList

        elif nzAll == 6:
            ### N-fixation
            id = 'nifH: nitrogenase Fe protein'
            idList = ['K02588']
            nzDict[id] = idList

            id = 'nifDK: nitrogenase Mo-Fe protein'
            idList = ['K02591', 'K02586']
            nzDict[id] = idList

            ### Ammonification
            id = 'amiE: amidase'
            idList = ['K01426']
            nzDict[id] = idList

            id = 'ureC: urease alpha subunit'
            idList = ['K01428']
            nzDict[id] = idList

            ### Nitrification
            id = 'pmoA-amoA: methane/ammonia monooxygenase'
            idList = ['K10944']
            nzDict[id] = idList

            id = 'hao: hydroxylamine dehydrogenase'
            idList = ['K10535']
            nzDict[id] = idList

            id = 'narGH: nitrate reductase'
            idList = ['K00370', 'K00371']
            nzDict[id] = idList

            ### Dissimilatory nitrate reduction
            id = 'nirBD: nitrite reductase (NADH)'
            idList = ['K00362', 'K00363']
            nzDict[id] = idList

            id = 'nrfA: nitrite reductase (cytochrome c-552)'
            idList = ['K03385']
            nzDict[id] = idList

            ### Assimilatory nitrate reduction
            id = 'nirA: ferrodoxin-nitrite reductase'
            idList = ['K00366']
            nzDict[id] = idList

            id = 'NIT-6: nitrite reductase (NAD(P)H)'
            idList = ['K17877']
            nzDict[id] = idList

            ### Denitrification
            id = 'nirK: nitrite reductase (NO forming)'
            idList = ['K00368']
            nzDict[id] = idList

            id = 'nirS:  nitrite reductase (NO forming)'
            idList = ['K15864']
            nzDict[id] = idList

            id = 'norBC: nitric oxide reductase'
            idList = ['K04561', 'K02305']
            nzDict[id] = idList

            id = 'nosZ: nitrous-oxide reductase'
            idList = ['K00376']
            nzDict[id] = idList

        elif nzAll == 10:
            ### PO4 Solubility
            id = 'pqqC: pyrroloquinoline-quinone synthase'
            idList = ['K06137']
            nzDict[id] = idList

            ### Biocontrol
            id = 'hcnABC: glycine dehydrogenase (cyanide-forming)'
            idList = ['K10814', 'K10815', 'K10816']
            nzDict[id] = idList

            id = 'budA: acetolactate decarboxylase'
            idList = ['K01575']
            nzDict[id] = idList

            id ='budC: (S,S)-butanediol dehydrogenase'
            idList = ['K18009']
            nzDict[id] = idList

            id = 'E3.2.1.6: endo-1,3(4)-beta-gulcanase'
            idList = ['K01180']
            nzDict[id] = idList

            id = 'E3.2.1.14: chitinase'
            idList = ['K01183']
            nzDict[id] = idList

            ### Root growth (drought/salt stress)
            id = 'ipdC: indolepyruvate decarboxylase'
            idList = ['K04103']
            nzDict[id] = idList

            id = 'acdS: 1-aminocyclopropane-1-carboxylate deaminase'
            idList = ['K01505']
            nzDict[id] = idList

            ### Siderophores
            id = 'iucC: aerobactin synthase'
            idList = ['K03895']
            nzDict[id] = idList

            id = 'entDF: enterobactin synthase'
            idList = ['K02362', 'K02364']
            nzDict[id] = idList

            id = 'mbtA: mycobactin salicyl-AMP ligase'
            idList = ['K04787']
            nzDict[id] = idList

            id = 'ybtE: yersiniabactin salicyl-AMP ligase'
            idList = ['K04783']
            nzDict[id] = idList

            ### C decomposition
            id = 'bglX: beta-glucosidase'
            idList = ['K05349']
            nzDict[id] = idList

            id = 'bglB: beta-glucosidase'
            idList = ['K05350']
            nzDict[id] = idList

            id = 'E3.2.1.21: beta-glucosidase'
            idList = ['K01188']
            nzDict[id] = idList

            ### P decomposition
            id = 'phoA: alkaline phosphatase'
            idList = ['K01077']
            nzDict[id] = idList

            id = 'phoD: alkaline phosphatase'
            idList = ['K01113']
            nzDict[id] = idList

            id = 'E3.1.3.2: acid phosphatase'
            idList = ['K01078']
            nzDict[id] = idList

            id = 'appA: acid phosphatase'
            idList = ['K01093']
            nzDict[id] = idList

            id = 'phoN: acid phosphatase'
            idList = ['K09474']
            nzDict[id] = idList

            ### S decomposition
            id = 'aslA: arylsulfatase'
            idList = ['K01130']
            nzDict[id] = idList

            ### N-fixation
            id = 'nifH: nitrogenase Fe protein'
            idList = ['K02588']
            nzDict[id] = idList

            id = 'nifDK: nitrogenase Mo-Fe protein'
            idList = ['K02591', 'K02586']
            nzDict[id] = idList

            ### Ammonification
            id = 'amiE: amidase'
            idList = ['K01426']
            nzDict[id] = idList

            id = 'ureC: urease alpha subunit'
            idList = ['K01428']
            nzDict[id] = idList

            ### Nitrification
            id = 'pmoA-amoA: methane/ammonia monooxygenase'
            idList = ['K10944']
            nzDict[id] = idList

            id = 'hao: hydroxylamine dehydrogenase'
            idList = ['K10535']
            nzDict[id] = idList

            id = 'narGH: nitrate reductase'
            idList = ['K00370', 'K00371']
            nzDict[id] = idList

            ### Dissimilatory nitrate reduction
            id = 'nirBD: nitrite reductase (NADH)'
            idList = ['K00362', 'K00363']
            nzDict[id] = idList

            id = 'nrfA: nitrite reductase (cytochrome c-552)'
            idList = ['K03385']
            nzDict[id] = idList

            ### Assimilatory nitrate reduction
            id = 'nirA: ferrodoxin-nitrite reductase'
            idList = ['K00366']
            nzDict[id] = idList

            id = 'NIT-6: nitrite reductase (NAD(P)H)'
            idList = ['K17877']
            nzDict[id] = idList

            ### Denitrification
            id = 'nirK: nitrite reductase (NO forming)'
            idList = ['K00368']
            nzDict[id] = idList

            id = 'nirS:  nitrite reductase (NO forming)'
            idList = ['K15864']
            nzDict[id] = idList

            # 1.7.2.5  nitric oxide reductase (cytochrome c)
            id = 'norBC: nitric oxide reductase'
            idList = ['K04561', 'K02305']
            nzDict[id] = idList

            # 1.7.2.4  nitrous-oxide reductase
            id = 'nosZ: nitrous-oxide reductase'
            idList = ['K00376']
            nzDict[id] = idList

        # create sample and otu lists based on meta data selection
        wanted = ['sampleid', 'otuid', 'abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']
        profileDF = savedDF.loc[:, wanted]
        profileDF.set_index('otuid', drop=True, inplace=True)

        # get PICRUSt data for otu
        otuList = pd.unique(profileDF.index.ravel().tolist())
        qs = PICRUSt.objects.using('picrust').filter(otuid__in=otuList)
        piList = []
        for item in qs:
            piList.append([item.otuid_id, item.geneList])
        picrustDF = pd.DataFrame(piList, columns=['otuid', 'geneList'])
        picrustDF.set_index('otuid', drop=True, inplace=True)

        curStep = functions.getBase(RID)
        sumKEGG(picrustDF, nzDict, nzAll, RID, PID, stops)
        functions.setBase(RID, curStep)

        picrustDF.drop('geneList', axis=1, inplace=True)
        picrustDF.fillna(value=0, inplace=True)
        levelList = picrustDF.columns.values.tolist()

        # convert profile to index (sampleid) and columns (keggid) and values (depvar)
        profileDF.reset_index(drop=False, inplace=True)
        if DepVar == 0:
            profileDF = profileDF.pivot(index='otuid', columns='sampleid', values='abund')
        elif DepVar == 1:
            profileDF = profileDF.pivot(index='otuid', columns='sampleid', values='rel_abund')
        elif DepVar == 2:
            profileDF = profileDF.pivot(index='otuid', columns='sampleid', values='rich')
        elif DepVar == 3:
            profileDF = profileDF.pivot(index='otuid', columns='sampleid', values='diversity')
        elif DepVar == 4:
            profileDF = profileDF.pivot(index='otuid', columns='sampleid', values='abund_16S')

        sampleList = profileDF.columns.values.tolist()
        taxaDF = pd.DataFrame()
        total = len(sampleList)
        counter = 1
        for i in sampleList:
            tempDF = pd.DataFrame(index=profileDF.index)
            for j in levelList:
                tempDF[j] = profileDF[i] * picrustDF[j]

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            tempDF = tempDF[levelList].sum()
            tempDF['sampleid'] = i
            taxaDF = taxaDF.append(tempDF, ignore_index=True)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            functions.setBase(RID, 'Calculating KEGG pathway abundances...sample ' + str(counter) + ' out of ' + str(total) + ' is finished!')
            counter += 1

        functions.setBase(RID, curStep)

        # df of mapped taxa to selected kegg orthologies
        if mapTaxa == 'yes':
            namesDict = {}
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
                else:
                    print str(level) + ' not found in database'
                namesDict[level] = name

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            allDF = picrustDF.reset_index(drop=False, inplace=False)
            allDF.rename(columns={'index': 'otuid'}, inplace=True)
            allDF.dropna(axis=1, how='any', inplace=False)
            allDF = allDF[(allDF.sum(axis=1) != 0)]

            for i in levelList:
                allDF[i] = allDF[i] / allDF[i]
            allDF.fillna(value=0, inplace=True)
            recordDict = {}
            for id in allDF.otuid.tolist():
                try:
                    qs = OTU_99.objects.all().filter(otuid=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyid_id__familyName', 'genusid_id__genusName', 'speciesid_id__speciesName', 'otuName')
                    record = '|'.join(qs[0])
                    recordDict[id] = record
                except:
                    recordDict[id] = 'No data'
            allDF['Taxonomy'] = allDF['otuid'].map(recordDict)
            order = ['otuid', 'Taxonomy'] + levelList
            allDF = allDF[order]
            allDF.rename(columns=namesDict, inplace=True)
        else:
            allDF = ''

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        #taxaDF.reset_index(drop=False, inplace=True)
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

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        if nzAll < 5:
            idList = functions.getFullNZ(list(finalDF.rank_id.unique()))
            finalDF['rank_name'] = finalDF['rank_id'].map(idList)
        else:
            finalDF['rank_name'] = finalDF['rank_id']

        # required to match getTaxaDF output
        metaDF.reset_index(drop=False, inplace=True)

        return finalDF, allDF

    except Exception:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "There was an error with your analysis!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = json.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def sumKEGG(picrustDF, keggDict, nzAll, RID, PID, stops):
    total = len(keggDict)
    counter = 0

    for key in keggDict:
        pathList = keggDict[key]

        functions.setBase(RID, 'Mapping phylotypes to KEGG pathways...pathway/enzyme ' +
                          str(counter) + ' out of ' + str(total) + ' is finished!')

        for row in zip(picrustDF.index.values, picrustDF['geneList']):
            if nzAll >= 5:
                if all(i in row[1] for i in pathList):
                    picrustDF.at[row[0], key] = 1.0
            else:
                if any(i in row[1] for i in pathList):
                    picrustDF.at[row[0], key] = 1.0

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

        counter += 1


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

        record = '|'.join(qs[0])
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
                    if ko_entry.objects.using('picrust').all().filter(ko_lvl4_id=id).exists():
                        qs = ko_entry.objects.using('picrust').all().filter(ko_lvl4_id=id).values_list('ko_lvl1_id_id__ko_lvl1_name', 'ko_lvl2_id_id__ko_lvl2_name', 'ko_lvl3_id_id__ko_lvl3_name', 'ko_desc')
                    else:
                        qs = ('Not found', )

        record = '|'.join(qs[0])
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
                        if nz_entry.objects.using('picrust').all().filter(nz_lvl5_id=id).exists():
                            qs = nz_entry.objects.using('picrust').all().filter(nz_lvl5_id=id).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_id_id__nz_lvl3_name', 'nz_lvl4_id_id__nz_lvl4_name', 'nz_desc')
                        else:
                            qs = ('Not found', )

        record = '|'.join(qs[0])
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
    if level == 1:
        myLevel = 'kingdomName'
        myID = 'kingdomid'
    elif level == 2:
        myLevel = 'phylaName'
        myID = 'phylaid'
    elif level == 3:
        myLevel = 'className'
        myID = 'classid'
    elif level == 4:
        myLevel = 'orderName'
        myID = 'orderid'
    elif level == 5:
        myLevel = 'familyName'
        myID = 'familyid'
    elif level == 6:
        myLevel = 'genusName'
        myID = 'genusid'
    elif level == 7:
        myLevel = 'speciesName'
        myID = 'speciesid'
    elif level == 9:
        myLevel = 'otuName'
        myID = 'otuid'

    if remUnclass == 'yes':
        # check if selecting based on level first, create else statement for tree usage
        savedDF = savedDF[~savedDF[myLevel].str.contains('unclassified')]

    if remZeroes == 'yes' and perZeroes > 0:
        threshold = int(perZeroes / 100.0 * numSamples)
        bytag = savedDF.groupby(myID).aggregate(np.count_nonzero)
        tags = bytag[bytag[myVar] >= threshold].index.tolist()
        savedDF = savedDF[savedDF[myID].isin(tags)]

    if filterData == 'yes' and not level == 8:
        threshold = int(filterPer)
        if filterMeth == 1:
            pass
        elif filterMeth == 2:
            sum = savedDF.groupby([myID, 'sampleid'])[myVar].sum()
            sumDF = sum.reset_index(drop=False)
            bytag_q3 = sumDF.groupby(myID)[myVar].quantile(0.75)
            bytag_q1 = sumDF.groupby(myID)[myVar].quantile(0.25)
            bytag = bytag_q3 - bytag_q1
            bytag.sort(axis=0, ascending=False, inplace=True)
            tags = bytag[:threshold].index.tolist()
            savedDF = savedDF[savedDF[myID].isin(tags)]
        elif filterMeth == 3:
            sum = savedDF.groupby([myID, 'sampleid'])[myVar].sum()
            sumDF = sum.reset_index(drop=False)
            bytag_sd = sumDF.groupby(myID)[myVar].std()
            bytag_mean = sumDF.groupby(myID)[myVar].mean()
            bytag = bytag_sd / bytag_mean
            bytag.sort(axis=0, ascending=False, inplace=True)
            tags = bytag[:threshold].index.tolist()
            savedDF = savedDF[savedDF[myID].isin(tags)]
        elif filterMeth == 4:
            sum = savedDF.groupby([myID, 'sampleid'])[myVar].sum()
            sumDF = sum.reset_index(drop=False)
            bytag = sumDF.groupby(myID)[myVar].std()
            bytag.sort(axis=0, ascending=False, inplace=True)
            tags = bytag[:threshold].index.tolist()
            savedDF = savedDF[savedDF[myID].isin(tags)]
        elif filterMeth == 5:
            sum = savedDF.groupby([myID, 'sampleid'])[myVar].sum()
            sumDF = sum.reset_index(drop=False)
            bytag = sumDF.groupby(myID)[myVar].mean()
            bytag.sort(axis=0, ascending=False, inplace=True)
            tags = bytag[:threshold].index.tolist()
            savedDF = savedDF[savedDF[myID].isin(tags)]
        elif filterMeth == 6:
            sum = savedDF.groupby([myID, 'sampleid'])[myVar].sum()
            sumDF = sum.reset_index(drop=False)
            bytag = sumDF.groupby(myID)[myVar].median()
            bytag.sort(axis=0, ascending=False, inplace=True)
            tags = bytag[:threshold].index.tolist()
            savedDF = savedDF[savedDF[myID].isin(tags)]

    return savedDF


