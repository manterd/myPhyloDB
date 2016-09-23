import ast
import datetime
from django import db
from django.http import HttpResponse
from django_pandas.io import read_frame
import logging
import math
import numpy as np
import pandas as pd
import simplejson

from database.models import Kingdom, Phyla, Class, Order, Family, Genus, Species
from database.models import PICRUSt
from database.models import ko_lvl1, ko_lvl2, ko_lvl3, ko_entry
from database.models import nz_lvl1, nz_lvl2, nz_lvl3, nz_lvl4, nz_entry


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getTaxaDF(finalSampleList, abund, selectAll, taxaDict, savedDF, metaDF, allFields, DepVar, RID, stops, PID):
    try:
        missingList = []
        taxaDF = pd.DataFrame(columns=['sampleid', 'rank', 'rank_id', 'rank_name', abund, 'abund_16S', 'rich', 'diversity'])
        if selectAll == 0:
            for key in taxaDict:
                taxaList = taxaDict[key]
                if isinstance(taxaList, str):
                    if key == 'Kingdom':
                        blankDF = pd.DataFrame(index=finalSampleList)
                        blankDF.loc[:, 'rank'] = 'kingdom'
                        blankDF.loc[:, 'rank_id'] = taxaList
                        blankDF.loc[:, 'rank_name'] = Kingdom.objects.get(kingdomid=taxaList).kingdomName
                        tempDF = savedDF.loc[savedDF['kingdomid'] == taxaList]
                        tempDF = tempDF[['sampleid', abund, 'abund_16S', 'rich', 'diversity']]
                        tempDF.set_index('sampleid', inplace=True)
                        tempDF = pd.merge(blankDF, tempDF, left_index=True, right_index=True, how='outer')
                        tempDF.fillna(0, inplace=True)
                        tempDF.reset_index(inplace=True)
                        tempDF.rename(columns={'index': 'sampleid'}, inplace=True)
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Phyla':
                        blankDF = pd.DataFrame(index=finalSampleList)
                        blankDF.loc[:, 'rank'] = 'phyla'
                        blankDF.loc[:, 'rank_id'] = taxaList
                        blankDF.loc[:, 'rank_name'] = Phyla.objects.get(phylaid=taxaList).phylaName
                        tempDF = savedDF.loc[savedDF['phylaid'] == taxaList]
                        tempDF = tempDF[['sampleid', abund, 'abund_16S', 'rich', 'diversity']]
                        tempDF.set_index('sampleid', inplace=True)
                        tempDF = pd.merge(blankDF, tempDF, left_index=True, right_index=True, how='outer')
                        tempDF.fillna(0, inplace=True)
                        tempDF.reset_index(inplace=True)
                        tempDF.rename(columns={'index': 'sampleid'}, inplace=True)
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Class':
                        blankDF = pd.DataFrame(index=finalSampleList)
                        blankDF.loc[:, 'rank'] = 'class'
                        blankDF.loc[:, 'rank_id'] = taxaList
                        blankDF.loc[:, 'rank_name'] = Class.objects.get(classid=taxaList).className
                        tempDF = savedDF.loc[savedDF['classid'] == taxaList]
                        tempDF = tempDF[['sampleid', abund, 'abund_16S', 'rich', 'diversity']]
                        tempDF.set_index('sampleid', inplace=True)
                        tempDF = pd.merge(blankDF, tempDF, left_index=True, right_index=True, how='outer')
                        tempDF.fillna(0, inplace=True)
                        tempDF.reset_index(inplace=True)
                        tempDF.rename(columns={'index': 'sampleid'}, inplace=True)
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Order':
                        blankDF = pd.DataFrame(index=finalSampleList)
                        blankDF.loc[:, 'rank'] = 'order'
                        blankDF.loc[:, 'rank_id'] = taxaList
                        blankDF.loc[:, 'rank_name'] = Order.objects.get(orderid=taxaList).orderName
                        tempDF = savedDF.loc[savedDF['orderid'] == taxaList]
                        tempDF = tempDF[['sampleid', abund, 'abund_16S', 'rich', 'diversity']]
                        tempDF.set_index('sampleid', inplace=True)
                        tempDF = pd.merge(blankDF, tempDF, left_index=True, right_index=True, how='outer')
                        tempDF.fillna(0, inplace=True)
                        tempDF.reset_index(inplace=True)
                        tempDF.rename(columns={'index': 'sampleid'}, inplace=True)
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Family':
                        blankDF = pd.DataFrame(index=finalSampleList)
                        blankDF.loc[:, 'rank'] = 'family'
                        blankDF.loc[:, 'rank_id'] = taxaList
                        blankDF.loc[:, 'rank_name'] = Family.objects.get(familyid=taxaList).familyName
                        tempDF = savedDF.loc[savedDF['familyid'] == taxaList]
                        tempDF = tempDF[['sampleid', abund, 'abund_16S', 'rich', 'diversity']]
                        tempDF.set_index('sampleid', inplace=True)
                        tempDF = pd.merge(blankDF, tempDF, left_index=True, right_index=True, how='outer')
                        tempDF.fillna(0, inplace=True)
                        tempDF.reset_index(inplace=True)
                        tempDF.rename(columns={'index': 'sampleid'}, inplace=True)
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Genus':
                        blankDF = pd.DataFrame(index=finalSampleList)
                        blankDF.loc[:, 'rank'] = 'genus'
                        blankDF.loc[:, 'rank_id'] = taxaList
                        blankDF.loc[:, 'rank_name'] = Genus.objects.get(genusid=taxaList).genusName
                        tempDF = savedDF.loc[savedDF['genusid'] == taxaList]
                        tempDF = tempDF[['sampleid', abund, 'abund_16S', 'rich', 'diversity']]
                        tempDF.set_index('sampleid', inplace=True)
                        tempDF = pd.merge(blankDF, tempDF, left_index=True, right_index=True, how='outer')
                        tempDF.fillna(0, inplace=True)
                        tempDF.reset_index(inplace=True)
                        tempDF.rename(columns={'index': 'sampleid'}, inplace=True)
                        taxaDF = taxaDF.append(tempDF)
                    elif key == 'Species':
                        blankDF = pd.DataFrame(index=finalSampleList)
                        blankDF.loc[:, 'rank'] = 'species'
                        blankDF.loc[:, 'rank_id'] = taxaList
                        blankDF.loc[:, 'rank_name'] = Species.objects.get(genusid=taxaList).speciesName
                        tempDF = savedDF.loc[savedDF['speciesid'] == taxaList]
                        tempDF = tempDF[['sampleid', abund, 'abund_16S', 'rich', 'diversity']]
                        tempDF.set_index('sampleid', inplace=True)
                        tempDF = pd.merge(blankDF, tempDF, left_index=True, right_index=True, how='outer')
                        tempDF.fillna(0, inplace=True)
                        tempDF.reset_index(inplace=True)
                        tempDF.rename(columns={'index': 'sampleid'}, inplace=True)
                        taxaDF = taxaDF.append(tempDF)
                else:
                    if key == 'Kingdom':
                        for item in taxaList:
                            blankDF = pd.DataFrame(index=finalSampleList)
                            blankDF.loc[:, 'rank'] = 'kingdom'
                            blankDF.loc[:, 'rank_id'] = item
                            blankDF.loc[:, 'rank_name'] = Kingdom.objects.get(kingdomid=item).kingdomName
                            tempDF = savedDF.loc[savedDF['kingdomid'] == item]
                            tempDF = tempDF[['sampleid', abund, 'abund_16S', 'rich', 'diversity']]
                            tempDF.set_index('sampleid', inplace=True)
                            tempDF = pd.merge(blankDF, tempDF, left_index=True, right_index=True, how='outer')
                            tempDF.fillna(0, inplace=True)
                            tempDF.reset_index(inplace=True)
                            tempDF.rename(columns={'index': 'sampleid'}, inplace=True)
                            taxaDF = taxaDF.append(tempDF)
                    elif key == 'Phyla':
                        for item in taxaList:
                            blankDF = pd.DataFrame(index=finalSampleList)
                            blankDF.loc[:, 'rank'] = 'phyla'
                            blankDF.loc[:, 'rank_id'] = item
                            blankDF.loc[:, 'rank_name'] = Phyla.objects.get(phylaid=item).phylaName
                            tempDF = savedDF.loc[savedDF['phylaid'] == item]
                            tempDF = tempDF[['sampleid', abund, 'abund_16S', 'rich', 'diversity']]
                            tempDF.set_index('sampleid', inplace=True)
                            tempDF = pd.merge(blankDF, tempDF, left_index=True, right_index=True, how='outer')
                            tempDF.fillna(0, inplace=True)
                            tempDF.reset_index(inplace=True)
                            tempDF.rename(columns={'index': 'sampleid'}, inplace=True)
                            taxaDF = taxaDF.append(tempDF)
                    elif key == 'Class':
                        for item in taxaList:
                            blankDF = pd.DataFrame(index=finalSampleList)
                            blankDF.loc[:, 'rank'] = 'class'
                            blankDF.loc[:, 'rank_id'] = item
                            blankDF.loc[:, 'rank_name'] = Class.objects.get(classid=item).className
                            tempDF = savedDF.loc[savedDF['classid'] == item]
                            tempDF = tempDF[['sampleid', abund, 'abund_16S', 'rich', 'diversity']]
                            tempDF.set_index('sampleid', inplace=True)
                            tempDF = pd.merge(blankDF, tempDF, left_index=True, right_index=True, how='outer')
                            tempDF.fillna(0, inplace=True)
                            tempDF.reset_index(inplace=True)
                            tempDF.rename(columns={'index': 'sampleid'}, inplace=True)
                            taxaDF = taxaDF.append(tempDF)
                    elif key == 'Order':
                        for item in taxaList:
                            blankDF = pd.DataFrame(index=finalSampleList)
                            blankDF.loc[:, 'rank'] = 'order'
                            blankDF.loc[:, 'rank_id'] = item
                            blankDF.loc[:, 'rank_name'] = Order.objects.get(orderid=item).orderName
                            tempDF = savedDF.loc[savedDF['orderid'] == item]
                            tempDF = tempDF[['sampleid', abund, 'abund_16S', 'rich', 'diversity']]
                            tempDF.set_index('sampleid', inplace=True)
                            tempDF = pd.merge(blankDF, tempDF, left_index=True, right_index=True, how='outer')
                            tempDF.fillna(0, inplace=True)
                            tempDF.reset_index(inplace=True)
                            tempDF.rename(columns={'index': 'sampleid'}, inplace=True)
                            taxaDF = taxaDF.append(tempDF)
                    elif key == 'Family':
                        for item in taxaList:
                            blankDF = pd.DataFrame(index=finalSampleList)
                            blankDF.loc[:, 'rank'] = 'family'
                            blankDF.loc[:, 'rank_id'] = item
                            blankDF.loc[:, 'rank_name'] = Family.objects.get(familyid=item).familyName
                            tempDF = savedDF.loc[savedDF['familyid'] == item]
                            tempDF = tempDF[['sampleid', abund, 'abund_16S', 'rich', 'diversity']]
                            tempDF.set_index('sampleid', inplace=True)
                            tempDF = pd.merge(blankDF, tempDF, left_index=True, right_index=True, how='outer')
                            tempDF.fillna(0, inplace=True)
                            tempDF.reset_index(inplace=True)
                            tempDF.rename(columns={'index': 'sampleid'}, inplace=True)
                            taxaDF = taxaDF.append(tempDF)
                    elif key == 'Genus':
                        for item in taxaList:
                            blankDF = pd.DataFrame(index=finalSampleList)
                            blankDF.loc[:, 'rank'] = 'genus'
                            blankDF.loc[:, 'rank_id'] = item
                            blankDF.loc[:, 'rank_name'] = Genus.objects.get(genusid=item).genusName
                            tempDF = savedDF.loc[savedDF['genusid'] == item]
                            tempDF = tempDF[['sampleid', abund, 'abund_16S', 'rich', 'diversity']]
                            tempDF.set_index('sampleid', inplace=True)
                            tempDF = pd.merge(blankDF, tempDF, left_index=True, right_index=True, how='outer')
                            tempDF.fillna(0, inplace=True)
                            tempDF.reset_index(inplace=True)
                            tempDF.rename(columns={'index': 'sampleid'}, inplace=True)
                            taxaDF = taxaDF.append(tempDF)
                    elif key == 'Species':
                        for item in taxaList:
                            blankDF = pd.DataFrame(index=finalSampleList)
                            blankDF.loc[:, 'rank'] = 'species'
                            blankDF.loc[:, 'rank_id'] = item
                            blankDF.loc[:, 'rank_name'] = Species.objects.get(speciesid=item).speciesName
                            tempDF = savedDF.loc[savedDF['speciesid'] == item]
                            tempDF = tempDF[['sampleid', abund, 'abund_16S', 'rich', 'diversity']]
                            tempDF.set_index('sampleid', inplace=True)
                            tempDF = pd.merge(blankDF, tempDF, left_index=True, right_index=True, how='outer')
                            tempDF.fillna(0, inplace=True)
                            tempDF.reset_index(inplace=True)
                            tempDF.rename(columns={'index': 'sampleid'}, inplace=True)
                            taxaDF = taxaDF.append(tempDF)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    res = ''
                    return HttpResponse(res, content_type='application/json')
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        #TODO get a list of taxa and then fix like above
        elif selectAll == 1:
            taxaDF = savedDF.loc[:, ['sampleid', 'kingdomid', 'kingdomName', abund, 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'kingdomid': 'rank_id', 'kingdomName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Kingdom'
        elif selectAll == 2:
            taxaDF = savedDF.loc[:, ['sampleid', 'phylaid', 'phylaName', abund, 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'phylaid': 'rank_id', 'phylaName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Phyla'
        elif selectAll == 3:
            taxaDF = savedDF.loc[:, ['sampleid', 'classid', 'className', abund, 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'classid': 'rank_id', 'className': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Class'
        elif selectAll == 4:
            taxaDF = savedDF.loc[:, ['sampleid', 'orderid', 'orderName', abund, 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'orderid': 'rank_id', 'orderName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Order'
        elif selectAll == 5:
            taxaDF = savedDF.loc[:, ['sampleid', 'familyid', 'familyName', abund, 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'familyid': 'rank_id', 'familyName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Family'
        elif selectAll == 6:
            taxaDF = savedDF.loc[:, ['sampleid', 'genusid', 'genusName', abund, 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'genusid': 'rank_id', 'genusName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Genus'
        elif selectAll == 7:
            taxaDF = savedDF.loc[:, ['sampleid', 'speciesid', 'speciesName', abund, 'abund_16S', 'rich', 'diversity']]
            taxaDF.rename(columns={'speciesid': 'rank_id', 'speciesName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Species'
        elif selectAll == 8:
            taxaDF = savedDF.loc[:, ['sampleid', 'genusid', 'genusName', abund, 'abund_16S', 'rich', 'diversity']]
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

        metaDF.reset_index(inplace=True, drop=False)
        finalDF = pd.merge(metaDF, taxaDF, left_on='sampleid', right_on='sampleid', how='outer')

        wantedList = allFields + ['sampleid', 'rank', 'rank_name', 'rank_id']
        if DepVar == 1:
            finalDF = finalDF.groupby(wantedList)[[abund]].sum()
        elif DepVar == 4:
            finalDF = finalDF.groupby(wantedList)[['abund_16S']].sum()
            finalDF['abund_16S'] = finalDF['abund_16S'].astype(int)
        elif DepVar == 2:
            finalDF = finalDF.groupby(wantedList)[['rich']].sum()
        elif DepVar == 3:
            finalDF = finalDF.groupby(wantedList)[['diversity']].sum()

        finalDF.reset_index(drop=False, inplace=True)
        return finalDF, missingList

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "There was an error with your analysis!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getKeggDF(abund, keggAll, keggDict, savedDF, tempDF, allFields, DepVar, RID, stops, PID):
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

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        # create sample and species lists based on meta data selection
        wanted = ['sampleid', 'speciesid', abund, 'abund_16S']
        profileDF = tempDF.loc[:, wanted]
        profileDF.set_index('speciesid', inplace=True)

        # get PICRUSt data for species
        speciesList = pd.unique(profileDF.index.ravel().tolist())
        qs = PICRUSt.objects.using('picrust').filter(speciesid__in=speciesList)
        picrustDF = read_frame(qs, fieldnames=['speciesid__speciesid', 'geneCount'])
        picrustDF.set_index('speciesid__speciesid', inplace=True)

        levelList = []
        for key in koDict:
            levelList.append(str(key))

        picrustDF = pd.concat([picrustDF, pd.DataFrame(columns=levelList)])
        picrustDF.fillna(0.0, inplace=True)
        picrustDF = sumKEGG(speciesList, picrustDF, koDict, RID, PID, stops)
        picrustDF.drop('geneCount', axis=1, inplace=True)
        picrustDF[picrustDF > 0.0] = 1.0

        # merge to get final gene counts for all selected samples
        taxaDF = pd.merge(profileDF, picrustDF, left_index=True, right_index=True, how='inner')

        for level in levelList:
            if DepVar == 1:
                taxaDF[level] = taxaDF[abund] * taxaDF[level]
            elif DepVar == 2:
                taxaDF[level] = np.where(taxaDF[abund] * taxaDF[level] > 0, 1, 0)
            elif DepVar == 3:
                taxaDF[level] = taxaDF[abund] * taxaDF[level]
                taxaDF[level] = taxaDF[level].div(taxaDF[level].sum(), axis=0)
                taxaDF[level] = taxaDF[level].apply(lambda x: -1 * x * math.log(x) if x > 0 else 0)
            elif DepVar == 4:
                taxaDF[level] = taxaDF['abund_16S'] * taxaDF[level]

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        taxaDF = taxaDF.groupby('sampleid')[levelList].agg('sum')
        taxaDF.reset_index(drop=False, inplace=True)

        if DepVar == 1:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name=abund)
        elif DepVar == 2:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rich')
        elif DepVar == 3:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='diversity')
        elif DepVar == 4:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund_16S')

        wanted = allFields + ['sampleid']
        metaDF = savedDF.loc[:, wanted]
        metaDF.set_index('sampleid', drop=True, inplace=True)
        grouped = metaDF.groupby(level=0)
        metaDF = grouped.last()

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

        return finalDF

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "There was an error with your analysis!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def getNZDF(abund, nzAll, myDict, savedDF, tempDF, allFields, DepVar, RID, stops, PID):
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
            # 1.18.6.1  nitrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.18.6.1  nitrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.3.3.11  pyrroloquinoline-quinone synthase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.3.3.11  pyrroloquinoline-quinone synthase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.4.99.5  glycine dehydrogenase (cyanide-forming)
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.4.99.5  glycine dehydrogenase (cyanide-forming)').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 1.1.1.76  (S,S)-butanediol dehydrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.1.1.76  (S,S)-butanediol dehydrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.2.1.14  chitinase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='3.2.1.14  chitinase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 4.1.1.74  indolepyruvate decarboxylase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='4.1.1.74  indolepyruvate decarboxylase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='3.5.99.7  1-aminocyclopropane-1-carboxylate deaminase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 6.3.2.39  aerobactin synthase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='6.3.2.39  aerobactin synthase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

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

            # 3.5.1.4  amidase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.5.1.4  amidase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.5.1.5  urease
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.5.1.5  urease').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.1.3.1  alkaline phosphatase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.1.3.1  alkaline phosphatase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.1.3.2  acid phosphatase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.1.3.2  acid phosphatase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.1.6.1  arylsulfatase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name__startswith='3.1.6.1  arylsulfatase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

        elif nzAll == 6:
            # 1.18.6.1  nitrogenase
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='1.18.6.1  nitrogenase').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

            # 3.5.1.5  urease
            id = nz_lvl4.objects.using('picrust').get(nz_lvl4_name='3.5.1.5  urease').nz_lvl4_id
            idList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=id).values_list('nz_orthology', flat=True)
            nzDict[id] = idList

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

        # create sample and species lists based on meta data selection
        wanted = ['sampleid', 'speciesid', abund, 'abund_16S']
        profileDF = tempDF.loc[:, wanted]
        profileDF.set_index('speciesid', inplace=True)

        # get PICRUSt data for species
        speciesList = pd.unique(profileDF.index.ravel().tolist())
        qs = PICRUSt.objects.using('picrust').filter(speciesid__in=speciesList)
        picrustDF = read_frame(qs, fieldnames=['speciesid__speciesid', 'geneCount'], verbose=False)
        picrustDF.set_index('speciesid__speciesid', inplace=True)

        levelList = []
        for key in nzDict:
            levelList.append(str(key))

        picrustDF = pd.concat([picrustDF, pd.DataFrame(columns=levelList)])
        picrustDF.fillna(0.0, inplace=True)
        picrustDF = sumKEGG(speciesList, picrustDF, nzDict, RID, PID, stops)
        picrustDF.drop('geneCount', axis=1, inplace=True)
        picrustDF[picrustDF > 0.0] = 1.0

        # merge to get final gene counts for all selected samples
        taxaDF = pd.merge(profileDF, picrustDF, left_index=True, right_index=True, how='inner')
        for level in levelList:
            if DepVar == 1:
                taxaDF[level] = taxaDF[abund] * taxaDF[level]
            elif DepVar == 2:
                taxaDF[level] = np.where(taxaDF[abund] * taxaDF[level] > 0, 1, 0)
            elif DepVar == 3:
                taxaDF[level] = taxaDF[abund] * taxaDF[level]
                taxaDF[level] = taxaDF[level].div(taxaDF[level].sum(), axis=0)
                taxaDF[level] = taxaDF[level].apply(lambda x: -1 * x * math.log(x) if x > 0 else 0)
            elif DepVar == 4:
                taxaDF[level] = taxaDF['abund_16S'] * taxaDF[level]

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        taxaDF = taxaDF.groupby('sampleid')[levelList].agg('sum')
        taxaDF.reset_index(drop=False, inplace=True)

        if DepVar == 1:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name=abund)
        elif DepVar == 2:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='rich')
        elif DepVar == 3:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='diversity')
        elif DepVar == 4:
            taxaDF = pd.melt(taxaDF, id_vars='sampleid', var_name='rank_id', value_name='abund_16S')

        wanted = allFields + ['sampleid']
        metaDF = savedDF.loc[:, wanted]
        metaDF.set_index('sampleid', drop=True, inplace=True)
        grouped = metaDF.groupby(level=0)
        metaDF = grouped.last()

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
                finalDF.loc[index, 'rank_name'] = nz_entry.objects.using('picrust').get(nz_lvl5_id=row['rank_id']).nz_name

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        return finalDF

    except:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "There was an error with your analysis!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = simplejson.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def sumKEGG(speciesList, picrustDF, keggDict, RID, PID, stops):
    db.close_old_connections()
    for species in speciesList:
        try:
            cell = picrustDF.at[species, 'geneCount']
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

                picrustDF.at[species, key] = sum
        except:
            pass

    return picrustDF


def getFullTaxonomy(idList):
    recordList = []
    for id in idList:
        if Phyla.objects.all().filter(phylaid=id).exists():
            qs = Phyla.objects.all().filter(phylaid=id).values_list('kingdomid_id__kingdomName', 'phylaName')
            record = qs[0] + ('N/A', 'N/A', 'N/A', 'N/A', 'N/A',)
            recordList.append(record)
        else:
            if Class.objects.all().filter(classid=id).exists():
                qs = Class.objects.all().filter(classid=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'className')
                record = qs[0] + ('N/A', 'N/A', 'N/A', 'N/A',)
                recordList.append(record)
            else:
                if Order.objects.all().filter(orderid=id).exists():
                    qs = Order.objects.all().filter(orderid=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderName')
                    record = qs[0] + ('N/A', 'N/A', 'N/A',)
                    recordList.append(record)
                else:
                    if Family.objects.all().filter(familyid=id).exists():
                        qs = Family.objects.all().filter(familyid=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyName')
                        record = qs[0] + ('N/A', 'N/A',)
                        recordList.append(record)
                    else:
                        if Genus.objects.all().filter(genusid=id).exists():
                            qs = Genus.objects.all().filter(genusid=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyid_id__familyName', 'genusName')
                            record = qs[0] + ('N/A',)
                            recordList.append(record)
                        else:
                            if Species.objects.all().filter(speciesid=id).exists():
                                qs = Species.objects.all().filter(speciesid=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyid_id__familyName', 'genusid_id__genusName', 'speciesName')
                                record = qs[0]
                                recordList.append(record)
                            else:
                                record = ('Not found', 'Not found', 'Not found', 'Not found', 'Not found', 'Not found', 'Not found',)
                                recordList.append(record)

    return recordList


def getFullKO(idList):
    recordList = []

    for id in idList:
            if ko_lvl1.objects.using('picrust').all().filter(ko_lvl1_id=id).exists():
                qs = ko_lvl1.objects.using('picrust').all().filter(ko_lvl1_id=id).values_list('ko_lvl1_name')
                record = qs[0] + ('N/A', 'N/A')
                recordList.append(record)
            else:
                if ko_lvl2.objects.using('picrust').all().filter(ko_lvl2_id=id).exists():
                    qs = ko_lvl2.objects.using('picrust').all().filter(ko_lvl2_id=id).values_list('ko_lvl1_id_id__ko_lvl1_name', 'ko_lvl2_name')
                    record = qs[0] + ('N/A',)
                    recordList.append(record)
                else:
                    if ko_lvl3.objects.using('picrust').all().filter(ko_lvl3_id=id).exists():
                        qs = ko_lvl3.objects.using('picrust').all().filter(ko_lvl3_id=id).values_list('ko_lvl1_id_id__ko_lvl1_name', 'ko_lvl2_id_id__ko_lvl2_name', 'ko_lvl3_name')
                        record = qs[0]
                        recordList.append(record)
                    else:
                        record = ('Not found', 'Not found', 'Not found',)
                        recordList.append(record)

    return recordList


def getFullNZ(idList):
    recordList = []

    for id in idList:
            if nz_lvl1.objects.using('picrust').all().filter(nz_lvl1_id=id).exists():
                qs = nz_lvl1.objects.using('picrust').all().filter(nz_lvl1_id=id).values_list('nz_lvl1_name')
                record = qs[0] + ('N/A', 'N/A', 'N/A',)
                recordList.append(record)
            else:
                if nz_lvl2.objects.using('picrust').all().filter(nz_lvl2_id=id).exists():
                    qs = nz_lvl2.objects.using('picrust').all().filter(nz_lvl2_id=id).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_name')
                    record = qs[0] + ('N/A', 'N/A',)
                    recordList.append(record)
                else:
                    if nz_lvl3.objects.using('picrust').all().filter(nz_lvl3_id=id).exists():
                        qs = nz_lvl3.objects.using('picrust').all().filter(nz_lvl3_id=id).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_name')
                        record = qs[0] + ('N/A',)
                        recordList.append(record)
                    else:
                        if nz_lvl4.objects.using('picrust').all().filter(nz_lvl4_id=id).exists():
                            qs = nz_lvl4.objects.using('picrust').all().filter(nz_lvl4_id=id).values_list('nz_lvl1_id_id__nz_lvl1_name', 'nz_lvl2_id_id__nz_lvl2_name', 'nz_lvl3_id_id__nz_lvl3_name', 'nz_lvl4_name')
                            record = qs[0]
                            recordList.append(record)
                        else:
                            record = ('Not found', 'Not found', 'Not found', 'Not found',)
                            recordList.append(record)

    return recordList


def insertTaxaInfo(button3, zipped, DF, pos=1):
    if button3 == 1:
        k, p, c, o, f, g, s = map(None, *zipped)
        DF.insert(pos, 'Kingdom', k)
        DF.insert(pos+1, 'Phyla', p)
        DF.insert(pos+2, 'Class', c)
        DF.insert(pos+3, 'Order', o)
        DF.insert(pos+4, 'Family', f)
        DF.insert(pos+5, 'Genus', g)
        DF.insert(pos+6, 'Species', s)
    elif button3 == 2:
        L1, L2, L3 = map(None, *zipped)
        DF.insert(pos, 'Level_1', L1)
        DF.insert(pos+1, 'Level_2', L2)
        DF.insert(pos+2, 'Level_3', L3)
    elif button3 == 3:
        L1, L2, L3, L4 = map(None, *zipped)
        DF.insert(pos, 'Level_1', L1)
        DF.insert(pos+1, 'Level_2', L2)
        DF.insert(pos+2, 'Level_3', L3)
        DF.insert(pos+3, 'Level_4', L4)

    DF.fillna(value='N/A', inplace=True)
