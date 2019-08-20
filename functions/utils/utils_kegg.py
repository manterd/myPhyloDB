import datetime
from django.http import HttpResponse
import logging
import numpy as np
import pandas as pd
import json
import time

from database.models import Kingdom, Phyla, Class, Order, Family, Genus, Species, OTU_99, \
    PICRUSt, \
    ko_lvl1, ko_lvl2, ko_lvl3, ko_entry, \
    nz_lvl1, nz_lvl2, nz_lvl3, nz_lvl4, nz_entry

import functions
from functions.utils.debug import debug


LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getTaxaDF(selectAll, taxaDict, savedDF, metaDF, allFields, DepVar, RID, stops, PID, soilHealth=False):
    try:
        debug("getTaxaDF")
        Cols = ['sampleid', 'rank', 'rank_id', 'rank_name']
        Dep = []
        if DepVar == 0:
            Dep = ['abund']
        elif DepVar == 1:
            Dep = ['rel_abund']
        elif DepVar == 2:
            Dep = ['rich']
        elif DepVar == 3:
            Dep = ['diversity']
        elif DepVar == 4:
            Dep = ['abund_16S']
        elif DepVar == 10:
            Dep = ['abund', 'rel_abund', 'abund_16S', 'rich', 'diversity']
        Cols.extend(Dep)
        debug("getTaxaDF: missingList")
        missingList = []
        taxaDF = pd.DataFrame(columns=Cols)
        if selectAll == 0:
            for key in taxaDict:
                taxaList = taxaDict[key]
                if key == 'Kingdom':
                    if isinstance(taxaList, unicode):
                        tempDF = savedDF.loc[savedDF['kingdomid'] == taxaList]
                    else:
                        tempDF = savedDF.loc[savedDF['kingdomid'].isin(taxaList)]
                    Cols = ['sampleid', 'kingdomid', 'kingdomName']
                    Cols.extend(Dep)
                    tempDF = tempDF[Cols]
                    tempDF.rename(columns={'kingdomid': 'rank_id', 'kingdomName': 'rank_name'}, inplace=True)
                    tempDF.loc[:, 'rank'] = 'Kingdom'
                    taxaDF = taxaDF.append(tempDF)
                elif key == 'Phyla':
                    if isinstance(taxaList, unicode):
                        tempDF = savedDF.loc[savedDF['phylaid'] == taxaList]
                    else:
                        tempDF = savedDF.loc[savedDF['phylaid'].isin(taxaList)]
                    Cols = ['sampleid', 'phylaid', 'phylaName']
                    Cols.extend(Dep)
                    tempDF = tempDF[Cols]
                    tempDF.rename(columns={'phylaid': 'rank_id', 'phylaName': 'rank_name'}, inplace=True)
                    tempDF.loc[:, 'rank'] = 'Phyla'
                    taxaDF = taxaDF.append(tempDF)
                elif key == 'Class':
                    if isinstance(taxaList, unicode):
                        tempDF = savedDF.loc[savedDF['classid'] == taxaList]
                    else:
                        tempDF = savedDF.loc[savedDF['classid'].isin(taxaList)]
                    Cols = ['sampleid', 'classid', 'className']
                    Cols.extend(Dep)
                    tempDF = tempDF[Cols]
                    tempDF.rename(columns={'classid': 'rank_id', 'className': 'rank_name'}, inplace=True)
                    tempDF.loc[:, 'rank'] = 'Class'
                    taxaDF = taxaDF.append(tempDF)
                elif key == 'Order':
                    if isinstance(taxaList, unicode):
                        tempDF = savedDF.loc[savedDF['orderid'] == taxaList]
                    else:
                        tempDF = savedDF.loc[savedDF['orderid'].isin(taxaList)]
                    Cols = ['sampleid', 'orderid', 'orderName']
                    Cols.extend(Dep)
                    tempDF = tempDF[Cols]
                    tempDF.rename(columns={'orderid': 'rank_id', 'orderName': 'rank_name'}, inplace=True)
                    tempDF.loc[:, 'rank'] = 'Order'
                    taxaDF = taxaDF.append(tempDF)
                elif key == 'Family':
                    if isinstance(taxaList, unicode):
                        tempDF = savedDF.loc[savedDF['familyid'] == taxaList]
                    else:
                        tempDF = savedDF.loc[savedDF['familyid'].isin(taxaList)]
                    Cols = ['sampleid', 'familyid', 'familyName']
                    Cols.extend(Dep)
                    tempDF = tempDF[Cols]
                    tempDF.rename(columns={'familyid': 'rank_id', 'familyName': 'rank_name'}, inplace=True)
                    tempDF.loc[:, 'rank'] = 'Family'
                    taxaDF = taxaDF.append(tempDF)
                elif key == 'Genus':
                    if isinstance(taxaList, unicode):
                        tempDF = savedDF.loc[savedDF['genusid'] == taxaList]
                    else:
                        tempDF = savedDF.loc[savedDF['genusid'].isin(taxaList)]
                    Cols = ['sampleid', 'genusid', 'genusName']
                    Cols.extend(Dep)
                    tempDF = tempDF[Cols]
                    tempDF.rename(columns={'genusid': 'rank_id', 'genusName': 'rank_name'}, inplace=True)
                    tempDF.loc[:, 'rank'] = 'Genus'
                    taxaDF = taxaDF.append(tempDF)
                elif key == 'Species':
                    if isinstance(taxaList, unicode):
                        tempDF = savedDF.loc[savedDF['speciesid'] == taxaList]
                    else:
                        tempDF = savedDF.loc[savedDF['speciesid'].isin(taxaList)]
                    Cols = ['sampleid', 'speciesid', 'speciesName']
                    Cols.extend(Dep)
                    tempDF = tempDF[Cols]
                    tempDF.rename(columns={'speciesid': 'rank_id', 'speciesName': 'rank_name'}, inplace=True)
                    tempDF.loc[:, 'rank'] = 'Species'
                    taxaDF = taxaDF.append(tempDF)
                elif key == 'OTU_99':
                    if isinstance(taxaList, unicode):
                        tempDF = savedDF.loc[savedDF['otuid'] == taxaList]
                    else:
                        tempDF = savedDF.loc[savedDF['otuid'].isin(taxaList)]
                    Cols = ['sampleid', 'otuid', 'otuName']
                    Cols.extend(Dep)
                    tempDF = tempDF[Cols]
                    tempDF.rename(columns={'otuid': 'rank_id', 'otuName': 'rank_name'}, inplace=True)
                    tempDF.loc[:, 'rank'] = 'OTU_99'
                    taxaDF = taxaDF.append(tempDF)

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    return '', None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif selectAll == 1:
            Cols = ['sampleid', 'kingdomid', 'kingdomName']
            Cols.extend(Dep)
            taxaDF = savedDF.loc[:, Cols]
            taxaDF.rename(columns={'kingdomid': 'rank_id', 'kingdomName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Kingdom'
        elif selectAll == 2:
            Cols = ['sampleid', 'phylaid', 'phylaName']
            Cols.extend(Dep)
            taxaDF = savedDF.loc[:, Cols]
            taxaDF.rename(columns={'phylaid': 'rank_id', 'phylaName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Phyla'
        elif selectAll == 3:
            Cols = ['sampleid', 'classid', 'className']
            Cols.extend(Dep)
            taxaDF = savedDF.loc[:, Cols]
            taxaDF.rename(columns={'classid': 'rank_id', 'className': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Class'
        elif selectAll == 4:
            Cols = ['sampleid', 'orderid', 'orderName']
            Cols.extend(Dep)
            taxaDF = savedDF.loc[:, Cols]
            taxaDF.rename(columns={'orderid': 'rank_id', 'orderName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Order'
        elif selectAll == 5:
            Cols = ['sampleid', 'familyid', 'familyName']
            Cols.extend(Dep)
            taxaDF = savedDF.loc[:, Cols]
            taxaDF.rename(columns={'familyid': 'rank_id', 'familyName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Family'
        elif selectAll == 6:
            Cols = ['sampleid', 'genusid', 'genusName']
            Cols.extend(Dep)
            taxaDF = savedDF.loc[:, Cols]
            taxaDF.rename(columns={'genusid': 'rank_id', 'genusName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Genus'
        elif selectAll == 7:
            Cols = ['sampleid', 'speciesid', 'speciesName']
            Cols.extend(Dep)
            taxaDF = savedDF.loc[:, Cols]
            taxaDF.rename(columns={'speciesid': 'rank_id', 'speciesName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'Species'
        elif selectAll == 9:
            Cols = ['sampleid', 'otuid', 'otuName']
            Cols.extend(Dep)
            taxaDF = savedDF.loc[:, Cols]
            taxaDF.rename(columns={'otuid': 'rank_id', 'otuName': 'rank_name'}, inplace=True)
            taxaDF.loc[:, 'rank'] = 'OTU_99'
        elif selectAll == 8:
            Cols = ['sampleid', 'genusid', 'genusName']
            Cols.extend(Dep)
            taxaDF = savedDF.loc[:, Cols]
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
            return '', None
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        debug("getTaxaDF: metaDF.reset_index")
        metaDF.reset_index(drop=False, inplace=True)
        #print "Type meta:", type(metaDF), "Type taxa:", type(taxaDF)
        finalDF = pd.merge(metaDF, taxaDF, left_on='sampleid', right_on='sampleid', how='inner')
        #print "Type final:", type(finalDF)
        finalDF.reset_index(drop=False, inplace=True)
        #print "Type final2:", type(finalDF)

        if soilHealth:
            cashList = ['sampleid', 'soil_active_C_rating', 'soil_organic_matter_rating', 'soil_texture_sand', 'soil_k_rating',
                        'soil_pH_rating', 'soil_texture_clay', 'soil_agg_stability_rating',
                        'soil_ACE_protein_index_rating', 'soil_texture_silt', 'soil_p_rating',
                        'soil_minor_elements_rating', 'soil_water_cap_rating', 'soil_respiration_four_day_rating', 'CASH_SHI_rating']
            allFields = list(set(allFields) - set(cashList))

        if 'sample_name' not in allFields:
            wantedList = allFields + ['sample_name', 'sampleid', 'rank', 'rank_name', 'rank_id']
        else:
            wantedList = allFields + ['sampleid', 'rank', 'rank_name', 'rank_id']

        if DepVar == 0:
            finalDF['abund'] = finalDF['abund'].astype("float64")
            finalDF = finalDF.groupby(wantedList)[['abund']].sum()
        elif DepVar == 1:
            finalDF['rel_abund'] = finalDF['rel_abund'].astype("float64")
            #print "Type final3:", type(finalDF)
            #print "Type groupBy:", type(finalDF.groupby(wantedList))
            #print "Type rel:", type(finalDF.groupby(wantedList)['rel_abund'])
            #print "Type [rel]:", type(finalDF.groupby(wantedList)[['rel_abund']])
            #print "Cols:", finalDF.columns
            #print "Wanted:", wantedList
            #print "rel_abund:", finalDF['rel_abund']
            # this groupby has too much data to handle properly, like 26 thousand rows of data here (errors with 7k too)
            # did not error with 1600 rows of data. Can we do this in chunks?
            # how to even split this up, if we want a column of sums. We'd need to sum each chunk, then sum the sums
            # "negative dimensions" has to deal with overflowing when performing the groupby, its fine with 1600 rows, but not 7k
            #print "finalDF preHead:", finalDF.head

            finalDF = finalDF.groupby(wantedList)[['rel_abund']].sum()      # TODO 1.3 error here: negative dimensions not allowed (?)

            # its fine to sum this as a series but not a dataframe (memory issue it seems, but only on dev build)
            # could try summing as a series, then updating the appropriate column in the dataframe
            #finalDF = finalDF.groupby(wantedList)
            # TODO 1.3 I cannot make sense of this bug, it doesn't occur using the same code on the same dataset on live
            # going to switch to qiime biom upload for now
            #finalDF['rel_abund'] = finalDF['rel_abund'].sum()
            #print "finalDF postHead:", finalDF.head
            #print "finalDF PostCols:", finalDF.columns  # so yeah this is only rel_abund left
            '''try:
                finalDF = finalDF.groupby(wantedList)[['rel_abund']]
                print "Grouped and selected", finalDF
                # /\ this broke: "Cannot access attribute 'columns' of 'DataFrameGroupBy' objects, try using the 'apply' method"
                finalDF = finalDF.sum()
                print "Summed", finalDF.columns
            except Exception as ex:
                print "Error during test:", ex'''
        elif DepVar == 4:
            finalDF['abund_16S'] = finalDF['abund_16S'].astype("float64")
            finalDF = finalDF.groupby(wantedList)[['abund_16S']].sum()
        elif DepVar == 2:
            finalDF['rich'] = finalDF['rich'].astype("float64")
            finalDF = finalDF.groupby(wantedList)[['rich']].sum()
        elif DepVar == 3:
            finalDF['diversity'] = finalDF['diversity'].astype("float64")
            finalDF = finalDF.groupby(wantedList)[['diversity']].sum()

        debug("getTaxaDF: finalDF.reset_index")
        finalDF.reset_index(drop=False, inplace=True)   # cannot reset_index on a series, cast to dataframe? what is going on?
        debug("TaxaDF: returning finalDF as", type(finalDF))
        return finalDF, missingList

    except Exception:
        if not stops[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {'error': "There was an error with your analysis!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
            res = json.dumps(myDict)
            return res, None


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
                    return '', None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif keggAll == 1:
            keys = ko_lvl1.objects.using('picrust').values_list('ko_lvl1_id', flat=True)
            for key in keys:
                koList = ko_entry.objects.using('picrust').filter(ko_lvl1_id_id=key).values_list('ko_orthology', flat=True)
                if koList:
                    koDict[key] = koList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    return '', None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif keggAll == 2:
            keys = ko_lvl2.objects.using('picrust').values_list('ko_lvl2_id', flat=True)
            for key in keys:
                koList = ko_entry.objects.using('picrust').filter(ko_lvl2_id_id=key).values_list('ko_orthology', flat=True)
                if koList:
                    koDict[key] = koList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    return '', None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif keggAll == 3:
            keys = ko_lvl3.objects.using('picrust').values_list('ko_lvl3_id', flat=True)
            for key in keys:
                koList = ko_entry.objects.using('picrust').filter(ko_lvl3_id_id=key).values_list('ko_orthology', flat=True)
                if koList:
                    koDict[key] = koList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    return '', None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif keggAll == 4:
            keys = ko_entry.objects.using('picrust').values_list('ko_lvl4_id', flat=True)
            for key in keys:
                koList = ko_entry.objects.using('picrust').filter(ko_lvl4_id=key).values_list('ko_orthology', flat=True)
                if koList:
                    koDict[key] = koList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    return '', None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            return '', None
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        # create sample and otu lists based on meta data selection
        wanted = ['sampleid', 'otuid']
        if DepVar == 0:
            wanted.append('abund')
        elif DepVar == 1:
            wanted.append('rel_abund')
        elif DepVar == 2:
            wanted.append('rich')
        elif DepVar == 3:
            wanted.append('diversity')
        elif DepVar == 4:
            wanted.append('abund_16S')

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
            profileDF['abund'] = profileDF['abund'].astype("float64")
            profileDF = profileDF.pivot(index='otuid', columns='sampleid', values='abund')
        elif DepVar == 1:
            profileDF['rel_abund'] = profileDF['rel_abund'].astype("float64")
            profileDF = profileDF.pivot(index='otuid', columns='sampleid', values='rel_abund')
        elif DepVar == 2:
            profileDF['rich'] = profileDF['rich'].astype("float64")
            profileDF = profileDF.pivot(index='otuid', columns='sampleid', values='rich')
        elif DepVar == 3:
            profileDF['diversity'] = profileDF['diversity'].astype("float64")
            profileDF = profileDF.pivot(index='otuid', columns='sampleid', values='diversity')
        elif DepVar == 4:
            profileDF['abund_16S'] = profileDF['abund_16S'].astype("float64")
            profileDF = profileDF.pivot(index='otuid', columns='sampleid', values='abund_16S')

        sampleList = profileDF.columns.values.tolist()
        taxaDF = pd.DataFrame()
        total = len(sampleList)
        counter = 1
        for i in sampleList:
            sampStartTime = time.time()
            tempDF = pd.DataFrame(index=profileDF.index)
            loopStartTime = time.time()
            for j in levelList:
                tempDF[j] = profileDF[i] * picrustDF[j]
            loopRunTime = time.time() - loopStartTime
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                return '', None
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            tempDF = tempDF[levelList].sum()
            tempDF['sampleid'] = i
            taxaDF = taxaDF.append(tempDF, ignore_index=True)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                return '', None
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            functions.setBase(RID, 'Calculating KEGG pathway abundances...sample ' + str(counter) + ' out of ' + str(total) + ' is finished!')
            sampRunTime = time.time() - sampStartTime
            percentTime = loopRunTime * 100.0 / sampRunTime
            debug("GAGE tracker: this sample took", sampRunTime, "and", percentTime, "% was spent in for loop")
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
                    str(level) + ' not found in database'
                namesDict[level] = name

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    return '', None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            allDF = picrustDF.reset_index(drop=False, inplace=False)
            allDF.rename(columns={'index': 'otuid'}, inplace=True)
            allDF.dropna(axis=1, how='any', inplace=False)
            allDF = allDF[(allDF.sum(axis=1) != 0)]

            for i in levelList:
                allDF[i] = allDF[i] / allDF[i]
            allDF.fillna(value=0, inplace=True)
            recordDict = {}
            otuList = allDF.otuid.tolist()
            for id in otuList:
                try:
                    qs = OTU_99.objects.all().filter(otuid=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyid_id__familyName', 'genusid_id__genusName', 'speciesid_id__speciesName', 'otuName')
                    record = '|'.join(qs[0])
                    recordDict[id] = record
                except:
                    recordDict[id] = 'No data'
            rnaList = PICRUSt.objects.using('picrust').filter(otuid_id__in=otuList).values_list('otuid', 'rRNACount')
            allDF['Taxonomy'] = allDF['otuid'].map(recordDict)
            allDF['rRNACount'] = allDF['otuid'].map(dict(rnaList))
            order = ['otuid', 'Taxonomy', 'rRNACount'] + levelList
            allDF = allDF[order]
            allDF.rename(columns=namesDict, inplace=True)
        else:
            allDF = ''

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            return '', None
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
            return '', None
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        idList = functions.getFullKO(list(finalDF.rank_id.unique()))
        finalDF['rank_name'] = finalDF['rank_id'].map(idList)

        # required to match getTaxaDF output
        metaDF.reset_index(drop=False, inplace=True)

        return finalDF, allDF

    except Exception as exc:
        logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
        myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
        logging.exception(myDate)
        myDict = {'error': "There was an error with your analysis:" + str(exc)}
        res = json.dumps(myDict)
        print "Error in keggDF:", exc
        return str(res), None


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
                    return '', None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 1:
            keys = nz_lvl1.objects.using('picrust').values_list('nz_lvl1_id', flat=True)
            for key in keys:
                nzList = nz_entry.objects.using('picrust').filter(nz_lvl1_id_id=key).values_list('nz_orthology', flat=True)
                if nzList:
                    nzDict[key] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    return '', None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 2:
            keys = nz_lvl2.objects.using('picrust').values_list('nz_lvl2_id', flat=True)
            for key in keys:
                nzList = nz_entry.objects.using('picrust').filter(nz_lvl2_id_id=key).values_list('nz_orthology', flat=True)
                if nzList:
                    nzDict[key] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    return '', None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 3:
            keys = nz_lvl3.objects.using('picrust').values_list('nz_lvl3_id', flat=True)
            for key in keys:
                nzList = nz_entry.objects.using('picrust').filter(nz_lvl3_id_id=key).values_list('nz_orthology', flat=True)
                if nzList:
                    nzDict[key] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    return '', None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif nzAll == 4:
            keys = nz_lvl4.objects.using('picrust').values_list('nz_lvl4_id', flat=True)
            for key in keys:
                nzList = nz_entry.objects.using('picrust').filter(nz_lvl4_id_id=key).values_list('nz_orthology', flat=True)
                if nzList:
                    nzDict[key] = nzList

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    return '', None
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
            id = 'hcnA: glycine dehydrogenase (cyanide-forming)'
            idList = ['K10814']
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

            id = 'prnD: aminopyrrolnitrin oxygenase'
            idList = ['K19982']
            nzDict[id] = idList

            id = 'phlD: phloroglucinol synthase'
            idList = ['K15431']
            nzDict[id] = idList

            id = 'ituA: iturin family lipopeptide synthetase A'
            idList = ['K15661']
            nzDict[id] = idList

            id = 'fenA: fengycin family lipopeptide synthetase D'
            idList = ['K15667']
            nzDict[id] = idList

            id = 'srfAA: surfactin family lipopeptide synthetase A'
            idList = ['K15654']
            nzDict[id] = idList

            id = 'rifM: AHBA synthesis associated protein'
            idList = ['K16017']
            nzDict[id] = idList

            id = 'phzE: 2-amino-4-deoxychorismate synthase'
            idList = ['K13063']
            nzDict[id] = idList

            ### Root growth (drought/salt stress)
            id = 'ipdC: indolepyruvate decarboxylase'
            idList = ['K04103']
            nzDict[id] = idList

            id = 'acdS: 1-aminocyclopropane-1-carboxylate deaminase'
            idList = ['K01505']
            nzDict[id] = idList

            ### Siderophores
            id = 'mbtI: salicylate synthetase'
            idList = ['K04781']
            nzDict[id] = idList

            id = 'entA: 2,3-dihydro-2,3-dihydroxybenzoate dehydrogenase'
            idList = ['K00216']
            nzDict[id] = idList

            id = 'pchB: isochorismate pyruvate lysase'
            idList = ['K02364']
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

        elif nzAll == 7:
            ### Phosphoesterase genes
            id = 'PHO: acid phosphatase'
            idList = ['K01078']
            nzDict[id] = idList

            id = 'phoN: acid phosphatase (class A)'
            idList = ['K09474']
            nzDict[id] = idList

            id = 'aphA: acid phosphatase (class B)'
            idList = ['K03788']
            nzDict[id] = idList

            id = 'phoA: E3.1.3.1, phoA, phoB'
            idList = ['K01077']
            nzDict[id] = idList

            id = ''  # TODO 1.3 more kegg
            # jump
            idList = ['K01113']
            nzDict[id] = idList

            id = ''
            idList = ['K01126']
            nzDict[id] = idList

            id = ''
            idList = ['K07048']
            nzDict[id] = idList

            ### Phytase genes
            id = ''
            idList = ['K01083']
            nzDict[id] = idList

            id = ''
            idList = ['K01093']
            nzDict[id] = idList

            ### Phosphonate degradation genes
            id = ''
            idList = ['K02043']
            nzDict[id] = idList

            id = ''
            idList = ['K06166']
            nzDict[id] = idList

            id = ''
            idList = ['K06165']
            nzDict[id] = idList

            id = ''
            idList = ['K06164']
            nzDict[id] = idList

            id = ''
            idList = ['K06163']
            nzDict[id] = idList

            id = ''
            idList = ['K05781']
            nzDict[id] = idList

            id = ''
            idList = ['K05780']
            nzDict[id] = idList

            id = ''
            idList = ['K06162']
            nzDict[id] = idList

            id = ''
            idList = ['K05774']
            nzDict[id] = idList

            id = ''
            idList = ['K09994']
            nzDict[id] = idList

            id = ''
            idList = ['K06167']
            nzDict[id] = idList

            id = ''
            idList = ['K03430']
            nzDict[id] = idList

            id = ''
            idList = ['K05306']
            nzDict[id] = idList

            id = ''
            idList = ['K06193']
            nzDict[id] = idList

            #Inorganic phosphate solubilizing genes
            id = ''
            idList = ['K01507']
            nzDict[id] = idList

            id = ''
            idList = ['K01524']
            nzDict[id] = idList

            id = ''
            idList = ['K00937']
            nzDict[id] = idList

            id = ''
            idList = ['K00117']
            nzDict[id] = idList

            # Phosphorus transporter genes
            id = ''
            idList = ['K03306']
            nzDict[id] = idList

            id = ''
            idList = ['K16322']
            nzDict[id] = idList

            id = ''
            idList = ['K02038']
            nzDict[id] = idList

            id = ''
            idList = ['K02036']
            nzDict[id] = idList

            id = ''
            idList = ['K02037']
            nzDict[id] = idList

            id = ''
            idList = ['K02040']
            nzDict[id] = idList

            id = ''
            idList = ['K02041']
            nzDict[id] = idList

            id = ''
            idList = ['K02044']
            nzDict[id] = idList

            id = ''
            idList = ['K02042']
            nzDict[id] = idList

            id = ''
            idList = ['K05814']
            nzDict[id] = idList

            id = ''
            idList = ['K05813']
            nzDict[id] = idList

            id = ''
            idList = ['K05816']
            nzDict[id] = idList

            id = ''
            idList = ['K05815']
            nzDict[id] = idList

            id = ''
            idList = ['K07657']
            nzDict[id] = idList

            id = ''
            idList = ['K07636']
            nzDict[id] = idList

            id = ''
            idList = ['K02039']
            nzDict[id] = idList

        elif nzAll == 10:
            ### PO4 Solubility
            id = 'pqqC: pyrroloquinoline-quinone synthase'
            idList = ['K06137']
            nzDict[id] = idList

            ### Biocontrol
            id = 'hcnA: glycine dehydrogenase (cyanide-forming)'
            idList = ['K10814']
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

            id = 'prnD: aminopyrrolnitrin oxygenase'
            idList = ['K19982']
            nzDict[id] = idList

            id = 'phlD: phloroglucinol synthase'
            idList = ['K15431']
            nzDict[id] = idList

            id = 'ituA: iturin family lipopeptide synthetase A'
            idList = ['K15661']
            nzDict[id] = idList

            id = 'fenA: fengycin family lipopeptide synthetase D'
            idList = ['K15667']
            nzDict[id] = idList

            id = 'srfAA: surfactin family lipopeptide synthetase A'
            idList = ['K15654']
            nzDict[id] = idList

            id = 'rifM: AHBA synthesis associated protein'
            idList = ['K16017']
            nzDict[id] = idList

            id = 'phzE: 2-amino-4-deoxychorismate synthase'
            idList = ['K13063']
            nzDict[id] = idList

            ### Root growth (drought/salt stress)
            id = 'ipdC: indolepyruvate decarboxylase'
            idList = ['K04103']
            nzDict[id] = idList

            id = 'acdS: 1-aminocyclopropane-1-carboxylate deaminase'
            idList = ['K01505']
            nzDict[id] = idList

            ### Siderophores
            id = 'mbtI: salicylate synthetase'
            idList = ['K04781']
            nzDict[id] = idList

            id = 'entA: 2,3-dihydro-2,3-dihydroxybenzoate dehydrogenase'
            idList = ['K00216']
            nzDict[id] = idList

            id = 'pchB: isochorismate pyruvate lysase'
            idList = ['K02364']
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
        wanted = ['sampleid', 'otuid']
        if DepVar == 0:
            wanted.append('abund')
        elif DepVar == 1:
            wanted.append('rel_abund')
        elif DepVar == 2:
            wanted.append('rich')
        elif DepVar == 3:
            wanted.append('diversity')
        elif DepVar == 4:
            wanted.append('abund_16S')

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
        debug("NZDF: Path")
        for i in sampleList:
            profileDF[i] = profileDF[i].astype("float64")
            tempDF = pd.DataFrame(index=profileDF.index)
            for j in levelList:
                tempDF[j] = profileDF[i] * picrustDF[j]

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    return '', None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            tempDF = tempDF[levelList].sum()
            tempDF['sampleid'] = i
            taxaDF = taxaDF.append(tempDF, ignore_index=True)

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stops[PID] == RID:
                return '', None
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            functions.setBase(RID, 'Calculating KEGG pathway abundances...sample ' + str(counter) + ' out of ' + str(total) + ' is finished!')
            counter += 1

        functions.setBase(RID, curStep)
        debug("NZDF: Map")
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
                    str(level) + ' not found in database'
                namesDict[level] = name

                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
                if stops[PID] == RID:
                    return '', None
                # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            allDF = picrustDF.reset_index(drop=False, inplace=False)
            allDF.rename(columns={'index': 'otuid'}, inplace=True)
            allDF.dropna(axis=1, how='any', inplace=False)
            allDF = allDF[(allDF.sum(axis=1) != 0)]

            for i in levelList:
                allDF[i] = allDF[i] / allDF[i]
            allDF.fillna(value=0, inplace=True)
            recordDict = {}
            otuList = allDF.otuid.tolist()
            for id in otuList:
                try:
                    qs = OTU_99.objects.all().filter(otuid=id).values_list('kingdomid_id__kingdomName', 'phylaid_id__phylaName', 'classid_id__className', 'orderid_id__orderName', 'familyid_id__familyName', 'genusid_id__genusName', 'speciesid_id__speciesName', 'otuName')
                    record = '|'.join(qs[0])
                    recordDict[id] = record
                except:
                    recordDict[id] = 'No data'
            rnaList = PICRUSt.objects.using('picrust').filter(otuid_id__in=otuList).values_list('otuid', 'rRNACount')
            allDF['Taxonomy'] = allDF['otuid'].map(recordDict)
            allDF['rRNACount'] = allDF['otuid'].map(dict(rnaList))
            order = ['otuid', 'Taxonomy', 'rRNACount'] + levelList
            allDF = allDF[order]
            if nzAll < 5:
                allDF.rename(columns=namesDict, inplace=True)
        else:
            allDF = ''  # soil health uses this

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            return '', None
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        #taxaDF.reset_index(drop=False, inplace=True)
        debug("NZDF: Melt")
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
        debug("NZDF: Merge")
        finalDF = pd.merge(metaDF, taxaDF, left_index=True, right_index=True, how='inner')
        finalDF.reset_index(drop=False, inplace=True)

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stops[PID] == RID:
            return '', None
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        if nzAll < 5:
            idList = functions.getFullNZ(list(finalDF.rank_id.unique()))
            finalDF['rank_name'] = finalDF['rank_id'].map(idList)
        else:
            finalDF['rank_name'] = finalDF['rank_id']

        # required to match getTaxaDF output
        metaDF.reset_index(drop=False, inplace=True)
        return finalDF, allDF

    except Exception as exc:
        print "Error in NZDF:", exc
        logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
        myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
        logging.exception(myDate)
        myDict = {'error': "There was an error with your analysis!\nMore info can be found in 'error_log.txt' located in your myPhyloDB dir."}
        res = json.dumps(myDict)
        return res, None


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
    # this makes no assumptions about the level of each id in list, each can have a different level
    # returns a dictionary with key : value of id : list of names of parent taxa (and self)
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


# Can assert that all ids given as args belong to the given level, does not support PGR atm TODO 1.3 PGR support?
def getFullTaxaFromID(id, level):
    results = []    # TODO 1.4 change this to a dictionary where k:v = taxID:taxName
    # given level (kingdom phyla class etc), and list of ids, return a dictionary with key:val of id:list_of_parent_ids (+ source id)
    if level == 1:
        # kingdom is a weird case since its the top, so an idlist of kingdoms will have no parent ids
        pass
    elif level == 2:
        if Phyla.objects.all().filter(phylaid=id).exists():
            myPhyla = Phyla.objects.get(phylaid=id)
            results.append(myPhyla.kingdomid.kingdomid)
    elif level == 3:
        if Class.objects.all().filter(classid=id).exists():
            myClass = Class.objects.get(classid=id)
            results.append(myClass.kingdomid.kingdomid)
            results.append(myClass.phylaid.phylaid)
    elif level == 4:
        if Order.objects.all().filter(orderid=id).exists():
            myOrder = Order.objects.get(orderid=id)
            results.append(myOrder.kingdomid.kingdomid)
            results.append(myOrder.phylaid.phylaid)
            results.append(myOrder.classid.classid)
    elif level == 5:
        if Family.objects.all().filter(familyid=id).exists():
            myFamily = Family.objects.get(familyid=id)
            results.append(myFamily.kingdomid.kingdomid)
            results.append(myFamily.phylaid.phylaid)
            results.append(myFamily.classid.classid)
            results.append(myFamily.orderid.orderid)
    elif level == 6:
        if Genus.objects.all().filter(genusid=id).exists():
            myGenus = Genus.objects.get(genusid=id)
            results.append(myGenus.kingdomid.kingdomid)
            results.append(myGenus.phylaid.phylaid)
            results.append(myGenus.classid.classid)
            results.append(myGenus.orderid.orderid)
            results.append(myGenus.familyid.familyid)
    elif level == 7:
        if Species.objects.all().filter(species=id).exists():
            mySpecies = Species.objects.get(species=id)
            results.append(mySpecies.kingdomid.kingdomid)
            results.append(mySpecies.phylaid.phylaid)
            results.append(mySpecies.classid.classid)
            results.append(mySpecies.orderid.orderid)
            results.append(mySpecies.familyid.familyid)
            results.append(mySpecies.genusid.genusid)
    elif level == 9:
        if OTU_99.objects.all().filter(otuid=id).exists():
            myOTU = OTU_99.objects.get(otuid=id)
            results.append(myOTU.kingdomid.kingdomid)
            results.append(myOTU.phylaid.phylaid)
            results.append(myOTU.classid.classid)
            results.append(myOTU.orderid.orderid)
            results.append(myOTU.familyid.familyid)
            results.append(myOTU.genusid.genusid)
            results.append(myOTU.speciesid.speciesid)
    results.append(id)
    return results


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
    debug("filterDF")
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
    debug("filterDF: myLevel")
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
    debug("filterDF: remUnclass")
    if remUnclass == 'yes':
        # check if selecting based on level first, create else statement for tree usage
        savedDF = savedDF[~savedDF[myLevel].str.contains('unclassified')]

    if remZeroes == 'yes' and perZeroes > 0:
        threshold = int(perZeroes / 100.0 * numSamples)
        bytag = savedDF.groupby(myID).aggregate(np.count_nonzero)
        tags = bytag[bytag[myVar] >= threshold].index.tolist()
        savedDF = savedDF[savedDF[myID].isin(tags)]
    debug("filterDF: filterData")
    if filterData == 'yes' and not level == 8:
        threshold = int(filterPer)
        if filterMeth == 1:
            debug("filterDF: filterMeth == 1")
            pass
        else:
            savedDF[myVar] = savedDF[myVar].astype("float64")
            sum = savedDF.groupby([myID, 'sampleid'])[myVar].sum()
            if filterMeth == 2:
                debug("filterDF: filterMeth == 2")
                sumDF = sum.reset_index(drop=False)
                bytag_q3 = sumDF.groupby(myID)[myVar].quantile(0.75)
                bytag_q1 = sumDF.groupby(myID)[myVar].quantile(0.25)
                bytag = bytag_q3 - bytag_q1
                bytag.sort_values(axis=0, ascending=False, inplace=True)
                tags = bytag[:threshold].index.tolist()
                savedDF = savedDF[savedDF[myID].isin(tags)]
            elif filterMeth == 3:
                debug("filterDF: filterMeth == 3")
                sumDF = sum.reset_index(drop=False)
                bytag_sd = sumDF.groupby(myID)[myVar].std()
                bytag_mean = sumDF.groupby(myID)[myVar].mean()
                bytag = bytag_sd / bytag_mean
                bytag.sort_values(axis=0, ascending=False, inplace=True)
                tags = bytag[:threshold].index.tolist()
                savedDF = savedDF[savedDF[myID].isin(tags)]
            elif filterMeth == 4:
                debug("filterDF: filterMeth == 4")
                sumDF = sum.reset_index(drop=False)
                bytag = sumDF.groupby(myID)[myVar].std()
                bytag.sort_values(axis=0, ascending=False, inplace=True)
                tags = bytag[:threshold].index.tolist()
                savedDF = savedDF[savedDF[myID].isin(tags)]
            elif filterMeth == 5:
                debug("filterDF: filterMeth == 5")
                sumDF = sum.reset_index(drop=False)
                bytag = sumDF.groupby(myID)[myVar].mean()
                bytag.sort_values(axis=0, ascending=False, inplace=True)
                tags = bytag[:threshold].index.tolist()
                savedDF = savedDF[savedDF[myID].isin(tags)]
            elif filterMeth == 6:
                debug("filterDF: filterMeth == 6")
                sumDF = sum.reset_index(drop=False)
                bytag = sumDF.groupby(myID)[myVar].median()
                bytag.sort_values(axis=0, ascending=False, inplace=True)
                tags = bytag[:threshold].index.tolist()
                savedDF = savedDF[savedDF[myID].isin(tags)]
    debug("filterDF: return")
    return savedDF


