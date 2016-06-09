from datetime import datetime
from django.http import HttpResponse
import math
import pandas as pd
import simplejson
import gzip
import re
from uuid import uuid4

from database.models import Kingdom, Phyla, Class, Order, Family, Genus, Species
from database.models import ko_lvl1, ko_lvl2, ko_lvl3, ko_entry
from database.models import nz_lvl1, nz_lvl2, nz_lvl3, nz_lvl4, nz_entry
from database.models import PICRUSt


stage = ''
perc = 0
pd.set_option('display.max_colwidth', -1)


def statusPyBake(request):
    global base, perc
    if request.is_ajax():
        myDict = {'stage': stage}
        json_data = simplejson.dumps(myDict, encoding="Latin-1")
        return HttpResponse(json_data, content_type='application/json')


def geneParse(file1, file2, file3):
    global stage
    time1 = datetime.now()

    # remove all PICRUSt data
    PICRUSt.objects.using('picrust').all().delete()

    #taxonomy file
    stage = 'Step 1 of 4: Parsing taxonomy file...'
    df1 = pd.read_csv(file1, header=0, sep='\t', index_col=0)
    file1.close()

    # rRNA count file
    stage = 'Step 2 of 4: Parsing rRNA count file...'
    dFile2 = gzip.GzipFile(fileobj=file2, mode='rb')
    df2 = pd.read_csv(dFile2, header=0, sep='\t', index_col=0)
    file2.close()
    dFile2.close()

    # merging dataframes
    stage = 'Step 3 of 4: Merging files...'
    mergedDF = pd.merge(df1, df2, left_index=True, right_index=True, how='inner')
    total, col = mergedDF.shape

    dFile3 = gzip.GzipFile(fileobj=file3, mode='rb')
    chunksize = 1000.0
    df3 = pd.read_csv(dFile3, header=0, sep='\t', index_col=0, iterator=True, chunksize=chunksize, error_bad_lines=False, low_memory=False)
    chunks = int(math.ceil(total/chunksize)) + 1

    # add to database
    counter = 0
    sumDF = pd.DataFrame()
    countDF = pd.DataFrame()
    colList = []
    for chunk in df3:
        temp = pd.merge(mergedDF, chunk, left_index=True, right_index=True, how='inner')
        colList = temp.columns.values.tolist()
        colList.remove('Taxonomy')
        colList.remove('metadata_NSTI')
        grouped = temp.groupby('Taxonomy')[colList]
        sumDF = pd.concat([sumDF, grouped.agg('sum')])
        countDF = pd.concat([countDF, grouped.agg('count')])
        counter += 1
        stage = 'Step 3 of 4: Merging files...chunk ' + str(counter) + ' out of ' + str(chunks) + ' is complete!'

    sumDF.reset_index(drop=False, inplace=True)
    countDF.reset_index(drop=False, inplace=True)

    sum = sumDF.groupby('Taxonomy')[colList].agg('sum')
    count = countDF.groupby('Taxonomy')[colList].agg('sum')

    finalDF_KO = pd.DataFrame(sum/count, columns=colList, index=sum.index)
    finalDF_KO.dropna(inplace=True)
    finalDF_KO.reset_index(drop=False, inplace=True)

    stage = 'Step 4 of 4: Adding PICRUSt data to your database...'
    total, col = finalDF_KO.shape
    counter = 0
    for index, row in finalDF_KO.iterrows():
        subbed = re.sub(r'(\(.*?\)|k__|p__|c__|o__|f__|g__|s__)', '', row['Taxonomy'])
        subbed = re.sub(r';;', 'unclassified', subbed)
        subbed = subbed[:-1]
        taxon = subbed.split(';')

        if not Kingdom.objects.using('default').filter(kingdomName=taxon[0]).exists():
            kid = uuid4().hex
            Kingdom.objects.using('default').create(kingdomid=kid, kingdomName=taxon[0])

        k = Kingdom.objects.using('default').get(kingdomName=taxon[0]).kingdomid
        if not Phyla.objects.using('default').filter(kingdomid_id=k, phylaName=taxon[1]).exists():
            pid = uuid4().hex
            Phyla.objects.using('default').create(kingdomid_id=k, phylaid=pid, phylaName=taxon[1])

        p = Phyla.objects.using('default').get(kingdomid_id=k, phylaName=taxon[1]).phylaid
        if not Class.objects.using('default').filter(kingdomid_id=k, phylaid_id=p, className=taxon[2]).exists():
            cid = uuid4().hex
            Class.objects.using('default').create(kingdomid_id=k, phylaid_id=p, classid=cid, className=taxon[2])

        c = Class.objects.using('default').get(kingdomid_id=k, phylaid_id=p, className=taxon[2]).classid
        if not Order.objects.using('default').filter(kingdomid_id=k, phylaid_id=p, classid_id=c, orderName=taxon[3]).exists():
            oid = uuid4().hex
            Order.objects.using('default').create(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid=oid, orderName=taxon[3])

        o = Order.objects.using('default').get(kingdomid_id=k, phylaid_id=p, classid_id=c, orderName=taxon[3]).orderid
        if not Family.objects.using('default').filter(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyName=taxon[4]).exists():
            fid = uuid4().hex
            Family.objects.using('default').create(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid=fid, familyName=taxon[4])

        f = Family.objects.using('default').get(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyName=taxon[4]).familyid
        if not Genus.objects.using('default').filter(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusName=taxon[5]).exists():
            gid = uuid4().hex
            Genus.objects.using('default').create(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusid=gid, genusName=taxon[5])

        g = Genus.objects.using('default').get(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusName=taxon[5]).genusid
        if not Species.objects.using('default').filter(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusid_id=g, speciesName=taxon[6]).exists():
            sid = uuid4().hex
            Species.objects.using('default').create(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusid_id=g, speciesid=sid, speciesName=taxon[6])

        s = Species.objects.using('default').get(kingdomid_id=k, phylaid_id=p, classid_id=c, orderid_id=o, familyid_id=f, genusid_id=g, speciesName=taxon[6]).speciesid

        rRNACount = row['16S_rRNA_Count']
        row.drop('16S_rRNA_Count', inplace=True)
        row.drop('Taxonomy', inplace=True)
        rowDict = row.to_dict()
        geneCount = {x:y for x,y in rowDict.items() if y != 0}
        for key in geneCount:
            geneCount[key] = round(geneCount[key] / rRNACount, 4)

        if not PICRUSt.objects.using('picrust').filter(speciesid_id=s).exists():
            PICRUSt.objects.using('picrust').create(speciesid_id=s, rRNACount=rRNACount, geneCount=geneCount)

        counter += 1
        stage = 'Step 4 of 4: Adding PICRUSt data to your database...OTU ' + str(counter) + ' out of ' + str(total) + ' is complete!'

    time2 = datetime.now()
    delta = time2 - time1
    print 'Total time to update PICRUSt database: ', delta
    return None


def koParse(file):
    global stage, time1
    stage = 'Parsing KEGG pathways...'
    time1 = datetime.now()

    # remove all KEGG data
    ko_lvl1.objects.using('picrust').all().delete()
    ko_lvl2.objects.using('picrust').all().delete()
    ko_lvl3.objects.using('picrust').all().delete()
    ko_entry.objects.using('picrust').all().delete()

    text = file.read().splitlines()
    total = 0
    for line in text:
        total += 1

    file.seek(0)
    level1, level2, level3, level4, KEGG_entry, KEGG_name, KEGG_desc = '', '', '', '', '', '', ''
    counter = 0
    for line in text:
        if len(line) > 2:
            if line.startswith('A'):
                level1 = line.split('<b>')[1].split('</b>')[0]
                continue
            if line.startswith('B'):
                level2 = line.split(' <b>')[1].split('</b>')[0]
                continue
            if line.startswith('C'):
                level3 = line.split('    ')[1].split(' ', 1)[1]
                continue
            if line.startswith('D'):
                level4 = line.split('      ')[1]

                name = level4.split('; ')
                KEGG_entry = name[0].split('  ')[0]
                if not KEGG_entry.startswith('K'):
                    print "Skipped line:", line
                    continue
                KEGG_name = name[0].split('  ')[1]
                KEGG_desc = name[1]

                if not ko_lvl1.objects.using('picrust').filter(ko_lvl1_name=level1).exists():
                    ko_lvl1_id = uuid4().hex
                    ko_lvl1.objects.using('picrust').create(ko_lvl1_id=ko_lvl1_id, ko_lvl1_name=level1)

                lv1 = ko_lvl1.objects.using('picrust').get(ko_lvl1_name=level1).ko_lvl1_id
                if not ko_lvl2.objects.using('picrust').filter(ko_lvl1_id_id=lv1, ko_lvl2_name=level2).exists():
                    ko_lvl2_id = uuid4().hex
                    ko_lvl2.objects.using('picrust').create(ko_lvl1_id_id=lv1, ko_lvl2_id=ko_lvl2_id, ko_lvl2_name=level2)

                lv2 = ko_lvl2.objects.using('picrust').get(ko_lvl2_name=level2).ko_lvl2_id
                if not ko_lvl3.objects.using('picrust').filter(ko_lvl1_id_id=lv1, ko_lvl2_id_id=lv2, ko_lvl3_name=level3).exists():
                    ko_lvl3_id = uuid4().hex
                    ko_lvl3.objects.using('picrust').create(ko_lvl1_id_id=lv1, ko_lvl2_id_id=lv2, ko_lvl3_id=ko_lvl3_id, ko_lvl3_name=level3)

                lv3 = ko_lvl3.objects.using('picrust').get(ko_lvl3_name=level3).ko_lvl3_id
                if not ko_entry.objects.using('picrust').filter(ko_lvl1_id_id=lv1, ko_lvl2_id_id=lv2, ko_lvl3_id_id=lv3, ko_orthology=KEGG_entry).exists():
                    ko_lvl4_id = uuid4().hex
                    ko_entry.objects.using('picrust').create(ko_lvl1_id_id=lv1, ko_lvl2_id_id=lv2, ko_lvl3_id_id=lv3, ko_lvl4_id=ko_lvl4_id, ko_orthology=KEGG_entry, ko_name=KEGG_name, ko_desc=KEGG_desc)

        stage = 'Parsing KEGG pathways...line ' + str(counter) + ' out of ' + str(total) + ' is complete!'
        counter += 1

    time2 = datetime.now()
    delta = time2 - time1
    print 'Total time to update KEGG database: ', delta

    return None


def nzParse(file):
    global stage, time1
    stage = 'Parsing KEGG enzyme...'
    time1 = datetime.now()

    # remove all KEGG data
    nz_lvl1.objects.using('picrust').all().delete()
    nz_lvl2.objects.using('picrust').all().delete()
    nz_lvl3.objects.using('picrust').all().delete()
    nz_entry.objects.using('picrust').all().delete()

    text = file.read().splitlines()
    total = 0
    for line in text:
        total += 1

    file.seek(0)
    level1, level2, level3, level4, level5, KEGG_entry, KEGG_name, KEGG_desc = '', '', '', '', '', '', '', ''
    counter = 0
    for line in text:
        if len(line) > 2:
            if line.startswith('A'):
                level1 = line.split('<b>')[1].split('</b>')[0]
                continue
            if line.startswith('B'):
                level2 = line.split('  ', 1)[1]
                continue
            if line.startswith('C'):
                level3 = line.split('    ')[1]
                continue
            if line.startswith('D'):
                level4 = line.split('      ')[1]
                continue
            if line.startswith('E'):
                level5 = line.split('        ')[1]
                if ';' in level5:
                    name = level5.split('; ')
                    KEGG_entry = name[0].split('  ')[0]
                    KEGG_name = name[0].split('  ')[1]
                    KEGG_desc = name[1]
                else:
                    name = level5.split()
                    KEGG_entry = name[0]
                    if not KEGG_entry.startswith('K'):
                        print "Skipped line:", line
                        continue
                    KEGG_name = name[1]
                    KEGG_desc = ''

                if not nz_lvl1.objects.using('picrust').filter(nz_lvl1_name=level1).exists():
                    nz_lvl1_id = uuid4().hex
                    nz_lvl1.objects.using('picrust').create(nz_lvl1_id=nz_lvl1_id, nz_lvl1_name=level1)

                lv1 = nz_lvl1.objects.using('picrust').get(nz_lvl1_name=level1).nz_lvl1_id
                if not nz_lvl2.objects.using('picrust').filter(nz_lvl1_id_id=lv1, nz_lvl2_name=level2).exists():
                    nz_lvl2_id = uuid4().hex
                    nz_lvl2.objects.using('picrust').create(nz_lvl1_id_id=lv1, nz_lvl2_id=nz_lvl2_id, nz_lvl2_name=level2)

                lv2 = nz_lvl2.objects.using('picrust').get(nz_lvl2_name=level2).nz_lvl2_id
                if not nz_lvl3.objects.using('picrust').filter(nz_lvl1_id_id=lv1, nz_lvl2_id_id=lv2, nz_lvl3_name=level3).exists():
                    nz_lvl3_id = uuid4().hex
                    nz_lvl3.objects.using('picrust').create(nz_lvl1_id_id=lv1, nz_lvl2_id_id=lv2, nz_lvl3_id=nz_lvl3_id, nz_lvl3_name=level3)

                lv3 = nz_lvl3.objects.using('picrust').get(nz_lvl3_name=level3).nz_lvl3_id
                if not nz_lvl4.objects.using('picrust').filter(nz_lvl1_id_id=lv1, nz_lvl2_id_id=lv2, nz_lvl3_id_id=lv3, nz_lvl4_name=level4).exists():
                    nz_lvl4_id = uuid4().hex
                    nz_lvl4.objects.using('picrust').create(nz_lvl1_id_id=lv1, nz_lvl2_id_id=lv2, nz_lvl3_id_id=lv3, nz_lvl4_id=nz_lvl4_id,  nz_lvl4_name=level4)

                lv4 = nz_lvl4.objects.using('picrust').get(nz_lvl4_name=level4).nz_lvl4_id
                if not nz_entry.objects.using('picrust').filter(nz_lvl1_id_id=lv1, nz_lvl2_id_id=lv2, nz_lvl3_id_id=lv3, nz_lvl4_id_id=lv4, nz_orthology=KEGG_entry).exists():
                    nz_lvl5_id = uuid4().hex
                    nz_entry.objects.using('picrust').create(nz_lvl1_id_id=lv1, nz_lvl2_id_id=lv2, nz_lvl3_id_id=lv3, nz_lvl4_id_id=lv4, nz_lvl5_id=nz_lvl5_id, nz_orthology=KEGG_entry, nz_name=KEGG_name, nz_desc=KEGG_desc)

        stage = 'Parsing KEGG enzymes...line ' + str(counter) + ' out of ' + str(total) + ' is complete!'
        counter += 1

    time2 = datetime.now()
    delta = time2 - time1
    print 'Total time to update KEGG database: ', delta

    return None

