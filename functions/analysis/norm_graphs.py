import datetime
from django.http import HttpResponse
from django_pandas.io import read_frame
import json
import logging
from numpy import *
import numpy as np
from numpy.random.mtrand import RandomState
import pandas as pd
import pickle
from pyper import *
import zipfile

from database.models import Sample, Air, Human_Associated, Microbial, Soil, Water, UserDefined, \
    OTU_99, Profile, DaymetData

import functions


curSamples = {}
totSamples = {}
LOG_FILENAME = 'error_log.txt'
pd.set_option('display.max_colwidth', -1)


def getNorm(request, RID, stopList, PID):
    try:
        if request.is_ajax():
            # Get variables from web page
            allJson = request.body.split('&')[0]
            all = json.loads(allJson)
            functions.setBase(RID, 'Step 1 of 6: Querying database...')

            NormMeth = int(all["NormMeth"])

            remove = int(all["Remove"])
            cutoff = int(all["Cutoff"])
            Iters = int(all["Iters"])
            Lambda = float(all["Lambda"])
            NormVal = all["NormVal"]
            size_on = int(all["MinSize"])

            # Get selected samples from user's folder and query database for sample info
            myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
            path = str(myDir) + 'usr_sel_samples.pkl'
            with open(path, 'rb') as f:
                selected = pickle.load(f)
            samples = Sample.objects.filter(sampleid__in=selected).values_list('sampleid')

            # Generate a list of sequence reads per sample and filter samples if minimum samplesize
            countList = []
            subList = []
            size = 1
            counts = Sample.objects.filter(sampleid__in=samples)
            if size_on == 1:
                size = int(all["MinVal"])
                for i in counts:
                    if int(i.reads) >= size:
                        countList.append(i.reads)
                        subList.append(i.sampleid)
            else:
                for i in counts:
                    if int(i.reads) is not None:
                        countList.append(i.reads)
                        subList.append(i.sampleid)

            if not countList:
                myDict = {}
                myDict['error'] = "Error with Normalization!\nYour minimum sample has caused all samples to be removed!"
                res = json.dumps(myDict)
                return HttpResponse(res, content_type='application/json')

            # Calculate min/median/max of sequence reads for rarefaction
            NormReads = 0
            if not NormMeth == 1:
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

            result = ''
            result += 'Data Normalization:\n'
            newList = []
            counts = counts.filter(sampleid__in=subList)
            if NormMeth == 2:
                for i in counts:
                    if int(i.reads) >= NormReads:
                        newList.append(i.sampleid)
            else:
                for i in counts:
                    if int(i.reads) > 0:
                        newList.append(i.sampleid)

            if not newList:
                myDict = {}
                myDict['error'] = "Error with Normalization!\nYour sub-sample size has caused all samples to be removed!"
                res = json.dumps(myDict)
                return HttpResponse(res, content_type='application/json')
            else:
                myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
                path = str(myDir) + 'usr_norm_samples.pkl'

                with open(path, 'wb') as f:
                    pickle.dump(newList, f)

            metaDF = UnivMetaDF(newList, RID, stopList, PID)

            # remove emptycols
            metaDF.replace('nan', np.nan, inplace=True)
            metaDF.dropna(axis=1, how='all', inplace=True)
            metaDF.dropna(axis=0, how='all', inplace=True)
            lenB, col = metaDF.shape
            normRem = len(selected) - lenB

            result += str(lenB) + ' selected sample(s) were included in the final analysis...\n'
            result += str(normRem) + ' sample(s) did not met the desired normalization criteria...\n'
            result += '\n'

            # Create unique list of samples in meta dataframe (may be different than selected samples)
            myList = metaDF.index.values.tolist()

            # Create dataframe with all taxa/count data by sample
            taxaDF = functions.taxaProfileDF(myList)

            # Select only the taxa of interest if user used the selectAll button
            taxaDict = {}
            qs3 = Profile.objects.filter(sampleid__in=myList).values_list('otuid', flat='True').distinct()
            taxaDict['OTU_99'] = qs3

            functions.setBase(RID, 'Step 1 of 6: Querying database...done!')

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stopList[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            functions.setBase(RID, 'Step 2 of 6: Sub-sampling data...')

            normDF, DESeq_error = normalizeUniv(taxaDF, taxaDict, myList, NormMeth, NormReads, metaDF, Iters, Lambda, RID, stopList, PID)

            functions.setBase(RID, 'Step 4 of 6: Calculating indices...')

            if remove == 1:
                grouped = normDF.groupby('otuid')
                goodIDs = []
                for name, group in grouped:
                    if group['abund'].sum() > cutoff:
                        goodIDs.append(name)
                normDF = normDF.loc[normDF['otuid'].isin(goodIDs)]

            finalDict = {}
            if NormMeth == 1:
                result += 'No normalization was performed...\n'
            elif NormMeth == 2 or NormMeth == 3:
                result += 'Data were rarefied to ' + str(NormReads) + ' sequence reads with ' + str(Iters) + ' iteration(s)...\n'
            elif NormMeth == 4:
                result += 'Data were normalized by the total number of sequence reads...\n'
            elif NormMeth == 5 and DESeq_error == 'no':
                result += 'Data were normalized by DESeq2...\n'
            elif NormMeth == 5 and DESeq_error == 'yes':
                result += 'DESeq2 cannot run estimateSizeFactors...\n'
                result += 'Analysis was run without normalization...\n'
                result += 'To try again, please select fewer samples or another normalization method...\n'

            if size_on == 1:
                result += "Samples with fewer than " + str(size) + " reads were removed from your analysis...\n"
            else:
                result += "No minimum samples size was applied...\n"

            if remove == 1:
                result += "Phylotypes with fewer than " + str(cutoff) + " read(s) were removed from your analysis\n"
            else:
                result += "No minimum otu size was applied...\n"

            result += '===============================================\n\n\n'
            finalDict['text'] = result

            normDF.set_index('sampleid', inplace=True)
            finalDF = pd.merge(metaDF, normDF, left_index=True, right_index=True, how='inner')
            finalDF.reset_index(drop=False, inplace=True)
            finalDF.rename(columns={'index': 'sampleid'}, inplace=True)

            #re-order finalDF
            metaDFList = list(metaDF.columns.values)
            removeList = ['projectid', 'refid', 'sample_name']
            for x in removeList:
                if x in metaDFList:
                    metaDFList.remove(x)
            metaDFList = ['projectid', 'refid', 'sampleid', 'sample_name'] + metaDFList
            metaDFList = metaDFList + ['kingdom', 'phyla', 'class', 'order', 'family', 'genus', 'species', 'otu', 'otuid', 'abund']
            finalDF = finalDF[metaDFList]

            # save location info to session and save in temp/norm
            myDir = 'myPhyloDB/media/temp/norm/'
            path = str(myDir) + str(RID) + '.pkl'
            request.session['savedDF'] = pickle.dumps(path)
            request.session['NormMeth'] = NormMeth

            functions.setBase(RID, 'Step 5 of 6: Writing data to disk...')

            if not os.path.exists(myDir):
                os.makedirs(myDir)

            functions.setBase(RID, 'Step 5 of 6: Writing data to disk...done')

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stopList[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            functions.setBase(RID, 'Step 6 of 6: Formatting biom data...')

            # regardless of daymet flag, delete old daymetData object for this user (or try to)
            try:
                DaymetData.objects.get(user=request.user).delete()   # isn't request.user a username?
            except Exception as daydelerr:
                #print "Daymet Deletion Error:", daydelerr
                pass
            # TODO turn daymet gathering section into a function in utils or somesuch in case the feature is moved
            # get data from dictionaries before saving to biom
            # check if daymet checkbox is selected
            daymetSuccess = False
            daymetData = None
            sampIDs = []
            daymetKeys = []
            if all['daymet']:
                # do daymet stuff
                minYear = all['minYear']
                maxYear = all['maxYear']
                # get month selection filters
                month_jan = all['month_jan']
                month_feb = all['month_feb']
                month_mar = all['month_mar']
                month_apr = all['month_apr']
                month_may = all['month_may']
                month_jun = all['month_jun']
                month_jul = all['month_jul']
                month_aug = all['month_aug']
                month_sep = all['month_sep']
                month_oct = all['month_oct']
                month_nov = all['month_nov']
                month_dec = all['month_dec']
                # print "User chose daymet, with range", minYear, "to", maxYear
                # use coordinate data to query daymet database AFTER checking 'daymet' flag is set
                coords = {}
                for samp in counts:
                    #print "Samp:", samp.latitude, ":", samp.longitude
                    coords[samp.sampleid] = str(samp.latitude)+":"+str(samp.longitude)
                # given coords dict now, pull up R and query daymet
                # get date range from input? should be from samples somewhere actually

                # get R ready, load daymetR package
                # pass coords and date range into R
                # run daymetR query with dataset provided
                # pull resulting data back into python
                # add section to biom file creation that uses python side data
                # add way to display daymet on norm page? too much data, just add to biom data
                # set up analyses to use daymet data from biom file

                # TODO need support for month selection via checkboxes in page, filter for selected months only DONE
                # we have the months now via 'month_jan' var set
                # need to figure out how to query more specifically
                # query like normal, results have rows by day, count out the months? (183 days per year iirc)

                if os.name == 'nt':
                    r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
                else:
                    r = R(RCMD="R/R-Linux/bin/R", use_pandas=True)

                installed = r('installed.packages()')
                if 'daymetr' not in installed:
                    print r("install.packages('daymetr', repos='http://cran.us.r-project.org', dependencies=T)")
                r("library(daymetr)")
                r.assign("coords", coords)
                r.assign("minYear", minYear)
                r.assign("maxYear", maxYear)
                # actual r code start
                '''
                R starts at 1, so Jan 1st is day 1, Feb 1st is 32, etc
                Daymet includes february 29th on leap years, but removed december 31st to keep the 365 total (WHY!!)
                Table of start and end dates for each month (ranges for days, parentheses are for leap years):
                January 1-31
                February 32-59(60)
                March 60(61)-90(91)
                April 91(92)-120(121)
                May 121(122)-151(152)
                June 152(153)-181(182)
                July 182(183)-212(213)
                August 213(214)-243(244)
                September 244(245)-273(274)
                October 274(275)-304(305)
                November 305(306)-334(335)
                December 335(336)-365(365)

                Leap years between 1980 and 2017: 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016
                '''
                r("daymetData <- list()")
                # making list of valid days in python, handing it to R to subset data with
                month_subset = []
                num_years = int(maxYear) - int(minYear)+1
                currentYear = int(minYear)
                leap = (currentYear % 4 == 0)
                for dayNum in range(1, num_years*365+1):
                    # going from 1 to end+1 instead of 0 to end because
                    # this list is for R, which starts index at 1 instead of 0
                    # and python ranges exclude the listed end element
                    curDay = dayNum % 365

                    if 1 <= curDay <= 31 and month_jan:
                        month_subset.append(dayNum)
                    elif leap:
                        if 32 <= curDay <= 60 and month_feb:
                            month_subset.append(dayNum)
                        elif 61 <= curDay <= 91 and month_mar:
                            month_subset.append(dayNum)
                        elif 92 <= curDay <= 121 and month_apr:
                            month_subset.append(dayNum)
                        elif 122 <= curDay <= 152 and month_may:
                            month_subset.append(dayNum)
                        elif 153 <= curDay <= 182 and month_jun:
                            month_subset.append(dayNum)
                        elif 183 <= curDay <= 213 and month_jul:
                            month_subset.append(dayNum)
                        elif 214 <= curDay <= 244 and month_aug:
                            month_subset.append(dayNum)
                        elif 245 <= curDay <= 274 and month_sep:
                            month_subset.append(dayNum)
                        elif 275 <= curDay <= 305 and month_oct:
                            month_subset.append(dayNum)
                        elif 306 <= curDay <= 335 and month_nov:
                            month_subset.append(dayNum)
                        elif (336 <= curDay <= 364 or curDay == 0) and month_dec:
                            month_subset.append(dayNum)
                    else:
                        if 32 <= curDay <= 59 and month_feb:
                            month_subset.append(dayNum)
                        elif 60 <= curDay <= 90 and month_mar:
                            month_subset.append(dayNum)
                        elif 91 <= curDay <= 120 and month_apr:
                            month_subset.append(dayNum)
                        elif 121 <= curDay <= 151 and month_may:
                            month_subset.append(dayNum)
                        elif 152 <= curDay <= 181 and month_jun:
                            month_subset.append(dayNum)
                        elif 182 <= curDay <= 212 and month_jul:
                            month_subset.append(dayNum)
                        elif 213 <= curDay <= 243 and month_aug:
                            month_subset.append(dayNum)
                        elif 244 <= curDay <= 273 and month_sep:
                            month_subset.append(dayNum)
                        elif 274 <= curDay <= 304 and month_oct:
                            month_subset.append(dayNum)
                        elif 305 <= curDay <= 334 and month_nov:
                            month_subset.append(dayNum)
                        elif (335 <= curDay <= 364 or curDay == 0) and month_dec:
                            month_subset.append(dayNum)
                    if curDay == 0:
                        currentYear += 1
                        leap = (currentYear % 4 == 0)
                r.assign("month_subset", month_subset)
                r('for (sampID in names(coords)){\n'
                    'thisSampData <- list()\n'
                    'lat <- strsplit(coords[[sampID]], ":")[[1]][1]\n'
                    'lon <- strsplit(coords[[sampID]], ":")[[1]][2]\n'
                    'try({data <- download_daymet(lat = lat, lon = lon, start = minYear, end = maxYear, silent = FALSE)})\n'
                    # before sum and mean, remove rows whose days are not in selected months (get specific months figured out)
                    'data <- data[["data"]][month_subset, ]\n'
                    #'data <- subset(data, select=-c("year","yday"))\n'
                    # sum and mean reduce our row count to one per sample, so subsetting must occur beforehand to work

                    # sum and mean of each data based on sample, do while in loop
                    # Sum precip and snow, mean otherwise
                    'dataCols <- names(x=data)\n'
                    'for (col in dataCols){\n'
                    'if (col == "prcp..mm.day." || col == "swe..kg.m.2."){\n'
                    'thisSampData[[col]] <- sum(data[[col]])\n'
                    '} else {\n'
                    'thisSampData[[col]] <- mean(data[[col]])\n'
                    '}\n'
                    '}\n'
                    # update final data with sum/mean of column
                    'daymetData[[sampID]] <- thisSampData\n'
                    '}')
                #print r('warnings()')
                daymetData = r.get("daymetData")
                first = True
                for sampID in daymetData:
                    sampIDs.append(sampID)
                    if first:
                        for dayKey in daymetData[sampID]:
                            daymetKeys.append(dayKey)
                        first = False


                # can confirm, data is good to here, system doesn't break if data fails (catch other errors though?)
                # now take data and get it into biom format (? vs database ?)
                # end goal is to have daymet data selectable in meta quant tree, some data is from biom, some from DB
                # BUT this stuff changes per normalization, so its database entries are inconsistent
                # The key is that the data must be visible from trees.py, ie when tree and its children are created
                # could make a "daymet" object that uses userID to key, data is most recent daymet set
                # would also need a flag to tell if daymet data SHOULD be used, ie if new norm happens without Daymet
                # Biom file hijack could work still, but it seems clunky since we aren't using biom for this step yet
                # so pulling up a file for only one section might be strange
                # make model for users most recent daymet data, its a nested dictionary so could be tricky
                # have ";" delimited strings for sampleID and all daymet columns
                # sync these strings on position, so sampleID[0][dayl..s] = dayl[0], etc

                #print "ALL THE DATA"
                sampID_str = ""
                year_str = ""
                yday_str = ""
                dayl_str = ""
                prcp_str = ""
                srad_str = ""
                swe_str = ""
                tmax_str = ""
                tmin_str = ""
                vp_str = ""
                # TODO drop year and yday columns
    # "year" "yday" "dayl..s." "prcp..mm.day." "srad..W.m.2."  "swe..kg.m.2."  "tmax..deg.c."  "tmin..deg.c."  "vp..Pa."
                for sampID in sampIDs:
                    sampID_str += str(sampID) + ";"
                    year_str += str(daymetData[sampID]["year"]) + ";"
                    yday_str += str(daymetData[sampID]["yday"]) + ";"
                    dayl_str += str(daymetData[sampID]["dayl..s."]) + ";"
                    prcp_str += str(daymetData[sampID]["prcp..mm.day."]) + ";"
                    srad_str += str(daymetData[sampID]["srad..W.m.2."]) + ";"
                    swe_str += str(daymetData[sampID]["swe..kg.m.2."]) + ";"
                    tmax_str += str(daymetData[sampID]["tmax..deg.c."]) + ";"
                    tmin_str += str(daymetData[sampID]["tmin..deg.c."]) + ";"
                    vp_str += str(daymetData[sampID]["vp..Pa."]) + ";"

                myData = DaymetData.objects.create(user=request.user, sampleIDs=sampID_str, year=year_str,
                                                   yday=yday_str, dayl=dayl_str, prcp=prcp_str, srad=srad_str,
                                                   swe=swe_str, tmax=tmax_str, tmin=tmin_str, vp=vp_str)
                myData.save()
                daymetSuccess = True
                #print "Done with daymet data saving"


            myBiom = {}
            nameList = []
            myList.sort()
            for i in myList:
                # sampleID
                nameDict = metaDF.loc[i].to_dict()
                # meta var name as key, value is actual value
                # note: nameDict is a mix of cat and quant vars; the biome file follows this idea
                for key in nameDict:
                    nameDict[key] = str(nameDict[key])
                    # metadata is where daymet should go, so append nameDict? er, add to metaDF.loc?
                if daymetSuccess:
                    if i in sampIDs:
                        for dayKey in daymetKeys:
                            nameDict[dayKey] = str(daymetData[i][dayKey])
                nameList.append({"id": str(i), "metadata": nameDict})

            # get list of lists with abundances
            taxaOnlyDF = finalDF.loc[:, ['sampleid', 'otuid', 'abund']]
            abundDF = taxaOnlyDF.pivot(index='otuid', columns='sampleid', values='abund')
            abundDF.sort_index(axis=0, inplace=True)
            abundList = abundDF.values.tolist()

            # get list of taxa
            namesDF = finalDF.loc[:, ['sampleid', 'otuid']]
            namesDF['taxa'] = finalDF.loc[:, ['kingdom', 'phyla', 'class', 'order', 'family', 'genus', 'species', 'otu']].values.tolist()
            namesDF = namesDF.pivot(index='otuid', columns='sampleid', values='taxa')
            namesDF.sort_index(axis=0, inplace=True)

            taxaList = []
            for index, row in namesDF.iterrows():
                metaDict = {'taxonomy':  row[0]}
                taxaList.append({"id": index, "metadata": metaDict})

            shape = [len(taxaList), len(nameList)]
            myBiom['id'] = 'None'
            myBiom['format'] = 'myPhyloDB v.1.2.0'
            myBiom['format_url'] = 'http://biom-format.org'
            myBiom['generated_by'] = 'myPhyloDB'
            myBiom['type'] = 'OTU table'
            myBiom['date'] = str(datetime.datetime.now())
            myBiom['matrix_type'] = 'dense'
            myBiom['matrix_element_type'] = 'float'
            myBiom["shape"] = shape
            myBiom['rows'] = taxaList
            myBiom['columns'] = nameList
            myBiom['data'] = abundList

            myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
            path = str(myDir) + 'myphylodb.biom'
            with open(path, 'w') as outfile:
                json.dump(myBiom, outfile, ensure_ascii=True, indent=4)

            myBiom = {}
            nameList = []
            myList.sort()
            for i in myList:
                tempPanda = metaDF.fillna("NaN")
                nameDict = tempPanda.loc[i].to_dict()
                for key in nameDict:
                    nameDict[key] = str(nameDict[key])
                nameList.append({"id": str(i), "metadata": nameDict})

            # get list of lists with abundances
            taxaOnlyDF = finalDF.loc[:, ['sampleid', 'otuid', 'abund']]
            abundDF = taxaOnlyDF.pivot(index='otuid', columns='sampleid', values='abund')
            abundDF.sort_index(axis=0, inplace=True)
            abundList = abundDF.values.tolist()

            # get list of taxa
            namesDF = finalDF.loc[:, ['sampleid', 'otuid']]
            namesDF['taxa'] = finalDF.loc[:, ['kingdom', 'phyla', 'class', 'order', 'family', 'genus', 'species', 'otu']].values.tolist()
            namesDF = namesDF.pivot(index='otuid', columns='sampleid', values='taxa')
            namesDF.sort_index(axis=0, inplace=True)

            taxaList = []
            for index, row in namesDF.iterrows():
                taxonomy = row[0][:-1]
                taxonomy = [i.split(': ', 1)[1] for i in taxonomy]
                taxonomy = filter(lambda a: a != 'unclassified', taxonomy)
                metaDict = {'taxonomy': taxonomy}
                taxaList.append({"id": row[0][-1].split(': ')[1], "metadata": metaDict})

            shape = [len(taxaList), len(nameList)]
            myBiom['id'] = 'None'
            myBiom['format'] = 'Biological Observation Matrix 1.0.0'
            myBiom['format_url'] = 'http://biom-format.org'
            myBiom['generated_by'] = 'myPhyloDB'
            myBiom['type'] = 'OTU table'
            myBiom['date'] = str(datetime.datetime.now())
            myBiom['matrix_type'] = 'dense'
            myBiom['matrix_element_type'] = 'float'
            myBiom["shape"] = shape
            myBiom['rows'] = taxaList
            myBiom['columns'] = nameList
            myBiom['data'] = abundList

            myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
            path = str(myDir) + 'phyloseq.biom'
            with open(path, 'w') as outfile:
                json.dump(myBiom, outfile, ensure_ascii=True, indent=4)

            functions.setBase(RID, 'Step 6 of 6: Formatting biome data...done!')

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stopList[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

            finalDict['error'] = 'none'
            res = json.dumps(finalDict)
            return HttpResponse(res, content_type='application/json')

    except Exception as e:
        if not stopList[PID] == RID:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,)
            myDate = "\nDate: " + str(datetime.datetime.now()) + "\n"
            logging.exception(myDate)
            myDict = {}
            myDict['error'] = "There was an error during your normalization:\nError: " + str(e.message) + "\nTimestamp: " + str(datetime.datetime.now())
            res = json.dumps(myDict)
            return HttpResponse(res, content_type='application/json')


def UnivMetaDF(sampleList, RID, stopList, PID):
    tableNames = ['project', 'reference', 'sample', 'air', 'human_associated', 'microbial', 'soil', 'water', 'userdefined', 'profile']
    idList = ['projectid', 'refid', u'id']

    your_fields = Sample._meta.local_fields
    sampleTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in sampleTableList:
            sampleTableList.remove(x)

    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
    if stopList[PID] == RID:
        res = ''
        return HttpResponse(res, content_type='application/json')
    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    your_fields = Air._meta.local_fields
    airTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in airTableList:
            airTableList.remove(x)
    for x in idList:
        if x in airTableList:
            airTableList.remove(x)

    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
    if stopList[PID] == RID:
        res = ''
        return HttpResponse(res, content_type='application/json')
    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    your_fields = Human_Associated._meta.local_fields
    humanAssocTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in humanAssocTableList:
            humanAssocTableList.remove(x)
    for x in idList:
        if x in humanAssocTableList:
            humanAssocTableList.remove(x)

    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
    if stopList[PID] == RID:
        res = ''
        return HttpResponse(res, content_type='application/json')
    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    your_fields = Microbial._meta.local_fields
    microbialTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in microbialTableList:
            microbialTableList.remove(x)
    for x in idList:
        if x in microbialTableList:
            microbialTableList.remove(x)

    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
    if stopList[PID] == RID:
        res = ''
        return HttpResponse(res, content_type='application/json')
    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    your_fields = Soil._meta.local_fields
    soilTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in soilTableList:
            soilTableList.remove(x)
    for x in idList:
        if x in soilTableList:
            soilTableList.remove(x)

    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
    if stopList[PID] == RID:
        res = ''
        return HttpResponse(res, content_type='application/json')
    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    your_fields = Water._meta.local_fields
    waterTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in waterTableList:
            waterTableList.remove(x)
    for x in idList:
        if x in waterTableList:
            waterTableList.remove(x)

    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
    if stopList[PID] == RID:
        res = ''
        return HttpResponse(res, content_type='application/json')
    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    your_fields = UserDefined._meta.local_fields
    usrTableList = [f.name for f in your_fields]
    for x in tableNames:
        if x in usrTableList:
            usrTableList.remove(x)
    for x in idList:
        if x in usrTableList:
            usrTableList.remove(x)

    metaDF = pd.DataFrame(list(Sample.objects.filter(sampleid__in=sampleList).values(*sampleTableList)))
    metaDF.set_index('sampleid', drop=True, inplace=True)

    tempDF = pd.DataFrame(list(Air.objects.filter(sampleid_id__in=sampleList).values(*airTableList)))
    if not tempDF.empty:
        tempDF.set_index('sampleid', drop=True, inplace=True)
        metaDF = pd.merge(metaDF, tempDF, left_index=True, right_index=True, how='inner')

    tempDF = pd.DataFrame(list(Human_Associated.objects.filter(sampleid_id__in=sampleList).values(*humanAssocTableList)))
    if not tempDF.empty:
        tempDF.set_index('sampleid', drop=True, inplace=True)
        metaDF = pd.merge(metaDF, tempDF, left_index=True, right_index=True, how='inner')

    tempDF = pd.DataFrame(list(Microbial.objects.filter(sampleid_id__in=sampleList).values(*microbialTableList)))
    if not tempDF.empty:
        tempDF.set_index('sampleid', drop=True, inplace=True)
        metaDF = pd.merge(metaDF, tempDF, left_index=True, right_index=True, how='inner')

    tempDF = pd.DataFrame(list(Soil.objects.filter(sampleid_id__in=sampleList).values(*soilTableList)))
    if not tempDF.empty:
        tempDF.set_index('sampleid', drop=True, inplace=True)
        metaDF = pd.merge(metaDF, tempDF, left_index=True, right_index=True, how='inner')

    tempDF = pd.DataFrame(list(Water.objects.filter(sampleid_id__in=sampleList).values(*waterTableList)))
    if not tempDF.empty:
        tempDF.set_index('sampleid', drop=True, inplace=True)
        metaDF = pd.merge(metaDF, tempDF, left_index=True, right_index=True, how='inner')

    tempDF = pd.DataFrame(list(UserDefined.objects.filter(sampleid_id__in=sampleList).values(*usrTableList)))
    if not tempDF.empty:
        tempDF.set_index('sampleid', drop=True, inplace=True)
        metaDF = pd.merge(metaDF, tempDF, left_index=True, right_index=True, how='inner')

    return metaDF


def normalizeUniv(df, taxaDict, mySet, meth, reads, metaDF, iters, Lambda, RID, stopList, PID):
    global curSamples, totSamples
    df2 = df.reset_index()
    taxaID = ['kingdomid', 'phylaid', 'classid', 'orderid', 'familyid', 'genusid', 'speciesid', 'otuid']

    countDF = pd.DataFrame()
    DESeq_error = 'no'
    if meth == 1 or meth == 4:
        countDF = df2.copy()

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stopList[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    elif meth == 2:
        if reads >= 0:
            countDF = df2[taxaID].reset_index(drop=True)

            curSamples[RID] = 0
            totSamples[RID] = len(mySet)

            myArr = df2.loc[:, mySet].as_matrix()
            probArr = myArr.T
            finalArr = rarefaction_remove(probArr, RID, reads=reads, iters=iters)
            for i in xrange(len(mySet)):
                countDF[mySet[i]] = finalArr[i]

            curSamples[RID] = 0

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stopList[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif reads < 0:
            countDF = df2.reset_index(drop=True)

    elif meth == 3:
        if reads >= 0:
            countDF = df2[taxaID].reset_index(drop=True)

            curSamples[RID] = 0
            totSamples[RID] = len(mySet)

            myArr = df2.loc[:, mySet].as_matrix()
            probArr = myArr.T
            finalArr = rarefaction_keep(probArr, RID, reads=reads, iters=iters, myLambda=Lambda)
            for i in xrange(len(mySet)):
                countDF[mySet[i]] = finalArr[i]

            curSamples[RID] = 0

            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
            if stopList[PID] == RID:
                res = ''
                return HttpResponse(res, content_type='application/json')
            # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        elif reads < 0:
            countDF = df2.reset_index(drop=True)

    elif meth == 5:
        countDF = df2[taxaID].reset_index(drop=True)
        if os.name == 'nt':
            r = R(RCMD="R/R-Portable/App/R-Portable/bin/R.exe", use_pandas=True)
        else:
            r = R(RCMD="R/R-Linux/bin/R")

        functions.setBase(RID, 'Verifying R packages...missing packages are being installed')

        r("list.of.packages <- c('DESeq2')")
        r("new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]")
        r("if (length(new.packages)) source('http://bioconductor.org/biocLite.R')")
        print r("if (length(new.packages)) biocLite(new.packages, type='source', suppressUpdate=T, dependencies=T)")

        functions.setBase(RID, 'Step 2 of 6: Sub-sampling data...')

        print r("library(DESeq2)")

        df3 = df2.drop(taxaID, axis=1)

        r.assign("count", df3)
        r.assign("metaDF", metaDF)

        r("rows <- rownames(metaDF)")
        r("colnames(count) <- rows")

        r("trt <- factor(metaDF$merge)")

        r("dds <- DESeqDataSetFromMatrix(countData = count, colData = metaDF, design = ~ sample_name)")

        r("dds <- estimateSizeFactors(dds)")
        pycds = r.get("sizeFactors(dds)")
        colList = df3.columns.tolist()

        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
        if stopList[PID] == RID:
            res = ''
            return HttpResponse(res, content_type='application/json')
        # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

        found = False
        if pycds is list:
            for thing in pycds:
                if str(thing) == "None":
                    found = True
        else:
            if pycds is None:
                found = True

        if not found:
            DESeq_error = 'no'
            cdsDF = pd.DataFrame(r.get("counts(dds, normalize=TRUE)"), columns=[colList])
            countDF[colList] = cdsDF[colList]
        else:
            DESeq_error = 'yes'
            countDF = df2.reset_index(drop=True)

    functions.setBase(RID, 'Step 2 of 6: Sub-sampling data...done!')
    functions.setBase(RID, 'Step 3 of 6: Tabulating data...')

    field = 'otuid'
    taxaList = taxaDict['OTU_99']
    qs = OTU_99.objects.filter(otuid__in=taxaList)
    namesDF = read_frame(qs, fieldnames=['kingdomid__kingdomName', 'phylaid__phylaName', 'classid__className', 'orderid__orderName', 'familyid__familyName', 'genusid__genusName', 'speciesid__speciesName', 'otuName', 'kingdomid__kingdomid', 'phylaid__phylaid', 'classid__classid', 'orderid__orderid', 'familyid__familyid', 'genusid__genusid', 'speciesid__speciesid', 'otuid'])
    namesDF.rename(columns={'kingdomid__kingdomid': 'kingdomid', 'phylaid__phylaid': 'phylaid', 'classid__classid' : 'classid', 'orderid__orderid' : 'orderid', 'familyid__familyid' : 'familyid', 'genusid__genusid' : 'genusid', 'speciesid__speciesid' : 'speciesid'}, inplace=True)
    namesDF.rename(columns={'kingdomid__kingdomName': 'kingdomName', 'phylaid__phylaName': 'phylaName', 'classid__className' : 'className', 'orderid__orderName' : 'orderName', 'familyid__familyName' : 'familyName', 'genusid__genusName' : 'genusName', 'speciesid__speciesName' : 'speciesName'}, inplace=True)
    namesDF['kingdom'] = namesDF['kingdomid'].astype(str) + ': ' + namesDF['kingdomName']
    namesDF.drop(['kingdomid', 'kingdomName'], axis=1, inplace=True)
    namesDF['phyla'] = namesDF['phylaid'].astype(str) + ': ' + namesDF['phylaName']
    namesDF.drop(['phylaid', 'phylaName'], axis=1, inplace=True)
    namesDF['class'] = namesDF['classid'].astype(str) + ': ' + namesDF['className']
    namesDF.drop(['classid', 'className'], axis=1, inplace=True)
    namesDF['order'] = namesDF['orderid'].astype(str) + ': ' + namesDF['orderName']
    namesDF.drop(['orderid', 'orderName'], axis=1, inplace=True)
    namesDF['family'] = namesDF['familyid'].astype(str) + ': ' + namesDF['familyName']
    namesDF.drop(['familyid', 'familyName'], axis=1, inplace=True)
    namesDF['genus'] = namesDF['genusid'].astype(str) + ': ' + namesDF['genusName']
    namesDF['species'] = namesDF['speciesid'].astype(str) + ': ' + namesDF['genusName'] + " " + namesDF['speciesName']
    namesDF.drop(['genusid', 'genusName'], axis=1, inplace=True)
    namesDF.drop(['speciesid', 'speciesName'], axis=1, inplace=True)
    namesDF['otu'] = namesDF['otuid'].astype(str) + ': ' + namesDF['otuName']
    namesDF.drop(['otuName'], axis=1, inplace=True)
    namesDF.replace('unclassified unclassified', 'unclassified', regex=True, inplace=True)

    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #
    if stopList[PID] == RID:
        res = ''
        return HttpResponse(res, content_type='application/json')
    # /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\ #

    curSamples[RID] = 0
    totSamples[RID] = len(mySet)

    ser1 = countDF.groupby(field)[mySet].sum()
    ser2 = ser1.stack()
    abundDF = pd.Series.to_frame(ser2, name='abund')
    abundDF.reset_index(drop=False, inplace=True)

    normDF = pd.merge(abundDF, namesDF, left_on='otuid', right_on='otuid', how='inner')

    return normDF, DESeq_error


def rarefaction_remove(M, RID, reads=0, iters=0):
    global curSamples, totSamples
    nsamp = M.shape[0]

    Mrarefied = np.empty_like(M)
    for i in range(nsamp):
        counts = M[i]
        nz = counts.nonzero()[0]
        unpacked = np.concatenate([np.repeat(np.array(j,), counts[j]) for j in nz])
        myArr = np.zeros(len(counts), dtype=int)
        for n in xrange(iters):
            permuted = np.random.permutation(unpacked)[:reads]
            binArr = np.zeros(len(counts), dtype=int)
            for p in permuted:
                binArr[p] += 1
            if n == 0:
                myArr = binArr
            else:
                myArr = np.vstack((myArr, binArr))

        if iters > 1:
            Mrarefied[i] = np.mean(myArr, axis=0)
        else:
            Mrarefied[i] = myArr
        curSamples[RID] += 1
        functions.setBase(RID, 'Step 2 of 6: Sub-sampling data...\nSub-sampling is complete for ' + str(curSamples[RID]) + ' out of ' + str(totSamples[RID]) + ' samples')
    return Mrarefied


def rarefaction_keep(M, RID, reads=0, iters=0, myLambda=0.1):
    global curSamples, totSamples
    noccur = np.sum(M, axis=1)  # number of occurrences for each sample
    nvar = M.shape[1]  # number of variables
    nsamp = M.shape[0]  # number of samples

    Mrarefied = np.empty_like(M)
    for i in range(nsamp):
        p = (M[i] + myLambda) / (float(noccur[i]) + nvar * myLambda)
        myArr = np.zeros(nvar)
        for n in xrange(iters):
            prng = RandomState()
            choice = prng.choice(nvar, size=reads, replace=True, p=p)
            binArr = np.bincount(choice, minlength=nvar)
            if n == 0:
                myArr = binArr
            else:
                myArr = np.vstack((myArr, binArr))

        if iters > 1:
            Mrarefied[i] = np.mean(myArr, axis=0)
        else:
            Mrarefied[i] = myArr
        curSamples[RID] += 1
        functions.setBase(RID, 'Step 2 of 6: Sub-sampling data...\nSub-sampling is complete for ' + str(curSamples[RID]) + ' out of ' + str(totSamples[RID]) + ' samples')
    return Mrarefied


def getTab(request):
    if request.is_ajax():
        myDict = {}
        myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
        path = str(myDir) + 'myphylodb.biom'

        df, metaDF, remCatFields = functions.exploding_panda(path)

        fileName2 = str(myDir) + 'usr_norm_data.csv'
        df.to_csv(fileName2)

        myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
        fileName3 = str(myDir) + 'usr_norm_data.csv.gz'
        zf = zipfile.ZipFile(fileName3, "w", zipfile.ZIP_DEFLATED, allowZip64=True)
        zf.write(fileName2, 'usr_norm_data.csv')
        zf.close()

        myDir = '../../myPhyloDB/media/usr_temp/' + str(request.user) + '/'
        fileName3 = str(myDir) + 'usr_norm_data.csv.gz'
        myDict['name'] = str(fileName3)
        res = json.dumps(myDict)
        return HttpResponse(res, content_type='application/json')


def getBiom(request):
    if request.is_ajax():
        myDict = {}
        myDir = 'myPhyloDB/media/usr_temp/' + str(request.user) + '/'
        fileName1 = str(myDir) + 'myphylodb.biom'
        fileName2 = str(myDir) + 'phyloseq.biom'

        zip_file = os.path.join(os.getcwd(), 'myPhyloDB', 'media', 'usr_temp', request.user.username, 'usr_norm_data.biom.gz')
        zf = zipfile.ZipFile(zip_file, "w", zipfile.ZIP_DEFLATED, allowZip64=True)

        zf.write(fileName1)
        zf.write(fileName2)
        zf.close()

        myDir = '../../myPhyloDB/media/usr_temp/' + str(request.user) + '/'
        fileName3 = str(myDir) + 'usr_norm_data.biom.gz'
        myDict['name'] = str(fileName3)
        res = json.dumps(myDict)
        return HttpResponse(res, content_type='application/json')
