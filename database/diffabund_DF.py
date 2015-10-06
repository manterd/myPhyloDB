import operator
import pandas as pd
from django.db.models import Q
from database.models import Sample, Human_Associated, Soil, UserDefined
import numpy as np
from numpy import *
from pyper import *


def catDiffAbundDF(idDict):
    sampleTableList = Sample._meta.get_all_field_names()
    human_associatedTableList = Human_Associated._meta.get_all_field_names()
    soilTableList = Soil._meta.get_all_field_names()
    usrTableList = UserDefined._meta.get_all_field_names()

    metaDF = pd.DataFrame()
    idList = []
    fieldList = []
    for key in idDict:
        fieldList.append(key)
        idList.extend(idDict[key])

    metaDF = pd.DataFrame()
    index = 0
    for key in idDict:
        tempDF = pd.DataFrame()
        mySet = idDict[key]

        if key in sampleTableList:
            fields = ['sampleid', 'sample_name', key]
            qs2 = Sample.objects.filter(sampleid__in=mySet).values(*fields)
            tempDF = pd.DataFrame.from_records(qs2, columns=fields).dropna()

        elif key in human_associatedTableList:
            fields = ['sampleid', 'sampleid__sample_name', key]
            qs2 = Human_Associated.objects.filter(sampleid__in=mySet).values(*fields)
            tempDF = pd.DataFrame.from_records(qs2, columns=fields).dropna()
            tempDF.rename(columns={'sampleid__sample_name': 'sample_name'}, inplace=True)

        elif key in soilTableList:
            fields = ['sampleid', 'sampleid__sample_name', key]
            qs2 = Soil.objects.filter(sampleid__in=mySet).values(*fields)
            tempDF = pd.DataFrame.from_records(qs2, columns=fields).dropna()
            tempDF.rename(columns={'sampleid__sample_name': 'sample_name'}, inplace=True)

        elif key in usrTableList:
            fields = ['sampleid', 'sampleid__sample_name', key]
            qs2 = UserDefined.objects.filter(sampleid__in=mySet).values(*fields)
            tempDF = pd.DataFrame.from_records(qs2, columns=fields).dropna()
            tempDF.rename(columns={'sampleid__sample_name': 'sample_name'}, inplace=True)

        tempDF.set_index('sampleid', inplace=True)

        if index == 0:
            metaDF = tempDF
        else:
            cols_to_use = tempDF.columns - metaDF.columns
            metaDF = pd.merge(metaDF, tempDF[cols_to_use], left_index=True, right_index=True)
        index += 1

    return metaDF


