import operator
import pandas as pd
from django.db.models import Q
from database.models import Sample, Human_Associated, Soil, UserDefined
import numpy as np
from numpy import *
from pyper import *


def catDiffAbundDF(qs1, metaDict):
    sampleTableList = Sample._meta.get_all_field_names()
    human_associatedTableList = Human_Associated._meta.get_all_field_names()
    soilTableList = Soil._meta.get_all_field_names()
    usrTableList = UserDefined._meta.get_all_field_names()

    metaDF = pd.DataFrame()
    for key in metaDict:
        value = metaDict[key]
        args_list = []
        field_list = []

        if key in sampleTableList:
            field_list.append('sampleid')
            if key != 'sample_name':
                field_list.append('sample_name')
            field = str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list)
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list).dropna()
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key in human_associatedTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'human_associated__' + str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list)
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list).dropna()
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key in soilTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'soil__' + str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list)
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list).dropna()
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key in usrTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'userdefined__' + str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list)
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list).dropna()
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

    return metaDF


