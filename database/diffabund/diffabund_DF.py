import pandas as pd
from database.models import Sample, Air, Human_Associated, Microbial, Soil, Water, UserDefined


def catDiffAbundDF(idDict):
    sampleTableList = Sample._meta.get_all_field_names()
    airTableList = Air._meta.get_all_field_names()
    human_associatedTableList = Human_Associated._meta.get_all_field_names()
    microbialTableList = Microbial._meta.get_all_field_names()
    soilTableList = Soil._meta.get_all_field_names()
    waterTableList = Water._meta.get_all_field_names()
    usrTableList = UserDefined._meta.get_all_field_names()

    metaDF = pd.DataFrame()
    index = 0
    for key in idDict:
        tempDF = pd.DataFrame()
        mySet = idDict[key]

        if key in sampleTableList:
            if key == 'sample_name':
                fields = ['sampleid', key]
            else:
                fields = ['sampleid', 'sample_name', key]
            qs2 = Sample.objects.filter(sampleid__in=mySet).values(*fields)
            tempDF = pd.DataFrame.from_records(qs2, columns=fields).dropna()

        elif key in airTableList:
            fields = ['sampleid', 'sampleid__sample_name', key]
            qs2 = Air.objects.filter(sampleid__in=mySet).values(*fields)
            tempDF = pd.DataFrame.from_records(qs2, columns=fields).dropna()
            tempDF.rename(columns={'sampleid__sample_name': 'sample_name'}, inplace=True)

        elif key in human_associatedTableList:
            fields = ['sampleid', 'sampleid__sample_name', key]
            qs2 = Human_Associated.objects.filter(sampleid__in=mySet).values(*fields)
            tempDF = pd.DataFrame.from_records(qs2, columns=fields).dropna()
            tempDF.rename(columns={'sampleid__sample_name': 'sample_name'}, inplace=True)

        elif key in microbialTableList:
            fields = ['sampleid', 'sampleid__sample_name', key]
            qs2 = Microbial.objects.filter(sampleid__in=mySet).values(*fields)
            tempDF = pd.DataFrame.from_records(qs2, columns=fields).dropna()
            tempDF.rename(columns={'sampleid__sample_name': 'sample_name'}, inplace=True)

        elif key in soilTableList:
            fields = ['sampleid', 'sampleid__sample_name', key]
            qs2 = Soil.objects.filter(sampleid__in=mySet).values(*fields)
            tempDF = pd.DataFrame.from_records(qs2, columns=fields).dropna()
            tempDF.rename(columns={'sampleid__sample_name': 'sample_name'}, inplace=True)

        elif key in waterTableList:
            fields = ['sampleid', 'sampleid__sample_name', key]
            qs2 = Water.objects.filter(sampleid__in=mySet).values(*fields)
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


