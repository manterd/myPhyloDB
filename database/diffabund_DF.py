import operator
import pandas as pd
from django.db.models import Q
from numpy import *
from pyper import *


def catDiffAbundDF(qs1, metaDict):
    sampleTableList = ['sample_name', 'organism', 'material', 'seq_platform', 'seq_gene', 'seq_gene_region', 'seq_for_primer', 'seq_rev_primer', 'collection_date', 'biome', 'feature', 'geo_loc_country', 'geo_loc_state', 'geo_loc_city', 'geo_loc_farm', 'geo_loc_plot']
    soilTableList = ['samp_collection_device', 'samp_size', 'samp_depth', 'sieve_size', 'storage_cond', 'samp_weight_dna_ext', 'pool_dna_extracts', 'fao_class', 'local_class', 'texture_class', 'porosity', 'profile_position', 'slope_aspect', 'slope_gradient', 'bulk_density', 'drainage_class', 'water_content_soil', 'cur_land_use', 'cur_vegetation', 'cur_crop', 'cur_cultivar', 'crop_rotation', 'cover_crop', 'fert_amendment_class', 'fert_placement', 'fert_type', 'fert_tot_amount', 'fert_N_tot_amount', 'fert_P_tot_amount', 'fert_K_tot_amount', 'irrigation_type', 'irrigation_tot_amount', 'residue_removal', 'residue_growth_stage', 'residue_removal_percent', 'tillage_event', 'tillage_event_depth', 'amend1_class', 'amend1_active_ingredient', 'amend1_tot_amount', 'amend2_class', 'amend2_active_ingredient', 'amend2_tot_amount', 'amend3_class', 'amend3_active_ingredient', 'amend3_tot_amount', 'rRNA_copies', 'microbial_biomass_C', 'microbial_biomass_N', 'microbial_respiration', 'soil_pH', 'soil_EC', 'soil_C', 'soil_OM', 'soil_N', 'soil_NO3_N', 'soil_NH4_N', 'soil_P', 'soil_K', 'soil_S', 'soil_Zn', 'soil_Fe', 'soil_Cu', 'soil_Mn', 'soil_Ca', 'soil_Mg', 'soil_Na', 'soil_B', 'plant_C', 'plant_N', 'plant_P', 'plant_K', 'plant_Ca', 'plant_Mg', 'plant_S', 'plant_Na', 'plant_Cl', 'plant_Al', 'plant_B', 'plant_Cu', 'plant_Fe', 'plant_Mn', 'plant_Zn', 'crop_tot_biomass_fw', 'crop_tot_biomass_dw', 'crop_tot_above_biomass_fw', 'crop_tot_above_biomass_dw', 'crop_tot_below_biomass_fw', 'crop_tot_below_biomass_dw', 'harv_fraction', 'harv_fresh_weight', 'harv_dry_weight', 'ghg_chamber_placement', 'ghg_N2O', 'ghg_CO2', 'ghg_NH4']
    microbialTableList = ['alkalinity', 'alkyl_diethers', 'altitude', 'aminopept_act', 'ammonium', 'bacteria_carb_prod', 'biomass', 'bishomohopanol', 'bromide', 'calcium', 'carb_nitro_ratio', 'chem_administration', 'chloride', 'chlorophyll', 'diether_lipids', 'diss_carb_dioxide', 'diss_hydrogen', 'diss_inorg_carb', 'diss_org_carb', 'diss_org_nitro', 'diss_oxygen', 'glucosidase_act', 'magnesium', 'mean_frict_vel', 'mean_peak_frict_vel', 'methane', 'n_alkanes', 'nitrate', 'nitrite', 'nitro', 'org_carb', 'org_matter', 'org_nitro', 'organism_count', 'oxy_stat_samp', 'part_org_carb', 'perturbation', 'petroleum_hydrocarb', 'ph', 'phaeopigments', 'phosphate', 'phosplipid_fatt_acid', 'potassium', 'pressure', 'redox_potential', 'rel_to_oxygen', 'salinity', 'samp_collect_device', 'samp_mat_process', 'samp_size', 'samp_store_dur', 'samp_store_loc', 'samp_store_temp', 'silicate', 'sodium', 'sulfate', 'sulfide', 'temp', 'tot_carb', 'tot_nitro', 'tot_org_carb', 'turbidity', 'water_content', 'user_defined']
    human_associatedTableList = ['samp_collect_device', 'samp_mat_process', 'samp_size', 'samp_store_temp', 'samp_store_dur', 'samp_type', 'samp_location', 'samp_temp', 'samp_ph', 'samp_oxy_stat', 'samp_salinity', 'host_subject_id', 'host_age', 'host_pulse', 'host_gender', 'host_ethnicity', 'host_height', 'host_weight', 'host_bmi', 'host_weight_loss_3_month', 'host_body_temp', 'host_occupation', 'pet_farm_animal', 'smoker', 'diet_type', 'diet_duration', 'diet_frequency', 'diet_last_six_month', 'last_meal', 'medic_hist_perform', 'disease_type', 'disease_location', 'disease_duration', 'organism_count', 'tumor_location', 'tumor_mass', 'tumor_stage', 'drug_usage', 'drug_type', 'drug_duration', 'drug_frequency', 'perturbation', 'pert_type', 'pert_duration', 'pert_frequency', 'fetal_health_stat', 'amniotic_fluid_color', 'gestation_stat', 'maternal_health_stat']
    waterTableList = ['alkalinity', 'alkyl_diethers', 'altitude', 'aminopept_act', 'ammonium', 'atmospheric_data', 'bac_prod', 'bac_resp', 'bacteria_carb_prod', 'biomass', 'bishomohopanol', 'bromide', 'calcium', 'carb_nitro_ratio', 'chem_administration', 'chloride', 'chlorophyll', 'conduc', 'density', 'diether_lipids', 'diss_carb_dioxide', 'diss_hydrogen', 'diss_inorg_carb', 'diss_inorg_nitro', 'diss_inorg_phosp', 'diss_org_carb', 'diss_org_nitro', 'diss_oxygen', 'down_par', 'elev', 'fluor', 'glucosidase_act', 'light_intensity', 'magnesium', 'mean_frict_vel', 'mean_peak_frict_vel', 'n_alkanes', 'nitrate', 'nitrite', 'nitro', 'org_carb', 'org_matter', 'org_nitro', 'organism_count', 'oxy_stat_samp', 'part_org_carb', 'part_org_nitro', 'perturbation', 'pretroleum_hydrocarb', 'ph', 'phaeopigments', 'phosphate', 'phosplipid_fatt_acid', 'photon_flux', 'potassium', 'pressure', 'primary_prod', 'redox_potential', 'rel_to_oxygen', 'samp_mat_process', 'samp_salinity', 'samp_size', 'samp_store_dur', 'samp_store_loc', 'samp_store_temp', 'samp_vol_we_dna_ext', 'silicate', 'sodium', 'soluble_react_phosp', 'source_material_id', 'sulfate', 'sulfide', 'suspen_part_matter', 'temp', 'tidal_stage', 'tot_depth_water_col', 'tot_diss_nitro', 'tot_inorg_nitro', 'tot_nitro', 'tot_part_carb', 'tot_phosp', 'water_current', 'user_defined']
    airTableList = ['barometric_press', 'carb_dioxide', 'carb_monoxide', 'chem_administration', 'elev', 'humidity', 'methane', 'organism_count', 'oxy_stat_samp', 'oxygen', 'perturbation', 'pollutants', 'rel_to_oxygen', 'resp_part_matter', 'samp_collect_device', 'samp_mat_process', 'samp_salinity', 'samp_size', 'samp_store_dur', 'samp_store_loc', 'samp_sotre_temp', 'solar_irradiance', 'temp', 'ventilation_rate', 'ventiliation_type', 'volatile_org_comp', 'wind_direction', 'wind_speed', 'user_defined']
    usrTableList = ['usr_cat1', 'usr_cat2', 'usr_cat3', 'usr_cat4', 'usr_cat5', 'usr_cat6']

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
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
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
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key in waterTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'water__' + str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key in airTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'air__' + str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: key}, inplace=True)
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
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key in microbialTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'microbial__' + str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

        elif key in usrTableList:
            field_list.append('sampleid')
            field_list.append('sample_name')
            field = 'user__' + str(key)
            field_list.append(field)
            if type(value) is unicode:
                args_list.append(Q(**{field: value}))
            else:
                for item in value:
                    args_list.append(Q(**{field: item}))
            exclude_list = []
            exclude_list.append(Q(**{field: 'null'}))
            qs2 = qs1.filter(reduce(operator.or_, args_list)).values(*field_list).exclude(reduce(operator.or_, exclude_list))
            tempDF = pd.DataFrame.from_records(qs2, columns=field_list)
            tempDF.rename(columns={field: key}, inplace=True)
            if metaDF.empty:
                metaDF = tempDF
            else:
                tempDF.drop('sample_name', axis=1, inplace=True)
                metaDF = metaDF.merge(tempDF, on='sampleid', how='outer')

    return metaDF


