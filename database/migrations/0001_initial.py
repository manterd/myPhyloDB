# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Project'
        db.create_table(u'database_project', (
            ('projectid', self.gf('django.db.models.fields.CharField')(max_length=36, primary_key=True)),
            ('status', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('projectType', self.gf('django.db.models.fields.CharField')(max_length=45)),
            ('project_name', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('project_desc', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('start_date', self.gf('django.db.models.fields.CharField')(max_length=15, blank=True)),
            ('end_date', self.gf('django.db.models.fields.CharField')(max_length=15, blank=True)),
            ('pi_last', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('pi_first', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('pi_affiliation', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('pi_email', self.gf('django.db.models.fields.EmailField')(max_length=75, blank=True)),
            ('pi_phone', self.gf('django.db.models.fields.CharField')(max_length=15, blank=True)),
        ))
        db.send_create_signal(u'database', ['Project'])

        # Adding model 'Reference'
        db.create_table(u'database_reference', (
            ('refid', self.gf('django.db.models.fields.CharField')(max_length=36, primary_key=True)),
            ('raw', self.gf('django.db.models.fields.BooleanField')()),
            ('projectid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Project'])),
            ('path', self.gf('django.db.models.fields.CharField')(max_length=90)),
            ('source', self.gf('django.db.models.fields.CharField')(max_length=90)),
            ('alignDB', self.gf('django.db.models.fields.CharField')(max_length=90, blank=True)),
            ('templateDB', self.gf('django.db.models.fields.CharField')(max_length=90, blank=True)),
            ('taxonomyDB', self.gf('django.db.models.fields.CharField')(max_length=90, blank=True)),
            ('author', self.gf('django.db.models.fields.related.ForeignKey')(related_name='entries', to=orm['auth.User'])),
        ))
        db.send_create_signal(u'database', ['Reference'])

        # Adding model 'Sample'
        db.create_table(u'database_sample', (
            ('projectid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Project'])),
            ('refid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Reference'])),
            ('sampleid', self.gf('django.db.models.fields.CharField')(max_length=36, primary_key=True)),
            ('sample_name', self.gf('django.db.models.fields.TextField')()),
            ('organism', self.gf('django.db.models.fields.CharField')(max_length=90, blank=True)),
            ('collection_date', self.gf('django.db.models.fields.CharField')(max_length=15, blank=True)),
            ('depth', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('elev', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('seq_platform', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('seq_gene', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('seq_gene_region', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('seq_for_primer', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('seq_rev_primer', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('env_biome', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('env_feature', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('env_material', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('geo_loc_country', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('geo_loc_state', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('geo_loc_city', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('geo_loc_farm', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('geo_loc_plot', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('latitude', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('longitude', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('annual_season_precpt', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('annual_season_temp', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'database', ['Sample'])

        # Adding model 'Human_Associated'
        db.create_table(u'database_human_associated', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sampleid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Sample'])),
            ('projectid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Project'])),
            ('refid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Reference'])),
            ('samp_collect_device', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('samp_mat_process', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('samp_size', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('samp_store_temp', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('samp_store_dur', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('samp_type', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('samp_location', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('samp_temp', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('samp_ph', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('samp_oxy_stat', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('samp_salinity', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('host_subject_id', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('host_age', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('host_pulse', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('host_gender', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('host_ethnicity', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('host_height', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('host_weight', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('host_bmi', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('host_weight_loss_3_month', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('host_body_temp', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('host_occupation', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('pet_farm_animal', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('smoker', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('diet_type', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('diet_duration', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('diet_frequency', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('diet_last_six_month', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('last_meal', self.gf('django.db.models.fields.CharField')(max_length=15, blank=True)),
            ('medic_hist_perform', self.gf('django.db.models.fields.CharField')(max_length=15, blank=True)),
            ('disease_type', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('disease_location', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('disease_duration', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('organism_count', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('tumor_location', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('tumor_mass', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('tumor_stage', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('drug_usage', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('durg_type', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('drug_duration', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('drug_frequency', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('perturbation', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('pert_type', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('pert_duration', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('pert_frequency', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('fetal_health_stat', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('amniotic_fluid_color', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('gestation_stat', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('maternal_health_stat', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
        ))
        db.send_create_signal(u'database', ['Human_Associated'])

        # Adding model 'Soil'
        db.create_table(u'database_soil', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sampleid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Sample'])),
            ('projectid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Project'])),
            ('refid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Reference'])),
            ('samp_collection_device', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('samp_size', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('samp_depth', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('samp_prep', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('sieve_size', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('storage_cond', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('samp_weight_dna_ext', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('pool_dna_extracts', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('fao_class', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('local_class', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('texture_class', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('porosity', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('profile_position', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('slope_aspect', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('slope_gradient', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('bulk_density', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('drainage_class', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('water_content_soil', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('cur_land_use', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('cur_vegetation', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('cur_crop', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('cur_cultivar', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('crop_rotation', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('cover_crop', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('fert_amendment_class', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('fert_placement', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('fert_type', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('fert_tot_amount', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('fert_N_tot_amount', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('fert_P_tot_amount', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('fert_K_tot_amount', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('irrigation_type', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('irrigation_tot_amount', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('residue_removal', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('residue_growth_stage', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('residue_removal_percent', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('tillage_event', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('tillage_event_depth', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('amend1_class', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('amend1_active_ingredient', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('amend1_tot_amount', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('amend2_class', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('amend2_active_ingredient', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('amend2_tot_amount', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('amend3_class', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('amend3_active_ingredient', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('amend3_tot_amount', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('rRNA_copies', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('microbial_biomass_C', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('microbial_biomass_N', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('microbial_respiration', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_pH', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_EC', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_C', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_OM', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_N', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_NO3_N', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_NH4_N', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_P', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_K', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_S', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_Zn', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_Fe', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_Cu', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_Mn', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_Ca', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_Mg', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_Na', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('soil_B', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('plant_C', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('plant_N', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('plant_P', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('plant_K', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('plant_Ca', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('plant_Mg', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('plant_S', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('plant_Na', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('plant_Cl', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('plant_Al', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('plant_B', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('plant_Cu', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('plant_Fe', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('plant_Mn', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('plant_Zn', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('crop_tot_biomass_fw', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('crop_tot_biomass_dw', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('crop_tot_above_biomass_fw', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('crop_tot_above_biomass_dw', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('crop_tot_below_biomass_fw', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('crop_tot_below_biomass_dw', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('harv_fraction', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('harv_fresh_weight', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('harv_dry_weight', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('ghg_chamber_placement', self.gf('django.db.models.fields.CharField')(max_length=45, blank=True)),
            ('ghg_N2O', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('ghg_CO2', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('ghg_NH4', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'database', ['Soil'])

        # Adding model 'UserDefined'
        db.create_table(u'database_userdefined', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sampleid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Sample'])),
            ('projectid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Project'])),
            ('refid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Reference'])),
            ('usr_cat1', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('usr_cat2', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('usr_cat3', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('usr_cat4', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('usr_cat5', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('usr_cat6', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('usr_quant1', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('usr_quant2', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('usr_quant3', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('usr_quant4', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('usr_quant5', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('usr_quant6', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'database', ['UserDefined'])

        # Adding model 'Kingdom'
        db.create_table(u'database_kingdom', (
            ('kingdomid', self.gf('django.db.models.fields.CharField')(max_length=36, primary_key=True)),
            ('kingdomName', self.gf('django.db.models.fields.CharField')(max_length=90, blank=True)),
        ))
        db.send_create_signal(u'database', ['Kingdom'])

        # Adding model 'Phyla'
        db.create_table(u'database_phyla', (
            ('kingdomid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Kingdom'])),
            ('phylaid', self.gf('django.db.models.fields.CharField')(max_length=36, primary_key=True)),
            ('phylaName', self.gf('django.db.models.fields.CharField')(max_length=90, blank=True)),
        ))
        db.send_create_signal(u'database', ['Phyla'])

        # Adding model 'Class'
        db.create_table(u'database_class', (
            ('kingdomid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Kingdom'])),
            ('phylaid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Phyla'])),
            ('classid', self.gf('django.db.models.fields.CharField')(max_length=36, primary_key=True)),
            ('className', self.gf('django.db.models.fields.CharField')(max_length=90, blank=True)),
        ))
        db.send_create_signal(u'database', ['Class'])

        # Adding model 'Order'
        db.create_table(u'database_order', (
            ('kingdomid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Kingdom'])),
            ('phylaid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Phyla'])),
            ('classid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Class'])),
            ('orderid', self.gf('django.db.models.fields.CharField')(max_length=36, primary_key=True)),
            ('orderName', self.gf('django.db.models.fields.CharField')(max_length=90, blank=True)),
        ))
        db.send_create_signal(u'database', ['Order'])

        # Adding model 'Family'
        db.create_table(u'database_family', (
            ('kingdomid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Kingdom'])),
            ('phylaid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Phyla'])),
            ('classid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Class'])),
            ('orderid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Order'])),
            ('familyid', self.gf('django.db.models.fields.CharField')(max_length=36, primary_key=True)),
            ('familyName', self.gf('django.db.models.fields.CharField')(max_length=90, blank=True)),
        ))
        db.send_create_signal(u'database', ['Family'])

        # Adding model 'Genus'
        db.create_table(u'database_genus', (
            ('kingdomid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Kingdom'])),
            ('phylaid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Phyla'])),
            ('classid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Class'])),
            ('orderid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Order'])),
            ('familyid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Family'])),
            ('genusid', self.gf('django.db.models.fields.CharField')(max_length=36, primary_key=True)),
            ('genusName', self.gf('django.db.models.fields.CharField')(max_length=90, blank=True)),
        ))
        db.send_create_signal(u'database', ['Genus'])

        # Adding model 'Species'
        db.create_table(u'database_species', (
            ('kingdomid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Kingdom'])),
            ('phylaid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Phyla'])),
            ('classid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Class'])),
            ('orderid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Order'])),
            ('familyid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Family'])),
            ('genusid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Genus'])),
            ('speciesid', self.gf('django.db.models.fields.CharField')(max_length=36, primary_key=True)),
            ('speciesName', self.gf('django.db.models.fields.CharField')(max_length=90, blank=True)),
        ))
        db.send_create_signal(u'database', ['Species'])

        # Adding model 'Profile'
        db.create_table(u'database_profile', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('projectid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Project'])),
            ('sampleid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Sample'])),
            ('refid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Reference'])),
            ('kingdomid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Kingdom'])),
            ('phylaid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Phyla'])),
            ('classid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Class'])),
            ('orderid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Order'])),
            ('familyid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Family'])),
            ('genusid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Genus'])),
            ('speciesid', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['database.Species'])),
            ('count', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'database', ['Profile'])


    def backwards(self, orm):
        # Deleting model 'Project'
        db.delete_table(u'database_project')

        # Deleting model 'Reference'
        db.delete_table(u'database_reference')

        # Deleting model 'Sample'
        db.delete_table(u'database_sample')

        # Deleting model 'Human_Associated'
        db.delete_table(u'database_human_associated')

        # Deleting model 'Soil'
        db.delete_table(u'database_soil')

        # Deleting model 'UserDefined'
        db.delete_table(u'database_userdefined')

        # Deleting model 'Kingdom'
        db.delete_table(u'database_kingdom')

        # Deleting model 'Phyla'
        db.delete_table(u'database_phyla')

        # Deleting model 'Class'
        db.delete_table(u'database_class')

        # Deleting model 'Order'
        db.delete_table(u'database_order')

        # Deleting model 'Family'
        db.delete_table(u'database_family')

        # Deleting model 'Genus'
        db.delete_table(u'database_genus')

        # Deleting model 'Species'
        db.delete_table(u'database_species')

        # Deleting model 'Profile'
        db.delete_table(u'database_profile')


    models = {
        u'auth.group': {
            'Meta': {'object_name': 'Group'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '80'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'})
        },
        u'auth.permission': {
            'Meta': {'ordering': "(u'content_type__app_label', u'content_type__model', u'codename')", 'unique_together': "((u'content_type', u'codename'),)", 'object_name': 'Permission'},
            'codename': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'auth.user': {
            'Meta': {'object_name': 'User'},
            'date_joined': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'first_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Group']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Permission']"}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        u'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'database.class': {
            'Meta': {'object_name': 'Class'},
            'className': ('django.db.models.fields.CharField', [], {'max_length': '90', 'blank': 'True'}),
            'classid': ('django.db.models.fields.CharField', [], {'max_length': '36', 'primary_key': 'True'}),
            'kingdomid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Kingdom']"}),
            'phylaid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Phyla']"})
        },
        u'database.family': {
            'Meta': {'object_name': 'Family'},
            'classid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Class']"}),
            'familyName': ('django.db.models.fields.CharField', [], {'max_length': '90', 'blank': 'True'}),
            'familyid': ('django.db.models.fields.CharField', [], {'max_length': '36', 'primary_key': 'True'}),
            'kingdomid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Kingdom']"}),
            'orderid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Order']"}),
            'phylaid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Phyla']"})
        },
        u'database.genus': {
            'Meta': {'object_name': 'Genus'},
            'classid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Class']"}),
            'familyid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Family']"}),
            'genusName': ('django.db.models.fields.CharField', [], {'max_length': '90', 'blank': 'True'}),
            'genusid': ('django.db.models.fields.CharField', [], {'max_length': '36', 'primary_key': 'True'}),
            'kingdomid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Kingdom']"}),
            'orderid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Order']"}),
            'phylaid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Phyla']"})
        },
        u'database.human_associated': {
            'Meta': {'object_name': 'Human_Associated'},
            'amniotic_fluid_color': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'diet_duration': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'diet_frequency': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'diet_last_six_month': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'diet_type': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'disease_duration': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'disease_location': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'disease_type': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'drug_duration': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'drug_frequency': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'drug_usage': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'durg_type': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'fetal_health_stat': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'gestation_stat': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'host_age': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'host_bmi': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'host_body_temp': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'host_ethnicity': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'host_gender': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'host_height': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'host_occupation': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'host_pulse': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'host_subject_id': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'host_weight': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'host_weight_loss_3_month': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'last_meal': ('django.db.models.fields.CharField', [], {'max_length': '15', 'blank': 'True'}),
            'maternal_health_stat': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'medic_hist_perform': ('django.db.models.fields.CharField', [], {'max_length': '15', 'blank': 'True'}),
            'organism_count': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'pert_duration': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'pert_frequency': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'pert_type': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'perturbation': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'pet_farm_animal': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'projectid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Project']"}),
            'refid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Reference']"}),
            'samp_collect_device': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'samp_location': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'samp_mat_process': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'samp_oxy_stat': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'samp_ph': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'samp_salinity': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'samp_size': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'samp_store_dur': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'samp_store_temp': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'samp_temp': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'samp_type': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'sampleid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Sample']"}),
            'smoker': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'tumor_location': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'tumor_mass': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'tumor_stage': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'})
        },
        u'database.kingdom': {
            'Meta': {'object_name': 'Kingdom'},
            'kingdomName': ('django.db.models.fields.CharField', [], {'max_length': '90', 'blank': 'True'}),
            'kingdomid': ('django.db.models.fields.CharField', [], {'max_length': '36', 'primary_key': 'True'})
        },
        u'database.order': {
            'Meta': {'object_name': 'Order'},
            'classid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Class']"}),
            'kingdomid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Kingdom']"}),
            'orderName': ('django.db.models.fields.CharField', [], {'max_length': '90', 'blank': 'True'}),
            'orderid': ('django.db.models.fields.CharField', [], {'max_length': '36', 'primary_key': 'True'}),
            'phylaid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Phyla']"})
        },
        u'database.phyla': {
            'Meta': {'object_name': 'Phyla'},
            'kingdomid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Kingdom']"}),
            'phylaName': ('django.db.models.fields.CharField', [], {'max_length': '90', 'blank': 'True'}),
            'phylaid': ('django.db.models.fields.CharField', [], {'max_length': '36', 'primary_key': 'True'})
        },
        u'database.profile': {
            'Meta': {'object_name': 'Profile'},
            'classid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Class']"}),
            'count': ('django.db.models.fields.IntegerField', [], {}),
            'familyid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Family']"}),
            'genusid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Genus']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'kingdomid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Kingdom']"}),
            'orderid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Order']"}),
            'phylaid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Phyla']"}),
            'projectid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Project']"}),
            'refid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Reference']"}),
            'sampleid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Sample']"}),
            'speciesid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Species']"})
        },
        u'database.project': {
            'Meta': {'object_name': 'Project'},
            'end_date': ('django.db.models.fields.CharField', [], {'max_length': '15', 'blank': 'True'}),
            'pi_affiliation': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'pi_email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'pi_first': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'pi_last': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'pi_phone': ('django.db.models.fields.CharField', [], {'max_length': '15', 'blank': 'True'}),
            'projectType': ('django.db.models.fields.CharField', [], {'max_length': '45'}),
            'project_desc': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'project_name': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'projectid': ('django.db.models.fields.CharField', [], {'max_length': '36', 'primary_key': 'True'}),
            'start_date': ('django.db.models.fields.CharField', [], {'max_length': '15', 'blank': 'True'}),
            'status': ('django.db.models.fields.CharField', [], {'max_length': '10'})
        },
        u'database.reference': {
            'Meta': {'object_name': 'Reference'},
            'alignDB': ('django.db.models.fields.CharField', [], {'max_length': '90', 'blank': 'True'}),
            'author': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'entries'", 'to': u"orm['auth.User']"}),
            'path': ('django.db.models.fields.CharField', [], {'max_length': '90'}),
            'projectid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Project']"}),
            'raw': ('django.db.models.fields.BooleanField', [], {}),
            'refid': ('django.db.models.fields.CharField', [], {'max_length': '36', 'primary_key': 'True'}),
            'source': ('django.db.models.fields.CharField', [], {'max_length': '90'}),
            'taxonomyDB': ('django.db.models.fields.CharField', [], {'max_length': '90', 'blank': 'True'}),
            'templateDB': ('django.db.models.fields.CharField', [], {'max_length': '90', 'blank': 'True'})
        },
        u'database.sample': {
            'Meta': {'object_name': 'Sample'},
            'annual_season_precpt': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'annual_season_temp': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'collection_date': ('django.db.models.fields.CharField', [], {'max_length': '15', 'blank': 'True'}),
            'depth': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'elev': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'env_biome': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'env_feature': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'env_material': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'geo_loc_city': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'geo_loc_country': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'geo_loc_farm': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'geo_loc_plot': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'geo_loc_state': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'latitude': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'longitude': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'organism': ('django.db.models.fields.CharField', [], {'max_length': '90', 'blank': 'True'}),
            'projectid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Project']"}),
            'refid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Reference']"}),
            'sample_name': ('django.db.models.fields.TextField', [], {}),
            'sampleid': ('django.db.models.fields.CharField', [], {'max_length': '36', 'primary_key': 'True'}),
            'seq_for_primer': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'seq_gene': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'seq_gene_region': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'seq_platform': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'seq_rev_primer': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'})
        },
        u'database.soil': {
            'Meta': {'object_name': 'Soil'},
            'amend1_active_ingredient': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'amend1_class': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'amend1_tot_amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'amend2_active_ingredient': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'amend2_class': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'amend2_tot_amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'amend3_active_ingredient': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'amend3_class': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'amend3_tot_amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'bulk_density': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'cover_crop': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'crop_rotation': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'crop_tot_above_biomass_dw': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'crop_tot_above_biomass_fw': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'crop_tot_below_biomass_dw': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'crop_tot_below_biomass_fw': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'crop_tot_biomass_dw': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'crop_tot_biomass_fw': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'cur_crop': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'cur_cultivar': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'cur_land_use': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'cur_vegetation': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'drainage_class': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'fao_class': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'fert_K_tot_amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'fert_N_tot_amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'fert_P_tot_amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'fert_amendment_class': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'fert_placement': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'fert_tot_amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'fert_type': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'ghg_CO2': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'ghg_N2O': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'ghg_NH4': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'ghg_chamber_placement': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'harv_dry_weight': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'harv_fraction': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'harv_fresh_weight': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'irrigation_tot_amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'irrigation_type': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'local_class': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'microbial_biomass_C': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'microbial_biomass_N': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'microbial_respiration': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plant_Al': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plant_B': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plant_C': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plant_Ca': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plant_Cl': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plant_Cu': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plant_Fe': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plant_K': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plant_Mg': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plant_Mn': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plant_N': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plant_Na': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plant_P': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plant_S': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plant_Zn': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'pool_dna_extracts': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'porosity': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'profile_position': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'projectid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Project']"}),
            'rRNA_copies': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'refid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Reference']"}),
            'residue_growth_stage': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'residue_removal': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'residue_removal_percent': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'samp_collection_device': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'samp_depth': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'samp_prep': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'samp_size': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'samp_weight_dna_ext': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'sampleid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Sample']"}),
            'sieve_size': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'slope_aspect': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'slope_gradient': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_B': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_C': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_Ca': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_Cu': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_EC': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_Fe': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_K': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_Mg': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_Mn': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_N': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_NH4_N': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_NO3_N': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_Na': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_OM': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_P': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_S': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_Zn': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'soil_pH': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'storage_cond': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'texture_class': ('django.db.models.fields.CharField', [], {'max_length': '45', 'blank': 'True'}),
            'tillage_event': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'tillage_event_depth': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'water_content_soil': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'})
        },
        u'database.species': {
            'Meta': {'object_name': 'Species'},
            'classid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Class']"}),
            'familyid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Family']"}),
            'genusid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Genus']"}),
            'kingdomid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Kingdom']"}),
            'orderid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Order']"}),
            'phylaid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Phyla']"}),
            'speciesName': ('django.db.models.fields.CharField', [], {'max_length': '90', 'blank': 'True'}),
            'speciesid': ('django.db.models.fields.CharField', [], {'max_length': '36', 'primary_key': 'True'})
        },
        u'database.userdefined': {
            'Meta': {'object_name': 'UserDefined'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'projectid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Project']"}),
            'refid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Reference']"}),
            'sampleid': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['database.Sample']"}),
            'usr_cat1': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'usr_cat2': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'usr_cat3': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'usr_cat4': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'usr_cat5': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'usr_cat6': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'usr_quant1': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'usr_quant2': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'usr_quant3': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'usr_quant4': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'usr_quant5': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'usr_quant6': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'})
        }
    }

    complete_apps = ['database']