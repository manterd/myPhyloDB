from django.db import models
from django_extensions.db.fields import UUIDField
from django.contrib.auth.models import User as Users
import time


class Project(models.Model):
    status = models.CharField(max_length=10, blank=False)
    projectType = models.CharField(max_length=45, blank=False)
    projectid = UUIDField(primary_key=True)
    project_name = models.TextField(blank=True)
    project_desc = models.TextField(blank=True)
    start_date = models.CharField(max_length=15, blank=True)
    end_date = models.CharField(max_length=15, blank=True)
    pi_last = models.TextField(blank=True)
    pi_first = models.TextField(blank=True)
    pi_affiliation = models.CharField(max_length=45, blank=True)
    pi_email = models.EmailField(blank=True)
    pi_phone = models.CharField(max_length=15, blank=True)

    def __unicode__(self):
        return unicode(self.project_name)


class Reference(models.Model):
    refid = UUIDField(primary_key=True)
    raw = models.BooleanField()
    projectid = models.ForeignKey(Project)
    path = models.CharField(max_length=90)
    source = models.CharField(max_length=90)
    alignDB = models.CharField(max_length=90, blank=True)
    templateDB = models.CharField(max_length=90, blank=True)
    taxonomyDB = models.CharField(max_length=90, blank=True)
    author = models.ForeignKey(Users, related_name='entries')

    class Meta:
        verbose_name_plural = 'entries'

    def __unicode__(self):
        return unicode(self.path)


class Sample(models.Model):
    projectid = models.ForeignKey(Project)
    refid = models.ForeignKey(Reference)
    sampleid = UUIDField(primary_key=True)
    sample_name = models.TextField(blank=False)
    organism = models.CharField(max_length=90, blank=True)
    collection_date = models.CharField(max_length=15, blank=True)
    depth = models.CharField(max_length=45, blank=True)
    elev = models.CharField(max_length=45, blank=True)
    seq_platform = models.CharField(max_length=45, blank=True)
    seq_gene = models.CharField(max_length=45, blank=True)
    seq_gene_region = models.CharField(max_length=45, blank=True)
    seq_for_primer = models.CharField(max_length=45, blank=True)
    seq_rev_primer = models.CharField(max_length=45, blank=True)
    env_biome = models.CharField(max_length=45, blank=True)
    env_feature = models.CharField(max_length=45, blank=True)
    env_material = models.CharField(max_length=45, blank=True)
    geo_loc_country = models.CharField(max_length=45, blank=True)
    geo_loc_state = models.CharField(max_length=45, blank=True)
    geo_loc_city = models.CharField(max_length=45, blank=True)
    geo_loc_farm = models.CharField(max_length=45, blank=True)
    geo_loc_plot = models.CharField(max_length=45, blank=True)
    latitude = models.FloatField(blank=True, null=True)
    longitude = models.FloatField(blank=True, null=True)
    annual_season_precpt = models.FloatField(blank=True, null=True)
    annual_season_temp = models.FloatField(blank=True, null=True)

    def natural_key(self):
        return self.sampleid

    def __unicode__(self):
        return unicode(self.sample_name)


class Human_Associated(models.Model):
    sampleid = models.ForeignKey(Sample)
    projectid = models.ForeignKey(Project)
    refid = models.ForeignKey(Reference)

    # sample collection
    samp_collect_device = models.CharField(max_length=45, blank=True)
    samp_mat_process = models.CharField(max_length=45, blank=True)
    samp_size = models.FloatField(blank=True, null=True)
    samp_store_temp = models.FloatField(blank=True, null=True)
    samp_store_dur = models.FloatField(blank=True, null=True)

    # sample classification
    samp_type = models.CharField(max_length=45, blank=True)
    samp_location = models.CharField(max_length=45, blank=True)
    samp_temp = models.FloatField(blank=True, null=True)
    samp_ph = models.FloatField(blank=True, null=True)
    samp_oxy_stat = models.CharField(max_length=45, blank=True)
    samp_salinity = models.FloatField(blank=True, null=True)

    host_subject_id = models.CharField(max_length=45, blank=True)
    host_age = models.FloatField(blank=True, null=True)
    host_pulse = models.FloatField(blank=True, null=True)
    host_gender = models.CharField(max_length=45, blank=True)
    host_ethnicity = models.CharField(max_length=45, blank=True)
    host_height = models.FloatField(blank=True, null=True)
    host_weight = models.FloatField(blank=True, null=True)
    host_bmi = models.FloatField(blank=True, null=True)
    host_weight_loss_3_month = models.FloatField(blank=True, null=True)
    host_body_temp = models.FloatField(blank=True, null=True)
    host_occupation = models.CharField(max_length=45, blank=True)
    pet_farm_animal = models.CharField(max_length=45, blank=True)
    smoker = models.CharField(max_length=45, blank=True)

    diet_type = models.CharField(max_length=45, blank=True)
    diet_duration = models.FloatField(blank=True, null=True)
    diet_frequency = models.CharField(max_length=45, blank=True)
    diet_last_six_month = models.CharField(max_length=45, blank=True)
    last_meal = models.CharField(max_length=15, blank=True)

    medic_hist_perform = models.CharField(max_length=15, blank=True)
    disease_type = models.CharField(max_length=45, blank=True)
    disease_location = models.CharField(max_length=45, blank=True)
    disease_duration = models.FloatField(blank=True, null=True)
    organism_count = models.FloatField(blank=True, null=True)
    tumor_location = models.CharField(max_length=45, blank=True)
    tumor_mass = models.FloatField(blank=True, null=True)
    tumor_stage = models.CharField(max_length=45, blank=True)

    drug_usage = models.CharField(max_length=45, blank=True)
    durg_type = models.CharField(max_length=45, blank=True)
    drug_duration = models.FloatField(blank=True, null=True)
    drug_frequency = models.CharField(max_length=45, blank=True)

    perturbation = models.CharField(max_length=45, blank=True)
    pert_type = models.CharField(max_length=45, blank=True)
    pert_duration = models.FloatField(blank=True, null=True)
    pert_frequency = models.CharField(max_length=45, blank=True)

    fetal_health_stat = models.CharField(max_length=45, blank=True)
    amniotic_fluid_color = models.CharField(max_length=45, blank=True)
    gestation_stat = models.CharField(max_length=45, blank=True)
    maternal_health_stat = models.CharField(max_length=45, blank=True)


class Soil(models.Model):
    sampleid = models.ForeignKey(Sample)
    projectid = models.ForeignKey(Project)
    refid = models.ForeignKey(Reference)

    samp_collection_device = models.CharField(max_length=45, blank=True)
    samp_size = models.FloatField(blank=True, null=True)
    samp_depth = models.CharField(max_length=45, blank=True)
    samp_prep = models.CharField(max_length=45, blank=True)
    sieve_size = models.FloatField(blank=True, null=True)
    storage_cond = models.FloatField(blank=True, null=True)
    samp_weight_dna_ext = models.FloatField(blank=True, null=True)
    pool_dna_extracts = models.IntegerField(blank=True, null=True)

    fao_class = models.CharField(max_length=45, blank=True)
    local_class = models.CharField(max_length=45, blank=True)
    texture_class = models.CharField(max_length=45, blank=True)
    porosity = models.FloatField(blank=True, null=True)
    profile_position = models.CharField(max_length=45, blank=True)
    slope_aspect = models.CharField(max_length=45, blank=True)
    slope_gradient = models.FloatField(blank=True, null=True)
    bulk_density = models.FloatField(blank=True, null=True)
    drainage_class = models.CharField(max_length=45, blank=True)
    water_content_soil = models.FloatField(blank=True, null=True)

    cur_land_use = models.TextField(blank=True)
    cur_vegetation = models.TextField(blank=True)
    cur_crop = models.TextField(blank=True)
    cur_cultivar = models.TextField(blank=True)
    crop_rotation = models.TextField(blank=True)
    cover_crop = models.TextField(blank=True)

    fert_amendment_class = models.TextField(blank=True)
    fert_placement = models.TextField(blank=True)
    fert_type = models.TextField(blank=True)
    fert_tot_amount = models.FloatField(blank=True, null=True)
    fert_N_tot_amount = models.FloatField(blank=True, null=True)
    fert_P_tot_amount = models.FloatField(blank=True, null=True)
    fert_K_tot_amount = models.FloatField(blank=True, null=True)

    irrigation_type = models.TextField(blank=True)
    irrigation_tot_amount = models.FloatField(blank=True, null=True)

    residue_removal = models.TextField(blank=True)
    residue_growth_stage = models.TextField(blank=True)
    residue_removal_percent = models.FloatField(blank=True, null=True)

    tillage_event = models.TextField(blank=True)
    tillage_event_depth = models.FloatField(blank=True, null=True)

    amend1_class = models.TextField(blank=True)
    amend1_active_ingredient = models.TextField(blank=True)
    amend1_tot_amount = models.FloatField(blank=True, null=True)
    amend2_class = models.TextField(blank=True)
    amend2_active_ingredient = models.TextField(blank=True)
    amend2_tot_amount = models.FloatField(blank=True, null=True)
    amend3_class = models.TextField(blank=True)
    amend3_active_ingredient = models.TextField(blank=True)
    amend3_tot_amount = models.FloatField(blank=True, null=True)

    rRNA_copies = models.FloatField(blank=True, null=True)
    microbial_biomass_C = models.FloatField(blank=True, null=True)
    microbial_biomass_N = models.FloatField(blank=True, null=True)
    microbial_respiration = models.FloatField(blank=True, null=True)

    soil_pH = models.FloatField(blank=True, null=True)
    soil_EC = models.FloatField(blank=True, null=True)
    soil_C = models.FloatField(blank=True, null=True)
    soil_OM = models.FloatField(blank=True, null=True)
    soil_N = models.FloatField(blank=True, null=True)
    soil_NO3_N = models.FloatField(blank=True, null=True)
    soil_NH4_N = models.FloatField(blank=True, null=True)
    soil_P = models.FloatField(blank=True, null=True)
    soil_K = models.FloatField(blank=True, null=True)
    soil_S = models.FloatField(blank=True, null=True)
    soil_Zn = models.FloatField(blank=True, null=True)
    soil_Fe = models.FloatField(blank=True, null=True)
    soil_Cu = models.FloatField(blank=True, null=True)
    soil_Mn = models.FloatField(blank=True, null=True)
    soil_Ca = models.FloatField(blank=True, null=True)
    soil_Mg = models.FloatField(blank=True, null=True)
    soil_Na = models.FloatField(blank=True, null=True)
    soil_B = models.FloatField(blank=True, null=True)

    plant_C = models.FloatField(blank=True, null=True)
    plant_N = models.FloatField(blank=True, null=True)
    plant_P = models.FloatField(blank=True, null=True)
    plant_K = models.FloatField(blank=True, null=True)
    plant_Ca = models.FloatField(blank=True, null=True)
    plant_Mg = models.FloatField(blank=True, null=True)
    plant_S = models.FloatField(blank=True, null=True)
    plant_Na = models.FloatField(blank=True, null=True)
    plant_Cl = models.FloatField(blank=True, null=True)
    plant_Al = models.FloatField(blank=True, null=True)
    plant_B = models.FloatField(blank=True, null=True)
    plant_Cu = models.FloatField(blank=True, null=True)
    plant_Fe = models.FloatField(blank=True, null=True)
    plant_Mn = models.FloatField(blank=True, null=True)
    plant_Zn = models.FloatField(blank=True, null=True)

    crop_tot_biomass_fw = models.FloatField(blank=True, null=True)
    crop_tot_biomass_dw = models.FloatField(blank=True, null=True)
    crop_tot_above_biomass_fw = models.FloatField(blank=True, null=True)
    crop_tot_above_biomass_dw = models.FloatField(blank=True, null=True)
    crop_tot_below_biomass_fw = models.FloatField(blank=True, null=True)
    crop_tot_below_biomass_dw = models.FloatField(blank=True, null=True)
    harv_fraction = models.CharField(max_length=45, blank=True)
    harv_fresh_weight = models.FloatField(blank=True, null=True)
    harv_dry_weight = models.FloatField(blank=True, null=True)

    ghg_chamber_placement = models.CharField(max_length=45, blank=True)
    ghg_N2O = models.FloatField(blank=True, null=True)
    ghg_CO2 = models.FloatField(blank=True, null=True)
    ghg_NH4 = models.FloatField(blank=True, null=True)


class UserDefined(models.Model):
    sampleid = models.ForeignKey(Sample)
    projectid = models.ForeignKey(Project)
    refid = models.ForeignKey(Reference)
    usr_cat1 = models.TextField(blank=True)
    usr_cat2 = models.TextField(blank=True)
    usr_cat3 = models.TextField(blank=True)
    usr_cat4 = models.TextField(blank=True)
    usr_cat5 = models.TextField(blank=True)
    usr_cat6 = models.TextField(blank=True)
    usr_quant1 = models.FloatField(blank=True, null=True)
    usr_quant2 = models.FloatField(blank=True, null=True)
    usr_quant3 = models.FloatField(blank=True, null=True)
    usr_quant4 = models.FloatField(blank=True, null=True)
    usr_quant5 = models.FloatField(blank=True, null=True)
    usr_quant6 = models.FloatField(blank=True, null=True)

    def natural_key(self):
        return self.sampleid


class Kingdom(models.Model):
    kingdomid = UUIDField(primary_key=True)
    kingdomName = models.CharField(max_length=90, blank=True)


class Phyla(models.Model):
    kingdomid = models.ForeignKey(Kingdom)
    phylaid = UUIDField(primary_key=True)
    phylaName = models.CharField(max_length=90, blank=True)


class Class(models.Model):
    kingdomid = models.ForeignKey(Kingdom)
    phylaid = models.ForeignKey(Phyla)
    classid = UUIDField(primary_key=True)
    className = models.CharField(max_length=90, blank=True)


class Order(models.Model):
    kingdomid = models.ForeignKey(Kingdom)
    phylaid = models.ForeignKey(Phyla)
    classid = models.ForeignKey(Class)
    orderid = UUIDField(primary_key=True)
    orderName = models.CharField(max_length=90, blank=True)


class Family(models.Model):
    kingdomid = models.ForeignKey(Kingdom)
    phylaid = models.ForeignKey(Phyla)
    classid = models.ForeignKey(Class)
    orderid = models.ForeignKey(Order)
    familyid = UUIDField(primary_key=True)
    familyName = models.CharField(max_length=90, blank=True)


class Genus(models.Model):
    kingdomid = models.ForeignKey(Kingdom)
    phylaid = models.ForeignKey(Phyla)
    classid = models.ForeignKey(Class)
    orderid = models.ForeignKey(Order)
    familyid = models.ForeignKey(Family)
    genusid = UUIDField(primary_key=True)
    genusName = models.CharField(max_length=90, blank=True)


class Species(models.Model):
    kingdomid = models.ForeignKey(Kingdom)
    phylaid = models.ForeignKey(Phyla)
    classid = models.ForeignKey(Class)
    orderid = models.ForeignKey(Order)
    familyid = models.ForeignKey(Family)
    genusid = models.ForeignKey(Genus)
    speciesid = UUIDField(primary_key=True)
    speciesName = models.CharField(max_length=90, blank=True)


class Profile(models.Model):
    projectid = models.ForeignKey(Project)
    sampleid = models.ForeignKey(Sample)
    refid = models.ForeignKey(Reference)
    kingdomid = models.ForeignKey(Kingdom)
    phylaid = models.ForeignKey(Phyla)
    classid = models.ForeignKey(Class)
    orderid = models.ForeignKey(Order)
    familyid = models.ForeignKey(Family)
    genusid = models.ForeignKey(Genus)
    speciesid = models.ForeignKey(Species)
    count = models.IntegerField()


