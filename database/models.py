from django.db import models
from django.contrib.auth.models import User as Users
from jsonfield import JSONField


queue = 0


def addQueue():
    global queue
    queue += 1


def getQueue():
    return queue


def subQueue():
    global queue
    queue -= 1


class Project(models.Model):
    projectid = models.TextField(primary_key=True)
    status = models.CharField(max_length=15, blank=False)
    projectType = models.CharField(max_length=90, blank=False)
    project_name = models.TextField(blank=True)
    project_desc = models.TextField(blank=True)
    start_date = models.CharField(max_length=15, blank=True)
    end_date = models.CharField(max_length=15, blank=True)
    pi_last = models.TextField(blank=True)
    pi_first = models.TextField(blank=True)
    pi_affiliation = models.CharField(max_length=90, blank=True)
    pi_email = models.EmailField(blank=True)
    pi_phone = models.CharField(max_length=15, blank=True)


class Reference(models.Model):
    refid = models.TextField(primary_key=True)
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


class Sample(models.Model):
    projectid = models.ForeignKey(Project)
    refid = models.ForeignKey(Reference)
    sampleid = models.TextField(primary_key=True)

    sample_name = models.TextField(blank=False)
    organism = models.CharField(max_length=90, blank=True)
    collection_date = models.CharField(max_length=15, blank=True)
    depth = models.CharField(max_length=90, blank=True)
    elev = models.CharField(max_length=90, blank=True)
    seq_platform = models.CharField(max_length=90, blank=True)
    seq_gene = models.CharField(max_length=90, blank=True)
    seq_gene_region = models.CharField(max_length=90, blank=True)
    seq_barcode = models.CharField(max_length=90, blank=True)
    seq_for_primer = models.CharField(max_length=90, blank=True)
    seq_rev_primer = models.CharField(max_length=90, blank=True)
    env_biome = models.CharField(max_length=90, blank=True)
    env_feature = models.CharField(max_length=90, blank=True)
    env_material = models.CharField(max_length=90, blank=True)
    geo_loc_country = models.CharField(max_length=90, blank=True)
    geo_loc_state = models.CharField(max_length=90, blank=True)
    geo_loc_city = models.CharField(max_length=90, blank=True)
    geo_loc_farm = models.CharField(max_length=90, blank=True)
    geo_loc_plot = models.CharField(max_length=90, blank=True)
    latitude = models.FloatField(blank=True, null=True)
    longitude = models.FloatField(blank=True, null=True)
    annual_season_precpt = models.FloatField(blank=True, null=True)
    annual_season_temp = models.FloatField(blank=True, null=True)

    def natural_key(self):
        return self.sampleid


class Human_Associated(models.Model):
    sampleid = models.ForeignKey(Sample)
    projectid = models.ForeignKey(Project)
    refid = models.ForeignKey(Reference)

    # sample collection
    samp_collect_device = models.CharField(max_length=90, blank=True)
    samp_mat_process = models.CharField(max_length=90, blank=True)
    samp_size = models.FloatField(blank=True, null=True)
    samp_store_dur = models.FloatField(blank=True, null=True)
    samp_store_loc = models.CharField(max_length=90, blank=True)
    samp_store_temp = models.FloatField(blank=True, null=True)

    # sample classification
    samp_type = models.CharField(max_length=90, blank=True)
    samp_location = models.CharField(max_length=90, blank=True)
    samp_temp = models.FloatField(blank=True, null=True)
    samp_ph = models.FloatField(blank=True, null=True)
    samp_oxy_stat = models.CharField(max_length=90, blank=True)
    samp_salinity = models.FloatField(blank=True, null=True)

    host_subject_id = models.CharField(max_length=90, blank=True)
    host_age = models.FloatField(blank=True, null=True)
    host_pulse = models.FloatField(blank=True, null=True)
    host_gender = models.CharField(max_length=90, blank=True)
    host_ethnicity = models.CharField(max_length=90, blank=True)
    host_height = models.FloatField(blank=True, null=True)
    host_weight = models.FloatField(blank=True, null=True)
    host_bmi = models.FloatField(blank=True, null=True)
    host_weight_loss_3_month = models.FloatField(blank=True, null=True)
    host_body_temp = models.FloatField(blank=True, null=True)
    host_occupation = models.CharField(max_length=90, blank=True)
    pet_farm_animal = models.CharField(max_length=90, blank=True)
    obesity = models.CharField(max_length=90, blank=True)
    smoker = models.CharField(max_length=90, blank=True)

    diet_type = models.CharField(max_length=90, blank=True)
    diet_duration = models.FloatField(blank=True, null=True)
    diet_frequency = models.CharField(max_length=90, blank=True)
    diet_last_six_month = models.CharField(max_length=90, blank=True)
    last_meal = models.CharField(max_length=90, blank=True)

    medic_hist_perform = models.CharField(max_length=90, blank=True)
    disease_type = models.CharField(max_length=90, blank=True)
    disease_location = models.CharField(max_length=90, blank=True)
    disease_duration = models.FloatField(blank=True, null=True)
    organism_count = models.FloatField(blank=True, null=True)
    tumor_location = models.CharField(max_length=90, blank=True)
    tumor_mass = models.FloatField(blank=True, null=True)
    tumor_stage = models.CharField(max_length=90, blank=True)

    drug_usage = models.CharField(max_length=90, blank=True)
    drug_type = models.CharField(max_length=90, blank=True)
    drug_duration = models.FloatField(blank=True, null=True)
    drug_frequency = models.CharField(max_length=90, blank=True)

    perturbation = models.CharField(max_length=90, blank=True)
    pert_type = models.CharField(max_length=90, blank=True)
    pert_duration = models.FloatField(blank=True, null=True)
    pert_frequency = models.CharField(max_length=90, blank=True)

    fetal_health_stat = models.CharField(max_length=90, blank=True)
    amniotic_fluid_color = models.CharField(max_length=90, blank=True)
    gestation_stat = models.CharField(max_length=90, blank=True)
    maternal_health_stat = models.CharField(max_length=90, blank=True)


class Soil(models.Model):
    sampleid = models.ForeignKey(Sample)
    projectid = models.ForeignKey(Project)
    refid = models.ForeignKey(Reference)

    samp_collection_device = models.CharField(max_length=90, blank=True)
    samp_size = models.FloatField(blank=True, null=True)
    samp_depth = models.CharField(max_length=90, blank=True)
    samp_prep = models.CharField(max_length=90, blank=True)
    samp_sieve_size = models.FloatField(blank=True, null=True)
    samp_store_dur = models.FloatField(blank=True, null=True)
    samp_store_loc = models.CharField(max_length=90, blank=True)
    samp_store_temp = models.FloatField(blank=True, null=True)
    samp_weight_dna_ext = models.FloatField(blank=True, null=True)
    pool_dna_extracts = models.FloatField(blank=True, null=True)

    fao_class = models.CharField(max_length=90, blank=True)
    local_class = models.CharField(max_length=90, blank=True)
    texture_class = models.CharField(max_length=90, blank=True)
    porosity = models.FloatField(blank=True, null=True)
    profile_position = models.CharField(max_length=90, blank=True)
    slope_aspect = models.CharField(max_length=90, blank=True)
    slope_gradient = models.FloatField(blank=True, null=True)
    bulk_density = models.FloatField(blank=True, null=True)
    drainage_class = models.CharField(max_length=90, blank=True)
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

    # new stuff from health assessment
    soil_water_cap = models.FloatField(blank=True, null=True)
    soil_surf_hard = models.FloatField(blank=True, null=True)
    soil_subsurf_hard = models.FloatField(blank=True, null=True)
    soil_agg_stability = models.FloatField(blank=True, null=True)
    soil_ACE_protein = models.FloatField(blank=True, null=True)
    soil_active_C = models.FloatField(blank=True, null=True)

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
    harv_fraction = models.CharField(max_length=90, blank=True)
    harv_fresh_weight = models.FloatField(blank=True, null=True)
    harv_dry_weight = models.FloatField(blank=True, null=True)

    ghg_chamber_placement = models.CharField(max_length=45, blank=True)
    ghg_N2O = models.FloatField(blank=True, null=True)
    ghg_CO2 = models.FloatField(blank=True, null=True)
    ghg_NH4 = models.FloatField(blank=True, null=True)


class Air(models.Model):
    sampleid = models.ForeignKey(Sample)
    projectid = models.ForeignKey(Project)
    refid = models.ForeignKey(Reference)

    barometric_press = models.FloatField(blank=True, null=True)
    carb_dioxide = models.FloatField(blank=True, null=True)
    carb_monoxide = models.FloatField(blank=True, null=True)
    chem_admin_term = models.CharField(max_length=90, blank=True)
    chem_admin_time = models.CharField(max_length=90, blank=True)
    elev = models.FloatField(blank=True, null=True)
    humidity = models.FloatField(blank=True, null=True)
    methane = models.FloatField(blank=True, null=True)
    organism_type = models.CharField(max_length=90, blank=True)
    organism_count = models.FloatField(blank=True, null=True)
    oxy_stat_samp = models.CharField(max_length=90, blank=True)
    oxygen = models.FloatField(blank=True, null=True)
    perturbation_type = models.CharField(max_length=90, blank=True)
    perturbation_interval = models.CharField(max_length=90, blank=True)
    pollutants_type = models.CharField(max_length=90, blank=True)
    pollutants_concentration = models.FloatField(blank=True, null=True)
    rel_to_oxygen = models.CharField(max_length=90, blank=True)
    resp_part_matter_substance = models.CharField(max_length=90, blank=True)
    resp_part_matter_concentration = models.FloatField(blank=True, null=True)
    samp_collect_device = models.CharField(max_length=90, blank=True)
    samp_mat_process = models.CharField(max_length=90, blank=True)
    samp_salinity = models.FloatField(blank=True, null=True)
    samp_size = models.FloatField(blank=True, null=True)
    samp_store_dur = models.FloatField(blank=True, null=True)
    samp_store_loc = models.CharField(max_length=90, blank=True)
    samp_store_temp = models.FloatField(blank=True, null=True)
    solar_irradiance = models.FloatField(blank=True, null=True)
    temp = models.FloatField(blank=True, null=True)
    ventilation_rate = models.FloatField(blank=True, null=True)
    ventilation_type = models.CharField(max_length=90, blank=True)
    volatile_org_comp_name = models.CharField(max_length=90, blank=True)
    volatile_org_comp_concentration = models.FloatField(blank=True, null=True)
    wind_direction = models.CharField(max_length=90, blank=True)
    wind_speed = models.FloatField(blank=True, null=True)


class Microbial(models.Model):
    sampleid = models.ForeignKey(Sample)
    projectid = models.ForeignKey(Project)
    refid = models.ForeignKey(Reference)

    alkalinity = models.FloatField(blank=True, null=True)
    alkyl_diethers = models.FloatField(blank=True, null=True)
    altitude = models.FloatField(blank=True, null=True)
    aminopept_act = models.FloatField(blank=True, null=True)
    ammonium = models.FloatField(blank=True, null=True)
    bacteria_carb_prod = models.FloatField(blank=True, null=True)
    biomass_part_name = models.CharField(max_length=90, blank=True)
    biomass_amount = models.FloatField(blank=True, null=True)
    bishomohopanol = models.FloatField(blank=True, null=True)
    bromide = models.FloatField(blank=True, null=True)
    calcium = models.FloatField(blank=True, null=True)
    carb_nitro_ratio = models.FloatField(blank=True, null=True)
    chem_administration_term = models.CharField(max_length=90, blank=True)
    chem_administration_time = models.CharField(max_length=90, blank=True)
    chloride = models.FloatField(blank=True, null=True)
    chlorophyll = models.FloatField(blank=True, null=True)
    diether_lipids_name = models.CharField(max_length=90, blank=True)
    diether_lipids_concentration = models.FloatField(blank=True, null=True)
    diss_carb_dioxide = models.FloatField(blank=True, null=True)
    diss_hydrogen = models.FloatField(blank=True, null=True)
    diss_inorg_carb = models.FloatField(blank=True, null=True)
    diss_org_carb = models.FloatField(blank=True, null=True)
    diss_org_nitro = models.FloatField(blank=True, null=True)
    diss_oxygen = models.FloatField(blank=True, null=True)
    glucosidase_act = models.FloatField(blank=True, null=True)
    magnesium = models.FloatField(blank=True, null=True)
    mean_frict_vel = models.FloatField(blank=True, null=True)
    mean_peak_frict_vel = models.FloatField(blank=True, null=True)
    methane = models.FloatField(blank=True, null=True)
    n_alkanes_name = models.CharField(max_length=90, blank=True)
    n_alkanes_concentration = models.FloatField(blank=True, null=True)
    nitrate = models.FloatField(blank=True, null=True)
    nitrite = models.FloatField(blank=True, null=True)
    nitro = models.FloatField(blank=True, null=True)
    org_carb = models.FloatField(blank=True, null=True)
    org_matter = models.FloatField(blank=True, null=True)
    org_nitro = models.FloatField(blank=True, null=True)
    organism_name = models.CharField(max_length=90, blank=True)
    organism_count = models.FloatField(blank=True, null=True)
    oxy_stat_samp = models.CharField(max_length=90, blank=True)
    part_org_carb = models.FloatField(blank=True, null=True)
    perturbation_type = models.CharField(max_length=90, blank=True)
    perturbation_interval = models.CharField(max_length=90, blank=True)
    petroleum_hydrocarb = models.FloatField(blank=True, null=True)
    ph = models.FloatField(blank=True, null=True)
    phaeopigments_type = models.CharField(max_length=90, blank=True)
    phaeopigments_concentration = models.FloatField(blank=True, null=True)
    phosphate = models.FloatField(blank=True, null=True)
    phosplipid_fatt_acid_name = models.CharField(max_length=90, blank=True)
    phosplipid_fatt_acid_concentration = models.FloatField(blank=True, null=True)
    potassium = models.FloatField(blank=True, null=True)
    pressure = models.FloatField(blank=True, null=True)
    redox_potential = models.FloatField(blank=True, null=True)
    rel_to_oxygen = models.CharField(max_length=90, blank=True)
    salinity = models.FloatField(blank=True, null=True)
    samp_collect_device = models.CharField(max_length=90, blank=True)
    samp_mat_process = models.CharField(max_length=90, blank=True)
    samp_size = models.FloatField(blank=True, null=True)
    samp_store_dur = models.FloatField(blank=True, null=True)
    samp_store_loc = models.CharField(max_length=90, blank=True)
    samp_store_temp = models.FloatField(blank=True, null=True)
    silicate = models.FloatField(blank=True, null=True)
    sodium = models.FloatField(blank=True, null=True)
    sulfate = models.FloatField(blank=True, null=True)
    sulfide = models.FloatField(blank=True, null=True)
    temp = models.FloatField(blank=True, null=True)
    tot_carb = models.FloatField(blank=True, null=True)
    tot_nitro = models.FloatField(blank=True, null=True)
    tot_org_carb = models.FloatField(blank=True, null=True)
    turbidity = models.FloatField(blank=True, null=True)
    water_content = models.FloatField(blank=True, null=True)


class Water(models.Model):
    sampleid = models.ForeignKey(Sample)
    projectid = models.ForeignKey(Project)
    refid = models.ForeignKey(Reference)

    alkalinity = models.FloatField(blank=True, null=True)
    alkyl_diethers = models.FloatField(blank=True, null=True)
    altitude = models.FloatField(blank=True, null=True)
    aminopept_act = models.FloatField(blank=True, null=True)
    ammonium = models.FloatField(blank=True, null=True)
    atmospheric_data = models.FloatField(blank=True, null=True)
    bac_prod = models.FloatField(blank=True, null=True)
    bac_resp = models.FloatField(blank=True, null=True)
    bacteria_carb_prod = models.FloatField(blank=True, null=True)
    biomass_part_name = models.CharField(max_length=90, blank=True)
    biomass_amount = models.FloatField(blank=True, null=True)
    bishomohopanol = models.FloatField(blank=True, null=True)
    bromide = models.FloatField(blank=True, null=True)
    calcium = models.FloatField(blank=True, null=True)
    carb_nitro_ratio = models.FloatField(blank=True, null=True)
    chem_administration_name = models.CharField(max_length=90, blank=True)
    chem_administration_time = models.CharField(max_length=90, blank=True)
    chloride = models.FloatField(blank=True, null=True)
    chlorophyll = models.FloatField(blank=True, null=True)
    conduc = models.FloatField(blank=True, null=True)
    density = models.FloatField(blank=True, null=True)
    diether_lipids = models.FloatField(blank=True, null=True)
    diss_carb_dioxide = models.FloatField(blank=True, null=True)
    diss_hydrogen = models.FloatField(blank=True, null=True)
    diss_inorg_carb = models.FloatField(blank=True, null=True)
    diss_inorg_nitro = models.FloatField(blank=True, null=True)
    diss_inorg_phosp = models.FloatField(blank=True, null=True)
    diss_org_carb = models.FloatField(blank=True, null=True)
    diss_org_nitro = models.FloatField(blank=True, null=True)
    diss_oxygen = models.FloatField(blank=True, null=True)
    down_par = models.FloatField(blank=True, null=True)
    elev = models.FloatField(blank=True, null=True)
    fluor = models.FloatField(blank=True, null=True)
    glucosidase_act = models.FloatField(blank=True, null=True)
    light_intensity = models.FloatField(blank=True, null=True)
    magnesium = models.FloatField(blank=True, null=True)
    mean_frict_vel = models.FloatField(blank=True, null=True)
    mean_peak_frict_vel = models.FloatField(blank=True, null=True)
    n_alkanes = models.FloatField(blank=True, null=True)
    nitrate = models.FloatField(blank=True, null=True)
    nitrite = models.FloatField(blank=True, null=True)
    nitro = models.FloatField(blank=True, null=True)
    org_carb = models.FloatField(blank=True, null=True)
    org_matter = models.FloatField(blank=True, null=True)
    org_nitro = models.FloatField(blank=True, null=True)
    organism_name = models.CharField(max_length=90, blank=True)
    organism_count = models.FloatField(blank=True, null=True)
    oxy_stat_samp = models.CharField(max_length=90, blank=True)
    part_org_carb = models.FloatField(blank=True, null=True)
    part_org_nitro = models.FloatField(blank=True, null=True)
    perturbation_type = models.CharField(max_length=90, blank=True)
    perturbation_interval = models.CharField(max_length=90, blank=True)
    pretroleum_hydrocarb = models.FloatField(blank=True, null=True)
    ph = models.FloatField(blank=True, null=True)
    phaeopigments = models.FloatField(blank=True, null=True)
    phosphate = models.FloatField(blank=True, null=True)
    phosplipid_fatt_acid = models.FloatField(blank=True, null=True)
    photon_flux = models.FloatField(blank=True, null=True)
    potassium = models.FloatField(blank=True, null=True)
    pressure = models.FloatField(blank=True, null=True)
    primary_prod = models.FloatField(blank=True, null=True)
    redox_potential = models.FloatField(blank=True, null=True)
    rel_to_oxygen = models.CharField(max_length=90, blank=True)
    samp_mat_process = models.CharField(max_length=90, blank=True)
    samp_salinity = models.FloatField(blank=True, null=True)
    samp_size = models.FloatField(blank=True, null=True)
    samp_store_dur = models.FloatField(blank=True, null=True)
    samp_store_loc = models.CharField(max_length=90, blank=True)
    samp_store_temp = models.FloatField(blank=True, null=True)
    samp_vol_we_dna_ext = models.FloatField(blank=True, null=True)
    silicate = models.FloatField(blank=True, null=True)
    sodium = models.FloatField(blank=True, null=True)
    soluble_react_phosp = models.FloatField(blank=True, null=True)
    source_material_id = models.CharField(max_length=90, blank=True)
    sulfate = models.FloatField(blank=True, null=True)
    sulfide = models.FloatField(blank=True, null=True)
    suspen_part_matter = models.FloatField(blank=True, null=True)
    temp = models.FloatField(blank=True, null=True)
    tidal_stage = models.CharField(max_length=90, blank=True)
    tot_depth_water_col = models.FloatField(blank=True, null=True)
    tot_diss_nitro = models.FloatField(blank=True, null=True)
    tot_inorg_nitro = models.FloatField(blank=True, null=True)
    tot_nitro = models.FloatField(blank=True, null=True)
    tot_part_carb = models.FloatField(blank=True, null=True)
    tot_phosp = models.FloatField(blank=True, null=True)
    water_current_direction = models.CharField(max_length=90, blank=True)
    water_current_magnitude = models.FloatField(blank=True, null=True)


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
    kingdomid = models.TextField(primary_key=True)
    kingdomName = models.CharField(max_length=90, blank=True)


class Phyla(models.Model):
    kingdomid = models.ForeignKey(Kingdom)
    phylaid = models.TextField(primary_key=True)
    phylaName = models.CharField(max_length=90, blank=True)


class Class(models.Model):
    kingdomid = models.ForeignKey(Kingdom)
    phylaid = models.ForeignKey(Phyla)
    classid = models.TextField(primary_key=True)
    className = models.CharField(max_length=90, blank=True)


class Order(models.Model):
    kingdomid = models.ForeignKey(Kingdom)
    phylaid = models.ForeignKey(Phyla)
    classid = models.ForeignKey(Class)
    orderid = models.TextField(primary_key=True)
    orderName = models.CharField(max_length=90, blank=True)


class Family(models.Model):
    kingdomid = models.ForeignKey(Kingdom)
    phylaid = models.ForeignKey(Phyla)
    classid = models.ForeignKey(Class)
    orderid = models.ForeignKey(Order)
    familyid = models.TextField(primary_key=True)
    familyName = models.CharField(max_length=90, blank=True)


class Genus(models.Model):
    kingdomid = models.ForeignKey(Kingdom)
    phylaid = models.ForeignKey(Phyla)
    classid = models.ForeignKey(Class)
    orderid = models.ForeignKey(Order)
    familyid = models.ForeignKey(Family)
    genusid = models.TextField(primary_key=True)
    genusName = models.CharField(max_length=90, blank=True)


class Species(models.Model):
    kingdomid = models.ForeignKey(Kingdom)
    phylaid = models.ForeignKey(Phyla)
    classid = models.ForeignKey(Class)
    orderid = models.ForeignKey(Order)
    familyid = models.ForeignKey(Family)
    genusid = models.ForeignKey(Genus)
    speciesid = models.TextField(primary_key=True)
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


class ko_lvl1(models.Model):
    ko_lvl1_id = models.TextField(primary_key=True)
    ko_lvl1_name = models.TextField()


class ko_lvl2(models.Model):
    ko_lvl1_id = models.ForeignKey(ko_lvl1)
    ko_lvl2_id = models.TextField(primary_key=True)
    ko_lvl2_name = models.TextField()


class ko_lvl3(models.Model):
    ko_lvl1_id = models.ForeignKey(ko_lvl1)
    ko_lvl2_id = models.ForeignKey(ko_lvl2)
    ko_lvl3_id = models.TextField(primary_key=True)
    ko_lvl3_name = models.TextField()


class ko_entry(models.Model):
    ko_lvl1_id = models.ForeignKey(ko_lvl1)
    ko_lvl2_id = models.ForeignKey(ko_lvl2)
    ko_lvl3_id = models.ForeignKey(ko_lvl3)
    ko_lvl4_id = models.TextField(primary_key=True)
    ko_orthology = models.TextField(blank=True)
    ko_name = models.TextField(blank=True)
    ko_desc = models.TextField(blank=True)


class nz_lvl1(models.Model):
    nz_lvl1_id = models.TextField(primary_key=True)
    nz_lvl1_name = models.TextField()


class nz_lvl2(models.Model):
    nz_lvl1_id = models.ForeignKey(nz_lvl1)
    nz_lvl2_id = models.TextField(primary_key=True)
    nz_lvl2_name = models.TextField()


class nz_lvl3(models.Model):
    nz_lvl1_id = models.ForeignKey(nz_lvl1)
    nz_lvl2_id = models.ForeignKey(nz_lvl2)
    nz_lvl3_id = models.TextField(primary_key=True)
    nz_lvl3_name = models.TextField()


class nz_lvl4(models.Model):
    nz_lvl1_id = models.ForeignKey(nz_lvl1)
    nz_lvl2_id = models.ForeignKey(nz_lvl2)
    nz_lvl3_id = models.ForeignKey(nz_lvl3)
    nz_lvl4_id = models.TextField(primary_key=True)
    nz_lvl4_name = models.TextField()


class nz_entry(models.Model):
    nz_lvl1_id = models.ForeignKey(nz_lvl1)
    nz_lvl2_id = models.ForeignKey(nz_lvl2)
    nz_lvl3_id = models.ForeignKey(nz_lvl3)
    nz_lvl4_id = models.ForeignKey(nz_lvl4)
    nz_lvl5_id = models.TextField(primary_key=True)
    nz_orthology = models.TextField(blank=True)
    nz_name = models.TextField(blank=True)
    nz_desc = models.TextField(blank=True)


class PICRUSt(models.Model):
    speciesid = models.ForeignKey(Species)
    rRNACount = models.FloatField(blank=True, null=True)
    geneCount = JSONField(default={})


class UserProfile(models.Model):
    user = models.OneToOneField(Users, on_delete=models.CASCADE)
    firstName = models.CharField(blank=True, max_length=200)
    lastName = models.CharField(blank=True, max_length=200)
    affiliation = models.CharField(blank=True, max_length=200)
    city = models.CharField(blank=True, max_length=200)
    state = models.CharField(blank=True, max_length=100)
    country = models.CharField(blank=True, max_length=100)
    zip = models.CharField(blank=True, max_length=50)
    phone = models.CharField(blank=True, max_length=100)
    reference = models.CharField(blank=True, max_length=100)
    purpose = models.CharField(blank=True, max_length=100)


'''
class PLFA(models.Model):
    content_type = models.ForeignKey(ContentType)
    object_id = models.PositiveIntegerField()
    taxaKey = GenericForeignKey('content_type', 'object_id')  # contains type and content of taxa data
    label = models.TextField(blank=True)  # label/flag/id for PLFA (not sure what this looks like)
'''




