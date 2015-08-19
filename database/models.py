from django.db import models
from django_extensions.db.fields import UUIDField


class Project(models.Model):
    projectType = models.CharField(max_length=45, blank=False)
    projectid = UUIDField(primary_key=True)
    path = models.CharField(max_length=90)
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
    projectid = models.ForeignKey(Project)
    path = models.CharField(max_length=90)
    alignDB = models.CharField(max_length=90, blank=True)
    templateDB = models.CharField(max_length=90, blank=True)
    taxonomyDB = models.CharField(max_length=90, blank=True)


class Sample(models.Model):
    sampleid = UUIDField(primary_key=True)
    projectid = models.ForeignKey(Project)
    sample_name = models.TextField(blank=False)
    organism = models.CharField(max_length=90, blank=True)
    title = models.TextField(blank=True)
    seq_method = models.TextField(blank=True)
    collection_date = models.CharField(max_length=15, blank=True)
    biome = models.CharField(max_length=45, blank=True)
    feature = models.CharField(max_length=45, blank=True)
    geo_loc_country = models.CharField(max_length=45, blank=True)
    geo_loc_state = models.CharField(max_length=45, blank=True)
    geo_loc_city = models.CharField(max_length=45, blank=True)
    geo_loc_farm = models.CharField(max_length=45, blank=True)
    geo_loc_plot = models.CharField(max_length=45, blank=True)
    latitude = models.CharField(max_length=45, blank=True)
    longitude = models.CharField(max_length=45, blank=True)
    material = models.CharField(max_length=45, blank=True)
    elevation = models.CharField(max_length=45, blank=True)

    def natural_key(self):
        return self.sampleid

    def __unicode__(self):
        return unicode(self.sample_name)


class Human_Gut(models.Model):
    sampleid = models.ForeignKey(Sample)
    projectid = models.ForeignKey(Project)

    age = models.CharField(max_length=15, blank=True)  # float unit
    body_mass_index = models.CharField(max_length=45, blank=True)  # float unit
    body_product = models.CharField(max_length=45, blank=True)  # text
    chem_administration = models.CharField(max_length=45, blank=True)  # term; timestamp
    diet = models.CharField(max_length=45, blank=True)  # text
    disease = models.CharField(max_length=45, blank=True)  # term
    ethnicity = models.CharField(max_length=45, blank=True)  # integer|text
    family_relationship = models.CharField(max_length=45, blank=True)  # text; text
    gastrointest_disord = models.CharField(max_length=45, blank=True)  # text
    genotype = models.CharField(max_length=45, blank=True)  # text
    height = models.CharField(max_length=45, blank=True)  # float unit
    host_body_temp = models.CharField(max_length=45, blank=True)  # float unit
    host_subject_id = models.CharField(max_length=45, blank=True)  # text
    ihmc_medication_code = models.CharField(max_length=45, blank=True)  # integer
    last_meal = models.CharField(max_length=45, blank=True)  # text; period
    liver_disord = models.CharField(max_length=45, blank=True)  # text
    medic_hist_perform = models.CharField(max_length=45, blank=True)  # boolean
    nose_throat_disord = models.CharField(max_length=45, blank=True)  # text
    occupation = models.CharField(max_length=45, blank=True)  # integer
    organism_count = models.CharField(max_length=45, blank=True)  # text; float unit
    oxy_stat_samp = models.CharField(max_length=45, blank=True)  # [, 'aerobic', 'anaerobic']
    perturbation = models.CharField(max_length=45, blank=True)  # text; interval
    phenotype = models.CharField(max_length=45, blank=True)  # term
    pulse = models.CharField(max_length=45, blank=True)  # float unit
    rel_to_oxygen = models.CharField(max_length=45, blank=True)  # [, 'aerobe', 'anaerobe', 'facultative', 'microaerophilic', 'microanerobe', 'obligate aerobe', 'obligate anaerobe']
    samp_collect_device = models.CharField(max_length=45, blank=True)  # text
    samp_mat_process = models.CharField(max_length=45, blank=True)  # text|term
    samp_salinity = models.CharField(max_length=45, blank=True)  # float unit
    samp_size = models.CharField(max_length=45, blank=True)  # float unit
    samp_store_dur = models.CharField(max_length=45, blank=True)  # interval
    samp_store_loc = models.CharField(max_length=45, blank=True)  # text
    samp_store_temp = models.CharField(max_length=45, blank=True)  # float unit
    sex = models.CharField(max_length=45, blank=True)  # [, 'male', 'female', 'neuter', 'hermaphrodite', 'not determined']
    special_diet = models.CharField(max_length=45, blank=True)  # [, 'low carb', 'reduced calorie', 'vegetarian', 'other(to be specified)']
    temp = models.CharField(max_length=45, blank=True)  # float unit
    tissue = models.CharField(max_length=45, blank=True)  # not specified, probably text
    tot_mass = models.CharField(max_length=45, blank=True)  # float unit
    user_defined = models.CharField(max_length=45, blank=True)  # Perhaps merge this with the other user defined section, put it under common


class Microbial(models.Model):
    sampleid = models.ForeignKey(Sample)
    projectid = models.ForeignKey(Project)

    alkalinity = models.CharField(max_length=15, blank=True)
    alkyl_diethers = models.CharField(max_length=15, blank=True)
    altitude = models.CharField(max_length=15, blank=True)
    aminopept_act = models.CharField(max_length=15, blank=True)
    ammonium = models.CharField(max_length=15, blank=True)
    bacteria_carb_prod = models.CharField(max_length=15, blank=True)
    biomass = models.CharField(max_length=15, blank=True)
    bishomohopanol = models.CharField(max_length=15, blank=True)
    bromide = models.CharField(max_length=15, blank=True)
    calcium = models.CharField(max_length=15, blank=True)
    carb_nitro_ratio = models.CharField(max_length=15, blank=True)
    chem_administration = models.CharField(max_length=15, blank=True)
    chloride = models.CharField(max_length=15, blank=True)
    chlorophyll = models.CharField(max_length=15, blank=True)
    diether_lipids = models.CharField(max_length=15, blank=True)
    diss_carb_dioxide = models.CharField(max_length=15, blank=True)
    diss_hydrogen = models.CharField(max_length=15, blank=True)
    diss_inorg_carb = models.CharField(max_length=15, blank=True)
    diss_org_carb = models.CharField(max_length=15, blank=True)
    diss_org_nitro = models.CharField(max_length=15, blank=True)
    diss_oxygen = models.CharField(max_length=15, blank=True)
    glucosidase_act = models.CharField(max_length=15, blank=True)
    magnesium = models.CharField(max_length=15, blank=True)
    mean_frict_vel = models.CharField(max_length=15, blank=True)
    mean_peak_frict_vel = models.CharField(max_length=15, blank=True)
    methane = models.CharField(max_length=15, blank=True)
    n_alkanes = models.CharField(max_length=15, blank=True)
    nitrate = models.CharField(max_length=15, blank=True)
    nitrite = models.CharField(max_length=15, blank=True)
    nitro = models.CharField(max_length=15, blank=True)
    org_carb = models.CharField(max_length=15, blank=True)
    org_matter = models.CharField(max_length=15, blank=True)
    org_nitro = models.CharField(max_length=15, blank=True)
    organism_count = models.CharField(max_length=15, blank=True)
    oxy_stat_samp = models.CharField(max_length=15, blank=True)
    part_org_carb = models.CharField(max_length=15, blank=True)
    perturbation = models.CharField(max_length=15, blank=True)
    petroleum_hydrocarb = models.CharField(max_length=15, blank=True)
    ph = models.CharField(max_length=15, blank=True)
    phaeopigments = models.CharField(max_length=15, blank=True)
    phosphate = models.CharField(max_length=15, blank=True)
    phosplipid_fatt_acid = models.CharField(max_length=15, blank=True)
    potassium = models.CharField(max_length=15, blank=True)
    pressure = models.CharField(max_length=15, blank=True)
    redox_potential = models.CharField(max_length=15, blank=True)
    rel_to_oxygen = models.CharField(max_length=15, blank=True)
    salinity = models.CharField(max_length=15, blank=True)
    samp_collect_device = models.CharField(max_length=15, blank=True)
    samp_mat_process = models.CharField(max_length=15, blank=True)
    samp_size = models.CharField(max_length=15, blank=True)
    samp_store_dur = models.CharField(max_length=15, blank=True)
    samp_store_loc = models.CharField(max_length=15, blank=True)
    samp_store_temp = models.CharField(max_length=15, blank=True)
    silicate = models.CharField(max_length=15, blank=True)
    sodium = models.CharField(max_length=15, blank=True)
    sulfate = models.CharField(max_length=15, blank=True)
    sulfide = models.CharField(max_length=15, blank=True)
    temp = models.CharField(max_length=15, blank=True)
    tot_carb = models.CharField(max_length=15, blank=True)
    tot_nitro = models.CharField(max_length=15, blank=True)
    tot_org_carb = models.CharField(max_length=15, blank=True)
    turbidity = models.CharField(max_length=15, blank=True)
    water_content = models.CharField(max_length=15, blank=True)
    user_defined = models.CharField(max_length=15, blank=True)


class Human_Associated(models.Model):
    sampleid = models.ForeignKey(Sample)
    projectid = models.ForeignKey(Project)

    age = models.CharField(max_length=15, blank=True)
    amniotic_fluid_color = models.CharField(max_length=15, blank=True)
    blood_blood_disord = models.CharField(max_length=15, blank=True)
    body_mass_index = models.CharField(max_length=15, blank=True)
    body_product = models.CharField(max_length=15, blank=True)
    chem_administration = models.CharField(max_length=15, blank=True)
    diet = models.CharField(max_length=15, blank=True)
    diet_last_six_month = models.CharField(max_length=15, blank=True)
    disease = models.CharField(max_length=15, blank=True)
    drug_usage = models.CharField(max_length=15, blank=True)
    ethnicity = models.CharField(max_length=15, blank=True)
    family_relationship = models.CharField(max_length=15, blank=True)
    fetal_health_stat = models.CharField(max_length=15, blank=True)
    genotype = models.CharField(max_length=15, blank=True)
    gestation_state = models.CharField(max_length=15, blank=True)
    height = models.CharField(max_length=15, blank=True)
    hiv_stat = models.CharField(max_length=15, blank=True)
    host_body_temp = models.CharField(max_length=15, blank=True)
    host_subject_id = models.CharField(max_length=15, blank=True)
    ihmc_medication_code = models.CharField(max_length=15, blank=True)
    kidney_disord = models.CharField(max_length=15, blank=True)
    last_meal = models.CharField(max_length=15, blank=True)
    maternal_health_stat= models.CharField(max_length=15, blank=True)
    medic_hist_perform = models.CharField(max_length=15, blank=True)
    nose_throat_disord = models.CharField(max_length=15, blank=True)
    occupation = models.CharField(max_length=15, blank=True)
    organism_count = models.CharField(max_length=15, blank=True)
    oxy_stat_samp = models.CharField(max_length=15, blank=True)
    perturbation = models.CharField(max_length=15, blank=True)
    pet_farm_animal = models.CharField(max_length=15, blank=True)
    phenotype = models.CharField(max_length=15, blank=True)
    pulmonary_disord = models.CharField(max_length=15, blank=True)
    pulse = models.CharField(max_length=15, blank=True)
    rel_to_oxygen = models.CharField(max_length=15, blank=True)
    samp_collect_device = models.CharField(max_length=15, blank=True)
    samp_mat_process = models.CharField(max_length=15, blank=True)
    samp_salinity = models.CharField(max_length=15, blank=True)
    samp_size = models.CharField(max_length=15, blank=True)
    samp_store_dur = models.CharField(max_length=15, blank=True)
    samp_store_loc = models.CharField(max_length=15, blank=True)
    samp_store_temp = models.CharField(max_length=15, blank=True)
    sex = models.CharField(max_length=15, blank=True)
    smoker = models.CharField(max_length=15, blank=True)
    study_complt_stat = models.CharField(max_length=15, blank=True)
    temp = models.CharField(max_length=15, blank=True)
    tissue = models.CharField(max_length=15, blank=True)
    tot_mass = models.CharField(max_length=15, blank=True)
    travel_out_six_month = models.CharField(max_length=15, blank=True)
    twin_sibling = models.CharField(max_length=15, blank=True)
    urine_collect_meth = models.CharField(max_length=15, blank=True)
    urogenit_tract_disor = models.CharField(max_length=15, blank=True)
    weight_loss_3_month = models.CharField(max_length=15, blank=True)
    user_defined = models.CharField(max_length=15, blank=True)


class Soil(models.Model):
    sampleid = models.ForeignKey(Sample)
    projectid = models.ForeignKey(Project)

    depth = models.CharField(max_length=45, blank=True)
    pool_dna_extracts = models.CharField(max_length=45, blank=True)
    samp_size = models.CharField(max_length=45, blank=True)
    samp_collection_device = models.TextField(blank=True)
    samp_weight_dna_ext = models.CharField(max_length=45, blank=True)
    sieving = models.CharField(max_length=45, blank=True)
    storage_cond = models.CharField(max_length=45, blank=True)

    annual_season_precpt = models.CharField(max_length=45, blank=True)
    annual_season_temp = models.CharField(max_length=45, blank=True)

    bulk_density = models.CharField(max_length=45, blank=True)
    drainage_class = models.CharField(max_length=45, blank=True)
    fao_class = models.CharField(max_length=45, blank=True)
    horizon = models.CharField(max_length=45, blank=True)
    local_class = models.CharField(max_length=45, blank=True)
    porosity = models.CharField(max_length=45, blank=True)
    profile_position = models.CharField(max_length=45, blank=True)
    slope_aspect = models.CharField(max_length=45, blank=True)
    slope_gradient = models.CharField(max_length=45, blank=True)
    soil_type = models.CharField(max_length=45, blank=True)
    texture_class = models.CharField(max_length=45, blank=True)
    water_content_soil = models.CharField(max_length=45, blank=True)

    pH = models.CharField(max_length=45, blank=True)
    EC = models.CharField(max_length=45, blank=True)
    tot_C = models.CharField(max_length=45, blank=True)
    tot_OM = models.CharField(max_length=45, blank=True)
    tot_N = models.CharField(max_length=45, blank=True)
    NO3_N = models.CharField(max_length=45, blank=True)
    NH4_N = models.CharField(max_length=45, blank=True)
    P = models.CharField(max_length=45, blank=True)
    K = models.CharField(max_length=45, blank=True)
    S = models.CharField(max_length=45, blank=True)
    Zn = models.CharField(max_length=45, blank=True)
    Fe = models.CharField(max_length=45, blank=True)
    Cu = models.CharField(max_length=45, blank=True)
    Mn = models.CharField(max_length=45, blank=True)
    Ca = models.CharField(max_length=45, blank=True)
    Mg = models.CharField(max_length=45, blank=True)
    Na = models.CharField(max_length=45, blank=True)
    B = models.CharField(max_length=45, blank=True)

    agrochem_amendments = models.CharField(max_length=45, blank=True)
    agrochem_amendments_desc = models.TextField(blank=True)
    biological_amendments = models.CharField(max_length=45, blank=True)
    biological_amendments_desc = models.TextField(blank=True)
    soil_amendments = models.CharField(max_length=45, blank=True)
    soil_amendments_desc = models.TextField(blank=True)
    cover_crop = models.TextField(blank=True)
    crop_rotation = models.TextField(blank=True)
    cur_land_use = models.TextField(blank=True)
    cur_vegetation = models.TextField(blank=True)
    cur_crop = models.TextField(blank=True)
    cur_cultivar = models.TextField(blank=True)
    organic = models.TextField(blank=True)
    previous_land_use = models.TextField(blank=True)
    tillage = models.TextField(blank=True)

    rRNA_copies = models.CharField(max_length=45, blank=True)
    microbial_biomass_C = models.CharField(max_length=45, blank=True)
    microbial_biomass_N = models.CharField(max_length=45, blank=True)
    microbial_respiration = models.CharField(max_length=45, blank=True)


class Air(models.Model):
    sampleid = models.ForeignKey(Sample)
    projectid = models.ForeignKey(Project)

    barometric_press = models.CharField(max_length=45, blank=True)
    carb_dioxide = models.CharField(max_length=45, blank=True)
    carb_monoxide = models.CharField(max_length=45, blank=True)
    chem_administration = models.CharField(max_length=45, blank=True)
    elev = models.CharField(max_length=45, blank=True)
    humidity = models.CharField(max_length=45, blank=True)
    methane = models.CharField(max_length=45, blank=True)
    organism_count = models.CharField(max_length=45, blank=True)
    oxy_stat_samp = models.CharField(max_length=45, blank=True)
    oxygen = models.CharField(max_length=45, blank=True)
    perturbation = models.CharField(max_length=45, blank=True)
    pollutants = models.CharField(max_length=45, blank=True)
    rel_to_oxygen = models.CharField(max_length=45, blank=True)
    resp_part_matter = models.CharField(max_length=45, blank=True)
    samp_collect_device = models.CharField(max_length=45, blank=True)
    samp_mat_process = models.CharField(max_length=45, blank=True)
    samp_salinity = models.CharField(max_length=45, blank=True)
    samp_size = models.CharField(max_length=45, blank=True)
    samp_store_dur = models.CharField(max_length=45, blank=True)
    samp_store_loc = models.CharField(max_length=45, blank=True)
    samp_store_temp = models.CharField(max_length=45, blank=True)
    solar_irradiance = models.CharField(max_length=45, blank=True)
    temp = models.CharField(max_length=45, blank=True)
    ventilation_rate = models.CharField(max_length=45, blank=True)
    ventilation_type = models.CharField(max_length=45, blank=True)
    volatile_org_comp = models.CharField(max_length=45, blank=True)
    wind_direction = models.CharField(max_length=45, blank=True)
    wind_speed = models.CharField(max_length=45, blank=True)
    user_defined = models.CharField(max_length=45, blank=True)


class Water(models.Model):
     sampleid = models.ForeignKey(Sample)
     projectid = models.ForeignKey(Project)

     alkalinity = models.CharField(max_length=45, blank=True)
     alkyl_diethers = models.CharField(max_length=45, blank=True)
     altitude = models.CharField(max_length=45, blank=True)
     aminopept_act = models.CharField(max_length=45, blank=True)
     ammonium = models.CharField(max_length=45, blank=True)
     atmospheric_data = models.CharField(max_length=45, blank=True)
     bac_prod = models.CharField(max_length=45, blank=True)
     bac_resp = models.CharField(max_length=45, blank=True)
     bacteria_carb_prod = models.CharField(max_length=45, blank=True)
     biomass = models.CharField(max_length=45, blank=True)
     bishomohopanol = models.CharField(max_length=45, blank=True)
     bromide = models.CharField(max_length=45, blank=True)
     calcium = models.CharField(max_length=45, blank=True)
     carb_nitro_ratio = models.CharField(max_length=45, blank=True)
     chem_administration = models.CharField(max_length=45, blank=True)
     chloride = models.CharField(max_length=45, blank=True)
     chlorophyll = models.CharField(max_length=45, blank=True)
     conduc = models.CharField(max_length=45, blank=True)
     density = models.CharField(max_length=45, blank=True)
     diether_lipids = models.CharField(max_length=45, blank=True)
     diss_carb_dioxide = models.CharField(max_length=45, blank=True)
     diss_hydrogen = models.CharField(max_length=45, blank=True)
     diss_inorg_carb = models.CharField(max_length=45, blank=True)
     diss_inorg_nitro = models.CharField(max_length=45, blank=True)
     diss_inorg_phosp = models.CharField(max_length=45, blank=True)
     diss_org_carb = models.CharField(max_length=45, blank=True)
     diss_org_nitro = models.CharField(max_length=45, blank=True)
     diss_oxygen = models.CharField(max_length=45, blank=True)
     down_par = models.CharField(max_length=45, blank=True)
     elev = models.CharField(max_length=45, blank=True)
     fluor = models.CharField(max_length=45, blank=True)
     glucosidase_act = models.CharField(max_length=45, blank=True)
     light_intensity = models.CharField(max_length=45, blank=True)
     magnesium = models.CharField(max_length=45, blank=True)
     mean_frict_vel = models.CharField(max_length=45, blank=True)
     mean_peak_frict_vel = models.CharField(max_length=45, blank=True)
     n_alkanes = models.CharField(max_length=45, blank=True)
     nitrate = models.CharField(max_length=45, blank=True)
     nitrite = models.CharField(max_length=45, blank=True)
     nitro = models.CharField(max_length=45, blank=True)
     org_carb = models.CharField(max_length=45, blank=True)
     org_matter = models.CharField(max_length=45, blank=True)
     org_nitro = models.CharField(max_length=45, blank=True)
     organism_count = models.CharField(max_length=45, blank=True)
     oxy_stat_samp = models.CharField(max_length=45, blank=True)
     part_org_carb = models.CharField(max_length=45, blank=True)
     part_org_nitro = models.CharField(max_length=45, blank=True)
     perturbation = models.CharField(max_length=45, blank=True)
     petroleum_hydrocarb = models.CharField(max_length=45, blank=True)
     ph = models.CharField(max_length=45, blank=True)
     phaeopigments = models.CharField(max_length=45, blank=True)
     phosphate = models.CharField(max_length=45, blank=True)
     phosplipid_fatt_acid = models.CharField(max_length=45, blank=True)
     photon_flux = models.CharField(max_length=45, blank=True)
     potassium = models.CharField(max_length=45, blank=True)
     pressure = models.CharField(max_length=45, blank=True)
     primary_prod = models.CharField(max_length=45, blank=True)
     redox_potential = models.CharField(max_length=45, blank=True)
     rel_to_oxygen = models.CharField(max_length=45, blank=True)
     samp_mat_process = models.CharField(max_length=45, blank=True)
     samp_salinity = models.CharField(max_length=45, blank=True)
     samp_size = models.CharField(max_length=45, blank=True)
     samp_store_dur = models.CharField(max_length=45, blank=True)
     samp_store_loc = models.CharField(max_length=45, blank=True)
     samp_store_temp = models.CharField(max_length=45, blank=True)
     samp_vol_we_dna_ext = models.CharField(max_length=45, blank=True)
     silicate = models.CharField(max_length=45, blank=True)
     sodium = models.CharField(max_length=45, blank=True)
     soluble_react_phosp = models.CharField(max_length=45, blank=True)
     source_material_id = models.CharField(max_length=45, blank=True)
     sulfate = models.CharField(max_length=45, blank=True)
     sulfide = models.CharField(max_length=45, blank=True)
     suspend_part_matter = models.CharField(max_length=45, blank=True)
     temp = models.CharField(max_length=45, blank=True)
     tidal_stage = models.CharField(max_length=45, blank=True)
     tot_depth_water_col = models.CharField(max_length=45, blank=True)
     tot_diss_nitro = models.CharField(max_length=45, blank=True)
     tot_inorg_nitro = models.CharField(max_length=45, blank=True)
     tot_nitro = models.CharField(max_length=45, blank=True)
     tot_part_carb = models.CharField(max_length=45, blank=True)
     tot_phosp = models.CharField(max_length=45, blank=True)
     water_current = models.CharField(max_length=45, blank=True)
     user_defined = models.CharField(max_length=45, blank=True)


class User(models.Model):
    sampleid = models.ForeignKey(Sample)
    projectid = models.ForeignKey(Project)
    usr_cat1 = models.TextField(blank=True)
    usr_cat2 = models.TextField(blank=True)
    usr_cat3 = models.TextField(blank=True)
    usr_cat4 = models.TextField(blank=True)
    usr_cat5 = models.TextField(blank=True)
    usr_cat6 = models.TextField(blank=True)
    usr_quant1 = models.CharField(max_length=45, blank=True)
    usr_quant2 = models.CharField(max_length=45, blank=True)
    usr_quant3 = models.CharField(max_length=45, blank=True)
    usr_quant4 = models.CharField(max_length=45, blank=True)
    usr_quant5 = models.CharField(max_length=45, blank=True)
    usr_quant6 = models.CharField(max_length=45, blank=True)

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
    kingdomid = models.ForeignKey(Kingdom)
    phylaid = models.ForeignKey(Phyla)
    classid = models.ForeignKey(Class)
    orderid = models.ForeignKey(Order)
    familyid = models.ForeignKey(Family)
    genusid = models.ForeignKey(Genus)
    speciesid = models.ForeignKey(Species)
    count = models.IntegerField()


