import ast
from django.db import models
from django.contrib.auth.models import User as Users
from django.db.models.signals import post_save
from django.dispatch import receiver

queue = 0


def addQueue():
    global queue
    queue += 1


def getQueue():
    return queue


def subQueue():
    global queue
    queue -= 1


# given a "list" in string form using given delimiter, add new element if it doesn't exist
def appendTextList(text, delimiter, addMe):
    if addMe in text.split(delimiter):
        return text
    if text == "":
        return text + addMe
    return text + delimiter + addMe


# given a "list" as a string using given delimiter, remove element removeMe if it exists
def redactTextList(text, delimiter, removeMe):
    mySplitText = text.split(delimiter)
    if removeMe in mySplitText:
        text = ""
        for part in mySplitText:
            if str(part) != str(removeMe) and str(part) != "":
                if text is not "":
                    text += delimiter
                text += part
    return text


def updateList(myList, delimiter, remove=None, add=None):
    # this helper function cannot deal with the double none case
    if remove is not None:
        myList = redactTextList(myList, delimiter, remove)
    if add is not None:
        myList = appendTextList(myList, delimiter, add)
    return myList


class Project(models.Model):
    projectid = models.CharField(max_length=50, primary_key=True)
    status = models.CharField(max_length=15, blank=False)  # visibility, "public" "private", error if other value
    wip = models.BooleanField(default=False)
    projectType = models.CharField(max_length=90, blank=False)
    project_name = models.CharField(max_length=100, blank=True)
    project_desc = models.TextField(blank=True)
    start_date = models.CharField(max_length=25, blank=True)
    end_date = models.CharField(max_length=25, blank=True)
    pi_last = models.CharField(max_length=75, blank=True)
    pi_first = models.CharField(max_length=50, blank=True)
    pi_affiliation = models.CharField(max_length=100, blank=True)
    pi_email = models.EmailField(blank=True)
    pi_phone = models.CharField(max_length=15, blank=True)

    owner = models.ForeignKey(Users, null=True, default=None)  # user who first uploaded files for this project
    # note that future additions to the project will require edit permissions to be granted
    # by the owner of the project. Retroactive solution is to fill this slot manually

    whitelist_view = models.TextField(blank=True)  # default empty, contains list of user ids
    # who can view the project if it is private
    # if project is public, this field is ignored (and gets reset if a private project is changed to public)
    # field is text, user ids separated by ;
    whitelist_edit = models.TextField(blank=True)  # default empty, contains list of user ids
    # who can edit the project (update, reprocess, or upload extra/new references)
    # edit perm does NOT grant removal or permissions management perms
    # field is text, user ids separated by ;

    def update_viewPerms(self, remove=None, add=None):
        if remove is None and add is None:
            self.whitelist_view = ""
        self.whitelist_view = updateList(self.whitelist_view, ";", remove, add)

    def update_editPerms(self, remove=None, add=None):
        if remove is None and add is None:
            self.whitelist_edit = ""
        self.whitelist_edit = updateList(self.whitelist_edit, ";", remove, add)

    def __unicode__(self):
        return self.project_name


class Reference(models.Model):
    # save the location and data associated with a reference file (for finding names and such)
    refid = models.CharField(max_length=50, primary_key=True)
    raw = models.BooleanField()
    projectid = models.ForeignKey(Project)
    path = models.CharField(max_length=90)
    source = models.CharField(max_length=90)
    author = models.ForeignKey(Users, related_name='entries')

    def __unicode__(self):
        return self.path


class Sample(models.Model):
    # core sample model, this stores most general information about a sample within a project
    # gets referenced by subtypes such as Soil, Air, Water, etc
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

    reads = models.IntegerField(default=0, null=True, editable=False)

    def natural_key(self):
        return self.sampleid

    def __unicode__(self):
        return self.sample_name


class Human_Associated(models.Model):
    # data from a sample specifically for human_associated projects, following the fields described by mothur
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
    # data from a sample specifically for soil related projects, following the fields described by mothur
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

    # entries added for Cornell's C.A.S.H.  (not duplicating values already stored)
    # these should be added to the standard myphylodb meta file, all are quantitative and some depend on existing soil vals
    # soil vals which appear in cornell's set and the base myphylodb set should have a redirect from cornell to main
    soil_texture_sand = models.FloatField(blank=True, null=True)
    soil_texture_silt = models.FloatField(blank=True, null=True)
    soil_texture_clay = models.FloatField(blank=True, null=True)
    soil_water_cap = models.FloatField(blank=True, null=True)
    soil_water_cap_rating = models.FloatField(blank=True, null=True)
    soil_surf_hardness = models.FloatField(blank=True, null=True)
    soil_surf_hardness_rating = models.FloatField(blank=True, null=True)
    soil_subsurf_hardness = models.FloatField(blank=True, null=True)
    soil_subsurf_hardness_rating = models.FloatField(blank=True, null=True)
    soil_agg_stability = models.FloatField(blank=True, null=True)
    soil_agg_stability_rating = models.FloatField(blank=True, null=True)
    soil_organic_matter = models.FloatField(blank=True, null=True)
    soil_organic_matter_rating = models.FloatField(blank=True, null=True)
    soil_ACE_protein_index = models.FloatField(blank=True, null=True)
    soil_ACE_protein_index_rating = models.FloatField(blank=True, null=True)
    soil_root_pathogen_pressure = models.FloatField(blank=True, null=True)
    soil_root_pathogen_pressure_rating = models.FloatField(blank=True, null=True)
    soil_respiration_four_day = models.FloatField(blank=True, null=True)
    soil_respiration_four_day_rating = models.FloatField(blank=True, null=True)
    soil_active_C = models.FloatField(blank=True, null=True)
    soil_active_C_rating = models.FloatField(blank=True, null=True)
    soil_pH_rating = models.FloatField(blank=True, null=True)
    soil_p_rating = models.FloatField(blank=True, null=True)
    soil_k_rating = models.FloatField(blank=True, null=True)
    soil_minor_elements_rating = models.FloatField(blank=True, null=True)
    CASH_SHI_rating = models.FloatField(blank=True, null=True)

    def __unicode__(self):
        return self.sampleid.sample_name


class Air(models.Model):
    # data from a sample specifically for air projects, following the fields described by mothur
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
    # data from a sample specifically for microbial projects, following the fields described by mothur
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
    phosplipid_fatt_acid_conc = models.FloatField(blank=True, null=True)
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
    # data from a sample specifically for water projects, following the fields described by mothur
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
    # data from a sample with fields likely specific to the project, up to the user to remember which fields are comparable
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

    def __unicode__(self):
        return self.sampleid.sample_name


class DaymetData(models.Model):
    # data from the daymet database, specific to regions and months the user has selected (to add to meta variables)
    user = models.ForeignKey(Users, unique=True)   # at most one per user (most recent norm or none)

    # delete the old daymet data at start of norm
    # then just check if daymet data exists, since it won't in most cases already

    # have ";" delimited strings for sampleID and all daymet columns
    # sync these strings on position, so sampleID[0][dayl..s] = dayl[0], etc
    sampleIDs = models.TextField(blank=True)
    # "year" "yday" "dayl..s." "prcp..mm.day." "srad..W.m.2."  "swe..kg.m.2."  "tmax..deg.c."  "tmin..deg.c."  "vp..Pa."
    #year = models.TextField(blank=True)
    #yday = models.TextField(blank=True)
    dayl = models.TextField(blank=True)
    prcp = models.TextField(blank=True)
    srad = models.TextField(blank=True)
    swe = models.TextField(blank=True)
    tmax = models.TextField(blank=True)
    tmin = models.TextField(blank=True)
    tmean = models.TextField(blank=True)
    vp = models.TextField(blank=True)


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


class OTU_99(models.Model):
    kingdomid = models.ForeignKey(Kingdom)
    phylaid = models.ForeignKey(Phyla)
    classid = models.ForeignKey(Class)
    orderid = models.ForeignKey(Order)
    familyid = models.ForeignKey(Family)
    genusid = models.ForeignKey(Genus)
    speciesid = models.ForeignKey(Species)
    otuid = models.TextField(primary_key=True)  # not specifying yet getting duplicates????
    otuName = models.CharField(max_length=90, blank=True)
    otuSeq = models.TextField(blank=True)


class Profile(models.Model):
    # taxonomic profile, not to be confused with UserProfile
    projectid = models.ForeignKey(Project)
    sampleid = models.ForeignKey(Sample)

    kingdomid = models.ForeignKey(Kingdom)
    phylaid = models.ForeignKey(Phyla)
    classid = models.ForeignKey(Class)
    orderid = models.ForeignKey(Order)
    familyid = models.ForeignKey(Family)
    genusid = models.ForeignKey(Genus)
    speciesid = models.ForeignKey(Species)
    otuid = models.ForeignKey(OTU_99)
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
    otuid = models.ForeignKey(OTU_99)
    rRNACount = models.FloatField(blank=True, null=True)
    geneList = models.TextField(blank=True)


class koOtuList(models.Model):
    # used as a means of storing kegg data per gene, for significantly faster querying on path pages
    koID = models.TextField(primary_key=True, db_index=True)
    otuList = models.TextField(blank=True)  # ';' separated list of otu names associated with this ko


class UserProfile(models.Model):
    # Profile of extra data for each user, such as permissions and personal info

    # actual user object from Django account manager
    user = models.OneToOneField(Users, related_name='profile')
    # personal data
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
    # a unique string attached to the user, changed whenever a new selection is made
    # sent from the browser to verify against, if they do not match, the user's browser is not using the correct dataset
    dataID = models.CharField(blank=True, max_length=200, default="TEMPID")

    # recent analysis rid storage field
    recentRIDs = models.CharField(blank=True, default="", max_length=32)   # semi-colon delimited list of RIDs, stored oldest to newest

    def popRecentRID(self):
        # get the oldest RID from the list and return it, removing it in the process
        splitRIDs = self.recentRIDs.split(";")
        if len(splitRIDs) > 0:
            oldestRID = splitRIDs[0]
            self.recentRIDs = redactTextList(self.recentRIDs, ";", oldestRID)
            return oldestRID
        return None

    def addRecentRID(self, RID):
        # given an RID, add to recentRIDs such that <= 3 RIDs remain afterwards (oldest gets returned if removed)
        # call popRecentRID for the oldest removal if len recentRIDs >= N
        oldRID = None
        if len(self.recentRIDs.split(";")) >= 3:    # currently limiting history to 3
            # but if the limit is to be changed, simply adjust the 3 here
            oldRID = self.popRecentRID()
        self.recentRIDs = appendTextList(self.recentRIDs, ";", RID)
        return oldRID

    # permissions fields and methods, using usernames (functionally interchangeable with ids, but names are readable)
    hasPermsFrom = models.TextField(blank=True)  # acts as a list via split with ';' delimiter
    # Note: all changes to these lists should occur within perms.py (anything using them as well)
    # (everything related to these fields should always call the appropriate functions)
    gavePermsTo = models.TextField(blank=True)

    def update_hasPerms(self, remove=None, add=None):
        if remove is None and add is None:
            # really needed one of those, can only do this with mode = "has"
            # recalc hasPermsFrom via all gavePermsTo fields
            allProfs = UserProfile.objects.all()
            wipHas = ""
            for prof in allProfs:
                if self.user.username in prof.gavePermsTo.split(";"):
                    wipHas = appendTextList(wipHas, ";", prof.user.username)
            self.hasPermsFrom = wipHas
        self.hasPermsFrom = updateList(self.hasPermsFrom, ";", remove, add)
        return

    def update_gavePerms(self, remove=None, add=None):
        self.gavePermsTo = updateList(self.gavePermsTo, ";", remove, add)


    # need to update visible list whenever this user uploads or deletes a project
    # and whenever another user gives or removes permission, which is stored per project? per user? Could do both
    # this list should NEVER be editable by the user, only updated by server, and user sees RESULT of list (projects)
    privateProjectList = models.TextField(blank=True)

    @receiver(post_save, sender=Users)
    def create_user_profile(sender, instance, created, **kwargs):
        if created:
            UserProfile.objects.create(user=instance)

    @receiver(post_save, sender=Users)
    def save_user_profile(sender, instance, **kwargs):
        instance.profile.save()

    def update_list(self, remove=None, add=None):
        # slightly more complicated than public list for main creation
        # during normal runtime use, only call remove and add
        # need to run the generic creation section once for each user
        if remove is None and add is None:
            # generic list creation function (makes private project list from scratch)
            projString = ""
            whitelist = self.hasPermsFrom.split(';')
            for proj in Project.objects.all():
                # get list of permitted users for this project
                checkList = proj.whitelist_view.split(';')
                # check if this user is in the list
                if self.user.id in checkList or proj.owner.id in whitelist:
                    # user has permission, add to list
                    projString = appendTextList(projString, ",", proj.projectid)
            self.privateProjectList = projString
        self.privateProjectList = updateList(self.privateProjectList, ",", remove, add)


class PublicProjects(models.Model):
    List = models.TextField(blank=True)  # CSV formatted string which contains the project ID's of all public projects
    # need only one of these to exist, must be updated whenever a public project is added, removed
    # or switched to private, also private projects switched to public
    # so this should get called by something in the end of an upload, project delete, and meta file update (if needed)
    # this does not need to be sorted, since this list gets used right before sorting by project name anyways

    def update_list(self, remove=None, add=None):
        # will recreate list if no individual ID is given
        # if both remove and add are given (don't do that), only remove will be respected

        if remove is None and add is None:
            # this is the generic one, technically works for updating a single project BUT
            # a more specific implementation would be significantly more efficient
            # make a list of all public projects
            projString = ""
            for proj in Project.objects.all():
                if proj.status == 'public':
                    projString = appendTextList(projString, ",", proj.projectid)
            self.List = projString
        # helper function call
        self.List = updateList(self.List, ",", remove, add)
