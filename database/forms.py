from django import forms
from django_countries import countries
from localflavor.us.us_states import STATE_CHOICES
#from registration.forms import RegistrationForm

from functions.utils.utils_df import MultiFileField


class UploadForm1(forms.Form):
    docfile1 = forms.FileField(label='Select meta.xlsx file:')


class UploadForm2(forms.Form):
    docfile3 = forms.FileField(label='Select conserved taxonomy file:')
    docfile4 = forms.FileField(label='Select .shared file:')
    sff_files = MultiFileField()
    oligo_files = MultiFileField()
    docfile5 = forms.FileField(label='Select filenames file:')
    fna_files = MultiFileField()
    qual_files = MultiFileField()
    docfile6 = forms.FileField(label='Select Oligos file:')
    docfile7 = forms.FileField(label='Select Mothur batch file:')
    docfile13 = forms.FileField(label='Select 3-column contig file:')
    fastq_files = MultiFileField()
    docfile15 = forms.FileField(label='Select Mothur batch file:')
    source = forms.ChoiceField(widget=forms.Select, choices=(('mothur', 'Pre-processed mothur files'), ('454_sff', 'sff files'), ('454_fastq', 'fna/qual files'), ('miseq', 'fastq files')))
    processors = forms.IntegerField(initial=2, min_value=1, max_value=100)


class UploadForm4(forms.Form):
    docfile8 = forms.FileField(label='Select alignment file (e.g., silva.seed_v119.align):')
    docfile9 = forms.FileField(label='Select template file (e.g., gg_13_5_99.fasta):')
    docfile10 = forms.FileField(label='Select taxonomy file (e.g., gg_13_5_99.pds.tax):')
    alignFile = forms.ChoiceField(widget=forms.Select, choices=('f1', 'File1'))
    templateFile = forms.ChoiceField(widget=forms.Select, choices=('f1', 'File1'))
    taxonomyFile = forms.ChoiceField(widget=forms.Select, choices=('f1', 'File1'))
    processors = forms.IntegerField(initial=2, min_value=1, max_value=100)


class UploadForm5(forms.Form):
    docfile11 = forms.FileField(label='Select meta.xlsx file:')


class UploadForm6(forms.Form):
    taxonomy = forms.FileField(label='Select file (e.g., gg_13_5_99.dkm.tax):')
    precalc_16S = forms.FileField(label='Select zip file (e.g., 16S_13_5_precalculated.tab.gz):')
    precalc_KEGG = forms.FileField(label='Select zip file (e.g., ko_13_5_precalculated.tab.gz):')


class UploadForm7(forms.Form):
    ko_htext = forms.FileField(label='Select file (e.g., ko00001.keg):')


class UploadForm8(forms.Form):
    nz_htext = forms.FileField(label='Select file (e.g., ko01000.keg):')


class UploadForm9(forms.Form):
    normFile = forms.FileField()

reference_choices = (
    ('No data', '(not selected)'),
    ('Publication', 'Publication'),
    ('Web Search', 'Web Search'),
    ('Colleague or Friend', 'Colleague or Friend'),
    ('Other', 'Other')
)

purpose_choices = (
    ('No data', '(not selected)'),
    ('Research', 'Research'),
    ('Teaching', 'Teaching'),
    ('Consulting', 'Consulting'),
    ('Software Development', 'Software Development'),
    ('Other', 'Other')
)

myStates = list(STATE_CHOICES)
myStates.insert(0, ('No data', '(not selected)'))

myCountries = list(countries)
myCountries.insert(0, ('No data', '(not selected)'))

'''
class UserRegForm(RegistrationForm):
    firstName = forms.CharField(label='First Name', required=False)
    lastName = forms.CharField(label='Last Name', required=False)
    affiliation = forms.CharField(label='Affiliation', required=False)
    city = forms.CharField(label='City', required=False)
    state = forms.ChoiceField(widget=forms.Select, choices=myStates)
    country = forms.ChoiceField(widget=forms.Select, choices=myCountries)
    zip = forms.CharField(label='Zip', required=False)
    phone = forms.CharField(label='Phone', required=False)
    reference = forms.ChoiceField(widget=forms.Select, choices=reference_choices)
    purpose = forms.ChoiceField(widget=forms.Select, choices=purpose_choices)
'''