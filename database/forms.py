from django import forms
from django_countries import countries
from localflavor.us.us_states import STATE_CHOICES
from uuid import uuid4
from multiupload.fields import MultiFileField


class UploadForm1(forms.Form):
    docfile1 = forms.FileField(label='Select meta.xlsx file:')


class UploadForm2(forms.Form):
    docfile4 = forms.FileField(label='Select shared file:', required=False)
    ref_data = forms.ChoiceField(label='Create new taxonomy file:', widget=forms.Select, choices=(('yes', 'yes'), ('no', 'no')))
    ref_fasta = forms.FileField(label='Select representative sequence file (e.g., dada.rep_seqs.fasta:', required=False)
    ref_var = forms.ChoiceField(label='Select the variable region:', widget=forms.Select, choices=(('V34', 'V34'), ('V4', 'V4'), ('V13', 'V13'), ('V19', 'V19')))
    docfile3 = forms.FileField(label='Select conserved taxonomy file:', required=False)
    sff_files = MultiFileField()
    oligo_files = MultiFileField()
    docfile5 = forms.FileField(label='Select filenames file:', required=False)
    fna_files = MultiFileField()
    qual_files = MultiFileField()
    docfile6 = forms.FileField(label='Select Oligos file:', required=False)
    docfile7 = forms.FileField(label='Select batch file:', required=False)
    docfile13 = forms.FileField(label='Select 3-column contig file:', required=False)
    platform = forms.ChoiceField(widget=forms.Select, choices=(('mothur', 'mothur'), ('dada2', 'dada2')), initial='dada2')
    fastq_files = MultiFileField()

    # biom upload files
    biom_file = forms.FileField(label='Select .biom file:', required=False)
    tsv_file = forms.FileField(label='Select .tsv file:', required=False)
    fasta_file = forms.FileField(label='Select .fasta file:', required=False)

    # make HTML version for file upload page
    source = forms.ChoiceField(widget=forms.Select, choices=(('mothur', 'Mothur: shared'), ('biom', 'Qiime: biom'), ('454_sff', 'sff files'), ('454_fastq', 'fna/qual files'), ('miseq', 'fastq files')))
    processors = forms.IntegerField(initial=2, min_value=1, max_value=100)


class UploadForm4(forms.Form):
    docfile8 = forms.FileField(label='Select alignment file (e.g., silva.seed_v119.align):', required=False)
    docfile9 = forms.FileField(label='Select template file (e.g., gg_13_5_99.fasta):', required=False)
    docfile10 = forms.FileField(label='Select taxonomy file (e.g., gg_13_5_99.pds.tax):', required=False)
    alignFile = forms.ChoiceField(widget=forms.Select, choices=('f1', 'File1'))
    templateFile = forms.ChoiceField(widget=forms.Select, choices=('f1', 'File1'))
    taxonomyFile = forms.ChoiceField(widget=forms.Select, choices=('f1', 'File1'))
    processors = forms.IntegerField(initial=2, min_value=1, max_value=100)


class UploadForm5(forms.Form):
    docfile11 = forms.FileField(label='Select meta.xlsx file:')


class UploadForm6(forms.Form):
    taxonomy = forms.FileField(label='Select file (e.g., gg_13_5_99.dkm.tax):', required=False)
    precalc_16S = forms.FileField(label='Select zip file (e.g., 16S_13_5_precalculated.tab.gz):', required=False)
    precalc_KEGG = forms.FileField(label='Select zip file (e.g., ko_13_5_precalculated.tab.gz):', required=False)


class UploadForm7(forms.Form):
    ko_htext = forms.FileField(label='Select file (e.g., ko00001.keg):', required=False)


class UploadForm8(forms.Form):
    nz_htext = forms.FileField(label='Select file (e.g., ko01000.keg):', required=False)


class UploadForm9(forms.Form):
    normFile = forms.FileField()


class UploadForm10(forms.Form):
    platform = forms.ChoiceField(widget=forms.Select, choices=(('mothur', 'mothur'), ('dada2', 'dada2')))
    mothurFile = forms.FileField(label='Select batch file:')


class BulkForm(forms.Form):
    metafiles = MultiFileField(required=False)
    sharedfiles = MultiFileField(required=False)
    taxafiles = MultiFileField(required=False)
    sequencefiles = MultiFileField(required=False)
    scriptfiles = MultiFileField(required=False)

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


class UserRegForm(forms.Form):
    firstName = forms.CharField(label='First Name', required=True)
    lastName = forms.CharField(label='Last Name', required=True)
    email = forms.EmailField(label='Email', required=True)
    affiliation = forms.CharField(label='Affiliation', required=False)
    city = forms.CharField(label='City', required=False)
    state = forms.ChoiceField(widget=forms.Select, choices=myStates)
    country = forms.ChoiceField(widget=forms.Select, choices=myCountries)
    zip = forms.CharField(label='Zip', required=False)
    phone = forms.CharField(label='Phone', required=False)
    reference = forms.ChoiceField(widget=forms.Select, choices=reference_choices)
    purpose = forms.ChoiceField(widget=forms.Select, choices=purpose_choices)

    def signup(self, request, user):
        form = UserRegForm(request.POST or None, request.FILES or None)
        if form.is_valid():
            user.first_name = self.cleaned_data['firstName']
            user.last_name = self.cleaned_data['lastName']
            user.email = self.cleaned_data['email']
            user.save()

            up = user.profile
            up.firstName = self.cleaned_data['firstName']
            up.lastName = self.cleaned_data['lastName']
            up.affiliation = self.cleaned_data['affiliation']
            up.city = self.cleaned_data['city']
            up.state = self.cleaned_data['state']
            up.country = self.cleaned_data['country']
            up.zip = self.cleaned_data['zip']
            up.phone = self.cleaned_data['phone']
            up.reference = self.cleaned_data['reference']
            up.purpose = self.cleaned_data['purpose']
            up.dataID = uuid4().hex
            up.save()


class UserUpdateForm(forms.Form):
    firstName = forms.CharField(label='First Name', required=True)
    lastName = forms.CharField(label='Last Name', required=True)
    email = forms.EmailField(label='Email', required=True)
    affiliation = forms.CharField(label='Affiliation', required=False)
    city = forms.CharField(label='City', required=False)
    state = forms.ChoiceField(widget=forms.Select, choices=myStates)
    country = forms.ChoiceField(widget=forms.Select, choices=myCountries)
    zip = forms.CharField(label='Zip', required=False)
    phone = forms.CharField(label='Phone', required=False)
    reference = forms.ChoiceField(widget=forms.Select, choices=reference_choices)
    purpose = forms.ChoiceField(widget=forms.Select, choices=purpose_choices)
    pword = forms.CharField(label="Password", widget=forms.PasswordInput)

    def update(self, request, user):
        form = UserRegForm(request.POST or None, request.FILES or None)
        if form.is_valid():
            user.first_name = self.cleaned_data['firstName']
            user.last_name = self.cleaned_data['lastName']
            user.email = self.cleaned_data['email']
            user.save()

            up = user.profile
            up.firstName = self.cleaned_data['firstName']
            up.lastName = self.cleaned_data['lastName']
            up.affiliation = self.cleaned_data['affiliation']
            up.city = self.cleaned_data['city']
            up.state = self.cleaned_data['state']
            up.country = self.cleaned_data['country']
            up.zip = self.cleaned_data['zip']
            up.phone = self.cleaned_data['phone']
            up.reference = self.cleaned_data['reference']
            up.purpose = self.cleaned_data['purpose']
            up.dataID = uuid4().hex
            up.save()
