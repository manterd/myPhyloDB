from django import forms
from utils import MultiFileField


class UploadForm1(forms.Form):
    docfile1 = forms.FileField(label='Select meta.csv file:')

class UploadForm2(forms.Form):
    docfile3 = forms.FileField(label='Select conserved taxonomy file:')
    docfile4 = forms.FileField(label='Select .shared file:')
    docfile5 = forms.FileField(label='Select sff file:')
    docfile6 = forms.FileField(label='Select Oligos file:')
    docfile7 = forms.FileField(label='Select Mothur batch file:')
    docfile13 = forms.FileField(label='Select 3-column contig file:')
    files = MultiFileField()
    docfile15 = forms.FileField(label='Select Mothur batch file:')
    source = forms.ChoiceField(widget=forms.Select, choices=(('mothur', 'Pre-processed Mothur Files'), ('454', 'Raw 454 Files'), ('miseq', 'Raw Illumina/MiSeq Files')))


class UploadForm4(forms.Form):
    docfile8 = forms.FileField(label='Select alignment file (e.g., silva.seed_v119.align):')
    docfile9 = forms.FileField(label='Select template file (e.g., gg_13_5_99.fasta):')
    docfile10 = forms.FileField(label='Select taxonomy file (e.g., gg_13_5_99.pds.tax):')
    alignFile = forms.ChoiceField(widget=forms.Select, choices=('f1', 'File1'))
    templateFile = forms.ChoiceField(widget=forms.Select, choices=('f1', 'File1'))
    taxonomyFile = forms.ChoiceField(widget=forms.Select, choices=('f1', 'File1'))


class UploadForm5(forms.Form):
    docfile11 = forms.FileField(label='Select meta.csv file:')

