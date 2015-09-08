from django import forms
from utils import MultiFileField


class UploadForm1(forms.Form):
    docfile1 = forms.FileField(label='Select meta_Project.csv file:')
    docfile2 = forms.FileField(label='Select meta_Sample.csv file:')
    type = forms.ChoiceField(widget=forms.Select, choices=(('soil', 'Soil'), ('human_associated', 'Human Associated'), ('human_gut', 'Human Gut'), ('air', 'Air'), ('water', 'Water'), ('microbial', 'Microbial')))


class UploadForm2(forms.Form):
    docfile3 = forms.FileField(label='Select conserved taxonomy file:')
    docfile4 = forms.FileField(label='Select .shared file:')


class UploadForm3(forms.Form):
    docfile5 = forms.FileField(label='Select sff file:')
    docfile6 = forms.FileField(label='Select Oligos file:')
    docfile7 = forms.FileField(label='Select Mothur batch file:')


class UploadForm6(forms.Form):
    docfile13 = forms.FileField(label='Select 3-column contig file:')
    files = MultiFileField()
    docfile15 = forms.FileField(label='Select Mothur batch file:')


class UploadForm4(forms.Form):
    docfile8 = forms.FileField(label='Select alignment file (e.g., silva.seed_v119.align):')
    docfile9 = forms.FileField(label='Select template file (e.g., gg_13_5_99.fasta):')
    docfile10 = forms.FileField(label='Select taxonomy file (e.g., gg_13_5_99.pds.tax):')
    alignFile = forms.ChoiceField(widget=forms.Select, choices=('f1', 'File1'))
    templateFile = forms.ChoiceField(widget=forms.Select, choices=('f1', 'File1'))
    taxonomyFile = forms.ChoiceField(widget=forms.Select, choices=('f1', 'File1'))


class UploadForm5(forms.Form):
    docfile11 = forms.FileField(label='Select meta_Project.csv file:')
    docfile12 = forms.FileField(label='Select meta_Sample.csv file:')
