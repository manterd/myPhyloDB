from django import forms



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

