from django import forms


class UploadForm1(forms.Form):
    docfile1 = forms.FileField(label='Select meta_Project.csv file:')
    docfile2 = forms.FileField(label='Select meta_Sample.csv file:')


class UploadForm2(forms.Form):
    docfile3 = forms.FileField(label='Select conserved taxonomy file:')
    docfile4 = forms.FileField(label='Select .shared file:')

'''
class UploadForm3(forms.Form):
    docfile5 = forms.FileField(label='Select QIIME BIOM file:')


class UploadForm4(forms.Form):
    docfile6 = forms.FileField(label='Select QIIME OTU Table file:')


class UploadForm5(forms.Form):
    docfile7 = forms.FileField(label='Select MG-RAST table file:')
'''

