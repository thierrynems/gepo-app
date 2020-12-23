from django import forms

class PredictionForm(forms.Form):
    modelSelect = forms.CharField(max_length=20,help_text="Select the model to use",required=True)
    actionSelect = forms.CharField(max_length=20,help_text="Select the action",required=True)
    geneSequence = forms.CharField(widget=forms.Textarea,help_text="insert the gene sequence",required=True)
    proteinSequence = forms.CharField(widget=forms.Textarea,help_text="insert the protein sequence",required=True)