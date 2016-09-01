# probes/urls.py
from django.conf.urls import url
from .views import cellline_microarray  # explicit relative import
from .views import similar_assessment
from .views import welcome,pca,help_similar_assessment,gene_signature,heatmap,upload
urlpatterns = [  
    url(r'^cellline_microarray$', cellline_microarray, name="cellline_microarray"),
    url(r'^similar_assessment$', similar_assessment, name="similar_assessment"),
    url(r'^$',welcome,name="welcome"),
    url(r'^pca$',pca,name="pca"),
    url(r'^heatmap$',heatmap,name="heatmap"),
    url(r'^gene_signature$',gene_signature,name="gene_signature"),
    url(r'^help_similar_assessment$',help_similar_assessment,name="help_similar_assessment"),
    url(r'^upload$',upload,name="upload"),
    
]