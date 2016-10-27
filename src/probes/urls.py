# probes/urls.py
from django.conf.urls import url
from .views import cellline_microarray  # explicit relative import
from .views import similar_assessment
from .views import welcome,pca,help_similar_assessment,gene_signature,heatmap,express_profiling,user_pca,sample_microarray,clinical_search
urlpatterns = [  
    url(r'^cellline_microarray$', cellline_microarray, name="cellline_microarray"),
    url(r'^similar_assessment$', similar_assessment, name="similar_assessment"),
    url(r'^$',welcome,name="welcome"),
    url(r'^pca$',pca,name="pca"),
    url(r'^heatmap$',heatmap,name="heatmap"),
    url(r'^gene_signature$',gene_signature,name="gene_signature"),
    url(r'^help_similar_assessment$',help_similar_assessment,name="help_similar_assessment"),
    url(r'^express_profiling$',express_profiling,name="express_profiling"),
    url(r'^user_pca$',user_pca,name="user_pca"),
    url(r'^clinical_search$',clinical_search,name="clinical_search"),
    url(r'^sample_microarray$', sample_microarray, name="sample_microarray"),
    
]