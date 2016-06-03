# probes/urls.py
from django.conf.urls import url
from .views import cellline_microarray  # explicit relative import
from .views import similar_assessment
from .views import welcome,pca,help_similar_assessment
urlpatterns = [  
    url(r'^cellline_microarray$', cellline_microarray, name="cellline_microarray"),
    url(r'^similar_assessment$', similar_assessment, name="similar_assessment"),
    url(r'^$',welcome,name="welcome"),
    url(r'^pca$',pca,name="pca"),
    url(r'^help_similar_assessment$',help_similar_assessment,name="help_similar_assessment"),
]