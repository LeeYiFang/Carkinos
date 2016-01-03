# probes/urls.py
from django.conf.urls import url
from .views import home  # explicit relative import

urlpatterns = [
    url(r'^$', home, name="home"),
]