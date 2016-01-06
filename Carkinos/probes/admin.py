from django.contrib import admin

# Register your models here.

from .models import ProbeID

@admin.register(ProbeID)
class ProbeIDAdmin(admin.ModelAdmin):
    list_display = ('Probe_id', 'Gene_symbol','Entrez_id','Gene_name')