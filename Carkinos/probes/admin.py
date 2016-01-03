from django.contrib import admin

# Register your models here.

from .models import ProbeID_GeneSymbol

@admin.register(ProbeID_GeneSymbol)
class ProbeID_GeneSymbolAdmin(admin.ModelAdmin):
    list_display = ('ProbeId', 'GeneSymbol',)