from django.db import models

# Create your models here.
#U133A ,U133B(in the future) ...platform        
'''class Platform(models.Model):

    level = models.CharField(max_length=20)  #RNA
    
    
    def __str__(self):
        return self.name'''
        
        
#cellline_short(only 30)
class CellLine(models.Model):

    name = models.CharField(max_length=20)
    primary_site = models.CharField(max_length=50)
    primary_hist = models.CharField(max_length=50)
    #platform = models.ForeignKey('Platform', related_name='plat_name')
    
    def __str__(self):
        return self.name


#Affy_U133A_probe_info_short(only 30)
class ProbeID(models.Model):

    Probe_id = models.CharField(max_length=20)
    Gene_symbol = models.CharField(max_length=20)
    Entrez_id = models.IntegerField()
    Gene_name = models.TextField(blank=True, default='')
    
    def __str__(self):
        return self.Probe_id



