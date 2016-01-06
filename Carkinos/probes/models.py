from django.db import models

# Create your models here.

#Affy_U133A_probe_info
class ProbeID(models.Model):

    Probe_id = models.CharField(max_length=20)
    Gene_symbol = models.CharField(max_length=20)
    Entrez_id = models.IntegerField()
    Gene_name = models.TextField(blank=True, default='')
    
    def __str__(self):
        return self.Probe_id,self.Gene_symbol


'''class MenuItem(models.Model):

    store = models.ForeignKey('Store', related_name='menu_items')
    name = models.CharField(max_length=20)
    price = models.IntegerField()

    def __str__(self):
        return self.name'''