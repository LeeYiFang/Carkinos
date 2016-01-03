from django.db import models

# Create your models here.

class ProbeID_GeneSymbol(models.Model):

    ProbeId = models.CharField(max_length=20)
    GeneSymbol = models.CharField(max_length=20)

    #def __str__(self):
    #    return self.ProbeId


'''class MenuItem(models.Model):

    store = models.ForeignKey('Store', related_name='menu_items')
    name = models.CharField(max_length=20)
    price = models.IntegerField()

    def __str__(self):
        return self.name'''