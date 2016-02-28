from django.db import models

class Dataset(models.Model):
    
    name=models.CharField(max_length=50)
    data_path=models.CharField(max_length=50)
    
    def __str__(self):
        return '%s of %s' % (self.name, self.data_path)
        


# Create your models here.
       
class Platform(models.Model):

    name=models.CharField(max_length=50)
    description=models.CharField(max_length=50)
    
    def __str__(self):
        return '%s of %s' % (self.name, self.description)
        
class Sample(models.Model):
    
    name=models.CharField(max_length=50)
    filename=models.CharField(max_length=50)
    cell_line_id=models.ForeignKey('CellLine', related_name='fcell_line_id')
    dataset_id=models.ForeignKey('Dataset',related_name='fdataset_id',default='')
    platform_id=models.ForeignKey('Platform', related_name='fplatform_id')
    
    def __str__(self):
        return 'name=%s filename=%s' % (self.name, self.filename)


        

class CellLine(models.Model):

    name = models.CharField(max_length=20)
    primary_site = models.CharField(max_length=50)
    primary_hist = models.CharField(max_length=50)
    
    
    def __str__(self):
        return '%s of %s of %s' % (self.name, self.primary_site, self.primary_hist)


#(still 30)
class ProbeID(models.Model):

    Probe_id = models.CharField(max_length=20)
    Gene_symbol = models.CharField(max_length=20)
    Entrez_id = models.IntegerField()
    Gene_name = models.TextField(blank=True, default='')
    platform = models.ForeignKey('Platform', related_name='fplatform',default='')
    
    def __str__(self):
        return self.Probe_id  #Error: not all arguments converted during string


