from django.db import models
from decimal import Decimal

#add:new dataset
class Dataset(models.Model):
    
    name=models.CharField(max_length=50)
    data_path=models.CharField(max_length=50)
    
    def __str__(self):
        return '%s of %s' % (self.name, self.data_path)

class Clinical_Dataset(models.Model):
    
    name=models.CharField(max_length=50)
    data_path=models.CharField(max_length=50)
    
    def __str__(self):
        return '%s of %s' % (self.name, self.data_path)
        

     
class Platform(models.Model):

    name=models.CharField(max_length=50)
    description=models.CharField(max_length=50)
    
    def __str__(self):
        return '%s of %s' % (self.name, self.description)
        
        
      
#need to update with new datasets  
class Sample(models.Model):
    
    name=models.CharField(max_length=50)
    filename=models.CharField(max_length=50)
    cell_line_id=models.ForeignKey('CellLine', related_name='fcell_line_id')
    dataset_id=models.ForeignKey('Dataset',related_name='fdataset_id',default='')
    platform_id=models.ForeignKey('Platform', related_name='fplatform_id')
    offset=models.IntegerField(default=0)
    def __str__(self):
        return 'name=%s filename=%s' % (self.name, self.filename)

class Clinical_sample(models.Model):
    
    name=models.CharField(max_length=20)
    filename=models.CharField(max_length=50)
    dataset_id=models.ForeignKey('Clinical_Dataset',related_name='clinic_dataset_id',default='')
    platform_id=models.ForeignKey('Platform', related_name='clinc_platform_id')
    offset=models.IntegerField(default=0)
    age=models.CharField(max_length=10,default='UNKNOWN')
    gender=models.CharField(max_length=10,default='UNKNOWN')
    ethnic=models.CharField(max_length=20,default='UNKNOWN')
    grade=models.CharField(max_length=10)
    stage=models.CharField(max_length=10)
    primary_site= models.CharField(max_length=50)
    primary_hist= models.CharField(max_length=100)
    stageT=models.CharField(max_length=10,default='UNKNOWN')
    stageN=models.CharField(max_length=10,default='UNKNOWN')
    stageM=models.CharField(max_length=10,default='UNKNOWN')
    metastatic=models.BooleanField(default=False)
    
    def __str__(self):
        return '%s of %s of %s of %s' % (self.name, self.primary_site, self.primary_hist,self.dataset_id)

        
#need to update with new datasets
class CellLine(models.Model):

    name = models.CharField(max_length=20)
    primary_site = models.CharField(max_length=50)
    primary_hist = models.CharField(max_length=100)
    
    
    def __str__(self):
        return '%s of %s of %s' % (self.name, self.primary_site, self.primary_hist)



#alread put in all three platform data
class ProbeID(models.Model):

    Probe_id = models.CharField(max_length=20)
    Gene_symbol = models.CharField(max_length=20)
    Entrez_id = models.IntegerField()
    Gene_name = models.TextField(blank=True, default='')
    platform = models.ForeignKey('Platform', related_name='fplatform',default='')
    offset=models.IntegerField(default=0)
    def __str__(self):
        return self.Probe_id
        




