from django.db import models

# Create your models here.
#U133A ,U133B(in the future) ...platform        
class Platform(models.Model):

    name=models.CharField(max_length=50)
    description=models.CharField(max_length=50)
    
    def __str__(self):
        return '%s of %s' % (self.name, self.description)
        
class Sample(models.Model):
    
    name=models.CharField(max_length=50)
    filename=models.CharField(max_length=50)
    cell_line_id=models.ForeignKey('CellLine', related_name='fcell_line_id')
    #dataset_id
    platform_id=models.ForeignKey('Platform', related_name='fplatform_id')
    
    def __str__(self):
        return '%s of %s' % (self.name, self.filename, self.cell_line_id.name,self.platform_id.name)


        
#Sanger project's cellline_short(only 30)
class CellLine(models.Model):

    name = models.CharField(max_length=20)
    primary_site = models.CharField(max_length=50)
    primary_hist = models.CharField(max_length=50)
    
    
    def __str__(self):
        return '%s of %s' % (self.name, self.primary_site, self.primary_hist)


#Affy_U133A_probe_info_short(only 30)
class ProbeID(models.Model):

    Probe_id = models.CharField(max_length=20)
    Gene_symbol = models.CharField(max_length=20)
    Entrez_id = models.IntegerField()
    Gene_name = models.TextField(blank=True, default='')
    platform = models.ForeignKey('Platform', related_name='fplatform',default='')
    
    def __str__(self):
        return self.Probe_id  #Error: not all arguments converted during string



