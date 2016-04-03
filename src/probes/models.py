import itertools
from django.db import models
from decimal import Decimal

#add:new dataset
class Dataset(models.Model):
    
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


class CellLineQuerySet(models.QuerySet):

    def annotate_dataset_name(self):
        '''Return cell line objects with dataset names

        Note:
            Cell lines can be duplicated if found in multiple datasets
        '''
        return self.raw('''
            SELECT DISTINCT
                c.id AS id,
                c.name AS name,
                c.primary_site AS primary_site,
                c.primary_hist AS primary_hist,
                sample_extended.dset_name AS dset_name
            FROM (
                -- Retrieve all dataset names for each sample
                SELECT
                    s.cell_line_id_id AS cl_id,
                    ds.id AS dset_id,
                    ds.name AS dset_name
                FROM probes_sample AS s
                JOIN probes_dataset AS ds
                ON s.dataset_id_id = ds.id
            ) AS sample_extended
            JOIN probes_cellline AS c
            ON sample_extended.cl_id = c.id
            ORDER BY id, dset_name
        ''')

    def collapsed(self):
        '''Collapse cell line dataset annotation.

        Note:
            Records for cell line will be distinct, and dataset names are
            collapsed into the attribute `dataset_names`
        '''
        for cell_line_id, grp in itertools.groupby(
            self.annotate_dataset_name(), lambda cl: cl.id
        ):
            records = list(grp)
            cell_line_obj = records[0]
            cell_line_obj.dataset_names = [r.dset_name for r in records]
            yield cell_line_obj


#need to update with new datasets
class CellLine(models.Model):

    name = models.CharField(max_length=20)
    primary_site = models.CharField(max_length=50)
    primary_hist = models.CharField(max_length=50)

    by_datasets = CellLineQuerySet.as_manager()

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
        




