from pathlib import Path
import pandas as pd
import numpy as np

import django
import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'Carkinos.settings.local'
django.setup()

from probes.models import Dataset,Platform,Sample,CellLine,ProbeID
root=Path('../').resolve()
sanger_path=root.joinpath('src','raw','Sanger','Sanger.xls')
nci60_path=root.joinpath('src','raw','nci60','nci60.xlsx')
gse_path=root.joinpath('src','raw','GSE36133','CCLE.xlsx')

sanger_sample = pd.read_excel(sanger_path.as_posix())
nci60_sample = pd.read_excel(nci60_path.as_posix())
gse_sample = pd.read_excel(gse_path.as_posix())

known_cell_lines = {}
for rowid, row in sanger_sample.iterrows():
    cl = str(row.cell_line)
    p_site=str(row.primary_site)
    p_hist=str(row.primary_hist)
    try:
        cell_line = known_cell_lines[cl+"/"+p_site+"/"+p_hist]
    except KeyError:
        cell_line = CellLine.objects.create(
            name=str(row.cell_line), 
            primary_site=row.primary_site, 
            primary_hist=row.primary_hist
        )
        known_cell_lines[cl+"/"+p_site+"/"+p_hist] = cell_line
        

for rowid, row in nci60_sample.iterrows():
    cl = str(row.cell_line)
    p_site=str(row.primary_site)
    p_hist=str(row.primary_hist)
    try:
        cell_line = known_cell_lines[cl+"/"+p_site+"/"+p_hist]
    except KeyError:
        cell_line = CellLine.objects.create(
            name=str(row.cell_line), 
            primary_site=row.primary_site, 
            primary_hist=row.primary_hist
        )
        known_cell_lines[cl+"/"+p_site+"/"+p_hist] = cell_line
        

for rowid, row in gse_sample.iterrows():
    cl = str(row.cell_line)
    p_site=str(row.primary_site)
    p_hist=str(row.primary_hist)
    try:
        cell_line = known_cell_lines[cl+"/"+p_site+"/"+p_hist]
    except KeyError:
        cell_line = CellLine.objects.create(
            name=str(row.cell_line), 
            primary_site=row.primary_site, 
            primary_hist=row.primary_hist
        )
        known_cell_lines[cl+"/"+p_site+"/"+p_hist] = cell_line
        
u133a=Platform.objects.get(id=1)#make sure again
plus2=Platform.objects.get(id=3) 
sanger=Dataset.objects.get(id=1)
nci60=Dataset.objects.get(id=2)
gse=Dataset.objects.get(id=3)      

i=0        
for rowid, row in sanger_sample.iterrows():
    name = row['name']
    # Create sample
    sample = Sample.objects.create(
        name=name,
        filename="sanger_cell_line_proj.npy",
        cell_line_id=known_cell_lines[str(row.cell_line)+"/"+str(row.primary_site)+"/"+str(row.primary_hist)],
        dataset_id=sanger,
        platform_id=u133a,
        offset=i
    )
    i=i+1
    
i=0        
for rowid, row in nci60_sample.iterrows():
    name = row['name']
    # Create sample
    sample = Sample.objects.create(
        name=name,
        filename="nci60.npy",
        cell_line_id=known_cell_lines[str(row.cell_line)+"/"+str(row.primary_site)+"/"+str(row.primary_hist)],
        dataset_id=nci60,
        platform_id=plus2,
        offset=i
    )
    i=i+1
    
i=0        
for rowid, row in gse_sample.iterrows():
    name = row['name']
    # Create sample
    sample = Sample.objects.create(
        name=name,
        filename="GSE36133.npy",
        cell_line_id=known_cell_lines[str(row.cell_line)+"/"+str(row.primary_site)+"/"+str(row.primary_hist)],
        dataset_id=gse,
        platform_id=plus2,
        offset=i
    )
    i=i+1