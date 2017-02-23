from pathlib import Path
import pandas as pd
import numpy as np

import django
import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'Carkinos.settings.local'
django.setup()

from probes.models import Dataset,Platform,Sample,CellLine,ProbeID
ProbeID.objects.all().delete()  #notice this!!!!!
root=Path('../').resolve()
probe_path=root.joinpath('src','Affy_U133A_probe_info.txt')
plat_info_df = pd.read_csv(probe_path.as_posix(),sep='\t')

i=-1
u133a=Platform.objects.get(id=1)#make sure again
plus2=Platform.objects.get(id=3) 
u133b=Platform.objects.get(id=2)

plat_info_df.SYMBOL.fillna('', inplace=True)
plat_info_df.GENENAME.fillna('', inplace=True)
plat_info_df.ENTREZID.fillna(0, inplace=True)

probes={}
for rowid, row in plat_info_df.iterrows():
    
    try:
        temp=probes[str(row.PROBEID)]
    except KeyError:
        i=i+1
        temp=i
        probes[str(row.PROBEID)]=i
    sample = ProbeID.objects.create(
    Probe_id =  str(row.PROBEID),
    Gene_symbol = str(row.SYMBOL),
    Entrez_id = row.ENTREZID,
    Gene_name = str(row.GENENAME),
    platform=u133a,
    offset=temp
    )
        
print(ProbeID.objects.all().count())       
    
probe_path=root.joinpath('src','Affy_U133plus2_probe_info.txt')
plat_info_df = pd.read_csv(probe_path.as_posix(),sep='\t')
plat_info_df.SYMBOL.fillna('', inplace=True)
plat_info_df.GENENAME.fillna('', inplace=True)
plat_info_df.ENTREZID.fillna(0, inplace=True)
i=-1
probes={}
for rowid, row in plat_info_df.iterrows():
    
    try:
        temp=probes[str(row.PROBEID)]
    except KeyError:
        i=i+1
        temp=i
        probes[str(row.PROBEID)]=i
    sample = ProbeID.objects.create(
    Probe_id =  str(row.PROBEID),
    Gene_symbol = str(row.SYMBOL),
    Entrez_id = row.ENTREZID,
    Gene_name = str(row.GENENAME),
    platform=plus2,
    offset=temp
    )
print(ProbeID.objects.all().count())     

from probes.models import Gene
Gene.objects.all().delete()
p2=pd.read_csv('uni_plus2.txt')
g=list(p2['SYMBOL'])
counter=0
plus2=Platform.objects.get(id=3)
for i in g:
    sample=Gene.objects.create(
    symbol=str(i),
    platform=plus2,
    offset=counter
    )
    counter=counter+1

print(Gene.objects.all().count())
    
p2=pd.read_csv('uni_u133a.txt')
g=list(p2['SYMBOL'])
counter=0
plus2=Platform.objects.get(id=1)
for i in g:
    sample=Gene.objects.create(
    symbol=str(i),
    platform=plus2,
    offset=counter
    )
    counter=counter+1
    
print(Gene.objects.all().count())