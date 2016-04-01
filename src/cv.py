from pathlib import Path
import pandas as pd
import numpy as np

import django
import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'Carkinos.settings.local'
django.setup()

from probes.models import Dataset,Platform,Sample,CellLine,ProbeID
root=Path('../').resolve()
u133a_path=root.joinpath('src','raw','Affy_U133A_probe_info.csv')
plus2_path=root.joinpath('src','raw','Affy_U133plus2_probe_info.csv')

u133a=pd.read_csv(u133a_path.as_posix())
plus2=pd.read_csv(plus2_path.as_posix())

sanger_val_pth=Path('../').resolve().joinpath('src','sanger_cell_line_proj.npy')
nci_val_pth=Path('../').resolve().joinpath('src','nci60.npy')
gse_val_pth=Path('../').resolve().joinpath('src','GSE36133.npy')
sanger_val=np.load(sanger_val_pth.as_posix(),mmap_mode='r')
nci_val=np.load(nci_val_pth.as_posix(),mmap_mode='r')
gse_val=np.load(gse_val_pth.as_posix(),mmap_mode='r')
plus2.SYMBOL.fillna('', inplace=True)
u133a.SYMBOL.fillna('', inplace=True)

#this is now for all platform that can find
#ugene=list(set(list(pd.unique(plus2.SYMBOL))+list(pd.unique(u133a.SYMBOL))))
ugene=list(pd.unique(u133a.SYMBOL))
ugene.remove('')
#ugene has all the gene symbols in U133A and U133PlUS2

#sanger=798,nci60=174,gse=917
sanger_offset=Sample.objects.filter(dataset_id=1).values_list('offset',flat=True)
nci60_offset=Sample.objects.filter(dataset_id__name__in=['NCI60']).values_list('offset',flat=True)
gse_offset=Sample.objects.filter(dataset_id__name__in=['GSE36133']).values_list('offset',flat=True)

min=10000000000
min_gene=[]
for gene in ugene:
    aprobe=ProbeID.objects.filter(platform=1,Gene_symbol=gene)
    pprobe=ProbeID.objects.filter(platform=3,Gene_symbol=gene)    
    aoffset=aprobe.values_list('offset',flat=True)
    aprobe_length=len(aoffset)
    poffset=pprobe.values_list('offset',flat=True)
    pprobe_length=len(poffset)    
    sanger_sample=sanger_val[np.ix_(aoffset,sanger_offset)]
    sanger_sum=np.sum(sanger_sample)    
    nci_sample=nci_val[np.ix_(poffset,nci60_offset)]
    nci_sum=np.sum(nci_sample)
    gse_sample=gse_val[np.ix_(poffset,gse_offset)]
    gse_sum=np.sum(gse_sample)
    mean=(sanger_sum+nci_sum+gse_sum)/(aprobe_length*798+pprobe_length*1091)
    sanger_square=np.sum(np.square(np.subtract(sanger_sample,mean)))
    nci_square=np.sum(np.square(np.subtract(nci_sample,mean)))
    gse_square=np.sum(np.square(np.subtract(gse_sample,mean)))
    std=((sanger_square+nci_square+gse_square)/(aprobe_length*798+pprobe_length*1091))**0.5
    temp=std/mean
    if min>temp:
        min=temp
        min_gene=gene
    
print(min)
print(min_gene)    
    
    
    
    



