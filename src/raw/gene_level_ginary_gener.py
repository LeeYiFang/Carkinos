from pathlib import Path
import pandas as pd
import numpy as np
import django
django.setup()
sanger_val_pth=Path('../').resolve().joinpath('src','sanger_cell_line_proj.npy')
nci60_path=Path('../').resolve().joinpath('src','nci60.npy')
gse_path=Path('../').resolve().joinpath('src','GSE36133.npy')



from probes.models import Sample,ProbeID,Platform
sanger_val=np.load(sanger_val_pth.as_posix(),mmap_mode='r')
root=Path('../').resolve()
probe_path=root.joinpath('src','uni_u133a.txt')
plat_info_df = pd.read_csv(probe_path.as_posix())
gene=plat_info_df['SYMBOL']
gene=list(gene)


samples=Sample.objects.filter(dataset_id__name__in=['Sanger Cell Line Project']).order_by('offset')        
offset=list(samples.values_list('offset',flat=True))
ps_id=str(Platform.objects.filter(name__in=["U133A"])[0].id)
t=[]
for norm_name in gene:    
    sanger_g=ProbeID.objects.filter(platform__in=ps_id).filter(Gene_symbol__in=[norm_name]) 
    sanger_probe_offset=list(sanger_g.values_list('offset',flat=True))
    temp=sanger_val[np.ix_(sanger_probe_offset,offset)]
    norm=np.mean(temp,axis=0, dtype=np.float64,keepdims=True)
    a=list(norm)
    t.append(a)

npy_path=root.joinpath('src','gene_sanger_cell_line_proj.npy')
with (npy_path).open('wb') as f:
    np.save(f, np.array(t).astype(np.float32))

    
    

probe_path=root.joinpath('src','uni_plus2.txt')
plat_info_df = pd.read_csv(probe_path.as_posix())
gene=plat_info_df['SYMBOL']
ps_id=str(Platform.objects.filter(name__in=["PLUS2"])[0].id)
    
nci_val=np.load(nci60_path.as_posix(),mmap_mode='r')
gse_val=np.load(gse_path.as_posix(),mmap_mode='r')

offset=[]
for i in range(0,174):
    offset.append(i)


t=[]
for norm_name in gene:    
    sanger_g=ProbeID.objects.filter(platform__in=ps_id).filter(Gene_symbol__in=[norm_name]) 
    sanger_probe_offset=list(sanger_g.values_list('offset',flat=True))
    temp=nci_val[np.ix_(sanger_probe_offset,offset)]
    norm=np.mean(temp,axis=0, dtype=np.float64,keepdims=True)
    a=list(norm)
    t.append(a)

npy_path=root.joinpath('src','gene_nci60.npy')
with (npy_path).open('wb') as f:
    np.save(f, np.array(t).astype(np.float32))

offset=[]
for i in range(0,917):
    offset.append(i)    

t=[]
for norm_name in gene:    
    sanger_g=ProbeID.objects.filter(platform__in=ps_id).filter(Gene_symbol__in=[norm_name]) 
    sanger_probe_offset=list(sanger_g.values_list('offset',flat=True))
    temp=gse_val[np.ix_(sanger_probe_offset,offset)]
    norm=np.mean(temp,axis=0, dtype=np.float64,keepdims=True)
    a=list(norm)
    t.append(a)

npy_path=root.joinpath('src','gene_GSE36133.npy')
with (npy_path).open('wb') as f:
    np.save(f, np.array(t).astype(np.float32))
    
expO_path=Path('../').resolve().joinpath('src','expO.npy')
roth_path=Path('../').resolve().joinpath('src','roth.npy')
expO_val=np.load(expO_path.as_posix(),mmap_mode='r')

offset=[]
for i in range(0,2152):
    offset.append(i)


t=[]
for norm_name in gene:    
    sanger_g=ProbeID.objects.filter(platform__in=ps_id).filter(Gene_symbol__in=[norm_name]) 
    sanger_probe_offset=list(sanger_g.values_list('offset',flat=True))
    temp=expO_val[np.ix_(sanger_probe_offset,offset)]
    norm=np.mean(temp,axis=0, dtype=np.float64,keepdims=True)
    a=list(norm)
    t.append(a)

npy_path=root.joinpath('src','gene_expO.npy')
with (npy_path).open('wb') as f:
    np.save(f, np.array(t).astype(np.float32))



roth_val=np.load(roth_path.as_posix(),mmap_mode='r')
offset=[]
for i in range(0,353):
    offset.append(i)


t=[]
for norm_name in gene:    
    sanger_g=ProbeID.objects.filter(platform__in=ps_id).filter(Gene_symbol__in=[norm_name]) 
    sanger_probe_offset=list(sanger_g.values_list('offset',flat=True))
    temp=roth_val[np.ix_(sanger_probe_offset,offset)]
    norm=np.mean(temp,axis=0, dtype=np.float64,keepdims=True)
    a=list(norm)
    t.append(a)

npy_path=root.joinpath('src','gene_roth.npy')
with (npy_path).open('wb') as f:
    np.save(f, np.array(t).astype(np.float32))
