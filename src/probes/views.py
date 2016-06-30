from django.shortcuts import render,render_to_response
from django.http import HttpResponse, Http404
from django.views.decorators.http import require_GET
from .models import Dataset, CellLine, ProbeID, Sample, Platform
from django.template import RequestContext
from django.utils.html import mark_safe
import json
import pandas as pd
import numpy as np
from pathlib import Path
import sklearn
from sklearn.decomposition import PCA



def welcome(request):
    return render_to_response('welcome.html',locals())

def help_similar_assessment(request):
    return render_to_response('help_similar_assessment.html',RequestContext(request))

def similar_assessment(request):    
    return render_to_response('similar_assessment.html',RequestContext(request))

def pca(request):
    
    if 'dataset' in request.POST:
        datas=request.POST.getlist('dataset')
    else:
        return render_to_response('help_similar_assessment.html',RequestContext(request))
    if 'cellline' in request.POST:
        ncell = request.POST['cellline']
        ncell = ncell.split()
    else:
        return render_to_response('help_similar_assessment.html',RequestContext(request))
    
    cell=CellLine.objects.filter(name__in=ncell)
    propotion=0
    sanger_cellline=[]
    output_cell=[]
    context={}
    colorX=[]
    colorY=[]
    colorZ=[]
    X=[]
    Y=[]
    Z=[]
    selected_name=[]
    all_name=[]
    total_offset=[]
    sanger_val_pth=Path('../').resolve().joinpath('src','sanger_cell_line_proj.npy')
    nci_val_pth=Path('../').resolve().joinpath('src','nci60.npy')
    gse_val_pth=Path('../').resolve().joinpath('src','GSE36133.npy')
    sanger_val=np.load(sanger_val_pth.as_posix(),mmap_mode='r')
    nci_val=np.load(nci_val_pth.as_posix(),mmap_mode='r')
    gse_val=np.load(gse_val_pth.as_posix(),mmap_mode='r')
    
    n=3  #the dimension for pca
    if(request.POST['show_type'] == 'd_center'):
        
    else:
        #if 'NCI60' or 'GSE36133' in datas:                          
        print(datas)
        if 'NCI60' in datas and 'GSE36133' in datas:
            dataset_name=['NCI60','GSE36133']
            nci_val=np.matrix(nci_val) #need fix
            tnci_val=np.transpose(nci_val)
            gse_val=np.matrix(gse_val) #need fix
            tgse_val=np.transpose(gse_val)
            val=np.concatenate((tnci_val, tgse_val))
            samples=Sample.objects.filter(dataset_id__name__in=dataset_name)  
            ss=samples.select_related('cell_line_id','dataset_id','cell_line_id__name')
            cell_line_name =list(ss.values_list('cell_line_id__name', flat=True))  #all name in datasets
            
            
        elif 'NCI60' in datas:
            dataset_name=['NCI60']
            nci_val=nci_val[~np.isnan(nci_val).any(axis=1)]
            nci_val=np.matrix(nci_val)
            val=np.transpose(nci_val)
            samples=Sample.objects.filter(dataset_id__name__in=dataset_name)  
            ss=samples.select_related('cell_line_id','dataset_id','cell_line_id__name')
            cell_line_name =list(ss.values_list('cell_line_id__name', flat=True))  #all name in datasets
        elif 'GSE36133' in datas:
            
            dataset_name=['GSE36133']
            gse_val=gse_val[~np.isnan(gse_val).any(axis=1)]
            gse_val=np.matrix(gse_val)
            val=np.transpose(gse_val)
            samples=Sample.objects.filter(dataset_id__name__in=dataset_name)  
            #print(samples)
            ss=samples.select_related('cell_line_id','dataset_id','cell_line_id__name')
            cell_line_name =list(ss.values_list('cell_line_id__name', flat=True))  #all name in datasets
        else:
            dataset_name=['Sanger Cell Line Project']
            sanger_val=sanger_val[~np.isnan(sanger_val).any(axis=1)]
            sanger_val=np.matrix(sanger_val)
            val=np.transpose(sanger_val)
            samples=Sample.objects.filter(dataset_id__name__in=dataset_name)  
            
            ss=samples.select_related('cell_line_id','dataset_id','cell_line_id__name')
            cell_line_name =list(ss.values_list('cell_line_id__name', flat=True))  #all name in datasets
        
        pca= PCA(n_components=n)
        Xplus2 = pca.fit_transform(val)
        tplus2=pca.explained_variance_ratio_
        propotion=sum(tplus2[0:n-1])
        
        all_name=cell_line_name
        output_cell=[]
        
    
        temp=0
             
        for k in cell:
            all_name=list(filter((k.name).__ne__, all_name))
            print(len(all_name))
            
            output_cell.append([k,[]])
            range_size=0
            count_set=0
            if 'NCI60' in datas and 'GSE36133' in datas:
                scell_nci=ss.filter(dataset_id__name__in=["NCI60"],cell_line_id__name=k.name) 
                scell_ccle=ss.filter(dataset_id__name__in=["GSE36133"],cell_line_id__name=k.name)
                offset_nci=list(scell_nci.values_list('offset',flat=True))
                offset_ccle=list(scell_ccle.values_list('offset',flat=True))
                new_index=Sample.objects.filter(dataset_id__name__in=["NCI60"]).count()
                offset_ccle=[x+new_index for x in offset_ccle]
                offset=offset_nci+offset_ccle
                add_size=0
            elif 'NCI60' in datas:
                scell_nci=ss.filter(dataset_id__name__in=["NCI60"],cell_line_id__name=k.name) 
                offset_nci=list(scell_nci.values_list('offset',flat=True))
                offset=offset_nci
                add_size=0
            elif 'GSE36133' in datas:
                scell_ccle=ss.filter(dataset_id__name__in=["GSE36133"],cell_line_id__name=k.name)
                offset_ccle=list(scell_ccle.values_list('offset',flat=True))
                offset=offset_ccle
                add_size=0
            else:    
                scell_sanger=ss.filter(dataset_id__name__in=['Sanger Cell Line Project'],cell_line_id__name=k.name)
                offset_sanger=list(scell_sanger.values_list('offset',flat=True))
                offset=offset_sanger
                add_size=0
                
            total_offset=total_offset+offset
            for d in dataset_name:
                
                range_size=Sample.objects.filter(dataset_id__name__in=[d]).count()
                
            
                counter=1
                print(add_size)
                print(range_size)
                
                    
                for i in offset:
                    for j in range(add_size,range_size+add_size):
                        if i!=j:
                            output_cell[temp][1].append([counter,cell_line_name[j],d,np.linalg.norm(Xplus2[j]-Xplus2[i])])
                    if count_set == 0:
                        selected_name.append(k.name+'('+str(counter)+')')
                        colorX.append(Xplus2[i][0])
                        colorY.append(Xplus2[i][1])
                        colorZ.append(Xplus2[i][2])
                    counter=counter+1
                count_set=count_set+1    
                add_size=add_size+range_size
                print(offset)
            
            temp=temp+1
            
            
        for j in range(0,len(Xplus2)):
            if j not in total_offset:
                X.append(Xplus2[j][0])
                Y.append(Xplus2[j][1])
                Z.append(Xplus2[j][2])     
        print("X:",len(X))    
            
    return render_to_response('pca.html',RequestContext(request,
    {
    'output_cell':output_cell,
    'propotion':propotion,
    'all_name':mark_safe(json.dumps(all_name)),
    'selected_name':mark_safe(json.dumps(selected_name)),
    'colorX':colorX,
    'colorY':colorY,
    'colorZ':colorZ,
    'X':X,
    'Y':Y,
    'Z':Z,
    }))

    

def cellline_microarray(request):
    # Pre-fetch the cell line field for all samples.
    # Reduce N query in to 1. N = number of samples
    samples = Sample.objects.filter(
        dataset_id__name__in=['Sanger Cell Line Project']
    ).select_related('cell_line_id')
    ncisamples = Sample.objects.filter(
        dataset_id__name__in=['NCI60']
    ).select_related('cell_line_id')
    CCsamples = Sample.objects.filter(
        dataset_id__name__in=['GSE36133']
    ).select_related('cell_line_id')
    # Get all distinct primary sites from selected samples
    primary_sites = sorted(
        samples.values_list('cell_line_id__primary_site', flat=True).distinct()
    )
    nciprimary_sites = sorted(
        ncisamples.values_list('cell_line_id__primary_site', flat=True).distinct()
    )
    CCprimary_sites = sorted(
        CCsamples.values_list('cell_line_id__primary_site', flat=True).distinct()
    )
    

    return render(request, 'cellline_microarray.html', {
        'samples': samples,
        'primary_sites': primary_sites,
        'ncisamples': ncisamples,
        'nciprimary_sites': nciprimary_sites,
        'CCsamples': CCsamples,
        'CCprimary_sites': CCprimary_sites,
        
    })


def cell_lines(request):
    #samples = Sample.objects.all().select_related('cell_line_id','dataset_id')    
    #lines=CellLine.objects.all().distinct()
    #val_pairs = (
    #            (l, l.fcell_line_id.prefetch_related('dataset_id__name').values_list('dataset_id__name',flat=True).distinct())                        
    #            for l in lines
    #        )
    #context['val_pairs']=val_pairs
    cell_line_dict={}
    context={}
    nr_samples=[]
    samples=Sample.objects.all().select_related('cell_line_id','dataset_id')
    for ss in samples:
        name=ss.cell_line_id.name
        primary_site=ss.cell_line_id.primary_site
        primary_hist=ss.cell_line_id.primary_hist
        comb=name+"/"+primary_site+"/"+primary_hist
        dataset=ss.dataset_id.name
        try: 
            sets=cell_line_dict[comb]
            if (dataset not in sets):
                cell_line_dict[comb]=dataset+"/"+sets
        except KeyError:
            cell_line_dict[comb]=dataset
            nr_samples.append(ss)
            
            
    val_pairs = (
                (ss,cell_line_dict[ss.cell_line_id.name+"/"+ss.cell_line_id.primary_site+"/"+ss.cell_line_id.primary_hist])                        
                for ss in nr_samples
            )                   
    context['val_pairs']=val_pairs
    return render_to_response('cell_line.html', RequestContext(request, context))
    

def data(request):
    SANGER=[]
    sanger_flag=0
    NCI=[]
    nci_flag=0
    GSE=[]
    gse_flag=0
    cell=[]
    ncicell=[]
    CCcell=[]
    ps_id='0'
    pn_id='0'
    if request.POST['cell_line_method'] == 'text':
        if request.POST['cellline'] =='':
            return HttpResponse("<p>please make sure to enter cell line name in Step3.</p>" )
        c = request.POST['cellline']
        c = c.split()
        sanger_flag=1
        samples=Sample.objects.filter(dataset_id__name__in=['Sanger Cell Line Project'])      
        cell=samples.select_related('cell_line_id','dataset_id').filter(cell_line_id__name__in=c)
        offset=cell.values_list('offset',flat=True)
        ps_id='1'
        
        nci_flag=1
        ncisamples=Sample.objects.filter(dataset_id__name__in=['NCI60']).select_related('cell_line_id','dataset_id')
        ncicell=ncisamples.filter(cell_line_id__name__in=c)
        ncioffset=ncicell.values_list('offset',flat=True)
        pn_id='3'
        
        gse_flag=1
        CCsamples=Sample.objects.filter(dataset_id__name__in=['GSE36133']).select_related('cell_line_id','dataset_id')
        CCcell=CCsamples.filter(cell_line_id__name__in=c)
        CCoffset=CCcell.values_list('offset',flat=True)
        pn_id='3'
    else:
        if 'dataset' in request.POST and request.POST['dataset'] != '':
            datas=request.POST.getlist('dataset')
            if 'Sanger Cell Line Project' in datas:
                sanger_flag=1
                SANGER=request.POST.getlist('select_sanger')
                samples=Sample.objects.filter(dataset_id__name__in=['Sanger Cell Line Project'])      
                cell=samples.select_related('cell_line_id','dataset_id').filter(cell_line_id__primary_site__in=SANGER)
                offset=cell.values_list('offset',flat=True)
                ps_id='1'
            if 'NCI60' in datas:
                nci_flag=1
                NCI=request.POST.getlist('select_nci')
                ncisamples=Sample.objects.filter(dataset_id__name__in=['NCI60']).select_related('cell_line_id','dataset_id')
                ncicell=ncisamples.filter(cell_line_id__primary_site__in=NCI)
                ncioffset=ncicell.values_list('offset',flat=True)
                pn_id='3'
            if 'GSE36133' in datas:
                gse_flag=1
                GSE=request.POST.getlist('select_gse')
                CCsamples=Sample.objects.filter(dataset_id__name__in=['GSE36133']).select_related('cell_line_id','dataset_id')
                CCcell=CCsamples.filter(cell_line_id__primary_site__in=GSE)
                CCoffset=CCcell.values_list('offset',flat=True)
                pn_id='3'
            if len(SANGER)==0 and len(NCI)==0 and len(GSE)==0:
                return HttpResponse("<p>please select primary sites.</p>" )
        else:
            return HttpResponse("<p>please check Step3 again.</p>" )
        
    
    
    if 'keyword' in request.POST and request.POST['keyword'] != '':
        words = request.POST['keyword']
        words = words.split()
    else:
        return HttpResponse("<p>where is your keyword?</p>")
      
    #open files
    sanger_val_pth=Path('../').resolve().joinpath('src','sanger_cell_line_proj.npy')
    nci_val_pth=Path('../').resolve().joinpath('src','nci60.npy')
    gse_val_pth=Path('../').resolve().joinpath('src','GSE36133.npy')
    sanger_val=np.load(sanger_val_pth.as_posix(),mmap_mode='r')
    nci_val=np.load(nci_val_pth.as_posix(),mmap_mode='r')
    gse_val=np.load(gse_val_pth.as_posix(),mmap_mode='r')
    
    gene = []
    ncigene = []
    CCgene = []
    context={}
    
    norm_name=[request.POST['normalize']]
    if sanger_flag==1:
        #if request.POST['normalize']!='NTRK3-AS1':
        sanger_g=ProbeID.objects.filter(platform__in=ps_id).filter(Gene_symbol__in=norm_name)
        sanger_probe_offset=sanger_g.values_list('offset',flat=True)
        temp=sanger_val[np.ix_(sanger_probe_offset,offset)]
        norm=np.mean(temp,axis=0, dtype=np.float64,keepdims=True)
        #else:
        #    norm=0.0
    else:
        norm=0.0  #if / should = 1   
    if nci_flag==1:
        nci_g=ProbeID.objects.filter(platform__in=pn_id).filter(Gene_symbol__in=norm_name)
        nci_probe_offset=nci_g.values_list('offset',flat=True)
        temp=nci_val[np.ix_(nci_probe_offset,ncioffset)]
        nci_norm=np.mean(temp,axis=0, dtype=np.float64,keepdims=True)
        #print(nci_norm)   
    else:
        nci_norm=0.0  #if / should = 1
    if gse_flag==1:
        CC_g=ProbeID.objects.filter(platform__in=pn_id).filter(Gene_symbol__in=norm_name)
        CC_probe_offset=CC_g.values_list('offset',flat=True)
        temp=gse_val[np.ix_(CC_probe_offset,CCoffset)]
        CC_norm=np.mean(temp,axis=0, dtype=np.float64,keepdims=True)
        #print(CC_norm)
    else:
        CC_norm=0.0  #if / should = 1             

        
            
    #dealing with probes    
    if 'gtype' in request.POST and request.POST['gtype'] == 'probeid':
        gene = ProbeID.objects.filter(platform__in=ps_id).filter(Probe_id__in=words)
        probe_offset=gene.values_list('offset',flat=True)
        
        ncigene = ProbeID.objects.filter(platform__in=pn_id).filter(Probe_id__in=words)
        nciprobe_offset=ncigene.values_list('offset',flat=True)
        #nci60 and ccle use same probe set(ncigene) and nicprobe
        
        # Make a generator to generate all (cell, probe, val) pairs
        if(len(gene)!=0 and len(cell)!=0):
            raw_test=sanger_val[np.ix_(probe_offset,offset)]
            normalize=np.subtract(raw_test,norm)#dimension different!!!!
            #normalize=np.around(normalize, decimals=1)
            cell_probe_val_pairs = (
                (c, p, raw_test[probe_ix, cell_ix],normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(gene)
                for cell_ix, c in enumerate(cell)
            )
            
        else:
            cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(ncicell)!=0):
            nci_raw_test=nci_val[np.ix_(nciprobe_offset,ncioffset)]
            nci_normalize=np.subtract(nci_raw_test,nci_norm)
            nci_cell_probe_val_pairs = (
                (c, p, nci_raw_test[probe_ix, cell_ix],nci_normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(ncigene)
                for cell_ix, c in enumerate(ncicell)
            )
            
        else:
            nci_cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(CCcell)!=0):
            CC_raw_test=gse_val[np.ix_(nciprobe_offset,CCoffset)]
            CC_normalize=np.subtract(CC_raw_test,CC_norm)
            CC_cell_probe_val_pairs = (
                (c, p, CC_raw_test[probe_ix, cell_ix],CC_normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(ncigene)
                for cell_ix, c in enumerate(CCcell)
            )
            
        else:
            CC_cell_probe_val_pairs =()
        context['cell_probe_val_pairs']=cell_probe_val_pairs
        context['nci_cell_probe_val_pairs']=nci_cell_probe_val_pairs
        context['CC_cell_probe_val_pairs']=CC_cell_probe_val_pairs
        return render_to_response('data.html', RequestContext(request,context))

    elif 'gtype' in request.POST and request.POST['gtype'] == 'symbol':
        gene = ProbeID.objects.filter(platform__in=ps_id).filter(Gene_symbol__in=words)
        probe_offset=gene.values_list('offset',flat=True)
        
        ncigene = ProbeID.objects.filter(platform__in=pn_id).filter(Gene_symbol__in=words)
        nciprobe_offset=ncigene.values_list('offset',flat=True)
        #nci60 and ccle use same probe set(ncigene) and nicprobe
        
        # Make a generator to generate all (cell, probe, val) pairs
        if(len(gene)!=0 and len(cell)!=0):
            raw_test=sanger_val[np.ix_(probe_offset,offset)]
            normalize=np.subtract(raw_test,norm)
            cell_probe_val_pairs = (
                (c, p, raw_test[probe_ix, cell_ix],normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(gene)
                for cell_ix, c in enumerate(cell)
            )
            
        else:
            cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(ncicell)!=0):
            nci_raw_test=nci_val[np.ix_(nciprobe_offset,ncioffset)]
            nci_normalize=np.subtract(nci_raw_test,nci_norm)
            nci_cell_probe_val_pairs = (
                (c, p, nci_raw_test[probe_ix, cell_ix],nci_normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(ncigene)
                for cell_ix, c in enumerate(ncicell)
            )
            
        else:
            nci_cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(CCcell)!=0):
            CC_raw_test=gse_val[np.ix_(nciprobe_offset,CCoffset)]
            CC_normalize=np.subtract(CC_raw_test,CC_norm)
            CC_cell_probe_val_pairs = (
                (c, p, CC_raw_test[probe_ix, cell_ix],CC_normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(ncigene)
                for cell_ix, c in enumerate(CCcell)
            )
            
        else:
            CC_cell_probe_val_pairs =()
        context['cell_probe_val_pairs']=cell_probe_val_pairs
        context['nci_cell_probe_val_pairs']=nci_cell_probe_val_pairs
        context['CC_cell_probe_val_pairs']=CC_cell_probe_val_pairs
        return render_to_response('data.html', RequestContext(request,context))

    elif 'gtype' in request.POST and request.POST['gtype'] == 'entrez':
        gene = ProbeID.objects.filter(platform__in=ps_id).filter(Entrez_id=words)
        probe_offset=gene.values_list('offset',flat=True)
        
        ncigene = ProbeID.objects.filter(platform__in=pn_id).filter(Entrez_id__in=words)
        nciprobe_offset=ncigene.values_list('offset',flat=True)
        #nci60 and ccle use same probe set(ncigene) and nicprobe
        
        # Make a generator to generate all (cell, probe, val) pairs
        if(len(gene)!=0 and len(cell)!=0):
            raw_test=sanger_val[np.ix_(probe_offset,offset)]
            normalize=np.subtract(raw_test,norm)
            cell_probe_val_pairs = (
                (c, p, raw_test[probe_ix, cell_ix],normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(gene)
                for cell_ix, c in enumerate(cell)
            )
            
        else:
            cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(ncicell)!=0):
            nci_raw_test=nci_val[np.ix_(nciprobe_offset,ncioffset)]
            nci_normalize=np.subtract(nci_raw_test,nci_norm)
            nci_cell_probe_val_pairs = (
                (c, p, nci_raw_test[probe_ix, cell_ix],nci_normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(ncigene)
                for cell_ix, c in enumerate(ncicell)
            )
            
        else:
            nci_cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(CCcell)!=0):
            CC_raw_test=gse_val[np.ix_(nciprobe_offset,CCoffset)]
            CC_normalize=np.subtract(CC_raw_test,CC_norm)
            CC_cell_probe_val_pairs = (
                (c, p, CC_raw_test[probe_ix, cell_ix],CC_normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(ncigene)
                for cell_ix, c in enumerate(CCcell)
            )
            
        else:
            CC_cell_probe_val_pairs =()
        context['cell_probe_val_pairs']=cell_probe_val_pairs
        context['nci_cell_probe_val_pairs']=nci_cell_probe_val_pairs
        context['CC_cell_probe_val_pairs']=CC_cell_probe_val_pairs
        return render_to_response('data.html', RequestContext(request,context))
    else:
        return HttpResponse(
            "<p>keyword type not match with your keyword input</p>"
        )
