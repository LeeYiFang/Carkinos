from django.shortcuts import render,render_to_response
from django.http import HttpResponse, Http404
from django.views.decorators.http import require_GET
from .models import Dataset, CellLine, ProbeID, Sample, Platform
from django.template import RequestContext
import pandas as pd
import numpy as np
from pathlib import Path


def home(request):
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
    

    return render(request, 'home.html', {
        'samples': samples,
        'primary_sites': primary_sites,
        'ncisamples': ncisamples,
        'nciprimary_sites': nciprimary_sites,
        'CCsamples': CCsamples,
        'CCprimary_sites': CCprimary_sites,
        
    })


def cell_lines(request):
    #lines = Sample.objects.select_related('cell_line_id','dataset_id')
    lines=CellLine.objects.all().distinct()
    
    val_pairs = (
                (l, l.fcell_line_id.prefetch_related('dataset_id__name').values_list('dataset_id__name',flat=True).distinct())                        
                for l in lines
            )
   
    return render_to_response('cell_line.html', RequestContext(request, locals()))
    

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
    
    #need to modify by correct number
    if request.POST['normalize'] == 'GAPDH' or request.POST['normalize'] == 'ACTB':
        CV_flag=0
        norm_name=[request.POST['normalize']]
        if sanger_flag==1:
            sanger_g=ProbeID.objects.filter(platform__in=ps_id).filter(Gene_symbol__in=norm_name)
            sanger_probe_offset=sanger_g.values_list('offset',flat=True)
            temp=sanger_val[np.ix_(sanger_probe_offset,offset)]
            norm=np.mean(temp,axis=0, dtype=np.float64,keepdims=True)
        else:
            norm=1.0  #if / should = 1   
        if nci_flag==1:
            nci_g=ProbeID.objects.filter(platform__in=pn_id).filter(Gene_symbol__in=norm_name)
            nci_probe_offset=nci_g.values_list('offset',flat=True)
            temp=nci_val[np.ix_(nci_probe_offset,ncioffset)]
            nci_norm=np.mean(temp,axis=0, dtype=np.float64,keepdims=True)
            #nci_norm=np.round(nci_norm, decimals=4)
        else:
            nci_norm=1.0  #if / should = 1
        if gse_flag==1:
            CC_g=ProbeID.objects.filter(platform__in=pn_id).filter(Gene_symbol__in=norm_name)
            CC_probe_offset=CC_g.values_list('offset',flat=True)
            temp=gse_val[np.ix_(CC_probe_offset,CCoffset)]
            CC_norm=np.mean(temp,axis=0, dtype=np.float64,keepdims=True)
            #CC_norm=np.around(CC_norm, decimals=1)
        else:
            CC_norm=1.0  #if / should = 1             
    else:
        norm=10
        nci_norm=10
        CC_norm=10
        
            
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
            normalize=np.divide(raw_test,norm)#dimension different!!!!
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
            nci_normalize=np.divide(nci_raw_test,nci_norm)
            nci_cell_probe_val_pairs = (
                (c, p, nci_raw_test[probe_ix, cell_ix],nci_normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(ncigene)
                for cell_ix, c in enumerate(ncicell)
            )
            
        else:
            nci_cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(CCcell)!=0):
            CC_raw_test=gse_val[np.ix_(nciprobe_offset,CCoffset)]
            CC_normalize=np.divide(CC_raw_test,CC_norm)
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
            normalize=np.divide(raw_test,norm)
            cell_probe_val_pairs = (
                (c, p, raw_test[probe_ix, cell_ix],normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(gene)
                for cell_ix, c in enumerate(cell)
            )
            
        else:
            cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(ncicell)!=0):
            nci_raw_test=nci_val[np.ix_(nciprobe_offset,ncioffset)]
            nci_normalize=np.divide(nci_raw_test,nci_norm)
            nci_cell_probe_val_pairs = (
                (c, p, nci_raw_test[probe_ix, cell_ix],nci_normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(ncigene)
                for cell_ix, c in enumerate(ncicell)
            )
            
        else:
            nci_cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(CCcell)!=0):
            CC_raw_test=gse_val[np.ix_(nciprobe_offset,CCoffset)]
            CC_normalize=np.divide(CC_raw_test,CC_norm)
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
            normalize=np.divide(raw_test,norm)
            cell_probe_val_pairs = (
                (c, p, raw_test[probe_ix, cell_ix],normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(gene)
                for cell_ix, c in enumerate(cell)
            )
            
        else:
            cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(ncicell)!=0):
            nci_raw_test=nci_val[np.ix_(nciprobe_offset,ncioffset)]
            nci_normalize=np.divide(nci_raw_test,nci_norm)
            nci_cell_probe_val_pairs = (
                (c, p, nci_raw_test[probe_ix, cell_ix],nci_normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(ncigene)
                for cell_ix, c in enumerate(ncicell)
            )
            
        else:
            nci_cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(CCcell)!=0):
            CC_raw_test=gse_val[np.ix_(nciprobe_offset,CCoffset)]
            CC_normalize=np.divide(CC_raw_test,CC_norm)
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
