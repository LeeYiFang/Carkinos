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
    lines = CellLine.objects.all()
    return render_to_response('cell_line.html', RequestContext(request, locals()))


def data(request):
    SANGER=[]
    NCI=[]
    GSE=[]
    cell=[]
    ncicell=[]
    CCcell=[]
    ps_id='0'
    pn_id='0'
    pc_id='0'
    if 'dataset' in request.POST and request.POST['dataset'] != '':
        datas=request.POST.getlist('dataset')
        if 'Sanger Cell Line Project' in datas:
            SANGER=request.POST.getlist('select_sanger')
            samples=Sample.objects.filter(dataset_id__name__in=['Sanger Cell Line Project']).select_related('cell_line_id')
            cell=samples.filter(cell_line_id__primary_site__in=SANGER)
            offset=cell.values_list('offset',flat=True)
            ps_id='1'
        if 'NCI60' in datas:
            NCI=request.POST.getlist('select_nci')
            ncisamples=Sample.objects.filter(dataset_id__name__in=['NCI60']).select_related('cell_line_id')
            ncicell=ncisamples.filter(cell_line_id__primary_site__in=NCI)
            pn_id='3'
        if 'GSE36133' in datas:
            GSE=request.POST.getlist('select_gse')
            CCsamples=Sample.objects.filter(dataset_id__name__in=['GSE36133']).select_related('cell_line_id')
            CCcell=CCsamples.filter(cell_line_id__primary_site__in=GSE)
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
        
    #if('dataset' in request.POST and request.POST['dataset']=='Sanger Cell Line Project'):
    #    p_id='1'
    #elif('dataset' in request.POST and request.POST['dataset']=='NCI60'):
    #    p_id='3'
    #else:
    #    return HttpResponse("<p>where is your cell line?</p>")
    dset_val_pth=Path('../').resolve().joinpath('src','sanger_cell_line_proj.npy')
    #Sfirst=Dataset.objects.first()
    #dset_val_pth=Path(dset_path,Sfirst.data_path)
    dset_val=np.load(dset_val_pth.as_posix(),mmap_mode='r')
    #raw_test=dset_val[np.ix_([1,2,3],[1,2])]
    
    
    gene = []
    ncigene = []
    CCgene = []
    if 'gtype' in request.POST and request.POST['gtype'] == 'probeid':
        #gene = ProbeID.objects.filter(platform__in=p_id).filter()
        gene = ProbeID.objects.filter(platform__in=ps_id).filter(Probe_id__in=words)
        probe_offset=gene.values_list('offset',flat=True)
        raw_test=dset_val[np.ix_(probe_offset,offset)]
        ncigene = ProbeID.objects.filter(platform__in=pn_id).filter(Probe_id__in=words)
        #CCgene = ProbeID.objects.filter(platform__in=pc_id).filter(Probe_id__in=words)
        #print(CCgene)

        # Make a generator to generate all (cell, probe, val) pairs
        cell_probe_val_pairs = (
            (c, p, raw_test[probe_ix, cell_ix])
            for probe_ix, p in enumerate(gene)
            for cell_ix, c in enumerate(cell)
        )
        return render_to_response('data.html', RequestContext(request,{
            'gene': gene,
            'cell': cell,
            'ncigene': ncigene,
            'ncicell': ncicell,
            'raw_test': raw_test,
            'cell_probe_val_pairs': cell_probe_val_pairs,
            'CCgene': CCgene,
            'CCcell': CCcell,
        }))

    elif 'gtype' in request.POST and request.POST['gtype'] == 'symbol':
        gene = ProbeID.objects.filter(platform__in=ps_id).filter(Gene_symbol__in=words)
        ncigene = ProbeID.objects.filter(platform__in=pn_id).filter(Gene_symbol__in=words)
        #CCgene = ProbeID.objects.filter(platform__in=pc_id).filter(Gene_symbol__in=words)
        return render_to_response('data.html', RequestContext(request,{
            'gene': gene,
            'cell': cell,
            'ncigene': ncigene,
            'ncicell': ncicell,
            'raw_test': raw_test,
            'CCgene': CCgene,
            'CCcell': CCcell,
        }))

    elif 'gtype' in request.POST and request.POST['gtype'] == 'entrez':
        gene = ProbeID.objects.filter(platform__in=ps_id).filter(Entrez_id=words)
        ncigene = ProbeID.objects.filter(platform__in=pn_id).filter(Entrez_id__in=words)
        #CCgene = ProbeID.objects.filter(platform__in=pc_id).filter(Entrez_id=words)
        return render_to_response('data.html', RequestContext(request,{
            'gene': gene,
            'cell': cell,
            'ncigene': ncigene,
            'ncicell': ncicell,
            'raw_test': raw_test,
            'CCgene': CCgene,
            'CCcell': CCcell,
        }))
    else:
        return HttpResponse(
            "<p>keyword type not match with your keyword input</p>"
        )
