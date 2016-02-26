from django.shortcuts import render,render_to_response
from django.http import HttpResponse, Http404
from django.views.decorators.http import require_GET
from .models import CellLine, ProbeID, Sample, Platform
from django.template import RequestContext


def home(request):
    # return render(request, 'home.html')
    #loading this part is very slow, how to fix it?
    set0 = Sample.objects.filter(dataset_id__in='1')
    #keys1=set0.values('cell_line_id__primary_site').distinct()
    known_cell_lines = {}
    
    #set0.cell_line_id.primary_site.distinct()
    for s in set0:
        cl=s.cell_line_id.primary_site
        try:
            test=known_cell_lines[cl]
        except KeyError:
            if cl=='nan':
                known_cell_lines[cl]='NAN' 
            else :
                known_cell_lines[cl]=cl
    keys1=known_cell_lines.keys()
    return render_to_response('home.html', RequestContext(request,locals()))


def cell_lines(request):
    lines = CellLine.objects.all()
    return render_to_response('cell_line.html', RequestContext(request, locals()))


def data(request):

    
    if 'cellline' in request.POST and request.POST['cellline'] != '':
        cell = CellLine.objects.filter(name__in=request.POST['cellline'].split())
    else:
        #return render_to_response('home.html', RequestContext(request,locals()))
        return HttpResponse("<p>where is the cell line? please check Step3 again.</p>")
    if 'keyword' in request.POST and request.POST['keyword'] != '':
        words = request.POST['keyword']
        words = words.split()
        # return HttpResponse(words)
    else:
        return HttpResponse("<p>where is your keyword?</p>")
        
    if(request.POST['dataset']!='' and request.POST['dataset']=='sanger'):
        p_id='1'
    elif(request.POST['dataset']!='' and request.POST['dataset']=='nci'):
        p_id='3'
    #else:
    #    p_id='2'
    gene = []
    if 'gtype' in request.POST and request.POST['gtype'] == 'probeid':
        gene = ProbeID.objects.filter(platform__in=p_id).filter(Probe_id__in=words)
        return render_to_response('data.html', RequestContext(request,{
            'gene': gene,
            'cell': cell,
        }))

    elif 'gtype' in request.POST and request.POST['gtype'] == 'symbol':
        gene = ProbeID.objects.filter(platform__in=p_id).filter(Gene_symbol__in=words)
        return render_to_response('data.html', RequestContext(request,{
            'gene': gene,
            'cell': cell,
        }))

    elif 'gtype' in request.POST and request.POST['gtype'] == 'entrez':
        gene = ProbeID.objects.filter(platform__in=p_id).filter(Entrez_id=words)
        return render_to_response('data.html', RequestContext(request,{
            'gene': gene,
            'cell': cell,
        }))
    else:
        return HttpResponse(
            "<p>keyword type not match with your keyword input</p>"
        )
