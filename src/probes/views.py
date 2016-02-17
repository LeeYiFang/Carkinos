from django.shortcuts import render,render_to_response
from django.http import HttpResponse, Http404
from django.views.decorators.http import require_GET
from .models import CellLine, ProbeID
from django.template import RequestContext


def home(request):
    # return render(request, 'home.html')
    return render_to_response('home.html', RequestContext(request, locals()))


def cell_lines(request):
    lines = CellLine.objects.all()
    return render_to_response('cell_line.html', RequestContext(request, locals()))


def data(request):

    if 'cellline' in request.POST and request.POST['cellline'] != '':
        cell = CellLine.objects.filter(name__in=request.POST['cellline'].split())
    else:
        return HttpResponse("<p>where is the cell line?</p>")
    if 'keyword' in request.POST and request.POST['keyword'] != '':
        words = request.POST['keyword']
        words = words.split()
        # return HttpResponse(words)
    else:
        return HttpResponse("<p>where is your keyword?</p>")

    gene = []
    if 'gtype' in request.POST and request.POST['gtype'] == 'probeid':
        gene = ProbeID.objects.filter(Probe_id__in=words)
        return render_to_response('data.html', RequestContext(request,{
            'gene': gene,
            'cell': cell,
        }))

    elif 'gtype' in request.POST and request.POST['gtype'] == 'symbol':
        gene = ProbeID.objects.filter(Gene_symbol__in=words)
        return render_to_response('data.html', RequestContext(request,{
            'gene': gene,
            'cell': cell,
        }))

    elif 'gtype' in request.POST and request.POST['gtype'] == 'entrez':
        gene = ProbeID.objects.filter(Entrez_id=words)
        return render_to_response('data.html', RequestContext(request,{
            'gene': gene,
            'cell': cell,
        }))
    else:
        return HttpResponse(
            "<p>keyword type not match with your keyword input</p>"
        )
