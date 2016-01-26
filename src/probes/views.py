from django.shortcuts import render,render_to_response
from django.http import HttpResponse, Http404
from django.views.decorators.http import require_GET
from .models import CellLine, ProbeID


def home(request):
    # return render(request, 'home.html')
    return render_to_response('home.html', locals())


def cell_lines(request):
    lines = CellLine.objects.all()
    return render_to_response('cell_line.html', locals())


def data(request):

    if 'cellline' in request.GET and request.GET['cellline'] != '':
        cell = CellLine.objects.filter(name__in=request.GET['cellline'].split())
    else:
        return HttpResponse("<p>where is the cell line?</p>")
    if 'keyword' in request.GET and request.GET['keyword'] != '':
        words = request.GET['keyword']
        words = words.split()
        # return HttpResponse(words)
    else:
        return HttpResponse("<p>where is your keyword?</p>")

    gene = []
    if 'gtype' in request.GET and request.GET['gtype'] == 'probeid':
        gene = ProbeID.objects.filter(Probe_id__in=words)
        return render_to_response('data.html', {
            'gene': gene,
            'cell': cell,
        })

    elif 'gtype' in request.GET and request.GET['gtype'] == 'symbol':
        for w in words:
            gene += ProbeID.objects.filter(Gene_symbol=w)
        return render_to_response('data.html', locals())

    elif 'gtype' in request.GET and request.GET['gtype'] == 'entrez':
        for w in words:
            gene += ProbeID.objects.filter(Entrez_id=w)
        return render_to_response('data.html', locals())
    else:
        return HttpResponse(
            "<p>keyword type not match with your keyword input</p>"
        )
