from django.shortcuts import render,render_to_response
from django.http import HttpResponse, Http404
from django.views.decorators.http import require_GET
from .models import CellLine, ProbeID
# Create your views here.


def home(request):
    
    #return render(request, 'home.html')
    return render_to_response('home.html',locals())

def data(request):
    if 'cellline' in request.GET and request.GET['cellline']!='':
        cell = CellLine.objects.get(name=request.GET['cellline'])
    if  'keyword' in request.GET and request.GET['keyword']!='':
        gene = ProbeID.objects.get(Gene_symbol=request.GET['keyword'])
        return render_to_response('data.html',locals())
    else:
        return HttpResponse("<p>datas</p>")