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
from scipy import stats


def upload(request):
    return render_to_response('upload.html',locals())

def welcome(request):
    return render_to_response('welcome.html',locals())

def help_similar_assessment(request):
    return render_to_response('help_similar_assessment.html',RequestContext(request))

def similar_assessment(request):   
    return render_to_response('similar_assessment.html',RequestContext(request))
    
def gene_signature(request):    
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
    celllines = sorted(
        samples.values_list('cell_line_id__name', flat=True).distinct()
    )
    ncicelllines = sorted(
        ncisamples.values_list('cell_line_id__name', flat=True).distinct()
    )
    CCcelllines = sorted(
        CCsamples.values_list('cell_line_id__name', flat=True).distinct()
    )
    
    check_celllines=list(celllines)
    plus2_celllines=list(ncicelllines)+list(CCcelllines)

    return render(request, 'gene_signature.html', {
        'check_celllines': mark_safe(json.dumps(check_celllines)),
        'plus2_celllines': mark_safe(json.dumps(plus2_celllines)),
        'celllines': celllines,
        'ncicelllines': ncicelllines,
        'CCcelllines': CCcelllines,
        
    })
    

def heatmap(request):   
    
    group1=[]
    group2=[]
    group_count=0
    presult={} #{probe object:p value}
    expression=[]
    probe_out=[]
    sample_out=[]
    
    ##open all the file
    sanger_val_pth=Path('../').resolve().joinpath('src','sanger_cell_line_proj.npy')
    nci_val_pth=Path('../').resolve().joinpath('src','nci60.npy')
    gse_val_pth=Path('../').resolve().joinpath('src','GSE36133.npy')
    sanger_val=np.load(sanger_val_pth.as_posix(),mmap_mode='r')
    nci_val=np.load(nci_val_pth.as_posix(),mmap_mode='r')
    gse_val=np.load(gse_val_pth.as_posix(),mmap_mode='r')
    
    
    if request.POST['cell_line_method'] == 'text':
        #this part is for input cell line

        pform=request.POST['data_platform']
        
        group_counter=1
        while True:
            temp_name='cellline_g'+str(group_counter)
            if temp_name in request.POST:
                #print(group_counter)
                group_counter=group_counter+1
            else:
                group_counter=group_counter-1
                break
        
        #text part with more than two group need to use one way ANOVA
        s_group_dict={}  #store sample
        offset_group_dict={} #store offset
        for i in range(1,group_counter+1):
            c='cellline_g'+str(i)
            if request.POST[c] !='':
                s=Sample.objects.filter(cell_line_id__name__in=(request.POST[c].split()),platform_id__name=pform).order_by('dataset_id').select_related('cell_line_id__name','dataset_id')
                s_group_dict['g'+str(i)]=s
                goffset=s.values_list('offset',flat=True)
                offset_group_dict['g'+str(i)]=goffset
        
        #get probe from different platform
        all_probe=ProbeID.objects.filter(platform__name=pform)
        probe_offset=list(all_probe.values_list('offset',flat=True))
                    
        #deal with "nan" in Sanger dataset
        if pform=="U133A":
            all_probe=list(all_probe)
            for x in range(22275,22281):
                probe_offset.index(x)
                probe_offset.remove(x)
                all_probe.pop(x)
            
        #for more than two datasets
        gseindex=-1
        
        val=[] #store value get from binary data 
        group_name=[]
        for i in range(1,group_counter+1):
            nci_data=[]
            gse_data=[]
            temp_name='g'+str(i)
            group_name.append(temp_name)
            if pform=="PLUS2":
                listg=list(s_group_dict[temp_name].values('dataset_id'))
                if {'dataset_id':3} in listg:
                    gseindex=listg.index({'dataset_id':3})
                    nci_data=nci_val[np.ix_(probe_offset,offset_group_dict[temp_name][:gseindex])]
                    gse_data=gse_val[np.ix_(probe_offset,offset_group_dict[temp_name][gseindex:])]
                    g_data=np.concatenate((nci_data,gse_data),axis=1)
                    val.append(g_data.tolist())
                else:
                    g_data=nci_val[np.ix_(probe_offset,offset_group_dict[temp_name])]
                    val.append(g_data.tolist())
            else:
                g_data=sanger_val[np.ix_(probe_offset,offset_group_dict[temp_name])]
                val.append(g_data.tolist())
                
        
          
    else:#this part is for select cell line
        pform=request.POST['data_platform']
        
        #get probe from different platform
        all_probe=ProbeID.objects.filter(platform__name=pform)
        probe_offset=list(all_probe.values_list('offset',flat=True))
        
        all_probe=list(all_probe)
        #deal with "nan" in Sanger dataset
        if pform=="U133A":
            #all_probe=list(all_probe)
            for x in range(22275,22281):
                probe_offset.index(x)
                probe_offset.remove(x)
                all_probe.pop(x)
        
        #count the number of group 
        group_counter=1
        while True:
            temp_name='dataset_g'+str(group_counter)
            if temp_name in request.POST:
                #print(group_counter)
                group_counter=group_counter+1
            else:
                group_counter=group_counter-1
                break
        
        #get binary data
        s_group_dict={}  #store sample
        val=[] #store value get from binary data 
        group_name=[]
        for i in range(1,group_counter+1):
            
            dname='dataset_g'+str(i)
            datasets=request.POST.getlist(dname)
            sanger_data=[]
            nci_data=[]
            gse_data=[]
            temp_name='g'+str(i)
            group_name.append(temp_name)
            if pform=="U133A":
                csanger='select_sanger_g'+str(i)
                if 'Sanger Cell Line Project' in datasets:
                    SANGER=request.POST.getlist(csanger)
                    s=Sample.objects.filter(cell_line_id__name__in=SANGER,platform_id__name=pform).select_related('cell_line_id__name','dataset_id')
                    goffset=s.values_list('offset',flat=True)
                    sanger_data=sanger_val[np.ix_(probe_offset,list(goffset))]
                    s_group_dict['g'+str(i)]=s
                    val.append(sanger_data.tolist())
            else:
                cnci='select_nci_g'+str(i)
                cgse='select_ccle_g'+str(i)
                s_nci=s_gse=[]
                #if 'NCI60' in datasets:
                NCI=request.POST.getlist(cnci)
                s_nci=Sample.objects.filter(cell_line_id__name__in=NCI,dataset_id__name__in=['NCI60']).select_related('cell_line_id__name','dataset_id')
                goffset=s_nci.values_list('offset',flat=True)
                nci_data=nci_val[np.ix_(probe_offset,list(goffset))]
                #if 'GSE36133' in datasets:
                GSE=request.POST.getlist(cgse)
                s_gse=Sample.objects.filter(cell_line_id__name__in=GSE,dataset_id__name__in=['GSE36133']).select_related('cell_line_id__name','dataset_id')
                goffset=s_gse.values_list('offset',flat=True)
                gse_data=gse_val[np.ix_(probe_offset,list(goffset))]
                
                #append nci60 and gse36133 as g_data
               
                s_group_dict['g'+str(i)]=list(s_nci)+list(s_gse)
                g_data=np.concatenate((nci_data,gse_data),axis=1)
                val.append(g_data.tolist())
            
            
            
    #run the one way ANOVA test or ttest for every probe base on the platform selected    
    express={}
    
    if group_counter<=2:
        if len(s_group_dict['g1'])==len(s_group_dict['g2']):
            print("use pair t test")
            for i in range(0,len(all_probe)):   #len(all_probe) need to fix if try to run on laptop
                presult[all_probe[i]]=stats.ttest_rel(list(val[0][i]),list(val[1][i]),nan_policy='omit')[1] 
                express[all_probe[i]]=np.append(val[0][i],val[1][i]).tolist()
        else:
            print("use unpair ttest")
            for i in range(0,len(all_probe)):    #need to fix if try to run on laptop
                presult[all_probe[i]]=stats.ttest_ind(list(val[0][i]),list(val[1][i]),equal_var=False,nan_policy='omit')[1]
                express[all_probe[i]]=np.append(val[0][i],val[1][i]).tolist()
    else:
        print("use one way anova")
        for i in range(0,len(all_probe)):   #need to fix if try to run on laptop
            to_anova=[]
            for n in range(0,group_counter):
                #val[n]=sum(val[n],[])
                to_anova.append(val[n][i])
              
            presult[all_probe[i]]=stats.f_oneway(*to_anova)[1]  
            express[all_probe[i]]=sum(to_anova,[])
       
        
    #sort the dictionary with p-value and need to get the expression data again (top20)  
    #presult[all_probe[0]]=float('nan')
    #presult[all_probe[11]]=float('nan')
    #how to deal with all "nan"?
    tempf=pd.DataFrame(list(presult.items()), columns=['probe', 'pvalue'])
    tempf=tempf.replace(to_replace=float('nan'),value=float('+inf'))
    presult=dict(zip(tempf.probe, tempf.pvalue))
    sortkey=sorted(presult,key=presult.get)
    
    counter=1
    pro_number=int(request.POST['probe_number'])
    stop_end=pro_number+1
    for w in sortkey:     
        #print(presult[w],":",w.Probe_id)
        expression.append(express[w])
        probe_out.append(w.Probe_id+"("+w.Gene_symbol+")")
        counter+=1
        if counter==stop_end:
            break
    
    n_counter=1
    for n in group_name:
        sample_counter=1
        for s in s_group_dict[n]:
            dataset_n=s.dataset_id.name
            if dataset_n=="Sanger Cell Line Project":
                sample_out.append(str(sample_counter)+"-SCLP("+s.cell_line_id.name+")"+"(group"+str(n_counter)+")")   
            else:
                sample_out.append(str(sample_counter)+"-"+s.dataset_id.name+"("+s.cell_line_id.name+")"+"(group"+str(n_counter)+")")   
            sample_counter+=1
        n_counter+=1
    

    return render_to_response('heatmap.html',RequestContext(request,
        {
        'pro_number':pro_number,
        'sample_out':mark_safe(json.dumps(sample_out)),
        'expression':expression,
        'probe_out':mark_safe(json.dumps(probe_out)),
        
        }))

    
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
        
    if 'show_type' in request.POST:
        display_method=request.POST.getlist('show_type')
    else:
        return render_to_response('help_similar_assessment.html',RequestContext(request))
        
    cell=CellLine.objects.filter(name__in=ncell)
    propotion=0
    output_cell=[]
    context={}
    cell_line_dict={}
    dataset_dict={}
    colorX=[]
    colorY=[]
    colorZ=[]
    X=[]
    Y=[]
    Z=[]
    selected_name=[]
    all_name=[]
    dataset_name=[]
    total_offset=[]
    ##open all the file
    sanger_val_pth=Path('../').resolve().joinpath('src','sanger_cell_line_proj.npy')
    nci_val_pth=Path('../').resolve().joinpath('src','nci60.npy')
    gse_val_pth=Path('../').resolve().joinpath('src','GSE36133.npy')
    sanger_val=np.load(sanger_val_pth.as_posix(),mmap_mode='r')
    nci_val=np.load(nci_val_pth.as_posix(),mmap_mode='r')
    gse_val=np.load(gse_val_pth.as_posix(),mmap_mode='r')
    
    if(request.POST['data_platform'] == 'U133A'):
        if 'Sanger Cell Line Project' in datas:
            dataset_name=['Sanger Cell Line Project']
            sanger_val=sanger_val[~np.isnan(sanger_val).any(axis=1)]
            sanger_val=np.matrix(sanger_val)
            val=np.transpose(sanger_val)
            dataset_size=[0,Sample.objects.filter(dataset_id__name__in=["Sanger Cell Line Project"]).count()]
        else:
            return render_to_response('help_similar_assessment.html',RequestContext(request))
    else:
        if 'NCI60' in datas and 'GSE36133' in datas:
            dataset_name=['NCI60','GSE36133']
            nci_val=np.matrix(nci_val)#[:30,:] #need fix
            nci_size=Sample.objects.filter(dataset_id__name__in=["NCI60"]).count()
            tnci_val=np.transpose(nci_val)
            gse_val=np.matrix(gse_val)#[:30,:] #need fix
            tgse_val=np.transpose(gse_val)
            val=np.concatenate((tnci_val, tgse_val))#combine together  ##notice that the offset number will change
            dataset_size=[0,nci_size,nci_size+Sample.objects.filter(dataset_id__name__in=["GSE36133"]).count()]
        elif 'NCI60' in datas:
            dataset_name=['NCI60']
            nci_val=nci_val[~np.isnan(nci_val).any(axis=1)]
            nci_val=np.matrix(nci_val)
            val=np.transpose(nci_val)
            dataset_size=[0,Sample.objects.filter(dataset_id__name__in=["NCI60"]).count()]
        else:            
            dataset_name=['GSE36133']
            gse_val=gse_val[~np.isnan(gse_val).any(axis=1)]
            gse_val=np.matrix(gse_val)
            val=np.transpose(gse_val)
            dataset_size=[0,Sample.objects.filter(dataset_id__name__in=["GSE36133"]).count()]
        

    samples=Sample.objects.filter(dataset_id__name__in=dataset_name).select_related('cell_line_id','dataset_id','cell_line_id__name')  
    sample_name=list(samples.values_list('name', flat=True))
    cell_line_name =list(samples.values_list('cell_line_id__name', flat=True))  #all name in datasets
    all_name=cell_line_name.copy()         #all_name need to delete the selected_name later
    all_name_distinct =list(samples.values_list('cell_line_id__name', flat=True).distinct()) 

    n=429  #need to fix to the best one #need to fix proportion  
    #429 0.900169151952
    #206 0.8
    #99 0.7
    pca= PCA(n_components=n)
    Xval = pca.fit_transform(val)
    ratio_temp=pca.explained_variance_ratio_
    propotion=sum(ratio_temp[0:3])
    table_propotion=sum(ratio_temp[0:n+1])
    
    if 'd_sample' in display_method:
        
        temp=0 #count the cell line number
        for k in cell:
            scell=samples.filter(cell_line_id__name=k.name)
            snames=list(scell.values_list('name', flat=True))
            all_name=list(filter((k.name).__ne__, all_name))
            if len(dataset_name)==1:
                
                offset=list(scell.values_list('offset',flat=True))
                total_offset=total_offset+offset  #combine all the cell line offset together
            else:  
                
                #if has more than 2 datasets: need to fix this part!!
                offset_nci=list(scell.filter(dataset_id__name__in=["NCI60"]).values_list('offset',flat=True))
                offset_ccle=list(scell.filter(dataset_id__name__in=["GSE36133"]).values_list('offset',flat=True))
               
                offset_ccle=[x+dataset_size[1] for x in offset_ccle]
                offset=offset_nci+offset_ccle
                total_offset=total_offset+offset  
            output_cell.append([k,[]])
            counter=1 #count the repeat sample with same cell line name
            
            x=0
            for i in offset:
                for j in range(0,dataset_size[-1]):
                    if j==dataset_size[x+1]:
                        x=x+1
                    if i!=j:  
                        output_cell[temp][1].append([k.name+'('+str(counter)+')',snames[counter-1],cell_line_name[j],sample_name[j],dataset_name[x],np.linalg.norm(Xval[j]-Xval[i])])
                selected_name.append(k.name+'('+str(counter)+')')
                colorX.append(Xval[i][0])
                colorY.append(Xval[i][1])
                colorZ.append(Xval[i][2])
                counter=counter+1
            temp=temp+1
        for j in range(0,dataset_size[-1]):
            if j not in total_offset:
                X.append(Xval[j][0])
                Y.append(Xval[j][1])
                Z.append(Xval[j][2])
                
        return render_to_response('pca.html',RequestContext(request,
        {
        'output_cell':output_cell,
        'propotion':propotion,
        'table_propotion':table_propotion,
        'all_name':mark_safe(json.dumps(all_name)),
        'selected_name':mark_safe(json.dumps(selected_name)),
        'colorX':colorX,
        'colorY':colorY,
        'colorZ':colorZ,
        'X':X,
        'Y':Y,
        'Z':Z,
        }))
    
    else:
        #count all centroid of cell lines in selected datasets
        for name in all_name_distinct:
            
            
            scell=samples.filter(cell_line_id__name=name) #find all the samples that has same cell line name in selected datasets
            if len(dataset_name)==1:
                dataset_dict[name]=dataset_name[0]
                
                offset=list(scell.values_list('offset',flat=True))
                total_offset=offset  #combine all the cell line offset together
            else:  
                
                for ss in scell:
                    d_name=ss.dataset_id.name
                    try: 
                        sets=dataset_dict[name]
                        if (d_name not in sets):
                            dataset_dict[name]=d_name+"/"+sets
                    except KeyError:
                        dataset_dict[name]=d_name
                        #print(name,",",d_name)
            
                
                #if has more than 2 datasets: need to fix this part!!
                offset_nci=list(scell.filter(dataset_id__name__in=["NCI60"]).values_list('offset',flat=True))
                offset_ccle=list(scell.filter(dataset_id__name__in=["GSE36133"]).values_list('offset',flat=True))
               
                offset_ccle=[x+dataset_size[1] for x in offset_ccle]
                offset=offset_nci+offset_ccle
                total_offset=offset 
                
            Xval_a=np.array(Xval)
            selected_Xval=Xval_a[total_offset]
            new_loca=(np.mean(selected_Xval,axis=0,dtype=np.float64,keepdims=True)).tolist()[0]

            cell_line_dict[name]=new_loca
        
        #count distance of selected cell lines and other cell lines base on centroid counted above
        temp=0 #count the cell line number
        temp_cell_line_dict=dict(cell_line_dict)
        
        for k in cell:
            
            
            selected_cell=k.name
            if selected_cell in cell_line_dict:
                del temp_cell_line_dict[selected_cell]
                
                output_cell.append([k,[]])
                for c_name in all_name_distinct:
                    if c_name!=selected_cell:
                        aX=np.array(cell_line_dict[selected_cell])
                        aY=np.array(cell_line_dict[c_name])
                        output_cell[temp][1].append([selected_cell,dataset_dict[selected_cell]
                                                    ,c_name,dataset_dict[c_name]
                                                    ,np.linalg.norm(aX-aY)])    
                
                selected_name.append(selected_cell)
                templocation=cell_line_dict[selected_cell]
                colorX.append(templocation[0])
                colorY.append(templocation[1])
                colorZ.append(templocation[2])
                temp=temp+1
            else:
                output_cell.append([k,[]])
        
        all_name=[]
        for key, value in temp_cell_line_dict.items():
            all_name.append(key)
            #print(value)
            X.append(value[0])    
            Y.append(value[1])
            Z.append(value[2])
            
        
        return render_to_response('pca_center.html',RequestContext(request,
        {
        'output_cell':output_cell,
        'propotion':propotion,
        'table_propotion':table_propotion,
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
