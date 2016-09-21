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

def generate_samples():
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
    nci=[]
    nciprime=Sample.objects.filter(dataset_id__name="NCI60").order_by('cell_line_id__primary_site').select_related('cell_line_id')
    ncisites=list(nciprime.values_list('cell_line_id__primary_site',flat=True))
    ncicells=list(nciprime.values_list('cell_line_id__name',flat=True))
    #ncidi=list(ncicells.values_list('cell_line_id__primary_site',flat=True))
    ncidi=list(nciprime.values_list('cell_line_id__primary_site',flat=True).distinct())
    ncicells=list(ncicells)
    id_counter=0
   
    for p in range(0,len(ncidi)):
        temp=ncisites.count(ncidi[p])
        nci.append((ncidi[p],list(set(ncicells[id_counter:id_counter+temp]))))
        id_counter+=temp
        
    sanger=[]
    sangerprime=Sample.objects.filter(dataset_id__name="Sanger Cell Line Project").order_by('cell_line_id__primary_site').select_related('cell_line_id')
    sangersites=list(sangerprime.values_list('cell_line_id__primary_site',flat=True))
    sangercells=list(sangerprime.values_list('cell_line_id__name',flat=True))
    #ncidi=list(ncicells.values_list('cell_line_id__primary_site',flat=True))
    sangerdi=list(sangerprime.values_list('cell_line_id__primary_site',flat=True).distinct())
    sangercells=list(sangercells)
    id_counter=0
   
    for p in range(0,len(sangerdi)):
        temp=sangersites.count(sangerdi[p])
        sanger.append((sangerdi[p],list(set(sangercells[id_counter:id_counter+temp]))))
        id_counter+=temp
        
    ccle=[]
    ccleprime=Sample.objects.filter(dataset_id__name="GSE36133").order_by('cell_line_id__primary_site').select_related('cell_line_id')
    cclesites=list(ccleprime.values_list('cell_line_id__primary_site',flat=True))
    cclecells=list(ccleprime.values_list('cell_line_id__name',flat=True))
    #ncidi=list(ncicells.values_list('cell_line_id__primary_site',flat=True))
    ccledi=list(ccleprime.values_list('cell_line_id__primary_site',flat=True).distinct())
    cclecells=list(cclecells)
    id_counter=0
   
    for p in range(0,len(ccledi)):
        temp=cclesites.count(ccledi[p])
        ccle.append((ccledi[p],list(set(cclecells[id_counter:id_counter+temp]))))
        id_counter+=temp
        
    check_celllines=list(celllines) #all U133A cell lines
    plus2_celllines=list(ncicelllines)+list(CCcelllines)

    return  {
        'check_celllines': mark_safe(json.dumps(check_celllines)),
        'plus2_celllines': mark_safe(json.dumps(plus2_celllines)),
        'nci':nci,
        'sanger':sanger,
        'ccle':ccle,
        
    }


def user_pca(request):
    text=request.FILES['user_file']
    data = pd.read_csv(text, sep=" ", header = None)
    print(data)
    #data.columns = ["a", "b", "c", "etc."]
    #print(text.read().decode('UTF-8').strip().split(' '))
    #print("===================")
    
    #print("===================")
#notice that we need to return a user_pca_center.html, too!!
    return render_to_response('welcome.html',locals())
def express_profiling(request):
    return render(request, 'express_profiling.html', generate_samples())
   

def welcome(request):
    return render_to_response('welcome.html',locals())

def help_similar_assessment(request):
    return render_to_response('help_similar_assessment.html',RequestContext(request))

def similar_assessment(request):   
    return render(request, 'similar_assessment.html', generate_samples())
    
def gene_signature(request):    
    return render(request, 'gene_signature.html', generate_samples())
    

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
                goffset=list(s.values_list('offset',flat=True))
                offset_group_dict['g'+str(i)]=goffset
        
        #get probe from different platform
        all_probe=ProbeID.objects.filter(platform__name=pform).order_by('offset')
        probe_offset=list(all_probe.values_list('offset',flat=True))
                    
        #deal with "nan" in Sanger dataset
        if pform=="U133A":
            all_probe=list(all_probe)
            for x in range(22275,22281):
                all_probe.pop(probe_offset.index(x))
                probe_offset.remove(x)
                
            
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
        all_probe=ProbeID.objects.filter(platform__name=pform).order_by('offset')
        probe_offset=list(all_probe.values_list('offset',flat=True))
        
        all_probe=list(all_probe)
        #deal with "nan" in Sanger dataset
        if pform=="U133A":
            #all_probe=list(all_probe)
            for x in range(22275,22281):
                all_probe.pop(probe_offset.index(x))
                probe_offset.remove(x)
        
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
                    s=Sample.objects.filter(cell_line_id__name__in=SANGER,platform_id__name=pform).order_by('dataset_id').select_related('cell_line_id__name','dataset_id')
                    goffset=list(s.values_list('offset',flat=True))
                    sanger_data=sanger_val[np.ix_(probe_offset,list(goffset))]
                    s_group_dict['g'+str(i)]=s
                    val.append(sanger_data.tolist())
            else:
                cnci='select_nci_g'+str(i)
                cgse='select_ccle_g'+str(i)
                s_nci=[]
                s_gse=[]
                nci_flag=0
                gse_flag=0
                if 'NCI60' in datasets:
                    nci_flag=1
                    NCI=request.POST.getlist(cnci)
                    s_nci=Sample.objects.filter(cell_line_id__name__in=NCI,dataset_id__name__in=['NCI60']).order_by('dataset_id').select_related('cell_line_id__name','dataset_id')
                    goffset=list(s_nci.values_list('offset',flat=True))
                    nci_data=nci_val[np.ix_(probe_offset,list(goffset))]
                if 'GSE36133' in datasets:
                    gse_flag=1
                    GSE=request.POST.getlist(cgse)
                    s_gse=Sample.objects.filter(cell_line_id__name__in=GSE,dataset_id__name__in=['GSE36133']).order_by('dataset_id').select_related('cell_line_id__name','dataset_id')
                    goffset=list(s_gse.values_list('offset',flat=True))
                    gse_data=gse_val[np.ix_(probe_offset,list(goffset))]
                
                #append nci60 and gse36133 as g_data
               
                s_group_dict['g'+str(i)]=list(s_nci)+list(s_gse)
                if(nci_flag==1 and gse_flag==1):
                    g_data=np.concatenate((nci_data,gse_data),axis=1)
                elif(nci_flag==1):
                    g_data=nci_data
                else:
                    g_data=gse_data
                val.append(g_data.tolist())
            
            
            
    #run the one way ANOVA test or ttest for every probe base on the platform selected    
    express={}
    
    if group_counter<=2:
        for i in range(0,len(all_probe)):    #need to fix if try to run on laptop
            presult[all_probe[i]]=stats.ttest_ind(list(val[0][i]),list(val[1][i]),equal_var=False,nan_policy='omit')[1]
            express[all_probe[i]]=np.append(val[0][i],val[1][i]).tolist()
    else:
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
                sample_out.append(s.cell_line_id.name+"(SCLP)(group"+str(n_counter)+"-"+str(sample_counter)+")")   
            else:
                sample_out.append(s.cell_line_id.name+"("+s.dataset_id.name+")"+"(group"+str(n_counter)+"-"+str(sample_counter)+")")   
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
    
    ##open all the file
    sanger_val_pth=Path('../').resolve().joinpath('src','sanger_cell_line_proj.npy')
    nci_val_pth=Path('../').resolve().joinpath('src','nci60.npy')
    gse_val_pth=Path('../').resolve().joinpath('src','GSE36133.npy')
    sanger_val=np.load(sanger_val_pth.as_posix(),mmap_mode='r')
    nci_val=np.load(nci_val_pth.as_posix(),mmap_mode='r')
    gse_val=np.load(gse_val_pth.as_posix(),mmap_mode='r')
    
    propotion=0
    table_propotion=0
    pform=request.POST['data_platform'] #get the platform
    show=request.POST['show_type']      #get the pca show type
    nci_size=Sample.objects.filter(dataset_id__name__in=["NCI60"]).count()
    group_counter=1
    s_group_dict={}  #store sample
    offset_group_dict={} #store offset
    cell_line_dict={}
    if request.POST['cell_line_method'] == 'text':
        
        #count how many group
        group_counter=1
        while True:
            temp_name='cellline_g'+str(group_counter)
            if temp_name in request.POST:
                group_counter=group_counter+1
            else:
                group_counter=group_counter-1
                break
    
        group_name=[]
        s_group_dict={}  #store sample
        offset_group_dict={} #store offset
        gse_flag=0
        nci_flag=0
        all_c=[]
        for i in range(1,group_counter+1):
            c='cellline_g'+str(i)
            if request.POST[c] !='':
                temp_name='g'+str(i)
                group_name.append(temp_name)
                
                temp=list(set(request.POST[c].split()))
                if all_c==[]:
                    all_c=all_c+temp
                    uni=temp
                else:
                    uni=list(set(temp)-set(all_c))
                    all_c=all_c+uni
                    
                s=Sample.objects.filter(cell_line_id__name__in=(uni),platform_id__name=pform).order_by('dataset_id'
                ).select_related('cell_line_id__name','cell_line_id__primary_site','cell_line_id__primary_hist','dataset_id','dataset_id__name')
                s_group_dict['g'+str(i)]=s
                goffset=list(s.values_list('offset',flat=True))
                offset_group_dict['g'+str(i)]=goffset
                
                cell_line_dict['g'+str(i)]=list(s.values_list('cell_line_id__name',flat=True))
                
                #deal with offset, because we have to combine u133plus2 data together PROBLEM!!!sample
                gseindex=-1
                if pform=="PLUS2":
                    listg=list(s_group_dict[temp_name].values('dataset_id'))
                    if {'dataset_id':3} in listg:
                        nci_flag=1
                    if {'dataset_id':3} in listg:
                        gse_flag=1
                        gseindex=listg.index({'dataset_id':3})
                        offset_group_dict['g'+str(i)]=offset_group_dict['g'+str(i)][:gseindex]+list(np.add(offset_group_dict['g'+str(i)][gseindex:],nci_size))
                
                s_group_dict['g'+str(i)]=list(s)        
        
        
    else:
        #this part is for selecting cell lines base on dataset
        #count how many group
        group_counter=1
        while True:
            temp_name='dataset_g'+str(group_counter)
            if temp_name in request.POST:
                group_counter=group_counter+1
            else:
                group_counter=group_counter-1
                break
    
        s_group_dict={}  #store sample
        group_name=[]
        offset_group_dict={} #store offset
        gse_flag=0
        nci_flag=0
        all_c=[]
        all_c_nci=[]
        all_c_gse=[]
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
                    temp=list(set(SANGER))
                    if all_c==[]:
                        all_c=all_c+temp
                        uni=temp
                    else:
                        uni=list(set(temp)-set(all_c))
                        all_c=all_c+uni
                    s=Sample.objects.filter(cell_line_id__name__in=uni,platform_id__name=pform).order_by('dataset_id'
                    ).select_related('cell_line_id__name','cell_line_id__primary_site','cell_line_id__primary_hist','dataset_id','dataset_id__name')
                    goffset=list(s.values_list('offset',flat=True))
                    s_group_dict['g'+str(i)]=list(s)
                    offset_group_dict['g'+str(i)]=goffset
                    cell_line_dict['g'+str(i)]=list(s.values_list('cell_line_id__name',flat=True))
            else:
                cnci='select_nci_g'+str(i)
                cgse='select_ccle_g'+str(i)
                s_nci=[]
                s_gse=[]
                cell_gse=[]
                cell_nci=[]
                goffset_nci=[]
                goffset_gse=[]
                if 'NCI60' in datasets:
                    nci_flag=1
                    NCI=request.POST.getlist(cnci)
                    
                    
                    temp_nci=list(set(NCI))
                    if 'd_sample' in show:
                        if all_c_nci==[]:
                            all_c_nci=all_c_nci+temp_nci
                            uni_nci=temp_nci
                        else:
                            uni_nci=list(set(temp_nci)-set(all_c_nci))
                            all_c_nci=all_c_nci+uni_nci
                    else:
                        uni_nci=list(temp_nci)          #do not filter duplicate input only when select+centroid
                    s_nci=Sample.objects.filter(cell_line_id__name__in=uni_nci,dataset_id__name__in=['NCI60']).order_by('dataset_id'
                    ).select_related('cell_line_id__name','cell_line_id__primary_site','cell_line_id__primary_hist','dataset_id','dataset_id__name')
                    goffset_nci=list(s_nci.values_list('offset',flat=True))
                    cell_nci=list(s_nci.values_list('cell_line_id__name',flat=True))
                if 'GSE36133' in datasets:
                    gse_flag=1
                    GSE=request.POST.getlist(cgse)
                    temp_gse=list(set(GSE))
                    
                    if 'd_sample' in show:
                        if all_c_gse==[]:
                            all_c_gse=all_c_gse+temp_gse
                            uni_gse=temp_gse
                        else:
                            uni_gse=list(set(temp_gse)-set(all_c_gse))
                            all_c_gse=all_c_gse+uni_gse
                    else:
                        uni_gse=list(temp_gse)
                    s_gse=Sample.objects.filter(cell_line_id__name__in=uni_gse,dataset_id__name__in=['GSE36133']).order_by('dataset_id'
                    ).select_related('cell_line_id__name','cell_line_id__primary_site','cell_line_id__primary_hist','dataset_id','dataset_id__name')
                    cell_gse=list(s_gse.values_list('cell_line_id__name',flat=True))
                    if(nci_flag==1):
                        goffset_gse=list(np.add(list(s_gse.values_list('offset',flat=True)),nci_size))
                    else:
                        goffset_gse=list(s_gse.values_list('offset',flat=True))
                
                
                #append nci60 and gse36133 as g_data
                s_group_dict['g'+str(i)]=list(s_nci)+list(s_gse)
                offset_group_dict['g'+str(i)]=goffset_nci+goffset_gse
                cell_line_dict['g'+str(i)]=cell_nci+cell_gse
                

    #delete nan, transpose matrix
    if pform=="U133A":
        sanger_val=sanger_val[~np.isnan(sanger_val).any(axis=1)]
        sanger_val=np.matrix(sanger_val)[:,:]#need fix
        val=np.transpose(sanger_val)
        #dataset_size=[0,Sample.objects.filter(dataset_id__name__in=["Sanger Cell Line Project"]).count()]
    elif((nci_flag==1) and (gse_flag==1)):
        ##PROBLEM:should we always use combined one to run pca?
        nci_val=np.matrix(nci_val)[:,:] #need fix
        tnci_val=np.transpose(nci_val)
        gse_val=np.matrix(gse_val)[:,:] #need fix
        tgse_val=np.transpose(gse_val)
        val=np.concatenate((tnci_val, tgse_val))#combine together  
        #dataset_size=[0,nci_size,nci_size+Sample.objects.filter(dataset_id__name__in=["GSE36133"]).count()]
    elif nci_flag==1:
        nci_val=np.matrix(nci_val)[:,:] #need fix
        val=np.transpose(nci_val)
    else:
        gse_val=np.matrix(gse_val)[:,:] #need fix
        val=np.transpose(gse_val)
        
    
    
    all_sample=[]
    all_cellline=[]
    cell_object=[]
    all_offset=[]
    sample_counter={}
    group_cell=[]
    g_s_counter=[0]
    
    for i in range(1,group_counter+1):
        all_sample=all_sample+s_group_dict['g'+str(i)] #will not exist duplicate sample if d_sample
        all_offset=all_offset+offset_group_dict['g'+str(i)]
        all_cellline=all_cellline+cell_line_dict['g'+str(i)]
        g_s_counter.append(g_s_counter[i-1]+len(s_group_dict['g'+str(i)]))
    
    for i in all_sample:
        sample_counter[i.name]=1
        cell_object.append(i.cell_line_id)
    
     
          
    pca_index=[]
    dis_offset=[]
      

    
    #PREMISE:same dataset same cell line will have only one type of primary site and primary histology
    name1=[]
    name2=[]
    name3=[]
    name4=[]
    name5=[]
    X1=[]
    Y1=[]
    Z1=[]
    X2=[]
    Y2=[]
    Z2=[]
    X3=[]
    Y3=[]
    Z3=[]
    X4=[]
    Y4=[]
    Z4=[]
    X5=[]
    Y5=[]
    Z5=[]
    n=3  #need to fix to the best one #need to fix proportion 
    
    if 'd_sample' in show:
        #count the pca first
        pca= PCA(n_components=n)
        Xval = pca.fit_transform(val[all_offset,:])  #cannot get Xval with original offset any more
        ratio_temp=pca.explained_variance_ratio_
        propotion=sum(ratio_temp[0:3])
        table_propotion=sum(ratio_temp[0:n+1])
        print(Xval)
        
        max=0
        min=10000000000
        out_group=[]
        exist_cell={}#cell line object:counter
        for g in range(1,group_counter+1):
              
            output_cell={}
            check={}
            for s in range(g_s_counter[g-1],g_s_counter[g]):
                  
                cell=all_sample[s].cell_line_id
                try: 
                    counter=exist_cell[cell]
                    exist_cell[cell]=counter+1
                    
                except KeyError:
                    exist_cell[cell]=1
                
                try:
                    t=output_cell[cell]
                except KeyError:
                    output_cell[cell]=[cell,[]]
                
                check[all_sample[s].name]=[]
                sample_counter[all_sample[s].name]=exist_cell[cell]    
                for i in range(0,len(all_sample)):
                    if i!=s:
                        try:
                            if(all_sample[s].name not in check[all_sample[i].name]):
                                distance=np.linalg.norm(Xval[i]-Xval[s])
                                if distance<min:
                                    min=distance
                                if distance>max:
                                    max=distance
                                output_cell[cell][1].append([all_cellline[s]+'('+str(exist_cell[cell])+')'
                                ,all_sample[s].name,all_sample[s].dataset_id.name,all_cellline[i],all_sample[i].name,all_sample[i].dataset_id.name,distance,cell_object[i]])
                                check[all_sample[s].name].append(all_sample[i].name)
                        except KeyError:
                            distance=np.linalg.norm(Xval[i]-Xval[s])
                            if distance<min:
                                min=distance
                            if distance>max:
                                max=distance
                            output_cell[cell][1].append([all_cellline[s]+'('+str(exist_cell[cell])+')'
                            ,all_sample[s].name,all_sample[s].dataset_id.name,all_cellline[i],all_sample[i].name,all_sample[i].dataset_id.name,distance,cell_object[i]])
                            check[all_sample[s].name].append(all_sample[i].name)
                        
                if(g==1):
                    name1.append(all_cellline[s]+'('+str(exist_cell[cell])+')'+'<br>'+all_sample[s].name)
                    X1.append(Xval[s][0])
                    Y1.append(Xval[s][1])
                    Z1.append(Xval[s][2])
                elif(g==2):
                    name2.append(all_cellline[s]+'('+str(exist_cell[cell])+')'+'<br>'+all_sample[s].name)
                    X2.append(Xval[s][0])
                    Y2.append(Xval[s][1])
                    Z2.append(Xval[s][2])
                elif(g==3):
                    name3.append(all_cellline[s]+'('+str(exist_cell[cell])+')'+'<br>'+all_sample[s].name)
                    X3.append(Xval[s][0])
                    Y3.append(Xval[s][1])
                    Z3.append(Xval[s][2])
                elif(g==4):
                    name4.append(all_cellline[s]+'('+str(exist_cell[cell])+')'+'<br>'+all_sample[s].name)
                    X4.append(Xval[s][0])
                    Y4.append(Xval[s][1])
                    Z4.append(Xval[s][2])
                elif(g==5):
                    name5.append(all_cellline[s]+'('+str(exist_cell[cell])+')'+'<br>'+all_sample[s].name)
                    X5.append(Xval[s][0])
                    Y5.append(Xval[s][1])
                    Z5.append(Xval[s][2])
            dictlist=[]
            for key, value in output_cell.items():
                temp = [value]
                dictlist+=temp
            output_cell=list(dictlist)
            out_group.append([g,output_cell])
        
        #[g,[[group_cell_line,[paired_cellline,......,]],[],[]]]
        for i in out_group:
            for temp_list in i[1]:
                for temp in temp_list[1]:
                    #print(temp)
                    temp[3]=temp[3]+'('+str(sample_counter[temp[4]])+')'
        
        return_html='pca.html'
    else:
        #This part is for centroid display
        return_html='pca_center.html'
        #deal with text part first, get all cell line base on platform instead of dataset--->different group need to filter same cell line name first
        #count the centroid--->use this new data to run pca--->new location to count distance
        if request.POST['cell_line_method'] == 'text':
            dis_cellline=list(set(cell_object))
            #location_dict={} #{cell object:new location}
            dataset_dict={}  #{cell object:dataset combined}
            a_cell_object=np.array(cell_object)
            X_val=[]
            
            for c in dis_cellline:
                total_offset=np.where(a_cell_object==c)[0]
                val_a=np.array(val)
                a_all_offset=np.array(all_offset)
                selected_val=val_a[a_all_offset[total_offset]]
                new_loca=(np.mean(selected_val,axis=0,dtype=np.float64,keepdims=True)).tolist()[0]
                #location_dict[c]=new_loca
                X_val.append(new_loca)   #in the order of dis_cellline
                
                a_sample=np.array(all_sample)
                selected_sample=a_sample[total_offset]
                
                for s in selected_sample:
                    
                    dataset=s.dataset_id.name
                    try:
                        sets=dataset_dict[c]
                        if (("NCI60" in sets) and ("GSE36133" in sets)):
                            break
                        if ("Sanger Cell Line Project" in sets):
                            break
                        if( dataset not in sets):
                            dataset_dict[c]=dataset+"/"+sets
                    except KeyError: 
                        dataset_dict[c]=dataset
            
            
            #run the pca again here and store it with new offset to get the new pca data
            X_val=np.matrix(X_val)
            pca= PCA(n_components=3)
            new_val = pca.fit_transform(X_val[:,:])  #cannot get Xval with original offset any more
            ratio_temp=pca.explained_variance_ratio_
            propotion=sum(ratio_temp[0:3])
            table_propotion=sum(ratio_temp[0:n+1])
            print(new_val)
            
            out_group=[]
            min=10000000000
            max=0
            #count distance base on X_val
            for g in range(1,group_counter+1):
                output_cell=[]
                exist_cell=[]
                check={} #to remove A-B and B-A
                for s in range(g_s_counter[g-1],g_s_counter[g]):
                    cell=all_sample[s].cell_line_id
                    index_cell=np.where(np.array(dis_cellline)==cell)[0][0]
                    if (cell not in exist_cell):
                        output_cell.append([cell,[]])
                        check[cell]=[]  
                        #count the distance
                        for c in dis_cellline:
                            if c != cell:
                                index_c=np.where(np.array(dis_cellline)==c)[0][0]
                                
                                try:
                                    if(cell not in check[c]):
                                        distance=np.linalg.norm(np.array(new_val[index_cell])-np.array(new_val[index_c]))
                                        if distance<min:
                                            min=distance
                                        if distance>max:
                                            max=distance
                                        output_cell[len(output_cell)-1][1].append([cell,dataset_dict[cell],c,dataset_dict[c],distance])
                                        check[cell].append(c)
                                except KeyError:
                                    
                                    distance=np.linalg.norm(np.array(new_val[index_cell])-np.array(new_val[index_c]))
                                    if distance<min:
                                        min=distance
                                    if distance>max:
                                        max=distance
                                    output_cell[len(output_cell)-1][1].append([cell,dataset_dict[cell],c,dataset_dict[c],distance])
                                    check[cell].append(c)
                                                    
                        exist_cell.append(cell) 
                        
                        
                        if(g==1):
                            name1.append(cell.name+'<br>'+dataset_dict[cell])
                            X1.append(new_val[index_cell][0])
                            Y1.append(new_val[index_cell][1])
                            Z1.append(new_val[index_cell][2])
                        elif(g==2):
                            name2.append(cell.name+'<br>'+dataset_dict[cell])
                            X2.append(new_val[index_cell][0])
                            Y2.append(new_val[index_cell][1])
                            Z2.append(new_val[index_cell][2])
                        elif(g==3):
                            name3.append(cell.name+'<br>'+dataset_dict[cell])
                            X3.append(new_val[index_cell][0])
                            Y3.append(new_val[index_cell][1])
                            Z3.append(new_val[index_cell][2])
                        elif(g==4):
                            name4.append(cell.name+'<br>'+dataset_dict[cell])
                            X4.append(new_val[index_cell][0])
                            Y4.append(new_val[index_cell][1])
                            Z4.append(new_val[index_cell][2])
                        elif(g==5):
                            name5.append(cell.name+'<br>'+dataset_dict[cell])
                            X5.append(new_val[index_cell][0])
                            Y5.append(new_val[index_cell][1])
                            Z5.append(new_val[index_cell][2]) 
                out_group.append([g,output_cell]) 
                           

        else:
        #This part is for select cell line base on dataset,count centroid base on the dataset
        #group中的cell line為單位來算重心
        
            for x in range(0,len(all_sample)):  #delete duplicate offset to prevent pca error
            
                if all_offset[x] not in dis_offset:
                    dis_offset.append(all_offset[x])
            
            for x in all_offset:
                pca_index.append(dis_offset.index(x))
        
        
            location_dict={} #{group number:[[cell object,dataset,new location]]}
            combined=[]
            sample_list=[]
            pca_index=np.array(pca_index)
            X_val=[]
            val_a=np.array(val)
            a_all_offset=np.array(all_offset)
            for i in range(1,group_counter+1):
                dis_cellline=list(set(cell_object[g_s_counter[i-1]:g_s_counter[i]]))  #cell object may have duplicate cell line since:NCI A + CCLE A===>[A,A]
                location_dict['g'+str(i)]=[]
                dataset_dict={}
                a_cell_object=np.array(cell_object)
                
                for c in dis_cellline:    #dis_cellline may not have the same order as cell_object
                    
                    temp1=np.where((a_cell_object==c))[0]
                    
                    temp2=np.where((temp1>=g_s_counter[i-1])&(temp1<g_s_counter[i]))
                    total_offset=temp1[temp2]
                    selected_val=val_a[a_all_offset[total_offset]]
                    new_loca=(np.mean(selected_val,axis=0,dtype=np.float64,keepdims=True)).tolist()[0]
                    
                    
                    a_sample=np.array(all_sample)
                    selected_sample=a_sample[total_offset]
                    
                    if list(selected_sample) in sample_list:   #to prevent two different colors in different group
                        continue
                    else:
                        sample_list.append(list(selected_sample))
                    
                    
                    #print(selected_sample)
                    for s in selected_sample:
                        
                        dataset=s.dataset_id.name
                        try:
                            sets=dataset_dict[c]
                            if(nci_flag==1 and gse_flag==1):
                                if (("NCI60" in sets) and ("GSE36133" in sets)):
                                    break
                            elif(nci_flag==1 and gse_flag==0):
                                if("NCI60" in sets):
                                    break
                            elif(gse_flag==1 and nci_flag==0):
                                if("GSE36133" in sets):
                                    break
                            elif ("Sanger Cell Line Project" in sets):
                                break
                            if( dataset not in sets):
                                dataset_dict[c]=dataset+"/"+sets
                        except KeyError: 
                            dataset_dict[c]=dataset
                    X_val.append(new_loca)
                    location_dict['g'+str(i)].append([c,dataset_dict[c],len(X_val)-1])  #the last part is the index to get pca result from new_val
                    combined.append([c,dataset_dict[c],len(X_val)-1])  #all cell line, do not matter order
            
            #run the pca
            X_val=np.matrix(X_val)
            pca= PCA(n_components=3)
            new_val = pca.fit_transform(X_val[:,:])  #cannot get Xval with original offset any more
            ratio_temp=pca.explained_variance_ratio_
            propotion=sum(ratio_temp[0:3])
            table_propotion=sum(ratio_temp[0:n+1])
            print(new_val)




            
            out_group=[]
            min=10000000000
            max=0
            for g in range(1,group_counter+1):
                output_cell=[]
                exist_cell={}
                
                for group_c in location_dict['g'+str(g)]:  #a list of [c,dataset_dict[c],new_val index] in group one
                    cell=group_c[0]
                    key_string=cell.name+'/'+cell.primary_site+'/'+cell.primary_hist+'/'+group_c[1]
                    exist_cell[key_string]=[]
                    output_cell.append([cell,[]])
                    
                    #count the distance
                    for temp_list in combined:
                        c=temp_list[0]
                        temp_string=c.name+'/'+c.primary_site+'/'+c.primary_hist+'/'+temp_list[1]
                        try:
                            if(key_string not in exist_cell[temp_string]):
                                distance=np.linalg.norm(np.array(new_val[group_c[2]])-np.array(new_val[temp_list[2]]))
                                if distance==0:
                                    continue
                                if distance<min:
                                    min=distance
                                if distance>max:
                                    max=distance
                                output_cell[len(output_cell)-1][1].append([cell,group_c[1],temp_list[0],temp_list[1],distance])
                                exist_cell[key_string].append(temp_string)
                        except KeyError:
                            distance=np.linalg.norm(np.array(new_val[group_c[2]])-np.array(new_val[temp_list[2]]))
                            if distance==0:
                                continue
                            if distance<min:
                                min=distance
                            if distance>max:
                                max=distance
                            output_cell[len(output_cell)-1][1].append([cell,group_c[1],temp_list[0],temp_list[1],distance])
                            exist_cell[key_string].append(temp_string)
                    if(g==1):
                        name1.append(cell.name+'<br>'+group_c[1])
                        X1.append(new_val[group_c[2]][0])
                        Y1.append(new_val[group_c[2]][1])
                        Z1.append(new_val[group_c[2]][2])
                    elif(g==2):
                        name2.append(cell.name+'<br>'+group_c[1])
                        X2.append(new_val[group_c[2]][0])
                        Y2.append(new_val[group_c[2]][1])
                        Z2.append(new_val[group_c[2]][2])
                    elif(g==3):
                        name3.append(cell.name+'<br>'+group_c[1])
                        X3.append(new_val[group_c[2]][0])
                        Y3.append(new_val[group_c[2]][1])
                        Z3.append(new_val[group_c[2]][2])
                    elif(g==4):
                        name4.append(cell.name+'<br>'+group_c[1])
                        X4.append(new_val[group_c[2]][0])
                        Y4.append(new_val[group_c[2]][1])
                        Z4.append(new_val[group_c[2]][2])
                    elif(g==5):
                        name5.append(cell.name+'<br>'+group_c[1])
                        X5.append(new_val[group_c[2]][0])
                        Y5.append(new_val[group_c[2]][1])
                        Z5.append(new_val[group_c[2]][2]) 
                out_group.append([g,output_cell])
                    
    return render_to_response(return_html,RequestContext(request,
    {
    'min':min,'max':max,
    'out_group':out_group,
    'propotion':propotion,
    'table_propotion':table_propotion,
    'X1':X1,'name1':mark_safe(json.dumps(name1)),
    'Y1':Y1,'name2':mark_safe(json.dumps(name2)),
    'Z1':Z1,'name3':mark_safe(json.dumps(name3)),
    'X2':X2,'name4':mark_safe(json.dumps(name4)),
    'Y2':Y2,'name5':mark_safe(json.dumps(name5)),
    'Z2':Z2,
    'X3':X3,
    'Y3':Y3,
    'Z3':Z3,
    'X4':X4,
    'Y4':Y4,
    'Z4':Z4,
    'X5':X5,
    'Y5':Y5,
    'Z5':Z5,
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
    '''
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
    '''
    nci=[]
    nciprime=Sample.objects.filter(dataset_id__name="NCI60").order_by('cell_line_id__primary_site').select_related('cell_line_id')
    ncisites=list(nciprime.values_list('cell_line_id__primary_site',flat=True))
    ncicells=list(nciprime.values_list('cell_line_id__name',flat=True))
    #ncidi=list(ncicells.values_list('cell_line_id__primary_site',flat=True))
    ncidi=list(nciprime.values_list('cell_line_id__primary_site',flat=True).distinct())
    ncicells=list(ncicells)
    id_counter=0
   
    for p in range(0,len(ncidi)):
        temp=ncisites.count(ncidi[p])
        nci.append((ncidi[p],list(set(ncicells[id_counter:id_counter+temp]))))
        id_counter+=temp
        
    sanger=[]
    sangerprime=Sample.objects.filter(dataset_id__name="Sanger Cell Line Project").order_by('cell_line_id__primary_site').select_related('cell_line_id')
    sangersites=list(sangerprime.values_list('cell_line_id__primary_site',flat=True))
    sangercells=list(sangerprime.values_list('cell_line_id__name',flat=True))
    #ncidi=list(ncicells.values_list('cell_line_id__primary_site',flat=True))
    sangerdi=list(sangerprime.values_list('cell_line_id__primary_site',flat=True).distinct())
    sangercells=list(sangercells)
    id_counter=0
   
    for p in range(0,len(sangerdi)):
        temp=sangersites.count(sangerdi[p])
        sanger.append((sangerdi[p],list(set(sangercells[id_counter:id_counter+temp]))))
        id_counter+=temp
        
    ccle=[]
    ccleprime=Sample.objects.filter(dataset_id__name="GSE36133").order_by('cell_line_id__primary_site').select_related('cell_line_id')
    cclesites=list(ccleprime.values_list('cell_line_id__primary_site',flat=True))
    cclecells=list(ccleprime.values_list('cell_line_id__name',flat=True))
    #ncidi=list(ncicells.values_list('cell_line_id__primary_site',flat=True))
    ccledi=list(ccleprime.values_list('cell_line_id__primary_site',flat=True).distinct())
    cclecells=list(cclecells)
    id_counter=0
   
    for p in range(0,len(ccledi)):
        temp=cclesites.count(ccledi[p])
        ccle.append((ccledi[p],list(set(cclecells[id_counter:id_counter+temp]))))
        id_counter+=temp
    
    
    

    return render(request, 'cellline_microarray.html', {
        'nci': nci,
        'sanger':sanger,
        'ccle':ccle,
        'samples': samples,
        #'primary_sites': primary_sites,
        'ncisamples': ncisamples,
        #'nciprimary_sites': nciprimary_sites,
        'CCsamples': CCsamples,
        #'CCprimary_sites': CCprimary_sites,
        
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
    samples=Sample.objects.all().select_related('cell_line_id','dataset_id').order_by('id') 
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
        c = list(set(c.split()))
        sanger_flag=1
        samples=Sample.objects.filter(dataset_id__name__in=['Sanger Cell Line Project']).order_by('id')      
        cell=samples.select_related('cell_line_id','dataset_id').filter(cell_line_id__name__in=c).order_by('id')
        offset=list(cell.values_list('offset',flat=True))
        ps_id='1'
        
        nci_flag=1
        ncisamples=Sample.objects.filter(dataset_id__name__in=['NCI60']).select_related('cell_line_id','dataset_id').order_by('id') 
        ncicell=ncisamples.filter(cell_line_id__name__in=c).order_by('id') 
        ncioffset=list(ncicell.values_list('offset',flat=True))
        pn_id='3'
        
        gse_flag=1
        CCsamples=Sample.objects.filter(dataset_id__name__in=['GSE36133']).select_related('cell_line_id','dataset_id').order_by('id') 
        CCcell=CCsamples.filter(cell_line_id__name__in=c).order_by('id') 
        CCoffset=list(CCcell.values_list('offset',flat=True))
        pn_id='3'
    else:
        if 'dataset' in request.POST and request.POST['dataset'] != '':
            datas=request.POST.getlist('dataset')
            if 'Sanger Cell Line Project' in datas:
                sanger_flag=1
                SANGER=list(set(request.POST.getlist('select_sanger')))
                samples=Sample.objects.filter(dataset_id__name__in=['Sanger Cell Line Project']).order_by('id')       
                cell=samples.select_related('cell_line_id','dataset_id').filter(cell_line_id__name__in=SANGER).order_by('id') 
                offset=list(cell.values_list('offset',flat=True))
                ps_id='1'
            if 'NCI60' in datas:
                nci_flag=1
                NCI=list(set(request.POST.getlist('select_nci')))
                ncisamples=Sample.objects.filter(dataset_id__name__in=['NCI60']).select_related('cell_line_id','dataset_id').order_by('id') 
                ncicell=ncisamples.filter(cell_line_id__name__in=NCI).order_by('id') 
                ncioffset=list(ncicell.values_list('offset',flat=True))
                pn_id='3'
            if 'GSE36133' in datas:
                gse_flag=1
                GSE=list(set(request.POST.getlist('select_gse')))
                CCsamples=Sample.objects.filter(dataset_id__name__in=['GSE36133']).select_related('cell_line_id','dataset_id').order_by('id') 
                CCcell=CCsamples.filter(cell_line_id__name__in=GSE).order_by('id') 
                CCoffset=list(CCcell.values_list('offset',flat=True))
                pn_id='3'
            if len(SANGER)==0 and len(NCI)==0 and len(GSE)==0:
                return HttpResponse("<p>please select primary sites.</p>" )
        else:
            return HttpResponse("<p>please check Step3 again.</p>" )
        
    
    
    if 'keyword' in request.POST and request.POST['keyword'] != '':
        words = request.POST['keyword']
        words = list(set(words.split()))
    else:
        return HttpResponse("<p>where is your keyword?</p>")
      
    #open files
    sanger_val_pth=Path('../').resolve().joinpath('src','sanger_cell_line_proj.npy')
    nci_val_pth=Path('../').resolve().joinpath('src','nci60.npy')
    gse_val_pth=Path('../').resolve().joinpath('src','GSE36133.npy')
    sanger_val=np.load(sanger_val_pth.as_posix(),mmap_mode='r')
    nci_val=np.load(nci_val_pth.as_posix(),mmap_mode='r')
    gse_val=np.load(gse_val_pth.as_posix(),mmap_mode='r')
    
    u133a_rank=np.load('ranking_u133a.npy')
    plus2_rank=np.load('ranking_u133plus2.npy')
    
    gene = []
    ncigene = []
    CCgene = []
    context={}
    
    norm_name=[request.POST['normalize']]
    if sanger_flag==1:
        #if request.POST['normalize']!='NTRK3-AS1':
        sanger_g=ProbeID.objects.filter(platform__in=ps_id).filter(Gene_symbol__in=norm_name).order_by('id') 
        sanger_probe_offset=list(sanger_g.values_list('offset',flat=True))
        temp=sanger_val[np.ix_(sanger_probe_offset,offset)]
        norm=np.mean(temp,axis=0, dtype=np.float64,keepdims=True)
        #else:
        #    norm=0.0
    else:
        norm=0.0  #if / should = 1   
    if nci_flag==1:
        nci_g=ProbeID.objects.filter(platform__in=pn_id).filter(Gene_symbol__in=norm_name).order_by('id') 
        nci_probe_offset=list(nci_g.values_list('offset',flat=True))
        temp=nci_val[np.ix_(nci_probe_offset,ncioffset)]
        nci_norm=np.mean(temp,axis=0, dtype=np.float64,keepdims=True)
        #print(nci_norm)   
    else:
        nci_norm=0.0  #if / should = 1
    if gse_flag==1:
        CC_g=ProbeID.objects.filter(platform__in=pn_id).filter(Gene_symbol__in=norm_name).order_by('id') 
        CC_probe_offset=list(CC_g.values_list('offset',flat=True))
        temp=gse_val[np.ix_(CC_probe_offset,CCoffset)]
        CC_norm=np.mean(temp,axis=0, dtype=np.float64,keepdims=True)
        #print(CC_norm)
    else:
        CC_norm=0.0  #if / should = 1             

        
            
    #dealing with probes    
    if 'gtype' in request.POST and request.POST['gtype'] == 'probeid':
        gene = ProbeID.objects.filter(platform__in=ps_id).filter(Probe_id__in=words).order_by('id') 
        probe_offset=list(gene.values_list('offset',flat=True))
        
        ncigene = ProbeID.objects.filter(platform__in=pn_id).filter(Probe_id__in=words).order_by('id') 
        nciprobe_offset=list(ncigene.values_list('offset',flat=True))
        #nci60 and ccle use same probe set(ncigene) and nicprobe
        
        # Make a generator to generate all (cell, probe, val) pairs
        if(len(gene)!=0 and len(cell)!=0):
            raw_test=sanger_val[np.ix_(probe_offset,offset)]
            normalize=np.subtract(raw_test,norm)#dimension different!!!!
            #normalize=np.around(normalize, decimals=1)
            cell_probe_val_pairs = (
                (c, p, raw_test[probe_ix, cell_ix],22277-np.where(u133a_rank==raw_test[probe_ix, cell_ix])[0],normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(gene)
                for cell_ix, c in enumerate(cell)
            )
            #because of nan, we only use 22277
        else:
            cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(ncicell)!=0):
            nci_raw_test=nci_val[np.ix_(nciprobe_offset,ncioffset)]
            nci_normalize=np.subtract(nci_raw_test,nci_norm)
            nci_cell_probe_val_pairs = (
                (c, p, nci_raw_test[probe_ix, cell_ix],54676-np.where(plus2_rank==nci_raw_test[probe_ix, cell_ix])[0],nci_normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(ncigene)
                for cell_ix, c in enumerate(ncicell)
            )
            
        else:
            nci_cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(CCcell)!=0):
            CC_raw_test=gse_val[np.ix_(nciprobe_offset,CCoffset)]
            CC_normalize=np.subtract(CC_raw_test,CC_norm)
            CC_cell_probe_val_pairs = (
                (c, p, CC_raw_test[probe_ix, cell_ix],54676-np.where(plus2_rank==CC_raw_test[probe_ix, cell_ix])[0],CC_normalize[probe_ix, cell_ix])                        
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
        gene = ProbeID.objects.filter(platform__in=ps_id).filter(Gene_symbol__in=words).order_by('id') 
        probe_offset=gene.values_list('offset',flat=True)
        
        ncigene = ProbeID.objects.filter(platform__in=pn_id).filter(Gene_symbol__in=words).order_by('id') 
        nciprobe_offset=ncigene.values_list('offset',flat=True)
        #nci60 and ccle use same probe set(ncigene) and nicprobe
        
        # Make a generator to generate all (cell, probe, val) pairs
        if(len(gene)!=0 and len(cell)!=0):
            raw_test=sanger_val[np.ix_(probe_offset,offset)]
            normalize=np.subtract(raw_test,norm)
            cell_probe_val_pairs = (
                (c, p, raw_test[probe_ix, cell_ix],22277-np.where(u133a_rank==raw_test[probe_ix, cell_ix])[0],normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(gene)
                for cell_ix, c in enumerate(cell)
            )
            
        else:
            cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(ncicell)!=0):
            nci_raw_test=nci_val[np.ix_(nciprobe_offset,ncioffset)]
            nci_normalize=np.subtract(nci_raw_test,nci_norm)
            nci_cell_probe_val_pairs = (
                (c, p, nci_raw_test[probe_ix, cell_ix],54676-np.where(plus2_rank==nci_raw_test[probe_ix, cell_ix])[0],nci_normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(ncigene)
                for cell_ix, c in enumerate(ncicell)
            )
            
        else:
            nci_cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(CCcell)!=0):
            CC_raw_test=gse_val[np.ix_(nciprobe_offset,CCoffset)]
            CC_normalize=np.subtract(CC_raw_test,CC_norm)
            CC_cell_probe_val_pairs = (
                (c, p, CC_raw_test[probe_ix, cell_ix],54676-np.where(plus2_rank==CC_raw_test[probe_ix, cell_ix])[0],CC_normalize[probe_ix, cell_ix])                        
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
        gene = ProbeID.objects.filter(platform__in=ps_id).filter(Entrez_id=words).order_by('id') 
        probe_offset=gene.values_list('offset',flat=True)
        
        ncigene = ProbeID.objects.filter(platform__in=pn_id).filter(Entrez_id__in=words).order_by('id') 
        nciprobe_offset=ncigene.values_list('offset',flat=True)
        #nci60 and ccle use same probe set(ncigene) and nicprobe
        
        # Make a generator to generate all (cell, probe, val) pairs
        if(len(gene)!=0 and len(cell)!=0):
            raw_test=sanger_val[np.ix_(probe_offset,offset)]
            normalize=np.subtract(raw_test,norm)
            cell_probe_val_pairs = (
                (c, p, raw_test[probe_ix, cell_ix],22277-np.where(u133a_rank==raw_test[probe_ix, cell_ix])[0],normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(gene)
                for cell_ix, c in enumerate(cell)
            )
            
        else:
            cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(ncicell)!=0):
            nci_raw_test=nci_val[np.ix_(nciprobe_offset,ncioffset)]
            nci_normalize=np.subtract(nci_raw_test,nci_norm)
            nci_cell_probe_val_pairs = (
                (c, p, nci_raw_test[probe_ix, cell_ix],54676-np.where(plus2_rank==nci_raw_test[probe_ix, cell_ix])[0],nci_normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(ncigene)
                for cell_ix, c in enumerate(ncicell)
            )
            
        else:
            nci_cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(CCcell)!=0):
            CC_raw_test=gse_val[np.ix_(nciprobe_offset,CCoffset)]
            CC_normalize=np.subtract(CC_raw_test,CC_norm)
            CC_cell_probe_val_pairs = (
                (c, p, CC_raw_test[probe_ix, cell_ix],54676-np.where(plus2_rank==CC_raw_test[probe_ix, cell_ix])[0],CC_normalize[probe_ix, cell_ix])                        
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
