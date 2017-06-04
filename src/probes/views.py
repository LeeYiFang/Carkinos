from django.shortcuts import render,render_to_response
from django.http import HttpResponse, Http404,JsonResponse
from django.views.decorators.http import require_GET
from .models import Dataset, CellLine, ProbeID, Sample, Platform, Clinical_Dataset,Clinical_sample,Gene
from django.template import RequestContext
from django.utils.html import mark_safe
import json
import pandas as pd
import numpy as np
from pathlib import Path
import sklearn
from sklearn.decomposition import PCA
from scipy import stats

import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

import uuid
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
r=ro.r
#lumi= importr('lumi')
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import csv
#import logging
#logger = logging.getLogger("__name__")

show_row=4000 #more than how many rows will become download file mode


def generate_samples():
    d=Dataset.objects.all()
    cell_d_name=list(d.values_list('name',flat=True))
    same_name=[]
    cell_datasets=[]  #[[dataset_name,[[primary_site,[cell line]]]]
    
    for i in cell_d_name:
        
        if i=="Sanger Cell Line Project":
            alias='sanger'
            same_name.append('sanger')
        elif i=="NCI60":
            alias='nci'
            same_name.append('nci')
        elif i=="GSE36133":
            alias='ccle'
            same_name.append('ccle')
        else:
            alias=i
            same_name.append(i)
        
        
        sample=Sample.objects.filter(dataset_id__name=i).order_by('cell_line_id__primary_site').select_related('cell_line_id')
        cell_datasets.append([i,alias,list(sample),[]])
        sites=list(sample.values_list('cell_line_id__primary_site',flat=True))
        hists=list(sample.values_list('cell_line_id__name',flat=True))
        
        dis_prim=list(sample.values_list('cell_line_id__primary_site',flat=True).distinct())
        hists=list(hists)
        id_counter=0
       
        for p in range(0,len(dis_prim)):
            temp=sites.count(dis_prim[p])
            cell_datasets[-1][3].append([dis_prim[p],list(set(hists[id_counter:id_counter+temp]))])
            id_counter+=temp
    
    
    d=Clinical_Dataset.objects.all()
    d_name=list(d.values_list('name',flat=True))
    datasets=[]  #[[dataset_name,[[primary_site,[primary_histology]]],[[filter_type,[filter_choice]]]]
    primarys=[]  #[[primary_site,[primary_hist]]]
    primh_filter=[] #[[filter_type,[filter_choice]]]
    f_type=['age','gender','ethnic','grade','stage','stageT','stageN','stageM','metastatic']
    for i in d_name:
        
        same_name.append(i)
        datasets.append([i,[],[]])
        sample=Clinical_sample.objects.filter(dataset_id__name=i).order_by('primary_site')
        sites=list(sample.values_list('primary_site',flat=True))
        hists=list(sample.values_list('primary_hist',flat=True))
        
        dis_prim=list(sample.values_list('primary_site',flat=True).distinct())
        hists=list(hists)
        id_counter=0
       
        for p in range(0,len(dis_prim)):
            temp=sites.count(dis_prim[p])
            datasets[-1][1].append([dis_prim[p],list(set(hists[id_counter:id_counter+temp]))])
            id_counter+=temp
        
        for f in f_type:
            temp=list(set(sample.values_list(f,flat=True)))
            datasets[-1][2].append([f,temp])
    
    sample=Clinical_sample.objects.all().order_by('primary_site')
    sites=list(sample.values_list('primary_site',flat=True))
    hists=list(sample.values_list('primary_hist',flat=True))
    dis_prim=list(sample.values_list('primary_site',flat=True).distinct())
    hists=list(hists)
    id_counter=0
       
    for p in range(0,len(dis_prim)):
        temp=sites.count(dis_prim[p])
        primarys.append([dis_prim[p],list(set(hists[id_counter:id_counter+temp]))])
        id_counter+=temp
    
    s=Clinical_sample.objects.all()
    for f in f_type:
        temp=list(set(s.values_list(f,flat=True)))
        primh_filter.append([f,temp])
    
    all_full_name=cell_d_name+d_name
    return  {
        'all_full_name':mark_safe(json.dumps(all_full_name)),  #full name of all datasets
        'same_name':mark_safe(json.dumps(same_name)),  #short name for all datasets
        'cell_d_name':mark_safe(json.dumps(cell_d_name)),  #cell line dataset name(full)
        'cell_datasets':cell_datasets, 
        'd_name': mark_safe(json.dumps(d_name)),  #clinical dataset name
        'datasets': datasets,
        'primarys': primarys,
        'primh_filter':primh_filter,
        
    }
    
def sample_microarray(request):
    
    d=Clinical_Dataset.objects.all()
    d_name=list(d.values_list('name',flat=True))
    datasets=[]  #[[dataset_name,[[primary_site,[primary_histology]]],[[filter_type,[filter_choice]]]]
    primarys=[]  #[[primary_site,[primary_hist]]]
    primh_filter=[] #[[filter_type,[filter_choice]]]
    f_type=['age','gender','ethnic','grade','stage','stageT','stageN','stageM','metastatic']
    for i in d_name:
        
        datasets.append([i,[],[]])
        sample=Clinical_sample.objects.filter(dataset_id__name=i).order_by('primary_site')
        sites=list(sample.values_list('primary_site',flat=True))
        hists=list(sample.values_list('primary_hist',flat=True))
        
        dis_prim=list(sample.values_list('primary_site',flat=True).distinct())
        hists=list(hists)
        id_counter=0
       
        for p in range(0,len(dis_prim)):
            temp=sites.count(dis_prim[p])
            datasets[-1][1].append([dis_prim[p],list(set(hists[id_counter:id_counter+temp]))])
            id_counter+=temp
        
        for f in f_type:
            temp=list(set(sample.values_list(f,flat=True)))
            datasets[-1][2].append([f,temp])
    
    sample=Clinical_sample.objects.all().order_by('primary_site')
    sites=list(sample.values_list('primary_site',flat=True))
    hists=list(sample.values_list('primary_hist',flat=True))
    dis_prim=list(sample.values_list('primary_site',flat=True).distinct())
    hists=list(hists)
    id_counter=0
       
    for p in range(0,len(dis_prim)):
        temp=sites.count(dis_prim[p])
        primarys.append([dis_prim[p],list(set(hists[id_counter:id_counter+temp]))])
        id_counter+=temp
    
    s=Clinical_sample.objects.all()
    for f in f_type:
        temp=list(set(s.values_list(f,flat=True)))
        primh_filter.append([f,temp])
    
    return render(request, 'sample_microarray.html', {
        'd_name': mark_safe(json.dumps(d_name)),
        'datasets': datasets,
        'primarys': primarys,
        'primh_filter':primh_filter,
        
    })
def user_pca(request):

    #load the ranking file and the table of probe first,open files
    pform=request.POST['data_platform']
    uni=[] #to store valid probe offset for getting the correct data
    uni_probe=[]
    gene_flag=0
    if(pform=="others"): #gene level
        gene_flag=1
        if(request.POST['ngs']=="ngs_u133a"):
            pform="U133A"
        else:
            pform="PLUS2"
        
    elif (pform=="U133A"):
        quantile=list(np.load('ranking_u133a.npy'))
        probe_path=Path('../').resolve().joinpath('src','Affy_U133A_probe_info.csv')
        probe_list = pd.read_csv(probe_path.as_posix())
        uni_probe=pd.unique(probe_list['PROBEID'])
        
    else:
        quantile=np.load('ranking_u133plus2.npy')
        probe_path=Path('../').resolve().joinpath('src','Affy_U133plus2_probe_info.csv')
        probe_list = pd.read_csv(probe_path.as_posix())
        uni_probe=pd.unique(probe_list['PROBEID'])
        
    
    
    
    
    propotion=0
    table_propotion=0
    show=request.POST['show_type']      #get the pca show type
    nci_size=Sample.objects.filter(dataset_id__name__in=["NCI60"]).count()
    gse_size=Sample.objects.filter(dataset_id__name__in=["GSE36133"]).count()
    group_counter=1
    user_out_group=[]
    s_group_dict={}  #store sample
    offset_group_dict={} #store offset
    cell_line_dict={}
    
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
    
    clinic=list(Clinical_Dataset.objects.all().values_list('name',flat=True))
    clline=list(Dataset.objects.all().values_list('name',flat=True))
    
    all_exist_dataset=[]
    for i in range(1,group_counter+1):
        dname='dataset_g'+str(i)
        all_exist_dataset=all_exist_dataset+request.POST.getlist(dname)
    all_exist_dataset=list(set(all_exist_dataset))
    
    all_base=[0]
    for i in range(0,len(all_exist_dataset)-1):
        if all_exist_dataset[i] in clline:
            all_base.append(all_base[i]+Sample.objects.filter(dataset_id__name__in=[all_exist_dataset[i]]).count())
        else:
            all_base.append(all_base[i]+Clinical_sample.objects.filter(dataset_id__name__in=[all_exist_dataset[i]]).count())
    
    
    all_c=[]
    
    for i in range(1,group_counter+1):
        s_group_dict['g'+str(i)]=[]
        offset_group_dict['g'+str(i)]=[]
        cell_line_dict['g'+str(i)]=[]
        dname='dataset_g'+str(i)
        datasets=request.POST.getlist(dname)
        group_name.append('g'+str(i))
        
        for dn in datasets:
            if dn=='Sanger Cell Line Project':
                c='select_sanger_g'+str(i)
            elif dn=='NCI60':
                c='select_nci_g'+str(i)
            elif dn=='GSE36133':
                c='select_ccle_g'+str(i)

            if dn in clline:
                temp=list(set(request.POST.getlist(c)))
                if 'd_sample' in show:
                    if all_c==[]:
                        all_c=all_c+temp
                        uni=temp
                    else:
                        uni=list(set(temp)-set(all_c))
                        all_c=all_c+uni
                else:
                    uni=list(temp)          #do not filter duplicate input only when select+centroid
                s=Sample.objects.filter(cell_line_id__name__in=uni,dataset_id__name__in=[dn]).order_by('dataset_id'
                ).select_related('cell_line_id__name','cell_line_id__primary_site','cell_line_id__primary_hist','dataset_id','dataset_id__name')
                
                cell_line_dict['g'+str(i)]=cell_line_dict['g'+str(i)]+list(s.values_list('cell_line_id__name',flat=True))
                s_group_dict['g'+str(i)]=s_group_dict['g'+str(i)]+list(s)
                offset_group_dict['g'+str(i)]=offset_group_dict['g'+str(i)]+list(np.add(list(s.values_list('offset',flat=True)),all_base[all_exist_dataset.index(dn)]))
                

            else: #dealing with clinical sample datasets
                com_hists=list(set(request.POST.getlist('primd_'+dn+'_g'+str(i))))    #can I get this by label to reduce number of queries?
                com_hists=[w1 for segments in com_hists for w1 in segments.split('/')]
                #print(com_hists)
                prims=com_hists[0::2]
                hists=com_hists[1::2]
                temp=request.POST.getlist('filter_'+dn+'_g'+str(i))
                age=[]
                gender=[]
                ethnic=[]
                grade=[]
                stage=[]
                T=[]
                N=[]
                M=[]
                metas=[]
                for t in temp:
                    if 'stage/' in t:
                        stage.append(t[6:])
                    elif 'gender/' in t:
                        gender.append(t[7:])
                    elif 'ethnic/' in t:
                        ethnic.append(t[7:])
                    elif 'grade/' in t:
                        grade.append(t[6:])
                    elif 'stageT/' in t:
                        T.append(t[7:])
                    elif 'stageN/' in t:
                        N.append(t[7:])
                    elif 'stageM/' in t:
                        M.append(t[7:])
                    elif 'metastatic/' in t:
                        if t[11:]=='False':
                            metas.append(0)
                        else:
                            metas.append(1)
                    else:  #"age/"
                        age.append(t[4:])
                #print(len(prims))
                #print(len(hists))
                
                for x in range(0,len(prims)):
                    s=Clinical_sample.objects.filter(dataset_id__name=dn,primary_site=prims[x],
                    primary_hist=hists[x],
                    age__in=age,
                    gender__in=gender,
                    ethnic__in=ethnic,
                    stage__in=stage,
                    grade__in=grade,
                    stageT__in=T,
                    stageN__in=N,
                    stageM__in=M,
                    metastatic__in=metas,
                    ).select_related('dataset_id').order_by('id')
                    s_group_dict['g'+str(i)]=s_group_dict['g'+str(i)]+list(s)
                    cell_line_dict['g'+str(i)]=cell_line_dict['g'+str(i)]+list(s.values_list('name',flat=True))
                    offset_group_dict['g'+str(i)]=offset_group_dict['g'+str(i)]+list(np.add(list(s.values_list('offset',flat=True)),all_base[all_exist_dataset.index(dn)]))
                
                #return render_to_response('welcome.html',locals())
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
        if str(type(i))=="<class 'probes.models.Sample'>":
            ##print("i am sample!!")
            cell_object.append(i.cell_line_id)
        else:
            ##print("i am clinical!!")
            cell_object.append(i)
            
    #read the user file   
    text=request.FILES.getlist('user_file')
    user_counter=len(text)
    if(gene_flag==1):
        user_counter=1
    user_dict={} #{user group number:user 2d array}
    samples=0
    nans=[] #to store the probe name that has nan 
    for x in range(1,user_counter+1):
        #check the file format and content here first        
        filetype=str(text[x-1]).split('.')
        if(filetype[-1]!="csv"):
            error_reason='You have the wrong file type. Please upload .csv files'
            return render_to_response('pca_error.html',RequestContext(request,
            {
            'error_reason':mark_safe(json.dumps(error_reason)),
            }))
        if(text[x-1].size>=80000000): #bytes
            error_reason='The file size is too big. Please upload .csv file with size less than 80MB.'
            return render_to_response('pca_error.html',RequestContext(request,
            {
            'error_reason':mark_safe(json.dumps(error_reason)),
            }))
        temp_data = pd.read_csv(text[x-1])
        col=list(temp_data.columns.values)
        samples=samples+len(col)-1
        if(samples==0):
           error_reason='The file does not have any samples.'
           return render_to_response('pca_error.html',RequestContext(request,
           {
            'error_reason':mark_safe(json.dumps(error_reason)),
           })) 
        if(gene_flag==0): #probe level check
            check_probe=[str(x) for x in list(temp_data.iloc[:,0]) if not str(x).lower().startswith('affx')]  
            #print(len(check_probe))
            if(len(check_probe)!=len(uni_probe)):
                error_reason='The probe number does not match with the platform you selected.'
                return render_to_response('pca_error.html',RequestContext(request,
                {
                'error_reason':mark_safe(json.dumps(error_reason)),
                }))
            if(set(check_probe)!=set(uni_probe)):
                error_reason='The probe number or probe name in your file does not match the platform you selected.'
                #error_reason+='</br>The probes that are not in the platform: '+str(set(check_probe)-set(uni_probe))[1:-1]
                #error_reason+='</br>The probes that are lacking: '+str(set(uni_probe)-set(check_probe))[1:-1]
                return render_to_response('pca_error.html',RequestContext(request,
                {
                'error_reason':mark_safe(json.dumps(error_reason)),
                }))
        col=list(temp_data.columns.values)
        n=pd.isnull(temp_data).any(1).nonzero()[0]
        nans += list(temp_data[col[0]][n])
        
        user_dict[x]=temp_data
    
    if 'd_sample' in show:
        if((len(all_sample)+samples)<4):
            error_reason='You should have at least 4 samples for PCA. The samples are not enough.<br />'\
            'The total number of samples in your uploaded file is '+str(samples)+'.<br />'\
            'The number of samples you selected is '+str(len(all_sample))+'.<br />'\
            'Total is '+str(len(all_sample)+samples)+'.'
            return render_to_response('pca_error.html',RequestContext(request,
            {
            'error_reason':mark_safe(json.dumps(error_reason)),
            }))
    else:
        s_count=0
        sample_list=[]
        a_sample=np.array(all_sample)
        for i in range(1,group_counter+1):
            dis_cellline=list(set(cell_object[g_s_counter[i-1]:g_s_counter[i]]))
            a_cell_object=np.array(cell_object)
            for c in dis_cellline:    
                
                temp1=np.where((a_cell_object==c))[0]
                
                temp2=np.where((temp1>=g_s_counter[i-1])&(temp1<g_s_counter[i]))
                total_offset=temp1[temp2]

                selected_sample=a_sample[total_offset]
                
                if list(selected_sample) in sample_list:   #to prevent two different colors in different group
                    continue
                else:
                    sample_list.append(list(selected_sample))
                    s_count=s_count+1
                    if(s_count>=4): #check this part
                        break
            if(s_count>=4): #check this part
                break
        if((s_count+samples)<4):
            error_reason='Since the display method is [centroid], you should have at least 4 dots for PCA. The total number is not enough.<br />'\
            'The total number of dots in you uploaded file is '+str(samples)+'.<br />'\
            'The number of centroid dots you selected is '+str(s_count)+'.<br />'\
            'Total is '+str(s_count+samples)+'.'
            return render_to_response('pca_error.html',RequestContext(request,
            {
            'error_reason':mark_safe(json.dumps(error_reason)),
            }))
    new_name=[]
    origin_name=[]
    com_gene=[] #for ngs select same gene
    nans=list(set(nans))
    for x in range(1,user_counter+1):
        
        #temp_data = pd.read_csv(text[x-1])
        temp_data=user_dict[x]
        col=list(temp_data.columns.values)
        col[0]='probe'
        temp_data.columns=col
        
        temp_data.index = temp_data['probe']
        temp_data.index.name = None
        temp_data=temp_data.iloc[:, 1:]
        
        #add "use_" to user's sample names
        col_name=list(temp_data.columns.values)   #have user's sample name list here
        origin_name=origin_name+list(temp_data.columns.values)
        col_name=[ "user_"+str(index)+"_"+s for index,s in enumerate(col_name)]
        temp_data.columns=col_name
        new_name=new_name+col_name
        
        if(gene_flag==0):
            try:
                temp_data=temp_data.reindex(uni_probe)
            except ValueError:
                return HttpResponse('The file has probes with the same names, please let them be unique.')
            #remove probe that has nan
            temp_data=temp_data.drop(nans)
            temp_data=temp_data.rank(method='dense')
    
            #this is for quantile
            for i in col_name:
                for j in range(0,len(temp_data[i])):
                    #if(not(np.isnan(temp_data[i][j]))):
                    temp_data[i][j]=quantile[int(temp_data[i][j]-1)]
            
            if x==1:
                data=temp_data
            else:
                data=np.concatenate((data,temp_data), axis=1)
            user_dict[x]=np.array(temp_data)
        else:
            temp_data=temp_data.drop(nans)  #drop nan
            temp_data=temp_data.loc[~(temp_data==0).all(axis=1)] #drop all rows with 0 here
            temp_data=temp_data.groupby(temp_data.index).first() #drop the duplicate gene row
            user_dict[x]=temp_data
        
        #print(temp_data)
        
        
    
    #delete nan, combine user data to the datasets,transpose matrix
    for x in range(0,len(all_exist_dataset)):
        if(gene_flag==0):
            if all_exist_dataset[x] in clline:
                pth=Path('../').resolve().joinpath('src',Dataset.objects.get(name=all_exist_dataset[x]).data_path)
            else:
                pth=Path('../').resolve().joinpath('src',Clinical_Dataset.objects.get(name=all_exist_dataset[x]).data_path)
        else:
            if all_exist_dataset[x] in clline:
                pth=Path('../').resolve().joinpath('src','gene_'+Dataset.objects.get(name=all_exist_dataset[x]).data_path)
            else:
                pth=Path('../').resolve().joinpath('src','gene_'+Clinical_Dataset.objects.get(name=all_exist_dataset[x]).data_path)
        if x==0:
            val=np.load(pth.as_posix())
        else:
            val=np.hstack((val, np.load(pth.as_posix())))#combine together
    
    #database dataset remove nan probes
    if(gene_flag==0):
        uni=[]
        p_offset=list(ProbeID.objects.filter(platform__name__in=[pform],Probe_id__in=nans).values_list('offset',flat=True))
        for n in range(0,len(uni_probe)):
            if(n not in p_offset):
                uni.append(n)
    else:
        #deal with ngs uploaded data here
        probe_path=Path('../').resolve().joinpath('src','new_human_gene_info.txt')
        #probe_list = pd.read_csv(probe_path.as_posix())
        #notice duplicate
        #get the match gene first, notice the size issue
        info=pd.read_csv(probe_path.as_posix(),sep='\t')
        col=list(info.columns.values)
        col[0]='symbol'
        info.columns=col
                
        info.index = info['symbol']
        info.index.name = None
        info=info.iloc[:, 1:]
        
        data=user_dict[1]
        #data=data.groupby(data.index).first() #drop the duplicate gene row
        com_gene=list(data.index)
        temp_data=data
        rloop=divmod(len(com_gene),990)
        if(rloop[1]==0):
            rloop=(rloop[0]-1,0)
        gg=[]
        for x in range(0,rloop[0]+1):
            gg+=list(Gene.objects.filter(platform__name__in=[pform],symbol__in=com_gene[x*990:(x+1)*990]).order_by('offset'))
        exist_gene=[]
        uni=[]
        for i in gg:
            exist_gene.append(i.symbol)
            uni.append(i.offset)
        info=info.drop(exist_gene,errors='ignore')
        new_data=temp_data.loc[data.index.isin(exist_gene)].reindex(exist_gene)
        ##print(exist_gene)
        ##print(new_data.index)
        
        #search remain symbol's alias and symbol
        search_alias=list(set(com_gene)-set(exist_gene))
        
        for i in search_alias:
            re_symbol=list(set(info.loc[info['alias'].isin([i])].index)) #find whether has alias first
            if(len(re_symbol)!=0):
                re_match=Gene.objects.filter(platform__name__in=[pform],symbol__in=re_symbol).order_by('offset') #check the symbol in database or not
                repeat=len(re_match)
                if(repeat!=0): #match gene symbol in database      
                    ##print(re_match)
                    for x in re_match:
                        to_copy=data.loc[i]
                        to_copy.name=x.symbol
                        new_data=new_data.append(to_copy)
                        uni.append(x.offset)
                        info=info.drop(x.symbol,errors='ignore')
        
        user_dict[1]=np.array(new_data)
        ##print("length of new data:"+str(len(new_data)))
        ##print("data:")
        ##print(data)
        ##print("new_data:")
        ##print(new_data)
        data=new_data
        
    if 'd_sample' in show:
        val=val[np.ix_(uni,all_offset)]
        #print(len(val))
        user_offset=len(val[0])
        if(gene_flag==1):
            #do the rank invariant here
            #print("sample with ngs data do rank invariant here")
            
            ref_path=Path('../').resolve().joinpath('src','cv_result.txt')
            ref=pd.read_csv(ref_path.as_posix())
            col=list(ref.columns.values)
            col[0]='symbol'
            ref.columns=col
                    
            ref.index = ref['symbol']
            ref.index.name = None
            ref=ref.iloc[:, 1:]
            ref=ref.iloc[:5000,:]  #rank invariant need 5000 genes
            same_gene=list(ref.index.intersection(data.index))
            
            #to lowess
            rref=pandas2ri.py2ri(ref.loc[same_gene])
            rngs=pandas2ri.py2ri(data.loc[same_gene].mean(axis=1))
            rall=pandas2ri.py2ri(data)
            ro.globalenv['x'] = rngs
            ro.globalenv['y'] = rref
            ro.globalenv['newx'] = rall

            r('x<-as.vector(as.matrix(x))')
            r('y<-as.vector(as.matrix(y))')
            r('newx<-as.matrix(newx)')
            try:
                if(request.POST['data_type']=='raw'):
                    r('y.loess<-loess(2**y~x,span=0.3)')
                    r('for(z in c(1:ncol(newx))) newx[,z]=log2(as.matrix(predict(y.loess,newx[,z])))')
                elif(request.POST['data_type']=='log2'):
                    r('y.loess<-loess(2**y~2**x,span=0.3)')
                    r('for(z in c(1:ncol(newx))) newx[,z]=log2(as.matrix(predict(y.loess,2**newx[,z])))')
                else:
                    r('y.loess<-loess(2**y~10**x,span=0.3)')
                    r('for(z in c(1:ncol(newx))) newx[,z]=log2(as.matrix(predict(y.loess,10**newx[,z])))')
            except RRuntimeError:
                error_reason='Match too less genes. Check your gene symbols again. We use NCBI standard gene symbol.'
                return render_to_response('pca_error.html',RequestContext(request,
                {
                'error_reason':mark_safe(json.dumps(error_reason)),
                }))
            
            #r('for(z in c(1:ncol(newx))) newx[,z]=log2(as.matrix(predict(y.loess,newx[,z])))')
                
            data=r('newx')
            #print(data[:10])
            #print(type(data))
            
        val=np.hstack((np.array(val), np.array(data)))
        val=val[~np.isnan(val).any(axis=1)]
        val=np.transpose(val)
    else:

        #val=np.array(val)
        val=val[np.ix_(uni)]
        user_offset=len(val[0])
        
        if(gene_flag==1):
            #print("sample with ngs data do rank invariant here")
            ref_path=Path('../').resolve().joinpath('src','cv_result.txt')
            ref=pd.read_csv(ref_path.as_posix())
            col=list(ref.columns.values)
            col[0]='symbol'
            ref.columns=col
                    
            ref.index = ref['symbol']
            ref.index.name = None
            ref=ref.iloc[:, 1:]
            ref=ref.iloc[:5000,:]  #rank invariant need 5000 genes
            same_gene=list(ref.index.intersection(data.index))
            
            #to lowess
            rref=pandas2ri.py2ri(ref.loc[same_gene])
            rngs=pandas2ri.py2ri(data.loc[same_gene].mean(axis=1))
            rall=pandas2ri.py2ri(data)
            ro.globalenv['x'] = rngs
            ro.globalenv['y'] = rref
            ro.globalenv['newx'] = rall

            r('x<-as.vector(as.matrix(x))')
            r('y<-as.vector(as.matrix(y))')
            r('newx<-as.matrix(newx)')
            try:
                if(request.POST['data_type']=='raw'):
                    r('y.loess<-loess(2**y~x,span=0.3)')
                    r('for(z in c(1:ncol(newx))) newx[,z]=log2(as.matrix(predict(y.loess,newx[,z])))')
                elif(request.POST['data_type']=='log2'):
                    r('y.loess<-loess(2**y~2**x,span=0.3)')
                    r('for(z in c(1:ncol(newx))) newx[,z]=log2(as.matrix(predict(y.loess,2**newx[,z])))')
                else:
                    r('y.loess<-loess(2**y~10**x,span=0.3)')
                    r('for(z in c(1:ncol(newx))) newx[,z]=log2(as.matrix(predict(y.loess,10**newx[,z])))')
            except RRuntimeError:
                error_reason='Match too less genes. Check your gene symbols again. We use NCBI standard gene symbol.'
                return render_to_response('pca_error.html',RequestContext(request,
                {
                'error_reason':mark_safe(json.dumps(error_reason)),
                }))

            #r('for(z in c(1:ncol(newx))) newx[,z]=log2(as.matrix(predict(y.loess,newx[,z])))')
                
            data=r('newx')
            #print(data[:10])
            #print(type(data))
            
        val=np.hstack((np.array(val), np.array(data)))
        val=val[~np.isnan(val).any(axis=1)]  
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
    n=4  #need to fix to the best one 
    
    if 'd_sample' in show:
        #count the pca first
        pca= PCA(n_components=n)
        
        #combine user sample's offset to all_offset in another variable
        
        Xval = pca.fit_transform(val[:,:])  #cannot get Xval with original offset any more
        ratio_temp=pca.explained_variance_ratio_
        propotion=sum(ratio_temp[1:n])
        table_propotion=sum(ratio_temp[0:n])
        user_new_offset=len(all_offset)
        ##print(Xval)
        
        max=0
        min=10000000000
        out_group=[]
        exist_cell={}#cell line object:counter
        for g in range(1,group_counter+1):
              
            output_cell={}
            check={}
            for s in range(g_s_counter[g-1],g_s_counter[g]):
                
                if str(type(all_sample[s]))=="<class 'probes.models.Sample'>":
                    cell=all_sample[s].cell_line_id
                else:
                    cell=all_sample[s]
               
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
                                distance=np.linalg.norm(Xval[i][n-3:n]-Xval[s][n-3:n])
                                if distance<min:
                                    min=distance
                                if distance>max:
                                    max=distance
                                output_cell[cell][1].append([all_cellline[s]+'('+str(exist_cell[cell])+')'
                                ,all_sample[s].name,cell.primary_site,cell.primary_hist,
                                all_sample[s].dataset_id.name,all_cellline[i],all_sample[i].name,cell_object[i].primary_site
                                ,cell_object[i].primary_hist,all_sample[i].dataset_id.name,distance])
                                
                                check[all_sample[s].name].append(all_sample[i].name)
                                
                        except KeyError:
                            distance=np.linalg.norm(Xval[i][n-3:n]-Xval[s][n-3:n])
                            if distance<min:
                                min=distance
                            if distance>max:
                                max=distance
                            output_cell[cell][1].append([all_cellline[s]+'('+str(exist_cell[cell])+')'
                                ,all_sample[s].name,cell.primary_site,cell.primary_hist,
                                all_sample[s].dataset_id.name,all_cellline[i],all_sample[i].name,cell_object[i].primary_site
                                ,cell_object[i].primary_hist,all_sample[i].dataset_id.name,distance])
                            check[all_sample[s].name].append(all_sample[i].name)
                g_count=1   
                u_count=len(user_dict[g_count][0])  #sample number in first user file
                for i in range(user_new_offset,user_new_offset+len(origin_name)):  #remember to prevent empty file uploaded
                    distance=np.linalg.norm(Xval[i][n-3:n]-Xval[s][n-3:n])
                    if distance<min:
                        min=distance
                    if distance>max:
                        max=distance
                    output_cell[cell][1].append([all_cellline[s]+'('+str(exist_cell[cell])+')'
                                ,all_sample[s].name,cell.primary_site,cell.primary_hist,
                                all_sample[s].dataset_id.name," ",origin_name[i-user_new_offset]," "," ","User Group"+str(g_count),distance])
                    if ((i-user_new_offset+1)==u_count):
                        g_count+=1
                        try:
                            u_count+=len(user_dict[g_count][0])
                        except KeyError:
                            u_count+=0
                    
                if(g==1):
                    name3.append(all_cellline[s]+'('+str(exist_cell[cell])+')'+'<br>'+all_sample[s].name)
                    X3.append(round(Xval[s][n-3],5))
                    Y3.append(round(Xval[s][n-2],5))
                    Z3.append(round(Xval[s][n-1],5))
                elif(g==2):
                    name4.append(all_cellline[s]+'('+str(exist_cell[cell])+')'+'<br>'+all_sample[s].name)
                    X4.append(round(Xval[s][n-3],5))
                    Y4.append(round(Xval[s][n-2],5))
                    Z4.append(round(Xval[s][n-1],5))
                elif(g==3):
                    name5.append(all_cellline[s]+'('+str(exist_cell[cell])+')'+'<br>'+all_sample[s].name)
                    X5.append(round(Xval[s][n-3],5))
                    Y5.append(round(Xval[s][n-2],5))
                    Z5.append(round(Xval[s][n-1],5))
                
            dictlist=[]
            for key, value in output_cell.items():
                temp = [value]
                dictlist+=temp
            output_cell=list(dictlist)
            out_group.append(["Dataset Group"+str(g),output_cell])
            
            
            if g==group_counter:
                output_cell={}
                g_count=1  
                output_cell[g_count]=[" ",[]]                
                u_count=len(user_dict[g_count][0])
                temp_count=u_count
                temp_g=1
                before=0
                for i in range(user_new_offset,user_new_offset+len(origin_name)):
                    for x in range(0,len(all_sample)):
                        distance=np.linalg.norm(Xval[x][n-3:n]-Xval[i][n-3:n])
                        if distance<min:
                            min=distance
                        if distance>max:
                            max=distance
                        output_cell[g_count][1].append([origin_name[i-user_new_offset],"User Group"+str(g_count),all_cellline[x]
                                    ,all_sample[x].name,cell_object[x].primary_site,cell_object[x].primary_hist,
                                    all_sample[x].dataset_id.name,distance])
                    
                    temp_g=1
                    temp_count=len(user_dict[temp_g][0])
                    for j in range(user_new_offset,user_new_offset+before):
                        if ((j-user_new_offset)==temp_count):
                            temp_g+=1
                            try:
                                temp_count+=len(user_dict[temp_g][0])
                            except KeyError:
                                temp_count+=0
                        distance=np.linalg.norm(Xval[j][n-3:n]-Xval[i][n-3:n])
                        if distance<min:
                            min=distance
                        if distance>max:
                            max=distance
                        output_cell[g_count][1].append([origin_name[i-user_new_offset],"User Group"+str(g_count)
                        ," ",origin_name[j-user_new_offset]," "," ","User Group"+str(temp_g),distance])
                        
                    temp_g=g_count  
                    temp_count=len(user_dict[g_count][0])
                    for j in range(i+1,user_new_offset+len(origin_name)):
                        if ((j-user_new_offset)==temp_count):
                            temp_g+=1
                            try:
                                temp_count+=len(user_dict[temp_g][0])
                            except KeyError:
                                temp_count+=0
                        distance=np.linalg.norm(Xval[j][n-3:n]-Xval[i][n-3:n])
                        if distance<min:
                            min=distance
                        if distance>max:
                            max=distance
                        output_cell[g_count][1].append([origin_name[i-user_new_offset],"User Group"+str(g_count)
                        ," ",origin_name[j-user_new_offset]," "," ","User Group"+str(temp_g),distance])                    

                    if g_count==1:
                        name1.append(origin_name[i-user_new_offset])
                        X1.append(round(Xval[i][n-3],5))
                        Y1.append(round(Xval[i][n-2],5))
                        Z1.append(round(Xval[i][n-1],5))
                    else:
                        name2.append(origin_name[i-user_new_offset])
                        X2.append(round(Xval[i][n-3],5))
                        Y2.append(round(Xval[i][n-2],5))
                        Z2.append(round(Xval[i][n-1],5))
                    if ((i-user_new_offset+1)==u_count):
                        dictlist=[]
                        for key, value in output_cell.items():
                            temp = [value]
                            dictlist+=temp
                        output_cell=list(dictlist)
                        user_out_group.append(["User Group"+str(g_count),output_cell])
                        
                        g_count+=1
                        before=u_count+before
                        #print("I am here!!")
                        try:
                            u_count+=len(user_dict[g_count][0])
                            output_cell={}
                            output_cell[g_count]=[" ",[]]
                        except KeyError:
                            u_count+=0
                    
        #[g,[group_cell_1 object,[[outputs paired1,......,],[paired2],[paired3]]],[group_cell_2 object,[[pair1],[pair2]]]]
        #for xx in origin_name:
            #sample_counter[xx]=1
        ##print(out_group)
        element_counter=0
        for i in out_group:
            for temp_list in i[1]:
                element_counter=element_counter+len(temp_list[1])
                for temp in temp_list[1]:
                    if(temp[5]!=" "):
                        temp[5]=temp[5]+'('+str(sample_counter[temp[6]])+')'
        for i in user_out_group:
            for temp_list in i[1]:
                for temp in temp_list[1]:
                    if(temp[2]!=" "):
                        temp[2]=temp[2]+'('+str(sample_counter[temp[3]])+')'

        return_html='user_pca.html'
    else:
        #This part is for centroid display
        return_html='user_pca_center.html'
        
        #This part is for select cell line base on dataset,count centroid base on the dataset
        #group中的cell line為單位來算重心

        location_dict={} #{group number:[[cell object,dataset,new location]]}
        combined=[]
        sample_list=[]
        pca_index=np.array(pca_index)
        X_val=[]
        val_a=np.array(val)
        a_all_offset=np.array(all_offset)
        a_sample=np.array(all_sample)
        for i in range(1,group_counter+1):
            dis_cellline=list(set(cell_object[g_s_counter[i-1]:g_s_counter[i]]))  #cell object may have duplicate cell line since:NCI A + CCLE A===>[A,A]
            location_dict['g'+str(i)]=[]
            dataset_dict={}
            a_cell_object=np.array(cell_object)
            
            for c in dis_cellline:    #dis_cellline may not have the same order as cell_object
                
                temp1=np.where((a_cell_object==c))[0]
                
                temp2=np.where((temp1>=g_s_counter[i-1])&(temp1<g_s_counter[i]))
                total_offset=temp1[temp2]
                selected_val=val_a[:,a_all_offset[total_offset]]
                selected_val=np.transpose(selected_val)
                new_loca=(np.mean(selected_val,axis=0,dtype=np.float64,keepdims=True)).tolist()[0]

                selected_sample=a_sample[total_offset]
                
                if list(selected_sample) in sample_list:   #to prevent two different colors in different group
                    continue
                else:
                    sample_list.append(list(selected_sample))
                
                
                d_temp=[]
                for s in selected_sample:
                    d_temp.append(s.dataset_id.name)
                dataset_dict[c]="/".join(list(set(d_temp)))  
                X_val.append(new_loca)
                location_dict['g'+str(i)].append([c,dataset_dict[c],len(X_val)-1])  #the last part is the index to get pca result from new_val
                combined.append([c,dataset_dict[c],len(X_val)-1])  #all cell line, do not matter order
        
        #run the pca
        user_new_offset=len(X_val)
        temp_val=np.transpose(val[:,user_offset:])
        for x in range(0,len(temp_val)):
            X_val.append(list(temp_val[x]))
        if(len(X_val)<4):
            error_reason='Since the display method is [centroid], you should have at least 4 dots for PCA. The total number is not enough.<br />'\
            'The total number of dots in you uploaded file is '+str(len(temp_val))+'.<br />'\
            'The number of centroid dots you selected is '+str(len(X_val)-len(temp_val))+'.<br />'\
            'Total is '+str(len(X_val))+'.'
            return render_to_response('pca_error.html',RequestContext(request,
            {
            'error_reason':mark_safe(json.dumps(error_reason)),
            }))
        X_val=np.matrix(X_val)
        pca= PCA(n_components=n)
        new_val = pca.fit_transform(X_val[:,:])  #cannot get Xval with original offset any more
        ratio_temp=pca.explained_variance_ratio_
        propotion=sum(ratio_temp[1:n])
        table_propotion=sum(ratio_temp[0:n])
        #print(new_val)

        out_group=[]
        min=10000000000
        max=0
        element_counter=0
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
                            distance=np.linalg.norm(np.array(new_val[group_c[2]][n-3:n])-np.array(new_val[temp_list[2]][n-3:n]))
                            if distance==0:
                                continue
                            if distance<min:
                                min=distance
                            if distance>max:
                                max=distance
                            output_cell[len(output_cell)-1][1].append([cell.name,cell.primary_site,cell.primary_hist
                        ,group_c[1],temp_list[0].name,temp_list[0].primary_site,temp_list[0].primary_hist,temp_list[1],distance])
                            element_counter=element_counter+1
                    except KeyError:
                        distance=np.linalg.norm(np.array(new_val[group_c[2]][n-3:n])-np.array(new_val[temp_list[2]][n-3:n]))
                        if distance==0:
                            continue
                        if distance<min:
                            min=distance
                        if distance>max:
                            max=distance
                        output_cell[len(output_cell)-1][1].append([cell.name,cell.primary_site,cell.primary_hist
                        ,group_c[1],temp_list[0].name,temp_list[0].primary_site,temp_list[0].primary_hist,temp_list[1],distance])
                        element_counter=element_counter+1
                        exist_cell[key_string].append(temp_string)
                g_count=1   
                u_count=len(user_dict[g_count][0])  #sample number in first user file
                for i in range(user_new_offset,user_new_offset+len(origin_name)):
                        distance=np.linalg.norm(np.array(new_val[group_c[2]][n-3:n])-np.array(new_val[i][n-3:n]))
                        if distance<min:
                            min=distance
                        if distance>max:
                            max=distance
                        output_cell[len(output_cell)-1][1].append([cell.name,cell.primary_site,cell.primary_hist,group_c[1]
                                    ,origin_name[i-user_new_offset]," "," ","User Group"+str(g_count),distance])
                        element_counter=element_counter+1
                        if ((i-user_new_offset+1)==u_count):
                            g_count+=1
                            try:
                                u_count+=len(user_dict[g_count][0])
                            except KeyError:
                                u_count+=0
                if(g==1):
                    name3.append(cell.name+'<br>'+group_c[1])
                    X3.append(round(new_val[group_c[2]][n-3],5))
                    Y3.append(round(new_val[group_c[2]][n-2],5))
                    Z3.append(round(new_val[group_c[2]][n-1],5))
                elif(g==2):
                    name4.append(cell.name+'<br>'+group_c[1])
                    X4.append(round(new_val[group_c[2]][n-3],5))
                    Y4.append(round(new_val[group_c[2]][n-2],5))
                    Z4.append(round(new_val[group_c[2]][n-1],5))
                elif(g==3):
                    name5.append(cell.name+'<br>'+group_c[1])
                    X5.append(round(new_val[group_c[2]][n-3],5))
                    Y5.append(round(new_val[group_c[2]][n-2],5))
                    Z5.append(round(new_val[group_c[2]][n-1],5))
                
            out_group.append(["Dataset Group"+str(g),output_cell])
            

            if g==group_counter:
                output_cell=[]
                g_count=1 
                output_cell.append([" ",[]])
                u_count=len(user_dict[g_count][0])
                temp_count=u_count
                temp_g=1
                before=0
                for i in range(user_new_offset,user_new_offset+len(origin_name)):
                    for temp_list in combined:
                        c=temp_list[0]
                        distance=np.linalg.norm(np.array(new_val[i][n-3:n])-np.array(new_val[temp_list[2]][n-3:n]))
                        if distance<min:
                            min=distance
                        if distance>max:
                            max=distance
                        output_cell[len(output_cell)-1][1].append([origin_name[i-user_new_offset],"User Group"+str(g_count)
                        ,c.name,c.primary_site,c.primary_hist,temp_list[1],distance])
                    
                    temp_g=1
                    temp_count=len(user_dict[temp_g][0])
                    for j in range(user_new_offset,user_new_offset+before):
                        if ((j-user_new_offset)==temp_count):
                            temp_g+=1
                            try:
                                temp_count+=len(user_dict[temp_g][0])
                            except KeyError:
                                temp_count+=0
                        distance=np.linalg.norm(np.array(new_val[i][n-3:n])-np.array(new_val[j][n-3:n]))
                        if distance<min:
                            min=distance
                        if distance>max:
                            max=distance
                        output_cell[len(output_cell)-1][1].append([origin_name[i-user_new_offset],"User Group"+str(g_count)
                        ,origin_name[j-user_new_offset]," "," ","User Group"+str(temp_g),distance])
                    
                    temp_g=g_count  
                    temp_count=len(user_dict[g_count][0])
                    for x in range(i+1,user_new_offset+len(origin_name)):
                        if ((x-user_new_offset)==temp_count):
                            temp_g+=1
                            try:
                                temp_count+=len(user_dict[temp_g][0])
                            except KeyError:
                                temp_count+=0
                        distance=np.linalg.norm(np.array(new_val[i][n-3:n])-np.array(new_val[x][n-3:n]))
                        if distance<min:
                            min=distance
                        if distance>max:
                            max=distance
                        output_cell[len(output_cell)-1][1].append([origin_name[i-user_new_offset],"User Group"+str(g_count)
                        ,origin_name[x-user_new_offset]," "," ","User Group"+str(temp_g),distance])
                    if g_count==1:        
                        name1.append(origin_name[i-user_new_offset])
                        X1.append(round(new_val[i][n-3],5))
                        Y1.append(round(new_val[i][n-2],5))
                        Z1.append(round(new_val[i][n-1],5))  
                    else:
                        name2.append(origin_name[i-user_new_offset])
                        X2.append(round(new_val[i][n-3],5))
                        Y2.append(round(new_val[i][n-2],5))
                        Z2.append(round(new_val[i][n-1],5)) 
                    if ((i-user_new_offset+1)==u_count):
                        user_out_group.append(["User Group"+str(g_count),output_cell])
                        g_count+=1
                        before=u_count+before
                        #print("I am here!!")
                        try:
                            u_count+=len(user_dict[g_count][0])
                            output_cell=[]
                            output_cell.append([" ",[]])
                        except KeyError:
                            u_count+=0
    
    
    
    #print(element_counter)
    #print(show_row)
    if(element_counter>show_row):
        big_flag=1    
        sid=str(uuid.uuid1())+".csv"
        if(return_html=='user_pca.html'):
            dataset_header=['Group Cell Line/Clinical Sample','Sample Name','Primary Site','Primary Histology'
                ,'Dataset','Paired Cell Line name/Clinical Sample','Sample Name','Primary Site','Primary Histology','Dataset','Distance']
            user_header=['User Sample Name','Dataset','Paired Cell Line name/Clinical Sample','Sample Name','Primary Site','Primary Histology','Dataset','Distance']
        else:
            dataset_header=['Group Cell Line/Clinical Sample','Primary Site','Primary Histology'
                ,'Dataset','Paired Cell Line name/Clinical Sample','Primary Site','Primary Histology','Dataset','Distance']
            user_header=['User Sample Name','Dataset','Paired Cell Line name/Clinical Sample','Primary Site','Primary Histology','Dataset','Distance']    
        P=Path('../').resolve().joinpath('src','static','csv',"dataset_"+sid)
        userP=Path('../').resolve().joinpath('src','static','csv',"user_"+sid)
        assP=Path('../').resolve().joinpath('src','assets','csv',"dataset_"+sid)
        assuserP=Path('../').resolve().joinpath('src','assets','csv',"user_"+sid)
        #print("start writing files")
        with open(str(P), "w", newline='') as f:
            writer = csv.writer(f)
            for index,output_cell in out_group:
                writer.writerows([[index]])
                writer.writerows([dataset_header])
                for cell_line,b in output_cell:
                    writer.writerows(b)
        #print("end writing first file")
        with open(str(assP), "w", newline='') as ff:
            writer = csv.writer(ff)
            for index,output_cell in out_group:
                writer.writerows([[index]])
                writer.writerows([dataset_header])
                for cell_line,b in output_cell:
                    writer.writerows(b)
        #print("end writing 2 file")
        with open(str(assuserP), "w", newline='') as ff:
            writer = csv.writer(ff)
            for index,output_cell in user_out_group:
                writer.writerows([[index]])
                writer.writerows([user_header])
                for cell_line,b in output_cell:
                    writer.writerows(b)
        #print("end writing 3 file")
        with open(str(userP), "w", newline='') as f:
            writer = csv.writer(f)
            for index,output_cell in user_out_group:
                writer.writerows([[index]])
                writer.writerows([user_header])
                for cell_line,b in output_cell:
                    writer.writerows(b)
        #print("end writing 4 file")
        data_file_name="dataset_"+sid
        user_file_name="user_"+sid
    else:
        big_flag=0
        data_file_name=0
        user_file_name=0
    return render_to_response(return_html,RequestContext(request,
    {
    'min':min,'max':max,
    'big_flag':big_flag,
    'out_group':out_group,'user_out_group':user_out_group,
    'propotion':propotion,
    'table_propotion':table_propotion,
    'data_file_name':data_file_name,
    'user_file_name':user_file_name,
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
    
    #notice that we need to return a user_pca_center.html, too!!
    #return render_to_response('welcome.html',locals())
    
def express_profiling(request):
    return render(request, 'express_profiling.html', generate_samples())
   

def welcome(request):
    return render_to_response('welcome.html',locals())

def help(request):
    example_name="CellExpress_Examples.pptx"
    tutorial_name="CellExpress_Tutorial.pptx"
    return render_to_response('help.html',locals())

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
    not_found=[]
    #get probe from different platform
    pform=request.POST['data_platform']
    stop_end=101  
    return_page_flag=0
    user_probe_flag=0
    if(request.POST['user_type']=="all"):
        all_probe=ProbeID.objects.filter(platform__name=pform).order_by('offset')
        probe_offset=list(all_probe.values_list('offset',flat=True))
        pro_number=float(request.POST['probe_number'])  #significant 0.05 or 0.01
        all_probe=list(all_probe)
    elif(request.POST['user_type']=="genes"): #for all genes
        return_page_flag=1
        if(pform=="U133A"):
            probe_path=Path('../').resolve().joinpath('src','uni_u133a.txt')
            gene_list = pd.read_csv(probe_path.as_posix())
            all_probe=list(gene_list['SYMBOL'])
        else:
            probe_path=Path('../').resolve().joinpath('src','uni_plus2.txt')
            gene_list = pd.read_csv(probe_path.as_posix())
            all_probe=list(gene_list['SYMBOL'])
        pro_number=float(request.POST['probe_number_gene'])
        probe_offset=[]
        for i in range(0,len(all_probe)):
            probe_offset.append(i)
    else:
        indata=request.POST['keyword']
        indata = list(set(indata.split()))
        
        if(request.POST['gtype']=="probeid"):
            user_probe_flag=1
            all_probe=ProbeID.objects.filter(platform__name=pform,Probe_id__in=indata).order_by('offset')
            probe_offset=list(all_probe.values_list('offset',flat=True))
            pro_number=10000
            not_found=list(set(set(indata) - set(all_probe.values_list('Probe_id',flat=True))))
            all_probe=list(all_probe)
            
        else:
            probe_offset=[]
            return_page_flag=1
            pro_number=10000
            if(pform=="U133A"):
                probe_path=Path('../').resolve().joinpath('src','uni_u133a.txt')
                gene_list = pd.read_csv(probe_path.as_posix())
                gene=list(gene_list['SYMBOL'])
            else:
                probe_path=Path('../').resolve().joinpath('src','uni_plus2.txt')
                gene_list = pd.read_csv(probe_path.as_posix())
                gene=list(gene_list['SYMBOL'])
            
            probe_path=Path('../').resolve().joinpath('src','new_human_gene_info.txt')
            info=pd.read_csv(probe_path.as_posix(),sep='\t')
            col=list(info.columns.values)
            col[0]='symbol'
            info.columns=col
                    
            info.index = info['symbol']
            info.index.name = None
            info=info.iloc[:, 1:]
            all_probe=[]
            
            
            for i in indata:
                try:
                    probe_offset.append(gene.index(i))
                    all_probe.append(i)
                    info=info.drop(i,errors='ignore')
                except ValueError:
                    re_symbol=list(set(info.loc[info['alias'].isin([i])].index)) #find whether has alias first
                    if(len(re_symbol)!=0):
                        re_match=Gene.objects.filter(platform__name__in=[pform],symbol__in=re_symbol).order_by('offset') #check the symbol in database or not
                        repeat=len(re_match)
                        if(repeat!=0): #match gene symbol in database      
                            ##print(re_match)
                            for x in re_match:
                                info=info.drop(x.symbol,errors='ignore')
                                probe_offset.append(x.offset)
                                all_probe.append(i+"("+x.symbol+")")
                        else:
                            not_found.append(i)
                    else:
                        not_found.append(i)
                    
            
    
    
    #count the number of group 
    group_counter=1
    check_set=[]
    while True:
        temp_name='dataset_g'+str(group_counter)
        if temp_name in request.POST:
            group_counter=group_counter+1
        else:
            group_counter=group_counter-1
            break
    
    #get binary data
    s_group_dict={}  #store sample
    val=[] #store value get from binary data 
    group_name=[]
    clinic=list(Clinical_Dataset.objects.all().values_list('name',flat=True))
    clline=list(Dataset.objects.all().values_list('name',flat=True))
    #print(clline)
    
    opened_name=[]
    opened_val=[]
    for i in range(1,group_counter+1):
        s_group_dict['g'+str(i)]=[]
        dname='dataset_g'+str(i)
        datasets=request.POST.getlist(dname)
        temp_name='g'+str(i)
        group_name.append(temp_name)  
        a_data=np.array([])
        for dn in datasets:
            if dn=='Sanger Cell Line Project':
                c='select_sanger_g'+str(i)
            elif dn=='NCI60':
                c='select_nci_g'+str(i)
            elif dn=='GSE36133':
                c='select_ccle_g'+str(i)

            if dn in clline:
                ACELL=request.POST.getlist(c)
                s=Sample.objects.filter(dataset_id__name__in=[dn],cell_line_id__name__in=ACELL).order_by('dataset_id').select_related('cell_line_id__name','dataset_id')
                s_group_dict['g'+str(i)]=list(s)+s_group_dict['g'+str(i)]
                goffset=list(s.values_list('offset',flat=True))
                #print(goffset)
                if dn not in opened_name:  #check if the file is opened
                    #print("opend file!!")
                    opened_name.append(dn)
                    if(return_page_flag==1):
                        pth=Path('../').resolve().joinpath('src','gene_'+Dataset.objects.get(name=dn).data_path)
                    else:
                        pth=Path('../').resolve().joinpath('src',Dataset.objects.get(name=dn).data_path)
                    raw_val=np.load(pth.as_posix(),mmap_mode='r')
                    opened_val.append(raw_val)
                    temp=raw_val[np.ix_(probe_offset,list(goffset))]
                    if (len(a_data)!=0 ) and (len(temp)!=0):
                        a_data=np.concatenate((a_data,temp),axis=1)
                    elif (len(temp)!=0):
                        a_data=raw_val[np.ix_(probe_offset,list(goffset))]
                    
                else:
                    temp=opened_val[opened_name.index(dn)][np.ix_(probe_offset,list(goffset))]
                    if (len(a_data)!=0 ) and (len(temp)!=0):
                        a_data=np.concatenate((a_data,temp),axis=1)
                    elif (len(temp)!=0):
                        a_data=opened_val[opened_name.index(dn)][np.ix_(probe_offset,list(goffset))]

                
            elif dn in clinic:
                #print("I am in clinical part")
                com_hists=list(set(request.POST.getlist('primd_'+dn+'_g'+str(i))))    #can I get this by label to reduce number of queries?
                com_hists=[w1 for segments in com_hists for w1 in segments.split('/')]
                prims=com_hists[0::2]
                hists=com_hists[1::2]
                temp=request.POST.getlist('filter_'+dn+'_g'+str(i))
                age=[]
                gender=[]
                ethnic=[]
                grade=[]
                stage=[]
                T=[]
                N=[]
                M=[]
                metas=[]
                for t in temp:
                    if 'stage/' in t:
                        stage.append(t[6:])
                    elif 'gender/' in t:
                        gender.append(t[7:])
                    elif 'ethnic/' in t:
                        ethnic.append(t[7:])
                    elif 'grade/' in t:
                        grade.append(t[6:])
                    elif 'stageT/' in t:
                        T.append(t[7:])
                    elif 'stageN/' in t:
                        N.append(t[7:])
                    elif 'stageM/' in t:
                        M.append(t[7:])
                    elif 'metastatic/' in t:
                        if t[11:]=='False':
                            metas.append(0)
                        else:
                            metas.append(1)
                    else:  #"age/"
                        age.append(t[4:])
                cgoffset=[]
                for x in range(0,len(prims)):
                    s=Clinical_sample.objects.filter(dataset_id__name=dn,primary_site=prims[x],
                    primary_hist=hists[x],
                    age__in=age,
                    gender__in=gender,
                    ethnic__in=ethnic,
                    stage__in=stage,
                    grade__in=grade,
                    stageT__in=T,
                    stageN__in=N,
                    stageM__in=M,
                    metastatic__in=metas,
                    ).select_related('dataset_id').order_by('id')
                    s_group_dict['g'+str(i)]=list(s)+s_group_dict['g'+str(i)]
                    cgoffset+=list(s.values_list('offset',flat=True))
                if dn not in opened_name:  #check if the file is opened
                    #print("opend file!!")
                    opened_name.append(dn)
                    if(return_page_flag==1):
                        pth=Path('../').resolve().joinpath('src','gene_'+Clinical_Dataset.objects.get(name=dn).data_path)
                    else:
                        pth=Path('../').resolve().joinpath('src',Clinical_Dataset.objects.get(name=dn).data_path)
                    raw_val=np.load(pth.as_posix(),mmap_mode='r')
                    opened_val.append(raw_val)
                    temp=raw_val[np.ix_(probe_offset,list(cgoffset))]
                    #print(temp)
                    if (len(a_data)!=0 ) and (len(temp)!=0):
                        a_data=np.concatenate((a_data,temp),axis=1)
                    elif (len(temp)!=0):
                        a_data=raw_val[np.ix_(probe_offset,list(cgoffset))]
                else:
                    temp=opened_val[opened_name.index(dn)][np.ix_(probe_offset,list(cgoffset))]
                    if (len(a_data)!=0 ) and (len(temp)!=0):
                        a_data=np.concatenate((a_data,temp),axis=1)
                    elif (len(temp)!=0):
                        a_data=opened_val[opened_name.index(dn)][np.ix_(probe_offset,list(cgoffset))]
        val.append(a_data.tolist())
        #print(len(val))
        #print(len(val[0]))
        ##print(val)
        
            
    #run the one way ANOVA test or ttest for every probe base on the platform selected    
    express={}
    #logger.info('run ttest or anova')
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
       
    #print("test done")    
    #sort the dictionary with p-value and need to get the expression data again (top20)  
    #presult[all_probe[0]]=float('nan')
    #presult[all_probe[11]]=float('nan')
    #how to deal with all "nan"?
    tempf=pd.DataFrame(list(presult.items()), columns=['probe', 'pvalue'])
    tempf=tempf.replace(to_replace=float('nan'),value=float('+inf'))
    presult=dict(zip(tempf.probe, tempf.pvalue))
    sortkey=sorted(presult,key=presult.get)   #can optimize here
    
    counter=1
    
    
    cell_probe_val=[]
    for w in sortkey: 
        ##print(presult[w],":",w.Probe_id)
        
        if (presult[w]<pro_number):
            cell_probe_val.append([w,presult[w]])
            
            express_mean=np.mean(np.array(express[w]))
            expression.append(list((np.array(express[w]))-express_mean))
            if(return_page_flag==1):
                probe_out.append(w)
            else:
                probe_out.append(w.Probe_id+"("+w.Gene_symbol+")")
            counter+=1
        else:
            break
        
        if counter>=stop_end:
            break

    n_counter=1
    for n in group_name:
        sample_counter=1
        for s in s_group_dict[n]:
            dataset_n=s.dataset_id.name
            if dataset_n=="Sanger Cell Line Project":
                sample_out.append(s.cell_line_id.name+"(SCLP)(group"+str(n_counter)+"-"+str(sample_counter)+")")   
            elif dataset_n in clline:
                #print(s.cell_line_id.name+"("+s.dataset_id.name+")"+"(group"+str(n_counter)+"-"+str(sample_counter)+")")
                sample_out.append(s.cell_line_id.name+"("+s.dataset_id.name+")"+"(group"+str(n_counter)+"-"+str(sample_counter)+")")   
            else:  #what to output for clinical part?
                #print(s.name+"("+s.dataset_id.name+")"+"(group"+str(n_counter)+"-"+str(sample_counter)+")")
                sample_out.append(s.name+"("+s.dataset_id.name+")"+"(group"+str(n_counter)+"-"+str(sample_counter)+")")
            sample_counter+=1
        n_counter+=1
    #logger.info('finish combine output samples')
    sns.set(font="monospace")
    test=pd.DataFrame(data=expression,index=probe_out,columns=sample_out)
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),
    
             'blue': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
    
             'green':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }
    
    my_cmap = LinearSegmentedColormap('my_colormap',cdict,256)
    try:
        g = sns.clustermap(test,cmap=my_cmap)
    
    except ValueError:
        if((return_page_flag==1) and (probe_out==[])):
            probe_out=indata  #probe_out is the rows of heatmap
        return render_to_response('noprobe.html',RequestContext(request,
        {
        'user_probe_flag':user_probe_flag,
        'return_page_flag':return_page_flag,
        'probe_out':probe_out,
        'not_found':not_found
        }))
    
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    if counter>=stop_end:
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=4)
    else:
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=7)
    
    
    
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=270,ha='center')
    
    sid=str(uuid.uuid1())+".png"
    #print(sid)
    P=Path('../').resolve().joinpath('src','static','image',sid)
    assP=Path('../').resolve().joinpath('src','assets','image',sid)
    g.savefig(str(P))
    g.savefig(str(assP))
    #g.savefig("heatmap_for_paper_600.svg")
    #test.to_csv('er_whole_gene.txt',sep='\t')
    file_name=sid
    
    return render_to_response('heatmap.html',RequestContext(request,
    {
    'user_probe_flag':user_probe_flag,
    'return_page_flag':return_page_flag,
    'cell_probe_val':cell_probe_val,
    'file_name':file_name,
    'pro_number':pro_number,
    'not_found':not_found
    }))
    

    
def pca(request):
    
    
    
    propotion=0
    table_propotion=0
    pform=request.POST['data_platform'] #get the platform
    show=request.POST['show_type']      #get the pca show type
    group_counter=1
    cell_line_dict={}
    
    
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

    clinic=list(Clinical_Dataset.objects.all().values_list('name',flat=True))
    clline=list(Dataset.objects.all().values_list('name',flat=True))
    
    all_exist_dataset=[]
    for i in range(1,group_counter+1):
        dname='dataset_g'+str(i)
        all_exist_dataset=all_exist_dataset+request.POST.getlist(dname)
    all_exist_dataset=list(set(all_exist_dataset))
    
    all_base=[0]
    for i in range(0,len(all_exist_dataset)-1):
        if all_exist_dataset[i] in clline:
            all_base.append(all_base[i]+Sample.objects.filter(dataset_id__name__in=[all_exist_dataset[i]]).count())
        else:
            all_base.append(all_base[i]+Clinical_sample.objects.filter(dataset_id__name__in=[all_exist_dataset[i]]).count())
    
    all_c=[]
    for i in range(1,group_counter+1):
        s_group_dict['g'+str(i)]=[]
        offset_group_dict['g'+str(i)]=[]
        cell_line_dict['g'+str(i)]=[]
        dname='dataset_g'+str(i)
        datasets=request.POST.getlist(dname)
        group_name.append('g'+str(i))

        goffset_nci=[]
        goffset_gse=[]
        
        for dn in datasets:
            if dn=='Sanger Cell Line Project':
                c='select_sanger_g'+str(i)
            elif dn=='NCI60':
                c='select_nci_g'+str(i)
            elif dn=='GSE36133':
                c='select_ccle_g'+str(i)

            if dn in clline:
                temp=list(set(request.POST.getlist(c)))
                if 'd_sample' in show:
                    if all_c==[]:
                        all_c=all_c+temp
                        uni=temp
                    else:
                        uni=list(set(temp)-set(all_c))
                        all_c=all_c+uni
                else:
                    uni=list(temp)          #do not filter duplicate input only when select+centroid
                s=Sample.objects.filter(cell_line_id__name__in=uni,dataset_id__name__in=[dn]).order_by('dataset_id'
                ).select_related('cell_line_id__name','cell_line_id__primary_site','cell_line_id__primary_hist','dataset_id','dataset_id__name')
                
                cell_line_dict['g'+str(i)]=cell_line_dict['g'+str(i)]+list(s.values_list('cell_line_id__name',flat=True))
                s_group_dict['g'+str(i)]=s_group_dict['g'+str(i)]+list(s)
                offset_group_dict['g'+str(i)]=offset_group_dict['g'+str(i)]+list(np.add(list(s.values_list('offset',flat=True)),all_base[all_exist_dataset.index(dn)]))
                

            else: #dealing with clinical sample datasets
                com_hists=list(set(request.POST.getlist('primd_'+dn+'_g'+str(i))))    #can I get this by label to reduce number of queries?
                com_hists=[w1 for segments in com_hists for w1 in segments.split('/')]
                prims=com_hists[0::2]
                hists=com_hists[1::2]
                temp=request.POST.getlist('filter_'+dn+'_g'+str(i))
                age=[]
                gender=[]
                ethnic=[]
                grade=[]
                stage=[]
                T=[]
                N=[]
                M=[]
                metas=[]
                for t in temp:
                    if 'stage/' in t:
                        stage.append(t[6:])
                    elif 'gender/' in t:
                        gender.append(t[7:])
                    elif 'ethnic/' in t:
                        ethnic.append(t[7:])
                    elif 'grade/' in t:
                        grade.append(t[6:])
                    elif 'stageT/' in t:
                        T.append(t[7:])
                    elif 'stageN/' in t:
                        N.append(t[7:])
                    elif 'stageM/' in t:
                        M.append(t[7:])
                    elif 'metastatic/' in t:
                        if t[11:]=='False':
                            metas.append(0)
                        else:
                            metas.append(1)
                    else:  #"age/"
                        age.append(t[4:])
                for x in range(0,len(prims)):
                    s=Clinical_sample.objects.filter(dataset_id__name=dn,primary_site=prims[x],
                    primary_hist=hists[x],
                    age__in=age,
                    gender__in=gender,
                    ethnic__in=ethnic,
                    stage__in=stage,
                    grade__in=grade,
                    stageT__in=T,
                    stageN__in=N,
                    stageM__in=M,
                    metastatic__in=metas,
                    ).select_related('dataset_id').order_by('id')
                    s_group_dict['g'+str(i)]=s_group_dict['g'+str(i)]+list(s)
                    cell_line_dict['g'+str(i)]=cell_line_dict['g'+str(i)]+list(s.values_list('name',flat=True))
                    offset_group_dict['g'+str(i)]=offset_group_dict['g'+str(i)]+list(np.add(list(s.values_list('offset',flat=True)),all_base[all_exist_dataset.index(dn)]))

    

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
    
    if 'd_sample' in show:
        if((len(all_sample))<4):
            error_reason='You should have at least 4 samples for PCA. The samples are not enough.<br />'\
            'The number of samples you selected is '+str(len(all_sample))+'.'
            return render_to_response('pca_error.html',RequestContext(request,
            {
            'error_reason':mark_safe(json.dumps(error_reason)),
            }))
    
    for i in all_sample:
        sample_counter[i.name]=1
        if str(type(i))=="<class 'probes.models.Sample'>":
            ##print("i am sample!!")
            cell_object.append(i.cell_line_id)
        else:
            ##print("i am clinical!!")
            cell_object.append(i)
    #delete nan, transpose matrix
    ##open file
    
    for x in range(0,len(all_exist_dataset)):
        if all_exist_dataset[x] in clline:
            pth=Path('../').resolve().joinpath('src',Dataset.objects.get(name=all_exist_dataset[x]).data_path)
        else:
            pth=Path('../').resolve().joinpath('src',Clinical_Dataset.objects.get(name=all_exist_dataset[x]).data_path)
        if x==0:
            val=np.load(pth.as_posix())
        else:
            val=np.hstack((val, np.load(pth.as_posix())))#combine together
    if 'd_sample' in show:
        val=val[:,all_offset]
        #val=val[~np.isnan(val).any(axis=1)]
        val=np.transpose(val)
   
        
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
    if(len(all_exist_dataset)==1):
        n=3  #need to fix to the best one #need to fix proportion 
    else:
        n=4
    #logger.info('pca show')
    if 'd_sample' in show:
        #count the pca first
        pca= PCA(n_components=n)
        Xval = pca.fit_transform(val[:,:])  #cannot get Xval with original all_offset any more
        ratio_temp=pca.explained_variance_ratio_
        propotion=sum(ratio_temp[n-3:n])
        table_propotion=sum(ratio_temp[0:n])
        ##print(Xval)
        ##print(all_cellline)
        ##print(all_sample)
        max=0
        min=10000000000
        out_group=[]
        exist_cell={}#cell line object:counter
        for g in range(1,group_counter+1):
              
            output_cell={}
            check={}
            for s in range(g_s_counter[g-1],g_s_counter[g]):
                
                if str(type(all_sample[s]))=="<class 'probes.models.Sample'>":
                    cell=all_sample[s].cell_line_id
                else:
                    cell=all_sample[s]
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
                                distance=np.linalg.norm(Xval[i][n-3:n]-Xval[s][n-3:n])
                                if distance<min:
                                    min=distance
                                if distance>max:
                                    max=distance
                                output_cell[cell][1].append([all_cellline[s]+'('+str(exist_cell[cell])+')'
                                ,all_sample[s].name,all_sample[s].dataset_id.name,all_cellline[i],all_sample[i].name,all_sample[i].dataset_id.name,distance,cell_object[i]])
                                check[all_sample[s].name].append(all_sample[i].name)
                        except KeyError:
                            distance=np.linalg.norm(Xval[i][n-3:n]-Xval[s][n-3:n])
                            if distance<min:
                                min=distance
                            if distance>max:
                                max=distance
                            output_cell[cell][1].append([all_cellline[s]+'('+str(exist_cell[cell])+')'
                            ,all_sample[s].name,all_sample[s].dataset_id.name,all_cellline[i],all_sample[i].name,all_sample[i].dataset_id.name,distance,cell_object[i]])
                            check[all_sample[s].name].append(all_sample[i].name)
                        
                if(g==1):
                    name1.append(all_cellline[s]+'('+str(exist_cell[cell])+')'+'<br>'+all_sample[s].name)
                    X1.append(round(Xval[s][n-3],5))
                    Y1.append(round(Xval[s][n-2],5))
                    Z1.append(round(Xval[s][n-1],5))
                elif(g==2):
                    name2.append(all_cellline[s]+'('+str(exist_cell[cell])+')'+'<br>'+all_sample[s].name)
                    X2.append(round(Xval[s][n-3],5))
                    Y2.append(round(Xval[s][n-2],5))
                    Z2.append(round(Xval[s][n-1],5))
                elif(g==3):
                    name3.append(all_cellline[s]+'('+str(exist_cell[cell])+')'+'<br>'+all_sample[s].name)
                    X3.append(round(Xval[s][n-3],5))
                    Y3.append(round(Xval[s][n-2],5))
                    Z3.append(round(Xval[s][n-1],5))
                elif(g==4):
                    name4.append(all_cellline[s]+'('+str(exist_cell[cell])+')'+'<br>'+all_sample[s].name)
                    X4.append(round(Xval[s][n-3],5))
                    Y4.append(round(Xval[s][n-2],5))
                    Z4.append(round(Xval[s][n-1],5))
                elif(g==5):
                    name5.append(all_cellline[s]+'('+str(exist_cell[cell])+')'+'<br>'+all_sample[s].name)
                    X5.append(round(Xval[s][n-3],5))
                    Y5.append(round(Xval[s][n-2],5))
                    Z5.append(round(Xval[s][n-1],5))
            dictlist=[]
            for key, value in output_cell.items():
                temp = [value]
                dictlist+=temp
            output_cell=list(dictlist)
            out_group.append([g,output_cell])
        element_counter=0
        #[g,[[group_cell_line,[paired_cellline,......,]],[],[]]]
        for i in out_group:
            for temp_list in i[1]:
                element_counter+=len(temp_list[1])
                for temp in temp_list[1]:
                    ##print(temp)
                    temp[3]=temp[3]+'('+str(sample_counter[temp[4]])+')'
        
        return_html='pca.html'
    else:
        #This part is for centroid display
        
        return_html='pca_center.html'
        element_counter=0
        #val=val[~np.isnan(val).any(axis=1)]  #bottle neck???
        
        #This part is for select cell line base on dataset,count centroid base on the dataset
        #group中的cell line為單位來算重心
        #logger.info('pca show centroid with selection')
        location_dict={} #{group number:[[cell object,dataset,new location]]}
        combined=[]
        sample_list=[]
        pca_index=np.array(pca_index)
        X_val=[]
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
                selected_val=val[:,a_all_offset[total_offset]]
                selected_val=np.transpose(selected_val)
                new_loca=(np.mean(selected_val,axis=0,dtype=np.float64,keepdims=True)).tolist()[0]
                
                
                a_sample=np.array(all_sample)
                selected_sample=a_sample[total_offset]
                
                if list(selected_sample) in sample_list:   #to prevent two different colors in different group
                    continue
                else:
                    sample_list.append(list(selected_sample))
                
                
                ##print(selected_sample)
                
                
                
                d_temp=[]
                for s in selected_sample:
                    d_temp.append(s.dataset_id.name)
                dataset_dict[c]="/".join(list(set(d_temp)))    
                ##print(dataset_dict[c])
                X_val.append(new_loca)
                location_dict['g'+str(i)].append([c,dataset_dict[c],len(X_val)-1])  #the last part is the index to get pca result from new_val
                combined.append([c,dataset_dict[c],len(X_val)-1])  #all cell line, do not matter order
        
        #run the pca
        ##print(len(X_val))
        if((len(X_val))<4):
            error_reason='Since the display method is [centroid], you should have at least 4 dots for PCA. The dots are not enough.<br />'\
            'The number of centroid dots you selected is '+str(len(X_val))+'.'
            return render_to_response('pca_error.html',RequestContext(request,
            {
            'error_reason':mark_safe(json.dumps(error_reason)),
            }))
        X_val=np.matrix(X_val)
        pca= PCA(n_components=n)
        new_val = pca.fit_transform(X_val[:,:])  #cannot get Xval with original offset any more
        ratio_temp=pca.explained_variance_ratio_
        propotion=sum(ratio_temp[n-3:n])
        table_propotion=sum(ratio_temp[0:n])
        ##print(new_val)


        
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
                            distance=np.linalg.norm(np.array(new_val[group_c[2]][n-3:n])-np.array(new_val[temp_list[2]][n-3:n]))
                            if distance==0:
                                continue
                            if distance<min:
                                min=distance
                            if distance>max:
                                max=distance
                            output_cell[len(output_cell)-1][1].append([cell,group_c[1],temp_list[0],temp_list[1],distance])
                            element_counter+=1
                            exist_cell[key_string].append(temp_string)
                    except KeyError:
                        distance=np.linalg.norm(np.array(new_val[group_c[2]][n-3:n])-np.array(new_val[temp_list[2]][n-3:n]))
                        if distance==0:
                            continue
                        if distance<min:
                            min=distance
                        if distance>max:
                            max=distance
                        output_cell[len(output_cell)-1][1].append([cell,group_c[1],temp_list[0],temp_list[1],distance])
                        element_counter+=1
                        exist_cell[key_string].append(temp_string)
                if(g==1):
                    name1.append(cell.name+'<br>'+group_c[1])
                    X1.append(round(new_val[group_c[2]][n-3],5))
                    Y1.append(round(new_val[group_c[2]][n-2],5))
                    Z1.append(round(new_val[group_c[2]][n-1],5))
                elif(g==2):
                    name2.append(cell.name+'<br>'+group_c[1])
                    X2.append(round(new_val[group_c[2]][n-3],5))
                    Y2.append(round(new_val[group_c[2]][n-2],5))
                    Z2.append(round(new_val[group_c[2]][n-1],5))
                elif(g==3):
                    name3.append(cell.name+'<br>'+group_c[1])
                    X3.append(round(new_val[group_c[2]][n-3],5))
                    Y3.append(round(new_val[group_c[2]][n-2],5))
                    Z3.append(round(new_val[group_c[2]][n-1],5))
                elif(g==4):
                    name4.append(cell.name+'<br>'+group_c[1])
                    X4.append(round(new_val[group_c[2]][n-3],5))
                    Y4.append(round(new_val[group_c[2]][n-2],5))
                    Z4.append(round(new_val[group_c[2]][n-1],5))
                elif(g==5):
                    name5.append(cell.name+'<br>'+group_c[1])
                    X5.append(round(new_val[group_c[2]][n-3],5))
                    Y5.append(round(new_val[group_c[2]][n-2],5))
                    Z5.append(round(new_val[group_c[2]][n-1],5)) 
            out_group.append([g,output_cell])
    #logger.info('end pca')    
    if(element_counter>show_row):
        big_flag=1    
        sid=str(uuid.uuid1())+".csv"
        if(return_html=='pca.html'):
            dataset_header=['Group Cell Line/Clinical Sample','Sample Name','Primary Site','Primary Histology'
                ,'Dataset','Paired Cell Line name/Clinical Sample','Sample Name','Primary Site','Primary Histology','Dataset','Distance']
        else:
            dataset_header=['Group Cell Line/Clinical Sample','Primary Site','Primary Histology'
                ,'Dataset','Paired Cell Line name/Clinical Sample','Primary Site','Primary Histology','Dataset','Distance']
        P=Path('../').resolve().joinpath('src','static','csv',sid)
        assP=Path('../').resolve().joinpath('src','assets','csv',sid)
        with open(str(P), "w", newline='') as f:
            writer = csv.writer(f)
            for index,output_cell in out_group:
                writer.writerows([["Group_"+str(index)]])
                writer.writerows([dataset_header])
                for cell_line,b in output_cell:
                    temp_b=[]
                    if(return_html=='pca.html'):
                        for group_cell,sn,dset,cname,sname,setname,dis,cell_object in b:
                            temp_b.append([group_cell,sn,cell_line.primary_site,cell_line.primary_hist,dset,cname
                            ,sname,cell_object.primary_site,cell_object.primary_hist,setname,dis])
                    else:
                        for group_cell,group_dataset,paired_cell,paired_dataset,dis in b:
                            temp_b.append([group_cell.name,group_cell.primary_site,group_cell.primary_hist,group_dataset
                            ,paired_cell.name,paired_cell.primary_site,paired_cell.primary_hist,paired_dataset,dis])
                    writer.writerows(temp_b)
        #print('write first file done')
        with open(str(assP), "w", newline='') as ff:
            writer = csv.writer(ff)
            for index,output_cell in out_group:
                writer.writerows([["Group_"+str(index)]])
                writer.writerows([dataset_header])
                for cell_line,b in output_cell:
                    temp_b=[]
                    if(return_html=='pca.html'):
                        for group_cell,sn,dset,cname,sname,setname,dis,cell_object in b:
                            temp_b.append([group_cell,sn,cell_line.primary_site,cell_line.primary_hist,dset,cname
                            ,sname,cell_object.primary_site,cell_object.primary_hist,setname,dis])
                    else:
                        for group_cell,group_dataset,paired_cell,paired_dataset,dis in b:
                            temp_b.append([group_cell.name,group_cell.primary_site,group_cell.primary_hist,group_dataset
                            ,paired_cell.name,paired_cell.primary_site,paired_cell.primary_hist,paired_dataset,dis])
                    writer.writerows(temp_b)
        #print('write second file done')
        data_file_name=sid
    else:
        big_flag=0
        data_file_name=0
    return render_to_response(return_html,RequestContext(request,
    {
    'min':min,'max':max,
    'out_group':out_group,
    'propotion':propotion,
    'big_flag':big_flag,
    'data_file_name':data_file_name,
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
    d=Dataset.objects.all()
    d_name=list(d.values_list('name',flat=True))
    datasets=[]  #[[dataset_name,[[primary_site,[cell line]]]]
    an=[]
    for i in d_name:
        
        if i=="Sanger Cell Line Project":
            alias='sanger'
        elif i=="NCI60":
            alias='nci'
        elif i=="GSE36133":
            alias='gse'
        else:
            alias=i
        
        an.append(alias)
        sample=Sample.objects.filter(dataset_id__name=i).order_by('cell_line_id__primary_site').select_related('cell_line_id')
        datasets.append([i,alias,list(sample),[]])
        sites=list(sample.values_list('cell_line_id__primary_site',flat=True))
        hists=list(sample.values_list('cell_line_id__name',flat=True))
        
        dis_prim=list(sample.values_list('cell_line_id__primary_site',flat=True).distinct())
        hists=list(hists)
        id_counter=0
       
        for p in range(0,len(dis_prim)):
            temp=sites.count(dis_prim[p])
            datasets[-1][3].append([dis_prim[p],list(set(hists[id_counter:id_counter+temp]))])
            id_counter+=temp

    return render(request, 'cellline_microarray.html', {
        'an':mark_safe(json.dumps(an)),
        'd_name':d_name,
        'datasets':datasets, 
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
    
def clinical_search(request):

    norm_name=[request.POST['normalize']]  #get the normalize gene name
    #f_type=['age','gender','ethnic','grade','stage','stageT','stageN','stageM','metastatic']
    age=[]
    gender=[]
    ethnic=[]
    grade=[]
    stage=[]
    T=[]
    N=[]
    M=[]
    metas=[]
    #get the probe/gene/id keywords
    if 'keyword' in request.POST and request.POST['keyword'] != '':
        words = request.POST['keyword']
        words = list(set(words.split()))
    else:
        return HttpResponse("<p>where is your keyword?</p>")
    
    plus2_rank=np.load('ranking_u133plus2.npy')   #open only plus2 platform rank
    sample_probe_val_pairs=[]  #for output
    
    if 'gtype' in request.POST and request.POST['gtype'] == 'probeid':
        gene = ProbeID.objects.filter(platform__name__in=["PLUS2"]).filter(Probe_id__in=words).order_by('id') 
        probe=list(gene.values_list('offset',flat=True))
        ##print(gene)
    elif 'gtype' in request.POST and request.POST['gtype'] == 'symbol':
        gene = ProbeID.objects.filter(platform__name__in=["PLUS2"]).filter(Gene_symbol__in=words).order_by('id') 
        probe=list(gene.values_list('offset',flat=True))
    else:
        gene = ProbeID.objects.filter(platform__name__in=["PLUS2"]).filter(Entrez_id__in=words).order_by('id') 
        probe=list(gene.values_list('offset',flat=True))
           
    
    if request.POST['clinical_method'] == 'prim_dataset':
        if 'dataset' in request.POST and request.POST['dataset'] != '':
            datas=request.POST.getlist('dataset')
            
    else:
        d=Clinical_Dataset.objects.all()
        datas=d.values_list('name',flat=True)
        com_hists=list(set(request.POST.getlist('primhist')))
        com_hists=[w1 for segments in com_hists for w1 in segments.split('/')]
        prims=com_hists[0::2]
        hists=com_hists[1::2]
        
        temp=request.POST.getlist('filter_primh')
        for i in temp:
            if 'stage/' in i:
                stage.append(i[6:])                
            elif 'gender/' in i:
                gender.append(i[7:])
            elif 'ethnic/' in i:
                ethnic.append(i[7:])
            elif 'grade/' in i:
                grade.append(i[6:])
            elif 'stageT/' in i:
                T.append(i[7:])
            elif 'stageN/' in i:
                N.append(i[7:])
            elif 'stageM/' in i:
                M.append(i[7:])
            elif 'metastatic/' in i:
                if i[11:]=='False':
                    metas.append(0)
                else:
                    metas.append(1)
            else:  #"age/"
                age.append(i[4:])
                                
    
    for sets in datas:
        samples=[]
        offset=[]        
        if request.POST['clinical_method'] == 'prim_dataset':
            com_hists=list(set(request.POST.getlist('primd_'+sets)))    #can I get this by label to reduce number of queries?
            com_hists=[w1 for segments in com_hists for w1 in segments.split('/')]
            prims=com_hists[0::2]
            hists=com_hists[1::2]
            temp=request.POST.getlist('filter_'+sets)
            age=[]
            gender=[]
            ethnic=[]
            grade=[]
            stage=[]
            T=[]
            N=[]
            M=[]
            metas=[]
            for i in temp:
                if 'stage/' in i:
                    stage.append(i[6:])
                elif 'gender/' in i:
                    gender.append(i[7:])
                elif 'ethnic/' in i:
                    ethnic.append(i[7:])
                elif 'grade/' in i:
                    grade.append(i[6:])
                elif 'stageT/' in i:
                    T.append(i[7:])
                elif 'stageN/' in i:
                    N.append(i[7:])
                elif 'stageM/' in i:
                    M.append(i[7:])
                elif 'metastatic/' in i:
                    if i[11:]=='False':
                        metas.append(0)
                    else:
                        metas.append(1)
                else:  #"age/"
                    age.append(i[4:])
        for i in range(0,len(prims)):
            #metas=[bool(x) for x in metas]
            s=Clinical_sample.objects.filter(dataset_id__name=sets,primary_site=prims[i],
            primary_hist=hists[i],
            age__in=age,
            gender__in=gender,
            ethnic__in=ethnic,
            stage__in=stage,
            grade__in=grade,
            stageT__in=T,
            stageN__in=N,
            stageM__in=M,
            metastatic__in=metas
            ).select_related('dataset_id').order_by('id')
            samples+=list(s) 
            offset+=list(s.values_list('offset',flat=True))
            ##print(s)
        pth=Path('../').resolve().joinpath('src',Clinical_Dataset.objects.get(name=sets).data_path)
        val=np.load(pth.as_posix(),mmap_mode='r')
        
        norm_probe=ProbeID.objects.filter(platform__name__in=["PLUS2"]).filter(Gene_symbol__in=norm_name).order_by('id') 
        probe_offset=list(norm_probe.values_list('offset',flat=True))
        temp=val[np.ix_(probe_offset,offset)]
        norm=np.mean(temp,axis=0, dtype=np.float64,keepdims=True)
        
        
        # Make a generator to generate all (cell, probe, val) pairs
        if(len(gene)!=0 and len(samples)!=0):
            raw_test=val[np.ix_(probe,offset)]
            normalize=np.subtract(raw_test,norm)#dimension different!!!!
           
            sample_probe_val_pairs += [
                (c, p, raw_test[probe_ix, cell_ix],54614-np.where(plus2_rank==raw_test[probe_ix, cell_ix])[0],normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(gene)
                for cell_ix, c in enumerate(samples)
            ]
            
               
    return render(request, 'clinical_search.html', {
        'sample_probe_val_pairs': sample_probe_val_pairs,
        
    })
    
    
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
                ps_id=str(Platform.objects.filter(name__in=["U133A"])[0].id)
            if 'NCI60' in datas:
                nci_flag=1
                NCI=list(set(request.POST.getlist('select_nci')))
                ncisamples=Sample.objects.filter(dataset_id__name__in=['NCI60']).select_related('cell_line_id','dataset_id').order_by('id') 
                ncicell=ncisamples.filter(cell_line_id__name__in=NCI).order_by('id') 
                ncioffset=list(ncicell.values_list('offset',flat=True))
                pn_id=str(Platform.objects.filter(name__in=["PLUS2"])[0].id)
            if 'GSE36133' in datas:
                gse_flag=1
                GSE=list(set(request.POST.getlist('select_gse')))
                CCsamples=Sample.objects.filter(dataset_id__name__in=['GSE36133']).select_related('cell_line_id','dataset_id').order_by('id') 
                CCcell=CCsamples.filter(cell_line_id__name__in=GSE).order_by('id') 
                CCoffset=list(CCcell.values_list('offset',flat=True))
                pn_id=str(Platform.objects.filter(name__in=["PLUS2"])[0].id)
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
        ##print(nci_norm)   
    else:
        nci_norm=0.0  #if / should = 1
    if gse_flag==1:
        CC_g=ProbeID.objects.filter(platform__in=pn_id).filter(Gene_symbol__in=norm_name).order_by('id') 
        CC_probe_offset=list(CC_g.values_list('offset',flat=True))
        temp=gse_val[np.ix_(CC_probe_offset,CCoffset)]
        CC_norm=np.mean(temp,axis=0, dtype=np.float64,keepdims=True)
        ##print(CC_norm)
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
                (c, p, raw_test[probe_ix, cell_ix],22216-np.where(u133a_rank==raw_test[probe_ix, cell_ix])[0],normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(gene)
                for cell_ix, c in enumerate(cell)
            )
            
        else:
            cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(ncicell)!=0):
            nci_raw_test=nci_val[np.ix_(nciprobe_offset,ncioffset)]
            nci_normalize=np.subtract(nci_raw_test,nci_norm)
            nci_cell_probe_val_pairs = (
                (c, p, nci_raw_test[probe_ix, cell_ix],54614-np.where(plus2_rank==nci_raw_test[probe_ix, cell_ix])[0],nci_normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(ncigene)
                for cell_ix, c in enumerate(ncicell)
            )
            
        else:
            nci_cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(CCcell)!=0):
            CC_raw_test=gse_val[np.ix_(nciprobe_offset,CCoffset)]
            CC_normalize=np.subtract(CC_raw_test,CC_norm)
            CC_cell_probe_val_pairs = (
                (c, p, CC_raw_test[probe_ix, cell_ix],54614-np.where(plus2_rank==CC_raw_test[probe_ix, cell_ix])[0],CC_normalize[probe_ix, cell_ix])                        
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
                (c, p, raw_test[probe_ix, cell_ix],22216-np.where(u133a_rank==raw_test[probe_ix, cell_ix])[0],normalize[probe_ix, cell_ix])                        
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
                (c, p, CC_raw_test[probe_ix, cell_ix],54614-np.where(plus2_rank==CC_raw_test[probe_ix, cell_ix])[0],CC_normalize[probe_ix, cell_ix])                        
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
                (c, p, raw_test[probe_ix, cell_ix],22216-np.where(u133a_rank==raw_test[probe_ix, cell_ix])[0],normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(gene)
                for cell_ix, c in enumerate(cell)
            )
            
        else:
            cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(ncicell)!=0):
            nci_raw_test=nci_val[np.ix_(nciprobe_offset,ncioffset)]
            nci_normalize=np.subtract(nci_raw_test,nci_norm)
            nci_cell_probe_val_pairs = (
                (c, p, nci_raw_test[probe_ix, cell_ix],54614-np.where(plus2_rank==nci_raw_test[probe_ix, cell_ix])[0],nci_normalize[probe_ix, cell_ix])                        
                for probe_ix, p in enumerate(ncigene)
                for cell_ix, c in enumerate(ncicell)
            )
            
        else:
            nci_cell_probe_val_pairs =()
            
        if(len(ncigene)!=0 and len(CCcell)!=0):
            CC_raw_test=gse_val[np.ix_(nciprobe_offset,CCoffset)]
            CC_normalize=np.subtract(CC_raw_test,CC_norm)
            CC_cell_probe_val_pairs = (
                (c, p, CC_raw_test[probe_ix, cell_ix],54614-np.where(plus2_rank==CC_raw_test[probe_ix, cell_ix])[0],CC_normalize[probe_ix, cell_ix])                        
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
