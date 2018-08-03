%process data with eeglab: mainly to remve line noise
%NOTE: this takes a while

for i=1:length(eeg)
    for j=1:length(eeg(i).sztype)
        datacube=eeg(i).szraw(j).datacube;
        datacube=permute(datacube,[2 1 3]);
        datacube=pop_importdata('dataformat','array','nbchan',size(datacube,1),'data',....
            'datacube','srate',512,'pnts',2560,'xmin',0);
        datacube=pop_cleanline(datacube, 'bandwidth',2,'chanlist',[1:size(datacube.data,1)],...
            'computepower',0,'linefreqs',[60 120] ,'normSpectrum',0,'p',0.01,...
            'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',...
            100,'verb',0,'winsize',5,'winstep',5);
        
        datacube=datacube.data;
        eeg(i).szraw(j).datacube=permute(datacube,[2 1 3]);
        
    end
end
%%
for i=1:length(eeg)
    datacube=eeg(i).bslraw.datacube;
    datacube=permute(datacube,[2 1 3]);
    datacube=pop_importdata('dataformat','array','nbchan',size(datacube,1),'data',....
        'datacube','srate',512,'pnts',2560,'xmin',0);
    datacube=pop_cleanline(datacube, 'bandwidth',2,'chanlist',[1:size(datacube.data,1)],...
        'computepower',0,'linefreqs',[60 120] ,'normSpectrum',0,'p',0.01,...
        'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',...
        100,'verb',0,'winsize',5,'winstep',5);
    
    datacube=datacube.data;
    eeg(i).bslraw.datacube=permute(datacube,[2 1 3]);
    
end
%% recalculate metrics based on cleaned data
for m=1:length(eeg)
    for k=1:size(eeg(m).szlabel,1)
        adj_cube=zeros(length(eeg(m).leads),length(eeg(m).leads),5);
        for j=1:5
            adj_cube(:,:,j)=corr(eeg(m).szraw(k).datacube(:,:,j));
        end
        eeg(m).szraw(k).adj_cube=adj_cube;
    end
    
end
%%
for m=1:length(eeg)
    adj_cube=zeros(length(eeg(m).leads),length(eeg(m).leads),size(eeg(m).bslraw.datacube,3));
    for j=1:size(eeg(m).bslraw.datacube,3)
        adj_cube(:,:,j)=corr(eeg(m).bslraw.datacube(:,:,j));
    end
    eeg(m).bslraw.adj_cube=adj_cube;
    
end
%%
totpt=num2str(length(eeg));
cutoff_method=1;
avg_edg=5;
for m=1:length(eeg)
    ptnum=num2str(m);
    disp(['Calculating metrics for ',ptnum,'/',totpt]);
    
    for k=1:length(eeg(m).szraw)
        adj_cube=eeg(m).szraw(k).adj_cube;
        if cutoff_method==1
            num_nodes=size(adj_cube,1);
            templist=reshape(adj_cube,[num_nodes^2,size(adj_cube,3)]); %makes each correlation matrix into a single row
            templist=sort(templist,1,'descend');
            min_wght=templist(num_nodes+2*avg_edg*num_nodes,:);
            
            gg=bsxfun(@ge,adj_cube,reshape(min_wght,[1 1 5])); %creates logical matrix containing 1 wherever the correlation is greater than the min_wght
            adj_cube(~gg)=0; %keeps only values greater than min_wght
            
        else
            
            num_nodes=size(adj_cube,1);
            gg=adj_cube>=min_wght;
            adj_cube(~gg)=0;
            
        end
        
        adj_cube(adj_cube>0)=1; %sets all nonzero values to 1
        eeg(m).szraw(k).adj_cube=adj_cube;
        
        %eeg(m).sz(k).metric=zeros(num_nodes,5,4);
    end
    metric=arrayfun(@(x) graphmetric(x.adj_cube),eeg(m).szraw,'UniformOutput',0);
    conn_cube=arrayfun(@(x) conn(x.adj_cube),eeg(m).szraw,'UniformOutput',0);
    
    for k=1:length(eeg(m).szraw)
        eeg(m).szraw(k).metric=metric{k};
        eeg(m).szraw(k).conn_cube=conn_cube{k};
    end
    [smlwrld,dataCL]=arrayfun(@(x) smi(x.conn_cube,100),eeg(m).szraw,'UniformOutput',0);
    snczb=arrayfun(@(x) synch(x.conn_cube)',eeg(m).szraw,'UniformOutput',0);
    glob_mat=arrayfun(@(x,y,z) horzcat(x{:},y{:},z{:}),smlwrld,dataCL,snczb,'UniformOutput',0);
    for k=1:length(eeg(m).szraw)
        eeg(m).szraw(k).glob_mat=glob_mat{k};
    end
    
    
end
%%
totpt=num2str(length(eeg));
cutoff_method=1;
avg_edg=5;
for m=1:length(eeg)
    ptnum=num2str(m);
    disp(['Calculating metrics for ',ptnum,'/',totpt]);
    
    adj_cube=eeg(m).bslraw.adj_cube;
    if cutoff_method==1
        num_nodes=size(adj_cube,1);
        templist=reshape(adj_cube,[num_nodes^2,size(adj_cube,3)]); %makes each correlation matrix into a single row
        templist=sort(templist,1,'descend');
        min_wght=templist(num_nodes+2*avg_edg*num_nodes,:);
        
        gg=bsxfun(@ge,adj_cube,reshape(min_wght,[1 1 length(min_wght)])); %creates logical matrix containing 1 wherever the correlation is greater than the min_wght
        adj_cube(~gg)=0; %keeps only values greater than min_wght
        
    else
        
        num_nodes=size(adj_cube,1);
        gg=adj_cube>=min_wght;
        adj_cube(~gg)=0;
        
    end
    
    adj_cube(adj_cube>0)=1; %sets all nonzero values to 1
    eeg(m).bslraw.adj_cube=adj_cube;
    
    %eeg(m).sz(k).metric=zeros(num_nodes,5,4);
    metric=arrayfun(@(x) graphmetric(x.adj_cube),eeg(m).bslraw,'UniformOutput',0);
    conn_cube=arrayfun(@(x) conn(x.adj_cube),eeg(m).bslraw,'UniformOutput',0);
    
    eeg(m).bslraw.metric=metric{:};
    eeg(m).bslraw.conn_cube=conn_cube{:};
    
    [smlwrld,dataCL]=arrayfun(@(x) smi(x.conn_cube,100),eeg(m).bslraw,'UniformOutput',0);
    snczb=arrayfun(@(x) synch(x.conn_cube)',eeg(m).bslraw,'UniformOutput',0);
    glob_mat=arrayfun(@(x,y,z) horzcat(x{:},y{:},z{:}),smlwrld,dataCL,snczb,'UniformOutput',0);
    eeg(m).bslraw.glob_mat=glob_mat;
    
    
end

%% calculate PDC and DTF
freq=exp(0.41:0.02:5.41);
for m=1:length(eeg)
    ptnum=num2str(m);
    for k=1:length(eeg(m).sztype)
        sznum=num2str(k);
        disp(['Calculating metrics for Sz ',sznum,' Pt ',ptnum]);
        datacube=eeg(m).szraw(k).datacube;
        nodes=size(datacube,2);
        ord=3;
        pdc=zeros(nodes,nodes,length(freq),size(datacube,3));
        dtf=pdc;
        for j=1:size(datacube,3)
            [~,A]=arfit(datacube(:,:,j),ord,ord);
            %[pdc(:,:,:,j),dtf(:,:,:,j)]=PDC_DTF_matrix(A,ord,512,256,256);
            [pdc(:,:,:,j),dtf(:,:,:,j)]=PDC_DTF_matrix2(A,ord,512,freq);
        end
        eeg(m).szraw(k).pdc=pdc;
        eeg(m).szraw(k).dtf=dtf;
        
    end
end
%%
freq=exp(0.41:0.02:5.41);
parfor m=1:length(eeg)
    ptnum=num2str(m);
    %for k=1:length(eeg(m).sztype)
        %sznum=num2str(k);
        disp(['Calculating metrics for Pt ',ptnum]);
        datacube=eeg(m).bslraw.datacube;
        nodes=size(datacube,2);
        ord=3;
        pdc=zeros(nodes,nodes,length(freq),size(datacube,3));
        dtf=pdc;
        for j=1:size(datacube,3)
            [~,A]=arfit(datacube(:,:,j),ord,ord);
            [pdc(:,:,:,j),dtf(:,:,:,j)]=PDC_DTF_matrix2(A,ord,512,freq);
        end
        eeg(m).bslraw.pdc=pdc;
        eeg(m).bslraw.dtf=dtf;
    %end
end


%%
cutoff_method=1;
avg_edg=5;
totpt=num2str(length(eeg));
for m=1:length(eeg);
    ptnum=num2str(m);
    disp(['Calculating metrics for ',ptnum,'/',totpt]);
    
    for k=1:length(eeg(m).szraw)
        adj_cube=squeeze(eeg(m).szraw(k).pdc(:,:,2,:));
        if cutoff_method==1
            num_nodes=size(adj_cube,1);
            templist=reshape(adj_cube,[num_nodes^2,size(adj_cube,3)]); %makes each correlation matrix into a single row
            templist=sort(templist,1,'descend');
            min_wght=templist(num_nodes+2*avg_edg*num_nodes,:);
            
            gg=bsxfun(@ge,adj_cube,reshape(min_wght,[1 1 5])); %creates logical matrix containing 1 wherever the correlation is greater than the min_wght
            adj_cube(~gg)=0; %keeps only values greater than min_wght
            
        else
            
            num_nodes=size(adj_cube,1);
            gg=adj_cube>=min_wght;
            adj_cube(~gg)=0;
            
        end
        
        adj_cube(adj_cube>0)=1; %sets all nonzero values to 1
        eeg(m).szraw(k).pdc_adj_cube=adj_cube;
        
        %eeg(m).sz(k).metric=zeros(num_nodes,5,4);
    end
    metric=arrayfun(@(x) graphmetric(x.pdc_adj_cube),eeg(m).szraw,'UniformOutput',0);
    conn_cube=arrayfun(@(x) conn(x.pdc_adj_cube),eeg(m).szraw,'UniformOutput',0);
    
    for k=1:length(eeg(m).szraw)
        eeg(m).szraw(k).metric_pdc=metric{k};
        eeg(m).szraw(k).conn_cube_pdc=conn_cube{k};
    end
    %[smlwrld,dataCL]=arrayfun(@(x) smi(x.pdc_adj_cube,100),eeg(m).szraw,'UniformOutput',0);
    %snczb=arrayfun(@(x) synch(x.pdc_adj_cube)',eeg(m).szraw,'UniformOutput',0);
    %glob_mat=arrayfun(@(x,y,z) horzcat(x{:},y{:},z{:}),smlwrld,dataCL,snczb,'UniformOutput',0);
    %for k=1:length(eeg(m).szraw)
    %    eeg(m).szraw(k).glob_mat_pdc=glob_mat{k};
    %end
    
    
end
%%
cutoff_method=1;
avg_edg=5;
totpt=num2str(length(eeg));

%lfreq=0:6:108;
%hfreq=20:6:128;

lfreq=[3,8,13,30];
hfreq=[7,12,29,50];

step=1;

lfreq=lfreq./step;
hfreq=hfreq./step;

for m=1:length(eeg);
    ptnum=num2str(m);
    disp(['Calculating metrics for ',ptnum,'/',totpt]);
    for k=1:length(eeg(m).szraw)
        eeg(m).szraw(k).dtf_adj_cube=zeros(length(eeg(m).leads),length(eeg(m).leads),length(lfreq),5);
        eeg(m).szraw(k).metric_dtf=zeros(length(eeg(m).leads),5,5,length(lfreq));
    end
    
    for z=1:length(lfreq)
        for k=1:length(eeg(m).szraw)
            adj_cube=squeeze(mean(eeg(m).szraw(k).dtf(:,:,lfreq(z):hfreq(z),:),3));
            if cutoff_method==1
                num_nodes=size(adj_cube,1);
                templist=reshape(adj_cube,[num_nodes^2,size(adj_cube,3)]); %makes each correlation matrix into a single row
                templist=sort(templist,1,'descend');
                %min_wght=templist(num_nodes+2*avg_edg*num_nodes,:);
                min_wght=templist(round(0.05*length(templist)),:);
                
                gg=bsxfun(@ge,adj_cube,reshape(min_wght,[1 1 5])); %creates logical matrix containing 1 wherever the correlation is greater than the min_wght
                adj_cube(~gg)=0; %keeps only values greater than min_wght
                
            else
                
                num_nodes=size(adj_cube,1);
                gg=adj_cube>=min_wght;
                adj_cube(~gg)=0;
                
            end
            
            adj_cube(adj_cube>0)=1; %sets all nonzero values to 1
            eeg(m).szraw(k).dtf_adj_cube(:,:,z,:)=adj_cube;
            
            %eeg(m).sz(k).metric=zeros(num_nodes,5,4);
        end
        metric=arrayfun(@(x) graphmetric(squeeze(x.dtf_adj_cube(:,:,z,:))),eeg(m).szraw,'UniformOutput',0);
        %conn_cube=arrayfun(@(x) conn(x.dtf_adj_cube),eeg(m).szraw,'UniformOutput',0);
        
        for k=1:length(eeg(m).szraw)
            eeg(m).szraw(k).metric_dtf(:,:,:,z)=metric{k};
            %eeg(m).szraw(k).conn_cube_dtf=conn_cube{k};
        end
    end
    %[smlwrld,dataCL]=arrayfun(@(x) smi(x.pdc_adj_cube,100),eeg(m).szraw,'UniformOutput',0);
    %snczb=arrayfun(@(x) synch(x.pdc_adj_cube)',eeg(m).szraw,'UniformOutput',0);
    %glob_mat=arrayfun(@(x,y,z) horzcat(x{:},y{:},z{:}),smlwrld,dataCL,snczb,'UniformOutput',0);
    %for k=1:length(eeg(m).szraw)
    %    eeg(m).szraw(k).glob_mat_pdc=glob_mat{k};
    %end
    
    
end
