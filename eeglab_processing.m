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
%%
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
%%
parfor m=1:length(eeg)
    ptnum=num2str(m);
    for k=1:length(eeg(m).sztype)
        sznum=num2str(k);
        disp(['Calculating metrics for Sz ',sznum,' Pt ',ptnum]);
        datacube=eeg(m).szraw(k).datacube;
        nodes=size(datacube,2);
        coharrdiag=zeros([size(datacube,2),nodes,33,5]);
        coharr=coharrdiag;
        subs=zeros(nodes^2,2);
        subs(:,1)=repmat(1:nodes,1,nodes);
        for j=1:5
            for i=1:size(datacube,2)
                coharrdiag(:,i,:,j)=mscohere(datacube(:,:,j),circshift(datacube(:,:,j),1-i,2),hamming(64),[],64,512)';
            end
            for i=1:33
                coharrtemp=spdiags(coharrdiag(:,:,i,j),[0:-1:-(nodes-1)],nodes,nodes);
                coharrtemp=tril(coharrtemp,-1)'+coharrtemp;
                coharr(:,:,i,j)=coharrtemp;
                %coharr(:,:,i,j)=spdiags(coharr(:,:,i,j),[0:-1:-(nodes-1)],nodes,nodes);
            end
        end
        eeg(m).szraw(k).coharr=coharr;
        
    end
end

parfor m=1:length(eeg)
    ptnum=num2str(m);
    disp(['Calculating metrics for Pt ',ptnum]);
    datacube=eeg(m).bslraw.datacube;
    nodes=size(datacube,2);
    coharrdiag=zeros([size(datacube,2),nodes,33,size(datacube,3)]);
    coharr=coharrdiag;
    subs=zeros(nodes^2,2);
    subs(:,1)=repmat(1:nodes,1,nodes);
    for j=1:size(datacube,3)
        for i=1:size(datacube,2)
            coharrdiag(:,i,:,j)=mscohere(datacube(:,:,j),circshift(datacube(:,:,j),1-i,2),hamming(64),[],64,512)';
        end
        for i=1:33
            coharrtemp=spdiags(coharrdiag(:,:,i,j),[0:-1:-(nodes-1)],nodes,nodes);
            coharrtemp=tril(coharrtemp,-1)'+coharrtemp;
            coharr(:,:,i,j)=coharrtemp;
            %coharr(:,:,i,j)=spdiags(coharr(:,:,i,j),[0:-1:-(nodes-1)],nodes,nodes);
        end
    end
    eeg(m).bslraw.coharr=coharr;
    
end
%%
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
figure
for i=1:6
    for j=1:10
        subplot(10,6,6*(j-1)+i)
        if i==1
            imagesc([eeg(1).res,eeg(1).bslraw.coharr(:,:,j,i)])
        else
            imagesc([eeg(1).res,eeg(1).szraw(1).coharr(:,:,j,i-1)])
        end
        set(gca','XTick',[],'YTick',[])
    end
end
%%
%tic
for m=4:4
    for k=1:length(eeg(m).sztype)
        clear C
        clear part
        clear wmd
        for i=1:6
            for j=1:17
                if i==1
                    AIJ=eeg(m).bslraw.coharr(:,:,j,i);
                else
                    AIJ=eeg(m).szraw(1).coharr(:,:,j,i-1);
                end
                AIJ(logical(eye(size(AIJ,1))))=0;
                [C(:,j,i),Q(m,i)] = community_louvain(AIJ);
                [X,Y,INDSORT] = grid_communities(C(:,j,i)); % call function
                %part(:,j,i)=participation_coef(eeg(m).szraw(1).coharr(:,:,j,i),C(:,j,i));
                %wmd(:,j,i)=module_degree_zscore(eeg(m).szraw(1).coharr(:,:,j,i),C(:,j,i));
                %subplot (10,5,5*(j-1)+i)
                %imagesc(AIJ(INDSORT,INDSORT));           % plot ordered adjacency matrix
                %hold on
                %plot(X,Y,'r','linewidth',2);
                %imagesc([C/max(C),C/max(C),C/max(C),AIJ])
            end
            
        end
        %eeg(m).szraw(k).part=part;
        %eeg(m).szraw(k).wmd=wmd;
        %toc
    end
    
    if mod(m-1,9)==0
        figure
    end
    subplot(3,3,mod(m-1,9)+1)
    colors=get(gca,'colororder');
    Cmat=squeeze(C(:,2,:));
    nres=nnz(logical(eeg(m).res));
    nnodes=length(eeg(m).res);
    %plot(Cmat'+0.2*rand([6 nnodes]),'b')
    for i=unique(Cmat(:,1))'
        %plot(Cmat((Cmat(:,1)==i),:)'+0.2*rand([6 nnz((Cmat(:,1)==i))])+i*max(Cmat(:,1))*ones(6,nnz((Cmat(:,1)==i))),'color',colors(i,:))
        plot(Cmat((Cmat(:,1)==i),:)'+0.2*rand([6 nnz((Cmat(:,1)==i))]),'color',colors(i,:))
        hold on
    end
    %plot(Cmat(~logical(eeg(m).res),:)'+0.2*rand([6 nnodes-nres]),'b')
    %hold on
    %if strcmp(eeg(m).engel,'1a')
    %    plot(Cmat(logical(eeg(m).res),:)'+0.2*rand([6 nres]),'b')
    %else
    %    plot(Cmat(logical(eeg(m).res),:)'+0.2*rand([6 nres]),'r')
    %end
end
%histogram(wmd(:,1,2),-3:3)
%hold on
%histogram(wmd(logical(eeg(m).res),1,2),-3:3)

%%
i=2;

if strcmp(eeg(m).engel,'1a')
    subplot 121
    title('unresected')
    scatter(part(~logical(eeg(m).res),1,i),wmd(~logical(eeg(m).res),1,i),20,'b','filled')
    hold on
    subplot 122
    title('resected')
    scatter(part(logical(eeg(m).res),1,i),wmd(logical(eeg(m).res),1,i),20,'b','filled')
    hold on
else
    subplot 121
    scatter(part(~logical(eeg(m).res),1,i),wmd(~logical(eeg(m).res),1,i),20,'r','filled')
    hold on
    subplot 122
    scatter(part(logical(eeg(m).res),1,i),wmd(logical(eeg(m).res),1,i),20,'r','filled')
    hold on
end


ylim([-3 3])
%title(eeg(m).engel)

end
end
%%
for j=1:1
    %subplot(1,10,j)
    Cmat=squeeze(C(:,j,:));
    Cmatnew=zeros(size(Cmat));
    Cmatnew(:,1)=Cmat(:,1);
    nclust=max(max(Cmat));
    for i=2:2
        ov=zeros(nclust,nclust);
        for t1=1:nclust;
            clust1=Cmat(:,i-1)==t1;
            for t2=1:nclust
                clust2=Cmat(:,i)==t2;
                ov(t1,t2)=nnz(clust1&clust2);
            end
            %[~,clust2new]=max(ov);
            %Cmatnew(Cmat(:,i)==clust2new,i)=t1;
        end
        %[ovsort,ind]=sort(ov,2,'descend');
        %for t1=1:nclust
        %    []
        %end
    end
    %imagesc(squeeze(C(:,j,:)));
end
%%
clust2new=zeros(1,nclust)

[ovsort,ind]=sort(ov,2,'descend')

[~,indi]=max(ovsort(:,1))
clust2new(indi)=ind(indi,1)
ovsort(indi,:)=zeros(1,5)
loser=ind(:,1)==ind(indi,1);
ind(loser,1:nclust-1)=ind(loser,2:nclust)
ind(loser,5)=0;
ind(indi,:)=zeros(1,5)

[~,indi]=max(ovsort(:,1))
clust2new(indi)=ind(indi,1)
ovsort(indi,:)=zeros(1,5)
loser=ind(:,1)==ind(indi,1);
ind(loser,1:nclust-1)=ind(loser,2:nclust)
ind(indi,:)=zeros(1,5)
%%

clust2new=zeros(1,nclust)

[ovsort,ind]=sort(ov,2,'descend')
indorig=ind
ovsortorig=ovsort
while nnz(clust2new)<nclust;
    [~,indi]=max(ovsort(:,1))
    indiorig=find(ovsortorig(indi,:)==ovsort(indi,1))
    loser=ind(:,1)==ind(indi,1)
    if nnz(clust2new==ind(indi,1))>0
        ind(loser,1:nclust-1)=ind(loser,2:nclust)
        ind(loser,nclust)=0
        ovsort(loser,1:nclust-1)=ovsort(loser,2:nclust)
        ovsort(loser,nclust)=0
    else
        clust2new(indiorig)=ind(indi,1)
        ovsort(indi,:)=zeros(1,nclust)
        ind(loser,1:nclust-1)=ind(loser,2:nclust)
        ind(loser,nclust)=0
        ind(indi,:)=zeros(1,nclust)
    end
end
for i=2:5
    for j=1:nclust
        Cmatnew(Cmat(:,i)==clust2new(j),i)=j;
    end
end
%%
figure
for j=1:17
    for i=1:5
        subplot(3,6,j)
        %[N,edges]=histcounts(log10(coharr(:,:,j,i)));
        %plot(edges(1:end-1),N-60);
        hold on
        connarr{j,i}=conn(coharr(:,:,j,i)>0.1);
    end
    %xlim([-5 0])
end
%%
cutoff_method=1;
avg_edg=5;
totpt=num2str(length(eeg));
for m=1:length(eeg);
    ptnum=num2str(m);
    disp(['Calculating metrics for ',ptnum,'/',totpt]);
    
    for k=1:length(eeg(m).szraw)
        adj_cube=squeeze(eeg(m).szraw(k).coharr(:,:,10,:));
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
        eeg(m).szraw(k).coh_adj_cube=adj_cube;
        
        %eeg(m).sz(k).metric=zeros(num_nodes,5,4);
    end
    metric=arrayfun(@(x) graphmetric(x.coh_adj_cube),eeg(m).szraw,'UniformOutput',0);
    conn_cube=arrayfun(@(x) conn(x.coh_adj_cube),eeg(m).szraw,'UniformOutput',0);
    
    for k=1:length(eeg(m).szraw)
        eeg(m).szraw(k).metric_coh=metric{k};
        eeg(m).szraw(k).conn_cube_coh=conn_cube{k};
    end
    [smlwrld,dataCL]=arrayfun(@(x) smi(x.conn_cube_coh,100),eeg(m).szraw,'UniformOutput',0);
    snczb=arrayfun(@(x) synch(x.conn_cube_coh)',eeg(m).szraw,'UniformOutput',0);
    glob_mat=arrayfun(@(x,y,z) horzcat(x{:},y{:},z{:}),smlwrld,dataCL,snczb,'UniformOutput',0);
    for k=1:length(eeg(m).szraw)
        eeg(m).szraw(k).glob_mat_coh=glob_mat{k};
    end
    
    
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
%%
pxx=zeros(17,length(eeg(14).leads),5);
for i=1:5
    pxx(:,:,i)=pwelch(eeg(14).szraw(1).datacube(:,:,i),hamming(32),[],32,512);
end
%%
for i=1:54
    figure
    imagesc(squeeze(log10(pxx(:,i,:))));
end