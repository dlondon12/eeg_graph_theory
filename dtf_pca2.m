p_account=zeros(45,4,4);
p=zeros(45,4);
for z=2:2
    nodemat=[];
    for m=1:length(eeg2)
        
        for k=1:length(eeg2(m).szraw)
            
            nodemat=cat(3,nodemat,permute(squeeze(eeg2(m).szraw(k).metric_dtf(:,z,:,:)),[3,2,1]));
            
        end
        
    end
    
    engel1a=arrayfun(@(x) repmat((strcmp(x.engel,'1a')),length(x.sz)*length(x.leads),1),eeg2,'UniformOutput',0)';
    engel1b=arrayfun(@(x) repmat((strcmp(x.engel,'1b')),length(x.sz)*length(x.leads),1),eeg2,'UniformOutput',0)';
    engel1c=arrayfun(@(x) repmat((strcmp(x.engel,'1c')),length(x.sz)*length(x.leads),1),eeg2,'UniformOutput',0)';
    engel2a=arrayfun(@(x) repmat((strcmp(x.engel,'2a')),length(x.sz)*length(x.leads),1),eeg2,'UniformOutput',0)';
    engel2b=arrayfun(@(x) repmat((strcmp(x.engel,'2b')),length(x.sz)*length(x.leads),1),eeg2,'UniformOutput',0)';
    engel3a=arrayfun(@(x) repmat((strcmp(x.engel,'3a')),length(x.sz)*length(x.leads),1),eeg2,'UniformOutput',0)';
    engel4b=arrayfun(@(x) repmat((strcmp(x.engel,'4b')),length(x.sz)*length(x.leads),1),eeg2,'UniformOutput',0)';
    engel1a=vertcat(engel1a{:});
    engel1b=vertcat(engel1b{:});
    engel1c=vertcat(engel1c{:});
    engel2a=vertcat(engel2a{:});
    engel2b=vertcat(engel2b{:});
    engel3a=vertcat(engel3a{:});
    engel4b=vertcat(engel4b{:});
    
    res=arrayfun(@(x) repmat(x.res,length(x.sz),1),eeg2,'UniformOutput',0)';
    res=logical(vertcat(res{:}));
    %subclass=isnan(szonset);
    %subclass=szonset==1|szonset==2|szonset==3;
    %subclass=isnan(szonset)|szonset==1|szonset==2|szonset==3;
    subclass=true(size(engel1a));
    res=res(subclass);
    
    nodematpca=squeeze(reshape(nodemat,[],1,size(nodemat,3)))';
    nodematpca=nodematpca(subclass,:);
    engel1a=engel1a(subclass);
    engel1b=engel1b(subclass);
    engel1c=engel1c(subclass);
    engel2a=engel2a(subclass);
    engel2b=engel2b(subclass);
    engel3a=engel3a(subclass);
    engel4b=engel4b(subclass);
    
    
    nodematpca_1a=nodematpca(engel1a,:);
    nodematpca_1b=nodematpca(engel1b|engel1c,:);
    nodematpca_2a=nodematpca(engel2a|engel2b,:);
    nodematpca_3=nodematpca(engel3a|engel4b,:);
    nodematpca_n=nodematpca(~engel1a,:);
    
    
    logmat=false(length(engel1a),5);
    logmat(:,1)=engel1a;
    logmat(:,2)=engel1b|engel1c|engel2a;
    logmat(:,3)=engel1b|engel1c|engel2a;
    logmat(:,4)=engel2b|engel3a|engel4b;
    logmat(:,5)=~engel1a;
    
    nodematpca_cell=cell(5,1);
    means=zeros(10000,size(nodematpca,2),5,4);
    diffact=zeros(size(nodematpca,2),4,4);
    diffmeans=zeros(10000,size(nodematpca,2),4,4);
    for j=1:5
        nodematpca_cell{j}=nodematpca(logmat(:,j),:);
        for i=1:10000
            means(i,:,j,z)=mean(nodematpca_cell{j}(randperm(nnz(logmat(:,j)),nnz(logmat(:,j)&res)),:));
        end
    end
    
    for j=1:3
        diffmeans(:,:,j,z)=means(:,:,j+1,z)-means(:,:,j,z);
        diffact(:,j,z)=mean(nodematpca(logmat(:,j+1)&res,:))-mean(nodematpca(logmat(:,j)&res,:));
    end
    diffmeans(:,:,4,z)=means(:,:,5,z)-means(:,:,1,z);
    diffact(:,4,z)=mean(nodematpca(logmat(:,5)&res,:))-mean(nodematpca(logmat(:,1)&res,:));
    
    
    %vals=-mean(nodematpca(logmat(:,1)&res,:))+mean(nodematpca(logmat(:,2)&res,:));
    
    %p_great=(sum(bsxfun(@ge,permute(diffact(:,:,z),[3,1,2]),diffmeans(:,:,:,z))))/10000;
    %p_less=(sum(bsxfun(@le,permute(diffact(:,:,z),[3,1,2]),diffmeans(:,:,:,z))))/10000;
    %p_account(:,:,z)=squeeze(min(cat(1,p_great,p_less))).*2;
    %%end
    
    %p=zeros(size(nodematpca,2),1);
    %tic
    %parfor i=1:size(nodematpca,2)
    %    p(i,z)=poisstestu(nodematpca(engel1a,i),nodematpca(~engel1a,i),10000);
    %    i
    %end
    %toc
end
%%
engelbsl=arrayfun(@(x) repmat(strcmp(x.engel,'1a'),length(x.leads)*size(x.bslraw.metric_dtf,3),1),eeg2,'UniformOutput',0)';
engelbsl=vertcat(engelbsl{:});
resbsl=arrayfun(@(x) repmat(x.res,size(x.bslraw.metric_dtf,3),1),eeg2,'UniformOutput',0)';
resbsl=logical(vertcat(resbsl{:}));
nodematbsl=[];
for m=1:length(eeg2)
    for k=1:size(eeg2(m).bslraw.metric_dtf,3)
        nodematbsl=cat(1,nodematbsl,squeeze(eeg2(m).bslraw.metric_dtf(:,z,k,:)));
    end
end

logmatbsl=false(length(engelbsl),2);
logmatbsl(:,1)=engelbsl;
logmatbsl(:,2)=~engelbsl;
nodematbsl_cell=cell(2,1);
meansbsl=zeros(10000,size(nodematbsl,2),2);
%diffactbsl=zeros(size(nodematbsl,2));
%diffmeansbsl=zeros(10000,size(nodematbsl,2));
for j=1:2
    nodematbsl_cell{j}=nodematbsl(logmatbsl(:,j),:);
    for i=1:10000
        meansbsl(i,:,j)=mean(nodematbsl_cell{j}(randperm(nnz(logmatbsl(:,j)),nnz(logmatbsl(:,j)&resbsl)),:));
    end
end

diffmeansbsl=meansbsl(:,:,2)-meansbsl(:,:,1);
diffactbsl=mean(nodematbsl(logmatbsl(:,2)&resbsl,:))-mean(nodematbsl(logmatbsl(:,1)&resbsl,:));

%valsbsl=-mean(nodematbsl(logmatbsl(:,1)&resbsl,:))+mean(nodematbsl(logmatbsl(:,2)&resbsl,:));
%p_greatbsl=(sum(bsxfun(@ge,diffactbsl,diffmeansbsl)))/10000;
%p_lessbsl=(sum(bsxfun(@le,diffactbsl,diffmeansbsl)))/10000;
%p_accountbsl=squeeze(min(cat(1,p_greatbsl,p_lessbsl))).*2;

%p_bsl=zeros(size(nodematbsl,2),1);
%tic
%parfor i=1:size(nodematbsl,2)
%    p_bsl(i)=poisstestu(nodematbsl(engelbsl,i),nodematbsl(~engelbsl,i),10000);
%    i
%end
%toc

%%
[sig,clusts]=clust_mass(nodemat(:,:,engel1a),nodemat(:,:,~engel1a),3,0.01,0.01,1000);
[sig_bsl,clusts_bsl]=clust_mass1d(nodematbsl(engelbsl,:),nodematbsl(~engelbsl,:),1,0.01,0.01,1000);
%map=colormap('jet');
%%
map=zeros(64,3);
map(:,1)=linspace(0,1,64);
map(:,2)=linspace(0,1,64);
map(:,3)=linspace(0.5,0,64);
map=colormap(map);
threshfreq=exp(0.91:0.5:4.91);
h1=subplot(24,2,5:2:15);
imagesc([mean(nodematbsl(engelbsl,:));reshape(mean(nodematpca(logmat(:,1),:)),9,5)'],[0.025,0.05])
set(gca,'XTick',1:9,'XTickLabel',num2str(threshfreq','%.0f'),'YTick',1:6,'YTickLabel',{'Int';'Pre';'Early';'Middle';'Late';'Post'},'FontSize',8,'box','off')

t1=title('Engel 1a');
%set(t1,'Units','normalized','FontSize',14,'FontWeight','bold');
%t1.Position(2)=1.28;
text(-0.24,1.2,'a','Units','Normalized','FontSize',16,'FontWeight','bold')
%text(1.16,1.02,{'Betweenness';'Centrality'},'Units','normalized',...
%    'horizontalAlignment','center','FontSize',7)
h1.Position(1:2)=[0.09,0.72];


h2=subplot(24,2,6:2:16);
imagesc([mean(nodematbsl(~engelbsl,:));reshape(mean(nodematpca(logmat(:,5),:)),9,5)'],[0.025,0.05])
set(gca,'XTick',1:9,'XTickLabel',num2str(threshfreq','%.0f'),'YTick',1:6,'YTickLabel',{'Int';'Pre';'Early';'Middle';'Late';'Post'},'FontSize',8,'box','off')
t1b=title('>Engel 1a');
%set(t1b,'Units','normalized','FontSize',14,'FontWeight','bold');
%t1b.Position(2)=1.28;
text(-0.15,1.2,'b','Units','Normalized','FontSize',16,'FontWeight','bold')
h2.Position(2)=0.72;
col=colorbar('eastoutside');
col.Position([1,3])=[0.92,0.01];
%set(col,'FontSize',10)
%subplot(3,4,9)
%imagesc(reshape(mean(nodematpca(logmat(:,4),:)),9,5)',[0.027,0.063])


h3=subplot(24,2,21:2:31);
vals=-[mean(nodematbsl(engelbsl,:))-mean(nodematbsl(~engelbsl,:));reshape(mean(nodematpca(engel1a,:))-mean(nodematpca(~engel1a,:)),9,5)'];

ind=round(interp1(linspace(-0.001,0.01,size(map,1)),1:size(map,1),vals));
cols(:,:,1)=reshape(map(ind,1),6,9);
cols(:,:,2)=reshape(map(ind,2),6,9);
cols(:,:,3)=reshape(map(ind,3),6,9);

%cols2=cols;
%cols(repmat([p_bsl'>=0.001;reshape(p(:,2)>=0.001,9,5)'],[1,1,3]))=0.7-(0.7-cols(repmat([p_bsl'>=0.001;reshape(p(:,2)>=0.001,9,5)'],[1,1,3])))*0.3;
cols(repmat([~sig_bsl;~sig'],[1,1,3]))=0.7-(0.7-cols(repmat([~sig_bsl;~sig'],[1,1,3])))*0.3;

imagesc(cols)
set(gca,'XTick',1:9,'XTickLabel',num2str(threshfreq','%.0f'),'YTick',1:6,'YTickLabel',{'Int';'Pre';'Early';'Middle';'Late';'Post'},'box','off','FontSize',8)
h3.Position(1:2)=[0.09,0.45];
col2=colorbar?('eastoutside');
col2.Position([1,3])=[0.44,0.01];
set(col2,'Ticks',[0:1/4:1],'TickLabels',num2str((-0.001:0.003:0.011)'))
t2=title('Difference');
%set(t2,'Units','normalized','FontSize',14,'FontWeight','bold');
%t2.Position(2)=1.4;
text(-0.24,1.2,'c','Units','Normalized','FontSize',16,'FontWeight','bold')
xlabel('Frequency (Hz)')

h4=subplot(24,2,22:2:32);

%colors=get(gca,'ColorOrder');
colors=[0,118,192;163,2,52]/255;
colorshade=[186,207,236;228,184,180]/255;
patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,1),37:45))+...
    std(nodematpca(logmat(:,1),37:45))/sqrt(nnz(logmat(:,1))),...
    fliplr(mean(nodematpca(logmat(:,1),37:45))-...
    std(nodematpca(logmat(:,1),37:45))/sqrt(nnz(logmat(:,1))))],...
    colorshade(1,:),'EdgeColor','none')

patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,5),37:45))+...
    std(nodematpca(logmat(:,5),37:45))/sqrt(nnz(logmat(:,5))),...
    fliplr(mean(nodematpca(logmat(:,5),37:45))-...
    std(nodematpca(logmat(:,5),37:45))/sqrt(nnz(logmat(:,5))))],...
    colorshade(2,:),'EdgeColor','none')

%patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,4),37:45))+...
%    std(nodematpca(logmat(:,4),37:45))/sqrt(nnz(logmat(:,4))),...
%    fliplr(mean(nodematpca(logmat(:,4),37:45))-...
%    std(nodematpca(logmat(:,4),37:45))/sqrt(nnz(logmat(:,4))))],...
%    mean([colors(3,:);1,1,1]),'EdgeColor','none')

hold on

e1=plot(1:9,mean(nodematpca(logmat(:,1),37:45)),'LineWidth',2,'color',colors(1,:));
ne1=plot(1:9,mean(nodematpca(logmat(:,5),37:45)),'LineWidth',2,'color',colors(2,:));
%plot(1:9,mean(nodematpca(logmat(:,4),37:45)),'LineWidth',2);

set(gca,'XTick',1:9,'XTickLabel',num2str(threshfreq','%.0f'),'YLim',[0.025,0.052],'YTick',0.025:0.01:0.045,'box','off','FontSize',8,'YAxisLocation','right')
text(-0.15,3.6,'d','Units','Normalized','FontSize',16,'FontWeight','bold')
text(-0.05,0.5,'Post','Units','Normalized','FontSize',12,'Rotation',90,'horizontalAlignment','center')
xlabel('Frequency (Hz)')
h4.Position(2)=0.45;

pos=h4.Position;
h4.Position=[pos(1:3),pos(4)/3];

h4_2=axes('position',[pos(1),pos(2)+pos(4)/3,pos(3),pos(4)/3]);

patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,1),28:36))+...
    std(nodematpca(logmat(:,1),28:36))/sqrt(nnz(logmat(:,1))),...
    fliplr(mean(nodematpca(logmat(:,1),28:36))-...
    std(nodematpca(logmat(:,1),28:36))/sqrt(nnz(logmat(:,1))))],...
    colorshade(1,:),'EdgeColor','none')

patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,5),28:36))+...
    std(nodematpca(logmat(:,5),28:36))/sqrt(nnz(logmat(:,5))),...
    fliplr(mean(nodematpca(logmat(:,5),28:36))-...
    std(nodematpca(logmat(:,5),28:36))/sqrt(nnz(logmat(:,5))))],...
    colorshade(2,:),'EdgeColor','none')

hold on

e1=plot(1:9,mean(nodematpca(logmat(:,1),28:36)),'LineWidth',2,'color',colors(1,:));
ne1=plot(1:9,mean(nodematpca(logmat(:,5),28:36)),'LineWidth',2,'color',colors(2,:));
%plot(1:9,mean(nodematpca(logmat(:,4),37:45)),'LineWidth',2);

set(gca,'XTick',1:9,'XTickLabel',[],'YLim',[0.025,0.044],'YTick',0.03:0.01:0.04,'box','off','FontSize',8,'YAxisLocation','right')
[leg,hobj]=legend([e1,ne1],'Engel 1a','>Engel 1a');
hobj(3).XData=[0.3,0.45];
hobj(5).XData=[0.3,0.45];
set(leg,'box','off','Position',[0.55,0.515,0.1,0.05])
text(-0.05,0.5,'Late','Units','Normalized','FontSize',12,'Rotation',90,'horizontalAlignment','center')
ylab=ylabel('BC');
%ylab.Position(2)=h4.YLim(2);

h4_3=axes('position',[pos(1),pos(2)+2*pos(4)/3,pos(3),pos(4)/3]);

patch([1:9,9:-1:1],[mean(nodematbsl(engelbsl,:))+...
    std(nodematbsl(engelbsl,:))/sqrt(nnz(engelbsl)),...
    fliplr(mean(nodematbsl(engelbsl,:))-...
    std(nodematbsl(engelbsl,:))/sqrt(nnz(engelbsl)))],...
    colorshade(1,:),'EdgeColor','none')

hold on

patch([1:9,9:-1:1],[mean(nodematbsl(~engelbsl,:))+...
    std(nodematbsl(~engelbsl,:))/sqrt(nnz(~engelbsl)),...
    fliplr(mean(nodematbsl(~engelbsl,:))-...
    std(nodematbsl(~engelbsl,:))/sqrt(nnz(~engelbsl)))],...
    colorshade(2,:),'EdgeColor','none')

e1_bsl=plot(1:9,mean(nodematbsl(engelbsl,:)),'LineWidth',2,'Color',colors(1,:));
ne1_bsl=plot(1:9,mean(nodematbsl(~engelbsl,:)),'LineWidth',2,'Color',colors(2,:));

set(gca,'XTick',1:9,'XTickLabel',{},'YLim',[0.02,0.037],'YTick',0.025:0.01:0.035,'box','off','FontSize',8,'YAxisLocation','right')
text(-0.05,0.5,'Int','Units','Normalized','FontSize',12,'Rotation',90,'horizontalAlignment','center')

h(1)=subplot(24,2,37:2:47);
nozer=mean(nodematpca(:,37:45),2);
nozer(nozer==0)=1e-5;
[count,edges]=histcounts(log10(nozer(engel1a)),'normalization','probability','BinWidth',0.125);
e2=plot(edges(1:end-1),count,'LineWidth',2,'color',colors(1,:));
hold on
[count,edges]=histcounts(log10(nozer(~engel1a)),'normalization','probability','BinWidth',0.125);
ne2=plot(edges(1:end-1),count,'LineWidth',2,'color',colors(2,:));
xlim([-5.5,0])

xlabel(texlabel('log10(BC)'),'interpreter','tex')
set(gca,'box','off','XTick',-5:1:0,'XTicklabel',{'-Inf';num2str((-4:1:0)')},'FontSize',8,'YTick',0:0.05:0.15,'YTickLabel',{'0';'';'0.1';''})
xlabel(texlabel('log10(BC)'),'interpreter','tex')
text(-0.24,2.4,'e','Units','Normalized','FontSize',16,'FontWeight','bold')
text(1.05,0.5,'Post','Units','Normalized','FontSize',12,'Rotation',90,'horizontalAlignment','center')
%title('Mean of All Frequencies, Post Seizure')
h(1).Position(1)=0.09;

ptengel1a=arrayfun(@(x) (strcmp(x.engel,'1a')),eeg2);
inds=1:length(eeg2);
ptind=arrayfun(@(x,y) repmat(y,1,length(x.sztype)*length(x.leads)),eeg2,inds,'UniformOutput',0);
ptind=horzcat(ptind{:})';
ptindbsl=arrayfun(@(x,y) repmat(y,1,size(x.bslraw.metric_dtf,3)*length(x.leads)),eeg2,inds,'UniformOutput',0);
ptindbsl=horzcat(ptindbsl{:})';

h(2)=subplot(24,2,38:2:44);
totalclass_post=logical(sum(bsxfun(@ge,nodematpca(:,37:45),mean(nodematpca(:,37:45))+3*std(nodematpca(:,37:45))),2));
totalclass_post=accumarray(ptind,totalclass_post,[],@mean);
cfier=0:0.01:0.6;
%[sorted,ind]=sort(totalclass);
fpr_post=zeros(length(cfier),1);
tpr_post=fpr_post;
for i=1:length(cfier)
tpr_post(i)=nnz(~ptengel1a(find(totalclass_post>=cfier(i))))./nnz(~ptengel1a);
fpr_post(i)=nnz(ptengel1a(find(totalclass_post>=cfier(i))))./nnz(ptengel1a);
end
post_roc=plot(fpr_post,tpr_post,'color',[0,0,0],'LineWidth',2,'LineStyle','-');
hold on
%[~,opt]=min((1-tpr).^2+fpr.^2); %distance from corner
[~,opt_post]=max(tpr_post-fpr_post); %J statistic

totalclass_late=logical(sum(bsxfun(@ge,nodematpca(:,33:36),mean(nodematpca(:,33:36))+3*std(nodematpca(:,33:36))),2));
totalclass_late=accumarray(ptind,totalclass_late,[],@mean);
fpr_late=zeros(length(cfier),1);
tpr_late=fpr_late;
for i=1:length(cfier)
tpr_late(i)=nnz(~ptengel1a(find(totalclass_late>=cfier(i))))./nnz(~ptengel1a);
fpr_late(i)=nnz(ptengel1a(find(totalclass_late>=cfier(i))))./nnz(ptengel1a);
end
late_roc=plot(fpr_late,tpr_late,'color',[0.35 0.35 0.35],'LineWidth',2,'LineStyle','-');
hold on
%[~,opt]=min((1-tpr).^2+fpr.^2); %distance from corner
[~,opt_late]=max(tpr_late-fpr_late); %J statistic

totalclass_bsl=logical(sum(bsxfun(@ge,nodematbsl(:,4:9),mean(nodematbsl(:,4:9))+3*std(nodematbsl(:,4:9))),2));
totalclass_bsl=accumarray(ptindbsl,totalclass_bsl,[],@mean);
%cfier=0:0.01:0.6;
%[sorted,ind]=sort(totalclass);
fpr_bsl=zeros(length(cfier),1);
tpr_bsl=fpr_bsl;
for i=1:length(cfier)
tpr_bsl(i)=nnz(~ptengel1a(find(totalclass_bsl>=cfier(i))))./nnz(~ptengel1a);
fpr_bsl(i)=nnz(ptengel1a(find(totalclass_bsl>=cfier(i))))./nnz(ptengel1a);
end
bsl_roc=plot(fpr_bsl,tpr_bsl,'color',[0.7 0.7 0.7],'LineWidth',2,'LineStyle','-');
hold on
%[~,opt]=min((1-tpr).^2+fpr.^2); %distance from corner
[~,opt_bsl]=max(tpr_bsl-fpr_bsl); %J statistic

plot(fpr_post(opt_post),tpr_post(opt_post),'Marker','d','MarkerFaceColor',...
    [0 0 0],'MarkerSize',10,'MarkerEdgeColor',[0 0 0]);

plot(fpr_late(opt_late),tpr_late(opt_late),'Marker','d','MarkerFaceColor',...
    [0.35 0.35 0.35],'MarkerSize',10,'MarkerEdgeColor',[0.35 0.35 0.35]);

plot(fpr_bsl(opt_bsl),tpr_bsl(opt_bsl),'Marker','d','MarkerFaceColor',...
    [0.7 0.7 0.7],'MarkerSize',10,'MarkerEdgeColor',[0.7 0.7 0.7]);

%q_post=annotation('arrow');
%set(q_post,'parent',h(2),'position',[0,1,0.95*fpr_post(opt_post),0.95*(tpr_post(opt_post)-1)],...
%    'color',[0,0,0],'LineStyle','-','HeadLength',5,'HeadWidth',5)
line([0,1],[0,1],'color',[0,0,0])

%q_late=annotation('arrow');
%set(q_late,'parent',h(2),'position',[0,1,0.95*fpr_late(opt_late),0.95*(tpr_late(opt_late)-1)],...
%    'color',[0.25 0.25 0.25],'LineStyle','-','HeadLength',5,'HeadWidth',5)

%q_bsl=annotation('arrow');
%set(q_bsl,'parent',h(2),'position',[0,1,0.95*fpr_bsl(opt_bsl),0.95*(tpr_bsl(opt_bsl)-1)],...
%    'color',[0.5 0.5 0.5],'LineStyle','-','HeadLength',5,'HeadWidth',5)

%qleg=annotation('arrow');
%set(qleg,'parent',h(2),'position',[0.4,0.2,0.17,0],...
%    'color','k','LineStyle','-','HeadLength',5,'HeadWidth',5)

[leg_roc,rocobj]=legend([bsl_roc,late_roc,post_roc],'Int','Late','Post');
set(leg_roc,'box','off','Position',[0.65,0.15,0.1,0.05]);

rocobj(4).XData(1)=0.4;
rocobj(6).XData(1)=0.4;
rocobj(8).XData(1)=0.4;
%text(0.67,0.18,['Optimal',char(10),'Discriminator'],'FontSize',8)

xlabel('FPR')
ylabel('TPR')
xlim([0 1])
ylim([0 1])
set(gca,'box','off','FontSize',8);
text(-0.25,1.2,'f','Units','Normalized','FontSize',16,'FontWeight','bold')

h(3)=subplot(24,2,46:2:48);
plotSpread(totalclass_bsl,'distributionIdx',~ptengel1a,'distributionColors',colors(1:2,:),'spreadWidth',1.5);
hold on
optline=line([0,3],[cfier(opt_bsl),cfier(opt_bsl)],'color','k','LineStyle','-');
plotSpread(totalclass_late,'distributionIdx',~ptengel1a,'distributionColors',colors(1:2,:),'spreadWidth',1.5,'xValues',[4,5]);
hold on
optline=line([3,6],[cfier(opt_late),cfier(opt_late)],'color','k','LineStyle','-');
plotSpread(totalclass_post,'distributionIdx',~ptengel1a,'distributionColors',colors(1:2,:),'spreadWidth',1.5,'xValues',[7,8]);
hold on
optline=line([6,9],[cfier(opt_post),cfier(opt_post)],'color','k','LineStyle','-');
ylabel('Fraction')
set(h(3),'YAxisLocation','right','XTick',[1,2,4,5,7,8],'XTickLabel',{'\color{white}>\color{black}1a','>1a','\color{white}>\color{black}1a','>1a'},'XTickLabelRotation',90,'FontSize',8)
xlim([0 9])
text(-0.35,1.2,'g','Units','Normalized','FontSize',16,'FontWeight','bold')
%legopt=legend(optline,['Optimal',char(10),'Discriminator']);
%set(legopt,'box','off')
%legopt.Position(1:2)=[0.63,0.1];
text(1/6,1.1,'Int','Units','Normalized','FontSize',10,'horizontalAlignment','center')
text(1/2,1.1,'Late','Units','Normalized','FontSize',10,'horizontalAlignment','center')
text(5/6,1.1,'Post','Units','Normalized','FontSize',10,'horizontalAlignment','center')

h(1).Position(2)=0.09;
h(2).Position([2,4])=h(1).Position([2,4]);
h(2).Position(3)=h(2).Position(4)-0.03;

h(3).Position(1)=h(2).Position(1)+h(2).Position(3)+0.05;
h(3).Position(2)=0.09;
h(3).Position(3)=h(1).Position(3)-h(2).Position(3)-0.05;
h(3).Position(4)=h(2).Position(4);

pos=h(1).Position;
h(1).Position=[pos(1:3),pos(4)/3];
ylab.Position(2)=h(1).YLim(2);

h_2(1)=axes('position',[pos(1),pos(2)+pos(4)/3,pos(3),pos(4)/3]);

nozer=mean(nodematpca(:,28:36),2);
nozer(nozer==0)=1e-5;
[count,edges]=histcounts(log10(nozer(engel1a)),'normalization','probability','BinWidth',0.125);
e2=plot(edges(1:end-1),count,'LineWidth',2,'color',colors(1,:));
hold on
[count,edges]=histcounts(log10(nozer(~engel1a)),'normalization','probability','BinWidth',0.125);
ne2=plot(edges(1:end-1),count,'LineWidth',2,'color',colors(2,:));
xlim([-5.5,0])
ylabel('Probability')

text(1.05,0.5,'Late','Units','Normalized','FontSize',12,'Rotation',90,'horizontalAlignment','center')
set(gca,'box','off','XTick',-5:1:0,'XTickLabel',{},'FontSize',8,'YTick',0:0.05:0.15,'YTickLabel',{'0';'';'0.1';''})

h_3(1)=axes('position',[pos(1),pos(2)+2*pos(4)/3,pos(3),pos(4)/3]);

nozer=mean(nodematbsl(:,4:9),2);
nozer(nozer==0)=1e-5;
[count,edges]=histcounts(log10(nozer(engelbsl)),'normalization','probability','BinWidth',0.125);
e3=plot(edges(1:end-1),count,'LineWidth',2,'color',colors(1,:));
hold on
[count,edges]=histcounts(log10(nozer(~engelbsl)),'normalization','probability','BinWidth',0.125);
ne3=plot(edges(1:end-1),count,'LineWidth',2,'color',colors(2,:));
xlim([-5.5,0])
ylab.Position(2)=h(1).YLim(2);
text(1.05,0.5,'Int','Units','Normalized','FontSize',12,'Rotation',90,'horizontalAlignment','center')
set(gca,'box','off','XTick',-5:1:0,'XTickLabel',{},'FontSize',8,'YTick',0:0.05:0.15,'YTickLabel',{'0';'';'0.1';''})

[leg2,hobj2]=legend([e3,ne3],'Engel 1a','>Engel 1a');
hobj2(3).XData=[0.3,0.45];
hobj2(5).XData=[0.3,0.45];
set(leg2,'box','off','Position',[0.12 0.23 0.13 0.049])

h(3).Children(6).MarkerSize=5;
h(3).Children(7).MarkerSize=5;
h(3).Children(9).MarkerSize=5;
h(3).Children(10).MarkerSize=5;
h(3).Children(12).MarkerSize=5;
h(3).Children(13).MarkerSize=5;
%%
subplot(1,3,1)
totalclass_bsl=logical(sum(bsxfun(@ge,nodematbsl(:,4:9),mean(nodematbsl(:,4:9))+3*std(nodematbsl(:,4:9))),2));
totalclass_bsl_ind=zeros(length(totalclass_bsl),6);
for i=1:6
totalclass_bsl_ind(:,i)=bsxfun(@ge,nodematbsl(:,3+i),mean(nodematbsl(:,3+i))+3*std(nodematbsl(:,3+i)));
end

ptotal_ind=sum(totalclass_bsl_ind)/length(totalclass_bsl);
pint=zeros(1,6);
for i=1:6
   
    pmat=nchoosek(ptotal_ind,i);
    nperms=size(pmat,1);
    invpmat=1-cell2mat(cellfun(@(x,y) setdiff(x,y),mat2cell(repmat(ptotal_ind,...
        nperms,1),ones(nperms,1)),mat2cell(pmat,ones(nperms,1)),'UniformOutput',0));
    pint(i)=sum(prod([pmat,invpmat],2));
    
    
end
dist_bsl=histcounts(sum(totalclass_bsl_ind,2),0.5:6.5);
dist_bsl_exp=pint.*size(nodematbsl,1);

bar(1:6,dist_bsl(1:6)./size(nodematbsl,1),0.4,'FaceColor',colors(1,:))
hold on
bar(1.3:6.3,dist_bsl_exp(1:6)./size(nodematbsl,1),0.4,'FaceColor',colors(2,:))
set(gca,'YScale','log','box','off')
xlabel('Number of Frequencies')
ylabel('Probability')
title('Int')
xlim([0,7])

subplot(1,3,2)
totalclass_late=logical(sum(bsxfun(@ge,nodematpca(:,33:36),mean(nodematpca(:,33:36))+3*std(nodematpca(:,33:36))),2));
totalclass_late_ind=zeros(length(totalclass_late),4);
for i=1:4
totalclass_late_ind(:,i)=bsxfun(@ge,nodematpca(:,32+i),mean(nodematpca(:,32+i))+3*std(nodematpca(:,32+i)));
end

ptotal_ind=sum(totalclass_late_ind)/length(totalclass_late);
pint=zeros(1,4);
for i=1:4
   
    pmat=nchoosek(ptotal_ind,i);
    nperms=size(pmat,1);
    invpmat=1-cell2mat(cellfun(@(x,y) setdiff(x,y),mat2cell(repmat(ptotal_ind,...
        nperms,1),ones(nperms,1)),mat2cell(pmat,ones(nperms,1)),'UniformOutput',0));
    pint(i)=sum(prod([pmat,invpmat],2));
    
    
end
dist_late=histcounts(sum(totalclass_late_ind,2),0.5:4.5);
dist_late_exp=pint.*size(nodematpca,1);

bar(1:4,dist_late(1:4)./size(nodematpca,1),0.4,'FaceColor',colors(1,:))
hold on
bar(1.3:4.3,dist_late_exp(1:4)./size(nodematpca,1),0.4,'FaceColor',colors(2,:))
set(gca,'YScale','log')
set(gca,'YScale','log','box','off')
xlabel('Number of Frequencies')
title('Late')
xlim([0,5])


subplot(1,3,3)
totalclass_post=logical(sum(bsxfun(@ge,nodematpca(:,37:45),mean(nodematpca(:,37:45))+3*std(nodematpca(:,37:45))),2));
totalclass_post_ind=zeros(length(totalclass_post),9);
for i=1:9
totalclass_post_ind(:,i)=bsxfun(@ge,nodematpca(:,36+i),mean(nodematpca(:,36+i))+3*std(nodematpca(:,36+i)));
end

ptotal_ind=sum(totalclass_post_ind)/length(totalclass_post);
pint=zeros(1,9);
for i=1:9
   
    pmat=nchoosek(ptotal_ind,i);
    nperms=size(pmat,1);
    invpmat=1-cell2mat(cellfun(@(x,y) setdiff(x,y),mat2cell(repmat(ptotal_ind,...
        nperms,1),ones(nperms,1)),mat2cell(pmat,ones(nperms,1)),'UniformOutput',0));
    pint(i)=sum(prod([pmat,invpmat],2));
    
    
end
dist_post=histcounts(sum(totalclass_post_ind,2),0.5:9.5);
dist_post_exp=pint.*size(nodematpca,1);

bar(1:9,dist_post(1:9)./size(nodematpca,1),0.4,'FaceColor',colors(1,:))
hold on
bar(1.3:9.3,dist_post_exp(1:9)./size(nodematpca,1),0.4,'FaceColor',colors(2,:))
set(gca,'YScale','log')
set(gca,'YScale','log','box','off')
xlabel('Number of Frequencies')
title('Post')
xlim([0,10])

[leg,hobj]=legend('Observed','Expected');
set(leg,'box','off','units','normalized','Position',[0.85,0.85,0.1,0.05])
hobj(3).Children.Vertices([1,2,5],1)=0.2
hobj(4).Children.Vertices([1,2,5],1)=0.2
%%
%map=colormap('jet');
tstat_rand=zeros(size(nodemat,1),size(nodemat,2),10000);
tstat_rand_bsl=zeros(size(nodematbsl,2),10000);

for i=1:10000
    tstat_rand(:,:,i)=tcalc(nodemat(:,:,randperm(nnz(logmat(:,1)),nnz(logmat(:,1)&res))),nodemat(:,:,randperm(nnz(logmat(:,5)),nnz(logmat(:,5)&res))),3);
    tstat_rand_bsl(:,i)=tcalc(nodematbsl(randperm(nnz(logmatbsl(:,1)),nnz(logmatbsl(:,1)&resbsl)),:),nodematbsl(randperm(nnz(logmatbsl(:,2)),nnz(logmatbsl(:,2)&resbsl)),:),1);
end
%%
[sig_account,clusts_account]=clust_mass(nodemat(:,:,logmat(:,1)&res),nodemat(:,:,logmat(:,5)&res),3,0.05,0.01,10000,tstat_rand);
[sig_account_bsl,clusts_account_bsl]=clust_mass1d(nodematbsl(logmatbsl(:,1)&resbsl,:),nodematbsl(logmatbsl(:,2)&resbsl,:),1,0.05,0.01,10000,tstat_rand_bsl);

%%
map=zeros(64,3);
map(:,1)=linspace(0,1,64);
map(:,2)=linspace(0,1,64);
map(:,3)=linspace(0.5,0,64);
map=colormap(map);
h1=subplot(32,2,5:2:15);
imagesc([mean(nodematbsl(engelbsl & resbsl,:));reshape(mean(nodematpca(logmat(:,1)&res,:)),9,5)'],[0.025,0.07])
set(gca,'XTick',1:9,'XTickLabel',num2str(threshfreq','%.0f'),'YTick',1:6,'YTickLabel',{'Int';'Pre';'Early';'Middle';'Late';'Post'},'FontSize',8,'box','off')
t1r=title('Engel 1a (Resected)');
%set(t1r,'Units','normalized','FontSize',14,'FontWeight','bold')
%t1r.Position(2)=1.4;
text(-0.24,1.2,'a','Units','Normalized','FontSize',16,'FontWeight','bold')
h1.Position([1,2])=[0.09,0.8];

h2=subplot(32,2,6:2:16);
imagesc([mean(nodematbsl(~engelbsl & resbsl,:));reshape(mean(nodematpca(logmat(:,5)&res,:)),9,5)'],[0.025,0.07])
set(gca,'XTick',1:9,'XTickLabel',num2str(threshfreq','%.0f'),'YTick',1:6,'YTickLabel',{'Int';'Pre';'Early';'Middle';'Late';'Post'},'FontSize',8,'box','off')
h2.Position(2)=0.8;
t1rb=title('>Engel 1a (Resected)');
%set(t1rb,'Units','normalized','FontSize',14,'FontWeight','bold')
text(-0.15,1.2,'b','Units','Normalized','FontSize',16,'FontWeight','bold')
col=colorbar('eastoutside');
col.Position([1,3])=[0.92,0.01];

h3=subplot(32,2,21:2:31);
ind=round(interp1(linspace(-0.015,0.015,size(map,1)),1:size(map,1),[diffactbsl';diffact(:,4,2)]));
cols(:,:,1)=reshape(map(ind,1),9,6)';
cols(:,:,2)=reshape(map(ind,2),9,6)';
cols(:,:,3)=reshape(map(ind,3),9,6)';

%cols(repmat(reshape([p_accountbsl';p_account(:,4,2)]>0.001,9,6)',[1,1,3]))=0.7-(0.7-cols(repmat(reshape([p_accountbsl';p_account(:,4,2)]>0.001,9,6)',[1,1,3])))*0.3;
cols(repmat([~sig_account_bsl;~sig_account'],[1,1,3]))=0.7-(0.7-cols(repmat([~sig_account_bsl;~sig_account'],[1,1,3])))*0.3;

imagesc(cols)
set(gca,'XTick',1:9,'XTickLabel',num2str(threshfreq','%.0f'),'YTick',1:6,'YTickLabel',{'Int';'Pre';'Early';'Middle';'Late';'Post'},'box','off','FontSize',8)
h3.Position([1,2])=[0.09,0.58];
t1r=title('Difference (Resected)');
%set(t1r,'Units','normalized','FontSize',14,'FontWeight','bold');
%set(t1r,'Units','normalized','FontSize',14,'FontWeight','bold')
%t1r.Position(2)=1.4;
text(-0.24,1.2,'c','Units','Normalized','FontSize',16,'FontWeight','bold')
col2=colorbar('eastoutside');
col2.Position([1,3])=[0.44,0.01];
set(col2,'Ticks',[1/26:6/26:25/18],'TickLabels',num2str((-0.012:0.006:0.012)'))


h4=subplot(32,2,22:2:32);

%colors=get(gca,'ColorOrder');
colors=[0,118,192;163,2,52]/255;
colorshade=[186,207,236;228,184,180]/255;

patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,1)&res,37:45))+...
    std(nodematpca(logmat(:,1)&res,37:45))/sqrt(nnz(logmat(:,1)&res)),...
    fliplr(mean(nodematpca(logmat(:,1)&res,37:45))-...
    std(nodematpca(logmat(:,1)&res,37:45))/sqrt(nnz(logmat(:,1)&res)))],...
    1-(1-colors(1,:))*1,'EdgeColor','none','FaceAlpha',0.5)
hold on

patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,5)&res,37:45))+...
    std(nodematpca(logmat(:,5)&res,37:45))/sqrt(nnz(logmat(:,5)&res)),...
    fliplr(mean(nodematpca(logmat(:,5)&res,37:45))-...
    std(nodematpca(logmat(:,5)&res,37:45))/sqrt(nnz(logmat(:,5)&res)))],...
    1-(1-colors(2,:))*1,'EdgeColor','none','FaceAlpha',0.5)

%patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,4)&res,37:45))+...
%    std(nodematpca(logmat(:,4)&res,37:45))/sqrt(nnz(logmat(:,4)&res)),...
%    fliplr(mean(nodematpca(logmat(:,4)&res,37:45))-...
%    std(nodematpca(logmat(:,4)&res,37:45))/sqrt(nnz(logmat(:,4)&res)))],...
%    1-(1-colors(3,:))*1,'EdgeColor','none','FaceAlpha',0.5)

e1=plot(1:9,mean(nodematpca(logmat(:,1)&res,37:45)),'Color',colors(1,:),'LineWidth',2);
ne1=plot(1:9,mean(nodematpca(logmat(:,5)&res,37:45)),'Color',colors(2,:),'LineWidth',2);
%plot(1:9,mean(nodematpca(logmat(:,4)&res,37:45)),'Color',0.5*colors(3,:),'LineWidth',2)

plot(1:9,mean(means(:,37:45,1,2)),'Color',colors(1,:),'LineWidth',2,'LineStyle',':');
%hold on
plot(1:9,mean(means(:,37:45,5,2)),'Color',colors(2,:),'LineWidth',2,'LineStyle',':');
%plot(1:9,mean(means(:,37:45,4)),'Color',colors(3,:),'LineWidth',2,'LineStyle',':');

set(gca,'XTick',1:9,'XTickLabel',num2str(threshfreq','%.0f'),'YLim',[0.025,0.07],'FontSize',8,'YAxisLocation','right')
text(-0.15,2.4,'d','Units','Normalized','FontSize',16,'FontWeight','bold')
t1rb=title('BC (Resected)');
%set(t1rb,'Units','normalized','FontSize',14,'FontWeight','bold')
ylab=ylabel('BC');
%[leg,hobj]=legend([e1,ne1],'Engel 1a','>Engel 1a');
%set(leg,'box','off','units','normalized');
%leg.Position(1:2)=[0.76,0.66];
%hobj(3).XData(1)=0.3;
%hobj(5).XData(1)=0.3;
h4.Position(2)=0.58;

pos=h4.Position;
h4.Position=[pos(1:3),pos(4)/2];
ylab.Position(2)=h4.YLim(2);
text(-0.05,0.5,'Post','Units','Normalized','FontSize',12,'Rotation',90,'horizontalAlignment','center')

h4_2=axes('position',[pos(1),pos(2)+pos(4)/2,pos(3),pos(4)/2]);


patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,1)&res,19:27))+...
    std(nodematpca(logmat(:,1)&res,19:27))/sqrt(nnz(logmat(:,1)&res)),...
    fliplr(mean(nodematpca(logmat(:,1)&res,19:27))-...
    std(nodematpca(logmat(:,1)&res,19:27))/sqrt(nnz(logmat(:,1)&res)))],...
    1-(1-colors(1,:))*1,'EdgeColor','none','FaceAlpha',0.5)
hold on

patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,5)&res,19:27))+...
    std(nodematpca(logmat(:,5)&res,19:27))/sqrt(nnz(logmat(:,5)&res)),...
    fliplr(mean(nodematpca(logmat(:,5)&res,19:27))-...
    std(nodematpca(logmat(:,5)&res,19:27))/sqrt(nnz(logmat(:,5)&res)))],...
    1-(1-colors(2,:))*1,'EdgeColor','none','FaceAlpha',0.5)

e1=plot(1:9,mean(nodematpca(logmat(:,1)&res,19:27)),'Color',colors(1,:),'LineWidth',2);
ne1=plot(1:9,mean(nodematpca(logmat(:,5)&res,19:27)),'Color',colors(2,:),'LineWidth',2);

plot(1:9,mean(means(:,19:27,1,2)),'Color',colors(1,:),'LineWidth',2,'LineStyle',':');
plot(1:9,mean(means(:,19:27,5,2)),'Color',colors(2,:),'LineWidth',2,'LineStyle',':');

set(gca,'XTick',1:9,'XTickLabel',[],'YLim',[0.025,0.07],'FontSize',8,'YAxisLocation','right')
t1rb=title('BC (Resected)');
%ylab=ylabel('BC');
text(-0.05,0.5,'Mid','Units','Normalized','FontSize',12,'Rotation',90,'horizontalAlignment','center')


h5=subplot(32,2,37:2:47);
m_all1a=mean(nodematpca(logmat(:,1),37:45));
sem_all1a=std(nodematpca(logmat(:,1),37:45))./sqrt(nnz(logmat(:,1)));
patch([1:9,9:-1:1],[m_all1a+sem_all1a,fliplr(m_all1a-sem_all1a)],colors(1,:),'FaceAlpha',0.5,'EdgeColor','none')
hold on

m_alln1a=mean(nodematpca(logmat(:,5),37:45));
sem_alln1a=std(nodematpca(logmat(:,5),37:45))./sqrt(nnz(logmat(:,5)));
patch([1:9,9:-1:1],[m_alln1a+sem_alln1a,fliplr(m_alln1a-sem_alln1a)],colors(2,:),'FaceAlpha',0.5,'EdgeColor','none')
plot(1:9,mean(nodematpca(logmat(:,1),37:45)),'Color',colors(1,:),'LineWidth',2)
hold on
plot(1:9,mean(nodematpca(logmat(:,5),37:45)),'Color',colors(2,:),'LineWidth',2)
ylim([0.025 0.07])
set(gca,'XTick',1:9,'XTickLabel',num2str(threshfreq','%.0f'),'box','off','FontSize',8,'YAxisLocation','right')
xlabel('Frequency (Hz)')
text(-0.24,2.4,'e','Units','Normalized','FontSize',16,'FontWeight','bold')
t1rb=title('BC (All)');
%set(t1rb,'Units','normalized','FontSize',14,'FontWeight','bold')
h5.Position(1:2)=[0.09,0.36];

pos=h5.Position;
h5.Position=[pos(1:3),pos(4)/2];
text(-0.05,0.5,'Post','Units','Normalized','FontSize',12,'Rotation',90,'horizontalAlignment','center')

h5_2=axes('position',[pos(1),pos(2)+pos(4)/2,pos(3),pos(4)/2]);

m_all1a=mean(nodematpca(logmat(:,1),19:27));
sem_all1a=std(nodematpca(logmat(:,1),19:27))./sqrt(nnz(logmat(:,1)));
patch([1:9,9:-1:1],[m_all1a+sem_all1a,fliplr(m_all1a-sem_all1a)],colors(1,:),'FaceAlpha',0.5,'EdgeColor','none')
hold on

m_alln1a=mean(nodematpca(logmat(:,5),19:27));
sem_alln1a=std(nodematpca(logmat(:,5),19:27))./sqrt(nnz(logmat(:,5)));
patch([1:9,9:-1:1],[m_alln1a+sem_alln1a,fliplr(m_alln1a-sem_alln1a)],colors(2,:),'FaceAlpha',0.5,'EdgeColor','none')
plot(1:9,mean(nodematpca(logmat(:,1),19:27)),'Color',colors(1,:),'LineWidth',2)
hold on
plot(1:9,mean(nodematpca(logmat(:,5),19:27)),'Color',colors(2,:),'LineWidth',2)
ylim([0.025 0.07])
set(gca,'XTick',1:9,'XTickLabel',[],'box','off','FontSize',8,'YAxisLocation','right')
t1rb=title('BC (All)');
text(-0.05,0.5,'Mid','Units','Normalized','FontSize',12,'Rotation',90,'horizontalAlignment','center')

h6=subplot(32,2,38:2:48);
m_unres1a=mean(nodematpca(logmat(:,1)&~res,37:45));
sem_unres1a=std(nodematpca(logmat(:,1)&~res,37:45))./sqrt(nnz(logmat(:,1)&~res));
e_patch=patch([1:9,9:-1:1],[m_unres1a+sem_unres1a,fliplr(m_unres1a-sem_unres1a)],mean([colors(1,:);1,1,1]),'FaceAlpha',0.5,'EdgeColor','none');
hold on

m_unresn1a=mean(nodematpca(logmat(:,5)&~res,37:45));
sem_unresn1a=std(nodematpca(logmat(:,5)&~res,37:45))./sqrt(nnz(logmat(:,5)&~res));
ne_patch=patch([1:9,9:-1:1],[m_unresn1a+sem_unresn1a,fliplr(m_unresn1a-sem_unresn1a)],mean([colors(2,:);1,1,1]),'FaceAlpha',0.5,'EdgeColor','none');
plot(1:9,mean(nodematpca(logmat(:,1)&~res,37:45)),'Color',colors(1,:),'LineWidth',2);
hold on
plot(1:9,mean(nodematpca(logmat(:,5)&~res,37:45)),'Color',colors(2,:),'LineWidth',2)
ylim([0.025 0.07])
set(gca,'XTick',1:9,'XTickLabel',num2str(threshfreq','%.0f'),'box','off','FontSize',8,'YAxisLocation','right')
%ylabel({'Betweenness';'Centrality'},'FontSize',14)
xlabel('Frequency (Hz)')
ylabel('BC')
text(-0.15,2.4,'f','Units','Normalized','FontSize',16,'FontWeight','bold')
t1rb=title('BC (Not Resected)');
%set(t1rb,'Units','normalized','FontSize',14,'FontWeight','bold')

h6.Position(2)=0.36;

pos=h6.Position;
h6.Position=[pos(1:3),pos(4)/2];
text(-0.05,0.5,'Post','Units','Normalized','FontSize',12,'Rotation',90,'horizontalAlignment','center')

h6_2=axes('position',[pos(1),pos(2)+pos(4)/2,pos(3),pos(4)/2]);

m_unres1a=mean(nodematpca(logmat(:,1)&~res,19:27));
sem_unres1a=std(nodematpca(logmat(:,1)&~res,19:27))./sqrt(nnz(logmat(:,1)&~res));
e_patch=patch([1:9,9:-1:1],[m_unres1a+sem_unres1a,fliplr(m_unres1a-sem_unres1a)],mean([colors(1,:);1,1,1]),'FaceAlpha',0.5,'EdgeColor','none');
hold on

m_unresn1a=mean(nodematpca(logmat(:,5)&~res,19:27));
sem_unresn1a=std(nodematpca(logmat(:,5)&~res,19:27))./sqrt(nnz(logmat(:,5)&~res));
ne_patch=patch([1:9,9:-1:1],[m_unresn1a+sem_unresn1a,fliplr(m_unresn1a-sem_unresn1a)],mean([colors(2,:);1,1,1]),'FaceAlpha',0.5,'EdgeColor','none');
plot(1:9,mean(nodematpca(logmat(:,1)&~res,19:27)),'Color',colors(1,:),'LineWidth',2);
hold on
plot(1:9,mean(nodematpca(logmat(:,5)&~res,19:27)),'Color',colors(2,:),'LineWidth',2)
ylim([0.025 0.07])
set(gca,'XTick',1:9,'XTickLabel',[],'box','off','FontSize',8,'YAxisLocation','right')
ylabel('BC')
t1rb=title('BC (Not Resected)');
text(-0.05,0.5,'Mid','Units','Normalized','FontSize',12,'Rotation',90,'horizontalAlignment','center')

h(1)=subplot(32,2,53:2:59);
resclass_mid=logical(sum(bsxfun(@ge,nodematpca(:,[23:26]),mean(nodematpca(:,[23:26]))+3*std(nodematpca(:,[23:26]))),2));
resclass_mid=resclass_mid(res);
ptindres=ptind(res);
resclass_mid=accumarray(ptindres,resclass_mid,[],@mean);

cfier=0:0.01:0.5;
%[sorted,ind]=sort(totalclass);
fpr_mid=zeros(length(cfier),1);
tpr_mid=fpr_mid;
for i=1:length(cfier)
tpr_mid(i)=nnz(~ptengel1a(find(resclass_mid>=cfier(i))))./nnz(~ptengel1a);
fpr_mid(i)=nnz(ptengel1a(find(resclass_mid>=cfier(i))))./nnz(ptengel1a);
end
%[~,opt]=min((1-tpr).^2+fpr.^2); %distance from corner
[~,opt_mid]=max(tpr_mid-fpr_mid); %J statistic
mid_roc=plot(fpr_mid,tpr_mid,'color',[0.5 0.5 0.5],'LineWidth',2);
hold on

resclass_post=logical(sum(bsxfun(@ge,nodematpca(:,[40:42]),mean(nodematpca(:,[40:42]))+3*std(nodematpca(:,[40:42]))),2));
resclass_post=resclass_post(res);
ptindres=ptind(res);
resclass_post=accumarray(ptindres,resclass_post,[],@mean);

cfier=0:0.01:0.5;
%[sorted,ind]=sort(totalclass);
fpr_post=zeros(length(cfier),1);
tpr_post=fpr_post;
for i=1:length(cfier)
tpr_post(i)=nnz(~ptengel1a(find(resclass_post>=cfier(i))))./nnz(~ptengel1a);
fpr_post(i)=nnz(ptengel1a(find(resclass_post>=cfier(i))))./nnz(ptengel1a);
end
%[~,opt]=min((1-tpr).^2+fpr.^2); %distance from corner
[~,opt_post]=max(tpr_post-fpr_post); %J statistic
post_roc=plot(fpr_post,tpr_post,'color',[0 0 0],'LineWidth',2);

plot(fpr_post(opt_post),tpr_post(opt_post),'Marker','d','MarkerFaceColor',...
    [0 0 0],'MarkerSize',10,'MarkerEdgeColor',[0 0 0]);
plot(fpr_mid(opt_mid),tpr_mid(opt_mid),'Marker','d','MarkerFaceColor',...
    [0.5 0.5 0.5],'MarkerSize',10,'MarkerEdgeColor',[0.5 0.5 0.5]);

hold on
%q=annotation('arrow');
%set(q,'parent',h(1),'position',[0,1,0.9*fpr(opt),0.9*(tpr(opt)-1)],...
%    'color','k','LineStyle','--','HeadLength',5,'HeadWidth',5)
hold on
line([0,1],[0,1],'color',[0.5,0.5,0.5])
xlabel('FPR')
ylabel('TPR')
xlim([0 1])
ylim([0 1])
set(gca,'box','off','FontSize',8);
text(-0.55,1.2,'g','Units','Normalized','FontSize',16,'FontWeight','bold')

[leg_roc,rocobj]=legend([mid_roc,post_roc],'Mid','Post');
set(leg_roc,'box','off','Position',[0.15,0.09,0.1,0.05]);

rocobj(3).XData(1)=0.4;
rocobj(5).XData(1)=0.4;

h(2)=subplot(32,2,61:2:63);
plotSpread(resclass_mid,'distributionIdx',~ptengel1a,'distributionColors',colors(1:2,:),'spreadWidth',1.5);
optline=line([0,3],[cfier(opt_mid),cfier(opt_mid)],'color','k','LineStyle','-');
plotSpread(resclass_post,'distributionIdx',~ptengel1a,'distributionColors',colors(1:2,:),'spreadWidth',1.5,'xValues',[4,5]);
optline=line([3,6],[cfier(opt_post),cfier(opt_post)],'color','k','LineStyle','-');

%legopt=legend(optline,'Optimal Discriminator');
%set(legopt,'box','off')
ylabel('Fraction')
set(h(2),'YAxisLocation','right','XTick',[1,2,4,5],'XTickLabel',{'\color{white}>\color{black}1a','>1a'},'XTickLabelRotation',90,'FontSize',8)
xlim([0 6])
text(-0.2,1.2,'h','Units','Normalized','FontSize',16,'FontWeight','bold')
text(0.25,1.1,'Mid','Units','Normalized','FontSize',10,'horizontalAlignment','center')
text(0.75,1.1,'Post','Units','Normalized','FontSize',10,'horizontalAlignment','center')

h(3)=subplot(32,2,54:2:64);
%scatter(totalclass(ptengel1a),resclass(ptengel1a),'filled');
%hold on
%scatter(totalclass(~ptengel1a),resclass(~ptengel1a),'filled');
%text(-0.15,1.2,'i','Units','Normalized','FontSize',16,'FontWeight','bold')
ylabel('Fraction Resected');
xlabel('Fraction Present');
set(h(3),'box','off','FontSize',8,'YAxisLocation','right');
legobj(1)=line([0,1],[0,1],'color',colors(1,:),'LineWidth',2,'Visible','off');
legobj(2)=line([0,1],[0,1],'color',colors(2,:),'LineWidth',2,'Visible','off');
legobj(3)=line([0,1],[0,1],'color',colors(1,:),'LineWidth',2,'LineStyle',':','Visible','off');
legobj(4)=line([0,1],[0,1],'color',colors(2,:),'LineWidth',2,'LineStyle',':','Visible','off');
%legobj(5)=line([0,1],[0,1],'color','k','LineWidth',2,'LineStyle','-','Visible','off');
[leg,lobj]=legend(legobj,'Engel 1a (Data)','Engel >1a (Data)','Engel 1a (Simulated Random Resection)','Engel >1a (Simulated Random Resection)');
set(leg,'box','off','color',[0,0,0])
leg.Position(1)=0.6;
for i=1:4
lobj(4+2*i).XData=[0.01,0.15];
end

for i=1:4
    lobj(i).Color=[0,0,0];
end
h(3).Position(2)=0.08;
h(1).Position(1)=0.09;
h(1).Position([2,4])=h(3).Position([2,4]);
h(2).Position([2,4])=h(3).Position([2,4]);
h(1).Position(3)=h1.Position(4);
h(2).Position(1)=h(1).Position(1)+h(1).Position(3)+0.05;
h(2).Position(3)=h(3).Position(3)-h(1).Position(3)-0.05;
h(3).Visible='off';
h(2).Children(3).MarkerSize=5;
h(2).Children(4).MarkerSize=5;
h(2).Children(6).MarkerSize=5;
h(2).Children(7).MarkerSize=5;

%%
subplot(1,3,1)
totalclass_mid=logical(sum(bsxfun(@ge,nodematpca(:,23:26),mean(nodematpca(:,23:26))+3*std(nodematpca(:,23:26))),2));
totalclass_mid_ind=zeros(length(totalclass_mid),4);
for i=1:4
totalclass_mid_ind(:,i)=bsxfun(@ge,nodematpca(:,22+i),mean(nodematpca(:,22+i))+3*std(nodematpca(:,22+i)));
end
totalclass_mid=totalclass_mid(res);
totalclass_mid_ind=totalclass_mid_ind(res,:);

ptotal_ind=sum(totalclass_mid_ind)/length(totalclass_mid);
pint=zeros(1,4);
for i=1:4
   
    pmat=nchoosek(ptotal_ind,i);
    nperms=size(pmat,1);
    invpmat=1-cell2mat(cellfun(@(x,y) setdiff(x,y),mat2cell(repmat(ptotal_ind,...
        nperms,1),ones(nperms,1)),mat2cell(pmat,ones(nperms,1)),'UniformOutput',0));
    pint(i)=sum(prod([pmat,invpmat],2));
    
    
end
dist_mid=histcounts(sum(totalclass_mid_ind,2),0.5:4.5);
dist_mid_exp=pint.*nnz(res);

bar(1:4,dist_mid(1:4)./size(nodematpca,1),0.4,'FaceColor',colors(1,:))
hold on
bar(1.3:4.3,dist_mid_exp(1:4)./size(nodematpca,1),0.4,'FaceColor',colors(2,:))
set(gca,'YScale','log','box','off')
xlabel('Number of Frequencies')
ylabel('Probability')
title('Mid')
xlim([0,5])

subplot(1,3,2)
totalclass_post=logical(sum(bsxfun(@ge,nodematpca(:,40:42),mean(nodematpca(:,40:42))+3*std(nodematpca(:,40:42))),2));
totalclass_post_ind=zeros(length(totalclass_post),3);
for i=1:3
totalclass_post_ind(:,i)=bsxfun(@ge,nodematpca(:,39+i),mean(nodematpca(:,39+i))+3*std(nodematpca(:,39+i)));
end
totalclass_post=totalclass_post(res);
totalclass_post_ind=totalclass_post_ind(res,:);

ptotal_ind=sum(totalclass_post_ind)/length(totalclass_post);
pint=zeros(1,3);
for i=1:3
   
    pmat=nchoosek(ptotal_ind,i);
    nperms=size(pmat,1);
    invpmat=1-cell2mat(cellfun(@(x,y) setdiff(x,y),mat2cell(repmat(ptotal_ind,...
        nperms,1),ones(nperms,1)),mat2cell(pmat,ones(nperms,1)),'UniformOutput',0));
    pint(i)=sum(prod([pmat,invpmat],2));
    
    
end
dist_post=histcounts(sum(totalclass_post_ind,2),0.5:3.5);
dist_post_exp=pint.*nnz(res);

bar(1:3,dist_post(1:3)./size(nodematpca,1),0.4,'FaceColor',colors(1,:))
hold on
bar(1.3:3.3,dist_post_exp(1:3)./size(nodematpca,1),0.4,'FaceColor',colors(2,:))
set(gca,'YScale','log','box','off')
xlabel('Number of Frequencies')
title('Post')
xlim([0,4])

subplot(1,3,3)
totalclass_ove=logical(sum(bsxfun(@ge,nodematpca(:,[24,33,42]),mean(nodematpca(:,[24,33,42]))+3*std(nodematpca(:,[24,33,42]))),2));
totalclass_ove_ind=zeros(length(totalclass_ove),4);
pts=[24,33,42];
for i=1:3
totalclass_ove_ind(:,i)=bsxfun(@ge,nodematpca(:,pts(i)),mean(nodematpca(:,pts(i)))+3*std(nodematpca(:,pts(i))));
end
totalclass_ove=totalclass_ove(res);
totalclass_ove_ind=totalclass_ove_ind(res,:);

ptotal_ind=sum(totalclass_ove_ind)/length(totalclass_ove);
pint=zeros(1,3);
for i=1:3
   
    pmat=nchoosek(ptotal_ind,i);
    nperms=size(pmat,1);
    invpmat=1-cell2mat(cellfun(@(x,y) setdiff(x,y),mat2cell(repmat(ptotal_ind,...
        nperms,1),ones(nperms,1)),mat2cell(pmat,ones(nperms,1)),'UniformOutput',0));
    pint(i)=sum(prod([pmat,invpmat],2));
    
    
end
dist_ove=histcounts(sum(totalclass_ove_ind,2),0.5:3.5);
dist_ove_exp=pint.*nnz(res);

bar(1:3,dist_ove(1:3)./size(nodematpca,1),0.4,'FaceColor',colors(1,:))
hold on
bar(1.3:3.3,dist_ove_exp(1:3)./size(nodematpca,1),0.4,'FaceColor',colors(2,:))
set(gca,'YScale','log','box','off')
xlabel('Number of Time Points')
title('30 Hz')
xlim([0,4])

[leg,hobj]=legend('Observed','Expected');
set(leg,'box','off','units','normalized','Position',[0.85,0.85,0.1,0.05])
hobj(3).Children.Vertices([1,2,5],1)=0.2
hobj(4).Children.Vertices([1,2,5],1)=0.2
%%
subplot 231
patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,1)&res,19:27))+...
    std(nodematpca(logmat(:,1)&res,19:27))/sqrt(nnz(logmat(:,1)&res)),...
    fliplr(mean(nodematpca(logmat(:,1)&res,19:27))-...
    std(nodematpca(logmat(:,1)&res,19:27))/sqrt(nnz(logmat(:,1)&res)))],...
    1-(1-colors(1,:))*1,'EdgeColor','none','FaceAlpha',0.5)
hold on
plot(1:9,mean(nodematpca(logmat(:,1)&res,19:27)),'Color',colors(1,:),'LineWidth',3)
plot(1:9,mean(means(:,19:27,1,2)),'Color',colors(1,:),'LineWidth',3,'LineStyle',':');

ylim([0.028 0.068])
ylabel('Betweenness Centrality')
text(-0.55,1.01,'A','Units','Normalized','FontSize',16,'FontWeight','bold')

set(gca,'XTick',1:9,'XTickLabel',[],'FontSize',12)

%subplot 232
%patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,2)&res,19:27))+...
%    std(nodematpca(logmat(:,2)&res,19:27))/sqrt(nnz(logmat(:,2)&res)),...
%    fliplr(mean(nodematpca(logmat(:,2)&res,19:27))-...
%    std(nodematpca(logmat(:,2)&res,19:27))/sqrt(nnz(logmat(:,2)&res)))],...
%    1-(1-colors(2,:))*1,'EdgeColor','none','FaceAlpha',0.5)
%hold on
%plot(1:9,mean(nodematpca(logmat(:,2)&res,19:27)),'Color',colors(2,:),'LineWidth',3)
%plot(1:9,mean(means(:,19:27,2)),'Color',colors(2,:),'LineWidth',3,'LineStyle',':');

%ylim([0.028 0.068])
%set(gca,'XTick',[],'FontSize',12)
%text(-0.3,1.01,'B','Units','Normalized','FontSize',16,'FontWeight','bold')


%patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,3)&res,19:27))+...
%    std(nodematpca(logmat(:,3)&res,19:27))/sqrt(nnz(logmat(:,3)&res)),...
%    fliplr(mean(nodematpca(logmat(:,3)&res,19:27))-...
%    std(nodematpca(logmat(:,3)&res,19:27))/sqrt(nnz(logmat(:,3)&res)))],...
%    1-(1-colors(3,:))*1,'EdgeColor','none','FaceAlpha',0.5)

subplot 232
patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,5)&res,19:27))+...
    std(nodematpca(logmat(:,5)&res,19:27))/sqrt(nnz(logmat(:,5)&res)),...
    fliplr(mean(nodematpca(logmat(:,5)&res,19:27))-...
    std(nodematpca(logmat(:,5)&res,19:27))/sqrt(nnz(logmat(:,5)&res)))],...
    1-(1-colors(2,:))*1,'EdgeColor','none','FaceAlpha',0.5)
hold on
plot(1:9,mean(nodematpca(logmat(:,5)&res,19:27)),'Color',colors(2,:),'LineWidth',3)
plot(1:9,mean(means(:,19:27,5,2)),'Color',colors(2,:),'LineWidth',3,'LineStyle',':');
text(-0.3,1.01,'B','Units','Normalized','FontSize',16,'FontWeight','bold')

ylim([0.028 0.068])
set(gca,'XTick',1:9,'XTickLabel',[],'FontSize',12)

h4=subplot(2,3,4);

patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,1)&res,19:27))+...
    std(nodematpca(logmat(:,1)&res,19:27))/sqrt(nnz(logmat(:,1)&res)),...
    fliplr(mean(nodematpca(logmat(:,1)&res,19:27))-...
    std(nodematpca(logmat(:,1)&res,19:27))/sqrt(nnz(logmat(:,1)&res)))],...
    1-(1-colors(1,:))*1,'EdgeColor','none','FaceAlpha',0.5)
hold on
l1=plot(1:9,mean(nodematpca(logmat(:,1)&res,19:27)),'Color',colors(1,:),'LineWidth',3);

%patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,2)&res,19:27))+...
%    std(nodematpca(logmat(:,2)&res,19:27))/sqrt(nnz(logmat(:,2)&res)),...
%    fliplr(mean(nodematpca(logmat(:,2)&res,19:27))-...
%    std(nodematpca(logmat(:,2)&res,19:27))/sqrt(nnz(logmat(:,2)&res)))],...
%    1-(1-colors(2,:))*1,'EdgeColor','none','FaceAlpha',0.5)
%l2=plot(1:9,mean(nodematpca(logmat(:,2)&res,19:27)),'Color',colors(2,:),'LineWidth',3);

patch([1:9,9:-1:1],[mean(nodematpca(logmat(:,5)&res,19:27))+...
    std(nodematpca(logmat(:,5)&res,19:27))/sqrt(nnz(logmat(:,5)&res)),...
    fliplr(mean(nodematpca(logmat(:,5)&res,19:27))-...
    std(nodematpca(logmat(:,5)&res,19:27))/sqrt(nnz(logmat(:,5)&res)))],...
    1-(1-colors(2,:))*1,'EdgeColor','none','FaceAlpha',0.5)
l3=plot(1:9,mean(nodematpca(logmat(:,5)&res,19:27)),'Color',colors(2,:),'LineWidth',3);
ylim([0.028 0.068])

set(gca,'XTick',1:9,'XTickLabel',num2str(threshfreq','%.0f'),'FontSize',12,'XTickLabelRotation',90)
h4.Position(2)=0.2;
leg=legend([l1,l3],'Engel 1a','>Engel 1a');
set(leg,'Orientation','Vertical','box','off','Units','Normalized','interpreter','tex')
leg.Position(1:2)=[0.65,0.8];
xlabel('Frequency (Hz)','FontSize',12)
ylabel('Betweenness Centrality');
text(-0.55,1.01,'C','Units','Normalized','FontSize',16,'FontWeight','bold')

h5=subplot(2,3,5);

plot(1:9,mean(means(:,19:27,1,2)),'Color',colors(1,:),'LineWidth',3,'LineStyle',':');
hold on
%plot(1:9,mean(means(:,19:27,2)),'Color',colors(2,:),'LineWidth',3,'LineStyle',':');
plot(1:9,mean(means(:,19:27,5,2)),'Color',colors(2,:),'LineWidth',3,'LineStyle',':');
set(gca,'box','off')

ylim([0.028 0.048])
set(gca,'XTick',1:9,'XTickLabel',num2str(threshfreq','%.0f'),'FontSize',12,'XTickLabelRotation',90)
solid=plot([2,4],[0.03,0.04],'Visible','off','color','k','LineStyle','-','LineWidth',3);
dotted=plot([2,4],[0.03,0.04],'Visible','off','color','k','LineStyle',':','LineWidth',3);
leg2=legend([solid,dotted],'Resection','Simulated Random Resection');
set(leg2,'Orientation','Vertical','box','off','Units','Normalized','interpreter','tex')
leg2.Position(1:2)=[0.65,0.7];
xlabel('Frequency (Hz)','FontSize',12)
h5.Position(2)=0.2;
text(-0.3,1.01,'D','Units','Normalized','FontSize',16,'FontWeight','bold')

h6=subplot(2,3,6);

relmean=zeros(3,9);
relmean(1,:)=mean(nodematpca(logmat(:,1)&res,19:27))./mean(means(:,19:27,1,2));
%relmean(2,:)=mean(nodematpca(logmat(:,2)&res,19:27))./mean(means(:,19:27,2,2));
relmean(3,:)=mean(nodematpca(logmat(:,5)&res,19:27))./mean(means(:,19:27,5,2));


err=zeros(3,9);
err(1,:)=std(nodematpca(logmat(:,1)&res,19:27))./(sqrt(nnz(logmat(:,1)))*mean(means(:,19:27,1,2)));
%err(2,:)=std(nodematpca(logmat(:,2)&res,19:27))./(sqrt(nnz(logmat(:,2)))*mean(means(:,19:27,2,2)));
err(3,:)=std(nodematpca(logmat(:,5)&res,19:27))./(sqrt(nnz(logmat(:,5)))*mean(means(:,19:27,5,2)));

patch([1:9,9:-1:1],[relmean(1,:)+err(1,:),fliplr(relmean(1,:)-err(1,:))],...
    1-(1-colors(1,:))*1,'EdgeColor','none','FaceAlpha',0.5)
hold on
%patch([1:9,9:-1:1],[relmean(2,:)+err(2,:),fliplr(relmean(2,:)-err(2,:))],...
%    1-(1-colors(2,:))*1,'EdgeColor','none','FaceAlpha',0.5)
patch([1:9,9:-1:1],[relmean(3,:)+err(3,:),fliplr(relmean(3,:)-err(3,:))],...
    1-(1-colors(2,:))*1,'EdgeColor','none','FaceAlpha',0.5)

plot(1:9,mean(nodematpca(logmat(:,1)&res,19:27))./mean(means(:,19:27,1,2)),'Color',colors(1,:),'LineWidth',3)
%plot(1:9,mean(nodematpca(logmat(:,2)&res,19:27))./mean(means(:,19:27,2,2)),'Color',colors(2,:),'LineWidth',3)
plot(1:9,mean(nodematpca(logmat(:,5)&res,19:27))./mean(means(:,19:27,5,2)),'Color',colors(2,:),'LineWidth',3)
set(gca,'XTick',1:9,'XTickLabel',num2str(threshfreq','%.0f'),'YAxisLocation','right','FontSize',12,'XTickLabelRotation',90)
xlabel('Frequency (Hz)','FontSize',12)
ylab=ylabel({'Relative BC'},'Rotation',-90);
set(ylab,'Units','normalized');
ylab.Position(1)=1.3;
h6.Position(2)=0.2;
text(-0.3,1.01,'E','Units','Normalized','FontSize',16,'FontWeight','bold')

%%
colors=colormap(jet(64));
close gcf
f=figure;
%for j=1:5
bin=round(interp1(linspace(min(eeg2(1).szraw(1).metric_dtf(:,2,1,1)),...
    max(eeg2(1).szraw(1).metric_dtf(:,2,1,1)),size(colors,1)),1:size(colors,1),...
    eeg2(1).szraw(1).metric_dtf(:,2,1,1)));
bg=biograph(eeg2(1).szraw(1).dtf_adj_cube(:,:,1,1));
for i=1:length(eeg2(1).leads)
    bg.Nodes(i).Color=colors(bin(i),:);
    bg.Nodes(i).ID=num2str(i);
    if eeg2(1).res(i)==0
        bg.Nodes(i).Shape='ellipse';
        bg.Nodes(i).LineColor=colors(bin(i),:);
    else
        bg.Nodes(i).Shape='trapezium';
        bg.Nodes(i).LineColor=[0 0 0];
        bg.Nodes(i).LineWidth=2;
    end
    bg.Nodes(i).Size=[30,10];
    bg.Nodes(i).FontSize=10;
    if bin(i)<10
        bg.Nodes(i).TextColor=[1,1,1];
    end
end
bg.NodeAutoSize='off';
bg.layouttype='radial';
bg.ArrowSize=3;
for i=1:length(bg.Edges)
    bg.Edges(i).LineWidth=1;
end

bg_view=biograph.bggui(bg);
%s(1)=subplot(1,2,1);
s=copyobj(bg_view.biograph.hgAxes,f);
s=subplot(6,1,[1:3],s);
set(s,'units','normalized')
%end
title('Sample Pre-Seizure Network at 2.5 Hz','FontSize',12)

h1=subplot(6,1,4);
colormap('jet');
clims=[min(eeg2(1).szraw(1).metric_dtf(:,2,1,1)),max(eeg2(1).szraw(1).metric_dtf(:,2,1,1))];
imagesc(nodemat(:,:,15)',clims)
title('Node 15')
set(gca,'XTick',[],'YTickLabel',{'Pre','Early','Mid','Late','Post'},'FontSize',8)
h2=subplot(6,1,5);
imagesc(nodemat(:,:,30)',clims)
title('Node 30')
set(gca,'XTick',1:9,'XTickLabel',num2str(threshfreq','%.0f'),'YTickLabel',{'Pre','Early','Mid','Late','Post'},'FontSize',8,'XTickLabelRotation',90)
xlabel('Frequency (Hz)')
h3=subplot(6,1,6);
imagesc(nodemat(:,:,7)',clims)
set(gca,'XTick',[],'YTick',1:5,'YTickLabel',{'Pre','Early','Mid','Late','Post'},'FontSize',8)
title('Node 7')

%subplot(3,3,[3,6,9])
%for i=10:20
%    surf(i/10*ones(10,6),squeeze(nodemat(:,:,i)),'CDataMapping','scaled')
%    %surf(zeros(9,5),squeeze(test(5,:,:)))
%    hold on
%end
%set(gca,'CLim',[0 0.3]);

s.Position=[-0.06 0.07,0.8,0.8];
h1.Position=[0.72,0.43,0.25,0.15];
h2.Position=[0.72,0.14,0.25,0.15];
h3.Position=[0.72,0.72,0.25,0.15];

axes(s);
text(0.1,1.07,'A','Units','Normalized','FontSize',16,'FontWeight','bold')
text(0.9,1.07,'B','Units','Normalized','FontSize',16,'FontWeight','bold')
axes(h3);
axes(h2);
axes(h1);
%%
h1=subplot(3,1,1);
m_all1a=mean(nodematpca(logmat(:,1),19:27));
sem_all1a=std(nodematpca(logmat(:,1),19:27))./sqrt(nnz(logmat(:,1)));
patch([1:9,9:-1:1],[m_all1a+sem_all1a,fliplr(m_all1a-sem_all1a)],colors(1,:),'FaceAlpha',0.5,'EdgeColor','none')
hold on

m_alln1a=mean(nodematpca(logmat(:,5),19:27));
sem_alln1a=std(nodematpca(logmat(:,5),19:27))./sqrt(nnz(logmat(:,5)));
patch([1:9,9:-1:1],[m_alln1a+sem_alln1a,fliplr(m_alln1a-sem_alln1a)],colors(2,:),'FaceAlpha',0.5,'EdgeColor','none')
plot(1:9,mean(nodematpca(logmat(:,1),19:27)),'Color',colors(1,:),'LineWidth',2)
hold on
plot(1:9,mean(nodematpca(logmat(:,5),19:27)),'Color',colors(2,:),'LineWidth',2)
ylim([0.02 0.07])
set(gca,'XTick',1:9,'XTickLabel',[],'box','off','FontSize',12)
ylabel({'Betweenness';'Centrality'},'FontSize',14)

h2=subplot(3,1,2);
m_unres1a=mean(nodematpca(logmat(:,1)&~res,19:27));
sem_unres1a=std(nodematpca(logmat(:,1)&~res,19:27))./sqrt(nnz(logmat(:,1)&~res));
patch([1:9,9:-1:1],[m_unres1a+sem_unres1a,fliplr(m_unres1a-sem_unres1a)],colors(1,:),'FaceAlpha',0.5,'EdgeColor','none')
hold on

m_unresn1a=mean(nodematpca(logmat(:,5)&~res,19:27));
sem_unresn1a=std(nodematpca(logmat(:,5)&~res,19:27))./sqrt(nnz(logmat(:,5)&~res));
patch([1:9,9:-1:1],[m_unresn1a+sem_unresn1a,fliplr(m_unresn1a-sem_unresn1a)],colors(2,:),'FaceAlpha',0.5,'EdgeColor','none')
plot(1:9,mean(nodematpca(logmat(:,1)&~res,19:27)),'Color',colors(1,:),'LineWidth',2)
hold on
plot(1:9,mean(nodematpca(logmat(:,5)&~res,19:27)),'Color',colors(2,:),'LineWidth',2)
ylim([0.02 0.07])
set(gca,'XTick',1:9,'XTickLabel',[],'box','off','FontSize',12)
ylabel({'Betweenness';'Centrality'},'FontSize',14)

h3=subplot(3,1,3);
m_res1a=mean(nodematpca(logmat(:,1)&res,19:27));
sem_res1a=std(nodematpca(logmat(:,1)&res,19:27))./sqrt(nnz(logmat(:,1)&res));
patch([1:9,9:-1:1],[m_res1a+sem_res1a,fliplr(m_res1a-sem_res1a)],colors(1,:),'FaceAlpha',0.5,'EdgeColor','none')
hold on

m_resn1a=mean(nodematpca(logmat(:,5)&res,19:27));
sem_resn1a=std(nodematpca(logmat(:,5)&res,19:27))./sqrt(nnz(logmat(:,5)&res));
patch([1:9,9:-1:1],[m_resn1a+sem_resn1a,fliplr(m_resn1a-sem_resn1a)],colors(2,:),'FaceAlpha',0.5,'EdgeColor','none')
plot(1:9,mean(nodematpca(logmat(:,1)&res,19:27)),'Color',colors(1,:),'LineWidth',2)
hold on
plot(1:9,mean(nodematpca(logmat(:,5)&res,19:27)),'Color',colors(2,:),'LineWidth',2)
ylim([0.02 0.07])
set(gca,'XTick',1:9,'XTickLabel',num2str(threshfreq','%.0f'),'box','off','FontSize',12)
xlabel('Frequency (Hz)','FontSize',14)
ylabel({'Betweenness';'Centrality'},'FontSize',14)

h1.Position(1)=0.16;
h2.Position(1)=0.16;
h3.Position(1)=0.16;

axes(h1);
text(-0.19,1.2,'A','units','normalized','FontSize',16,'FontWeight','bold')
leg=legend('Engel 1a','>Engel 1a');
set(leg,'box','off')
axes(h2);
text(-0.19,1.2,'B','units','normalized','FontSize',16,'FontWeight','bold')
axes(h3);
text(-0.19,1.2,'C','units','normalized','FontSize',16,'FontWeight','bold')
%%
colors=[0,118,192;163,2,52]/255;
for i=1:9
    subplot(3,3,i)
    nozer=nodematpca(:,36+i);
    nozer(nozer==0)=1e-6;
    [count,edges]=histcounts(log10(nozer(engel1a)),'normalization','probability','BinWidth',0.25);
    plot(edges(1:end-1),count,'LineWidth',1,'color',colors(1,:))
    hold on
    [count,edges]=histcounts(log10(nozer(~engel1a)),'normalization','probability','BinWidth',0.25);
    plot(edges(1:end-1),count,'LineWidth',1,'color',colors(2,:))
    xlim([-7,0])
    if mod(i,3)==1
        ylabel('Fraction')
    end
    set(gca,'box','off','XTick',-6:1:0,'FontSize',8)
    if i>6
        xlabel(texlabel('log10(BC)'),'interpreter','tex')
        set(gca,'XTick',-6:1:0,'XTicklabel',{'-Inf';num2str((-5:1:0)')})
    else
        set(gca,'XTickLabel',[])
    end
    title([num2str(threshfreq(i),'%.0f'),' Hz'])
    %poisstestu(nodematpca(engel1a,36+i),nodematpca(~engel1a,36+i),10000)
    ranksum(nodematpca(engel1a,36+i),nodematpca(~engel1a,36+i))
    if i==1
        [leg,hobj]=legend('Engel 1a','>Engel 1a');
        set(leg,'box','off')
        hobj(3).XData(1)=0.2;
        hobj(5).XData(1)=0.2;
    end
end
suptitle('All Nodes, Post Seizure')
%%
%%
colors=[0,118,192;163,2,52]/255;
for i=1:9
    subplot(3,3,i)
    nozer=nodematbsl(:,i);
    nozer(nozer==0)=1e-6;
    [count,edges]=histcounts(log10(nozer(engelbsl)),'normalization','probability','BinWidth',0.25);
    plot(edges(1:end-1),count,'LineWidth',1,'color',colors(1,:))
    hold on
    [count,edges]=histcounts(log10(nozer(~engelbsl)),'normalization','probability','BinWidth',0.25);
    plot(edges(1:end-1),count,'LineWidth',1,'color',colors(2,:))
    xlim([-7,0])
    if mod(i,3)==1
        ylabel('Fraction')
    end
    set(gca,'box','off','XTick',-6:1:0,'FontSize',8)
    if i>6
        xlabel(texlabel('log10(BC)'),'interpreter','tex')
        set(gca,'XTick',-6:1:0,'XTicklabel',{'-Inf';num2str((-5:1:0)')})
    else
        set(gca,'XTickLabel',[])
    end
    title([num2str(threshfreq(i),'%.0f'),' Hz'])
    %poisstestu(nodematpca(engel1a,36+i),nodematpca(~engel1a,36+i),10000)
    ranksum(nodematbsl(engelbsl,i),nodematpca(~engelbsl,i))
    if i==1
        [leg,hobj]=legend('Engel 1a','>Engel 1a');
        set(leg,'box','off')
        hobj(3).XData(1)=0.2;
        hobj(5).XData(1)=0.2;
    end
end
suptitle('All Nodes, Interictal')
%%
nozer=mean(nodematpca(:,37:45),2);
nozer(nozer==0)=1e-7;
[count,edges]=histcounts(log10(nozer(engel1a)),'normalization','probability','BinWidth',0.5);
plot(edges(1:end-1),count,'LineWidth',2)
hold on
[count,edges]=histcounts(log10(nozer(~engel1a)),'normalization','probability','BinWidth',0.5);
plot(edges(1:end-1),count,'LineWidth',2)
xlim([-8,0])
if mod(i,3)==1
    ylabel('Probability')
end
xlabel(texlabel('log10(BC)'),'interpreter','tex')
set(gca,'box','off','XTick',-7:1:0,'XTicklabel',{'-Inf';num2str((-6:1:0)')})
xlabel(texlabel('log10(BC)'),'interpreter','tex')
ylabel('Probability')
title('Mean of All Frequencies, Post Seizure')

%%
for i=1:9
    subplot(3,3,i)
    nozer=nodematpca(:,18+i);
    nozer(nozer==0)=1e-6;
    [count,edges]=histcounts(log10(nozer(engel1a&res)),'normalization','probability','BinWidth',0.5);
    plot(edges(1:end-1),count,'LineWidth',1,'color',colors(1,:))
    hold on
    [count,edges]=histcounts(log10(nozer(~engel1a&res)),'normalization','probability','BinWidth',0.5);
    plot(edges(1:end-1),count,'LineWidth',1,'color',colors(2,:))
    xlim([-7,0])
    if mod(i,3)==1
        ylabel('Probability')
    end
    set(gca,'box','off','XTick',-6:1:0,'FontSize',8)
    if i>6
        xlabel(texlabel('log10(BC)'),'interpreter','tex')
        set(gca,'XTick',-6:1:0,'XTicklabel',{'-Inf';num2str((-5:1:0)')})
    else
        set(gca,'XTickLabel',[])
    end
    title([num2str(threshfreq(i),'%.0f'),' Hz'])
    if i==1
        [leg,hobj]=legend('Engel 1a','>Engel 1a');
        set(leg,'box','off')
        leg.Position(1:2)=[0.135,0.805];
        hobj(3).XData(1)=0.2;
        hobj(5).XData(1)=0.2;
    end
end
suptitle('Resected Nodes, Mid Seizure')
%%
inds=1:length(eeg2);
ptind=arrayfun(@(x,y) repmat(y,1,length(x.sztype)*length(x.leads)),eeg2,inds,'UniformOutput',0);
ptind=horzcat(ptind{:})';
ptengel1a=arrayfun(@(x) (strcmp(x.engel,'1a')),eeg2);
cps=arrayfun(@(x) reshape(repmat((strcmp(x.sztype,'cps')),1,length(x.leads))',[],1),eeg2,'UniformOutput',0)';
cps=vertcat(cps{:});
figure
for i=1:3
    subplot(2,4,i)
    resgreat=nodematpca(:,22+i)>mean(nodematpca(:,22+i))+3*std(nodematpca(:,22+i));
    resgreat(~res)=false;
    %resgreat(~cps)=false;
    pt_resgreat=accumarray(ptind,resgreat,[],@sum);
    plotSpread(pt_resgreat,'distributionIdx',~ptengel1a,'xNames',{'1a';'>1a'},...
        'distributionColors',{colors(1,:),colors(2,:)})
    h=get(gca,'Children');
    h(1).MarkerSize=7;
    h(2).MarkerSize=7;
    set(gca,'box','off','FontSize',12);
    title([num2str(threshfreq(4+i),'%.0f'),' Hz'])
    if i==1
        ylabel({'Nodes Resected';'Above BC = \mu +3\sigma'})
    end
end

subplot(2,4,4)
resgreat_mult=logical(sum(bsxfun(@ge,nodematpca(:,23:25),mean(nodematpca(:,23:25))+3*std(nodematpca(:,23:25))),2));
resgreat_mult(~res)=false;
%resgreat_mult(~cps)=false;
pt_resgreat_mult=accumarray(ptind,resgreat_mult,[],@sum);
plotSpread(pt_resgreat_mult,'distributionIdx',~ptengel1a,'xNames',{'1a';'>1a'},...
    'distributionColors',{colors(1,:),colors(2,:)})
h=get(gca,'Children');
h(1).MarkerSize=7;
h(2).MarkerSize=7;
set(gca,'box','off','FontSize',12);
title('18-50 Hz')

tot_res=nodematpca(:,24)>-1;
%tot_res=logical(sum(bsxfun(@ge,nodematpca(:,23:25),mean(nodematpca(:,23:25))+3*std(nodematpca(:,23:25))),2));
tot_res(~res)=false;
fracres=pt_resgreat./accumarray(ptind,tot_res,[],@sum);

for i=1:3
    subplot(2,4,4+i)
    resgreat=nodematpca(:,22+i)>mean(nodematpca(:,22+i))+3*std(nodematpca(:,22+i));
    resgreat(~res)=false;
    %resgreat(~cps)=false;
    pt_resgreat=accumarray(ptind,resgreat,[],@sum);
    fracres=pt_resgreat./accumarray(ptind,tot_res,[],@sum);
    plotSpread(fracres,'distributionIdx',~ptengel1a,'xNames',{'1a';'>1a'},...
        'distributionColors',{colors(1,:),colors(2,:)})
    h=get(gca,'Children');
    h(1).MarkerSize=7;
    h(2).MarkerSize=7;
    set(gca,'box','off','FontSize',12);
    if i==1
        ylabel({'Fraction of Resected';'Nodes Above BC = \mu + 3\sigma'})
    end
end

subplot(2,4,8)
fracres_mult=pt_resgreat_mult./accumarray(ptind,tot_res,[],@sum);
plotSpread(fracres_mult,'distributionIdx',~ptengel1a,'xNames',{'1a';'>1a'},...
    'distributionColors',{colors(1,:),colors(2,:)})
h=get(gca,'Children');
h(1).MarkerSize=7;
h(2).MarkerSize=7;
set(gca,'box','off','FontSize',12);
%%
%import temporal only pts_2
name=string(name);
type=string(type);
%onset=string(onset);
%pt_onset=string(pt_onset);
for i=1:length(name)
    pt=find(cell2mat(cellfun(@(x) ~isempty(x),strfind({eeg2.name},lower(name{i})),'UniformOutput',0)));
    if ~isempty(pt)
        sztypes=cellfun(@(x,y) strcat(x,y),eeg2(pt).sztype,eeg2(pt).szlabel,'UniformOutput',0);
        sz=find(cell2mat(cellfun(@(x) ~isempty(x),strfind(lower(sztypes),lower(type{i})),'UniformOutput',0)));
        eeg2(pt).szonset(sz)=onset(i);
        eeg2(pt).ptonset=pt_onset(i);
    end
end
%%
ptind=arrayfun(@(x,y) repmat(y,1,length(x.sztype)*length(x.leads)),eeg2,inds,'UniformOutput',0);
ptind=horzcat(ptind{:})';
szo=cellfun(@(x) reshape(x,[],1),{eeg2.szo}','UniformOutput',0);
szo=vertcat(szo{:});
szo_overlap_mid=szo&logical(sum(bsxfun(@ge,nodematpca(:,23:26),mean(nodematpca(:,23:26))+3*std(nodematpca(:,23:26))),2));
pt_szo_overlap_mid=accumarray(ptind,szo_overlap_mid,[],@sum);

szo_overlap_post=szo&logical(sum(bsxfun(@ge,nodematpca(:,40:42),mean(nodematpca(:,40:42))+3*std(nodematpca(:,40:42))),2));
pt_szo_overlap_post=accumarray(ptind,szo_overlap_post,[],@sum);

subplot(1,2,1)
histogram(pt_szo_overlap_mid)
set(gca,'box','off','FontSize',8,'XTick',0:3)
xlabel('Number of Nodes')
ylabel('Number of Patients')
t=title('Overlap');
set(t,'units','normalized');
t.Position(2)=1;

subplot(1,2,2)
histogram(pt_szo_overlap_post)
set(gca,'box','off','FontSize',8,'XTick',0:3)
xlabel('Number of Nodes')
ylabel('Number of Patients')
t=title('Overlap');
set(t,'units','normalized');
t.Position(2)=1;
%%
leads=arrayfun(@(x) reshape(repmat(x.leads,length(x.sz),1),[],1),eeg2,'UniformOutput',0)';
leads=vertcat(leads{:});
overlap_table=table(leads(szo_overlap),{eeg2(ptind(szo_overlap)).name}','VariableNames',{'Lead','Patient'});

ptonset={eeg2.ptonset}';
empty=cellfun(@isempty,ptonset,'UniformOutput',0);
empty=vertcat(empty{:});
ptonset(empty)={NaN};
ptonset=vertcat(ptonset{:});

szonset=arrayfun(@(x) reshape(repmat(x.szonset,length(x.leads),1),[],1),eeg2,'UniformOutput',0)';
empty_ind=find(empty);
for i=1:length(empty_ind)
    szonset(empty_ind(i))={NaN*ones(length(eeg2(empty_ind(i)).szlabel)*length(eeg2(empty_ind(i)).leads),1)};
end
szonset=vertcat(szonset{:});

bad_nodes_mid=logical(sum(bsxfun(@ge,nodematpca(:,23:26),mean(nodematpca(:,23:26))+3*std(nodematpca(:,23:26))),2));
bad_nodes_mid_pt_onset={eeg2(ptind(find(bad_nodes_mid))).ptonset}';
empty_bad_nodes_mid=cell2mat(cellfun(@isempty,bad_nodes_mid_pt_onset,'UniformOutput',0));

bad_nodes_mid_pt_onset(empty_bad_nodes_mid)={NaN};
bad_nodes_mid_pt_onset=cell2mat(bad_nodes_mid_pt_onset);
onset_table_mid=table({eeg2(ptind(find(bad_nodes_mid))).name}',szonset(find(bad_nodes_mid)),...
    bad_nodes_mid_pt_onset,leads(find(bad_nodes_mid)),'VariableNames',{'Patient','Seizure_Onset','Patient_Onset','Lead'});

bad_nodes_post=logical(sum(bsxfun(@ge,nodematpca(:,40:42),mean(nodematpca(:,40:42))+3*std(nodematpca(:,40:42))),2));
bad_nodes_post_pt_onset={eeg2(ptind(find(bad_nodes_post))).ptonset}';
empty_bad_nodes_post=cell2mat(cellfun(@isempty,bad_nodes_post_pt_onset,'UniformOutput',0));

bad_nodes_post_pt_onset(empty_bad_nodes_post)={NaN};
bad_nodes_post_pt_onset=cell2mat(bad_nodes_post_pt_onset);
onset_table_post=table({eeg2(ptind(find(bad_nodes_post))).name}',szonset(find(bad_nodes_post)),...
    bad_nodes_post_pt_onset,leads(find(bad_nodes_post)),'VariableNames',{'Patient','Seizure_Onset','Patient_Onset','Lead'});
%%
sig=6;
szind=[];
sznum=1;
for i=1:length(eeg2)
    for j=1:length(eeg2(i).sztype)
        szind=[szind;sznum*ones(length(eeg2(i).leads),1)];
        sznum=sznum+1;
    end
end
szengel1a=arrayfun(@(x) repmat(strcmp(x.engel,'1a'),1,length(x.sztype)),eeg2,'UniformOutput',0);
szengel1a=horzcat(szengel1a{:})';

cps=arrayfun(@(x) reshape(repmat((strcmp(x.sztype,'cps')),1,length(x.leads))',[],1),eeg2,'UniformOutput',0)';
cps=vertcat(cps{:});
figure
maxval=5;
for i=1:3
    subplot(2,4,i)
    resgreat=nodematpca(:,22+i)>mean(nodematpca(:,22+i))+sig*std(nodematpca(:,22+i));
    resgreat(~res)=false;
    %resgreat(~cps)=false;
    sz_resgreat=accumarray(szind,resgreat,[],@sum);
    count=zeros(maxval,2);
    [count(:,1),edges]=histcounts(sz_resgreat(szengel1a),0:maxval,'normalization','probability');
    count(:,2)=histcounts(sz_resgreat(~szengel1a),0:maxval,'normalization','probability');
    %b=bar([-0.03:3.97;0.03:4.03]',count,1.1);
    b=bar(edges(1:end-1),count,1.25);
    b(1).FaceColor=colors(1,:);
    b(2).FaceColor=colors(2,:);
    xlim([-0.5,max(edges)])
    set(gca,'box','off','FontSize',12);
    title([num2str(threshfreq(4+i),'%.0f'),' Hz'])
    xlabel('Nodes Resected')
    if i==1
        ylabel('Fraction of Seizures')
    end
end

subplot(2,4,4)
resgreat_mult=logical(sum(bsxfun(@ge,nodematpca(:,23:25),mean(nodematpca(:,23:25))+sig*std(nodematpca(:,23:25))),2));
resgreat_mult(~res)=false;
%resgreat_mult(~cps)=false;
sz_resgreat_mult=accumarray(szind,resgreat_mult,[],@sum);
count=zeros(maxval,2);
[count(:,1),edges]=histcounts(sz_resgreat_mult(szengel1a),0:maxval,'normalization','probability');
count(:,2)=histcounts(sz_resgreat_mult(~szengel1a),0:maxval,'normalization','probability');
%b=bar([-0.03:3.97;0.03:4.03]',count,1.1);
b=bar(edges(1:end-1),count,1.25);
b(1).FaceColor=colors(1,:);
b(2).FaceColor=colors(2,:);
xlim([-0.5,max(edges)])
set(gca,'box','off','FontSize',12);
title('18-50 Hz')
xlabel('Nodes Resected')

tot_res=nodematpca(:,24)>-1;
%tot_res=logical(sum(bsxfun(@ge,nodematpca(:,23:25),mean(nodematpca(:,23:25))+3*std(nodematpca(:,23:25))),2));
tot_res(~res)=false;
fracres=sz_resgreat./accumarray(szind,tot_res,[],@sum);

maxval=0.3;
for i=1:3
    subplot(2,4,4+i)
    resgreat=nodematpca(:,22+i)>mean(nodematpca(:,22+i))+sig*std(nodematpca(:,22+i));
    resgreat(~res)=false;
    %resgreat(~cps)=false;
    sz_resgreat=accumarray(szind,resgreat,[],@sum);
    fracres=sz_resgreat./accumarray(szind,tot_res,[],@sum);
    count=zeros(round(maxval/0.05),2);
    [count(:,1),edges]=histcounts(fracres(szengel1a),0:0.05:maxval,'normalization','probability');
    count(:,2)=histcounts(fracres(~szengel1a),0:0.05:maxval,'normalization','probability');
    %b=bar([-0.0025:0.05:0.1475;0.0025:0.05:0.1525]',count,1.1);
    b=bar(edges(1:end-1),count,1.25);
    b(1).FaceColor=colors(1,:);
    b(2).FaceColor=colors(2,:);
    xlim([-0.04,max(edges)])
    set(gca,'box','off','FontSize',12);
    if i==1
        ylabel('Fraction of Seizures')
    end
    xlabel('Fraction Resected')
end

subplot(2,4,8)
fracres_mult=sz_resgreat_mult./accumarray(szind,tot_res,[],@sum);
count=zeros(round(maxval/0.05),2);
[count(:,1),edges]=histcounts(fracres(szengel1a),0:0.05:maxval,'normalization','probability');
count(:,2)=histcounts(fracres_mult(~szengel1a),0:0.05:maxval,'normalization','probability');
%b=bar([-0.0025:0.05:0.1475;0.0025:0.05:0.1525]',count,1.1);
b=bar(edges(1:end-1),count,1.25);
b(1).FaceColor=colors(1,:);
b(2).FaceColor=colors(2,:);
xlim([-0.04,max(edges)])
set(gca,'box','off','FontSize',12);
xlabel('Fraction Resected')
suptitle(['Nodes Resected','Above BC = \mu +',num2str(sig),'\sigma'])
l1=line([0,1],[0,1],'Color',colors(1,:),'LineWidth',5,'Visible','off');
l2=line([0,1],[0,1],'Color',colors(2,:),'LineWidth',5,'Visible','off');
[leg,legobj]=legend([l1,l2],'Engel 1a','Engel >1a');
set(leg,'box','off','Location','NorthEast')
leg.Position(1)=0.8;
legobj(3).XData=[0.26,0.36];
legobj(5).XData=[0.26,0.36];
%%
figure
subplot 121
totalclass=logical(sum(bsxfun(@ge,nodematpca(:,7:8),mean(nodematpca(:,7:8))+3*std(nodematpca(:,7:8))),2));
resclass=logical(sum(bsxfun(@ge,nodematpca(:,23:25),mean(nodematpca(:,23:25))+3*std(nodematpca(:,23:25))),2));
totalclass=accumarray(ptind,totalclass,[],@mean);
resclass=accumarray(ptind,resclass,[],@mean);
scatter(totalclass(ptengel1a&~isnan(ptonset)'),resclass(ptengel1a&~isnan(ptonset)'),30,'filled')
hold on
scatter(totalclass(~ptengel1a&~isnan(ptonset)'),resclass(~ptengel1a&~isnan(ptonset)'),30,'filled')

xlabel('Fraction of \mu + 3\sigma Pre Nodes')
ylabel('Fraction of Resected \mu + 3\sigma Mid Nodes')
title('Temporal Lobe Epilepsy Patients')
set(gca,'FontSize',12,'box','off')
xlim([0,0.2])
ylim([0,0.14])

subplot 122
totalclass=logical(sum(bsxfun(@ge,nodematpca(:,36:45),mean(nodematpca(:,36:45))+3*std(nodematpca(:,36:45))),2));
resclass=logical(sum(bsxfun(@ge,nodematpca(:,23:25),mean(nodematpca(:,23:25))+3*std(nodematpca(:,23:25))),2));
totalclass=accumarray(ptind,totalclass,[],@mean);
resclass=accumarray(ptind,resclass,[],@mean);
scatter(totalclass(ptengel1a&isnan(ptonset)'),resclass(ptengel1a&isnan(ptonset)'),30,'filled','Marker','o','MarkerFaceColor',colors(1,:))
hold on
scatter(totalclass(~ptengel1a&isnan(ptonset)'),resclass(~ptengel1a&isnan(ptonset)'),30,'filled','Marker','o','MarkerFaceColor',colors(2,:))

xlabel('Fraction of \mu + 3\sigma Post Nodes')
%ylabel('Fraction of Resected \mu + 3\sigma Mid Nodes')
title('Extratemporal Epilepsy Patients')
set(gca,'FontSize',12,'box','off')
xlim([0,0.4])
ylim([0,0.14])

leg=legend('Engel 1a','>Engel 1a');
set(leg,'box','on','location','NorthWest')
%%
f=figure;
%colors=get(gca,'colororder');
colors=[0,118,192;163,2,52;227,124,29]/255;
colorshade=[186,207,236;228,184,180]/255;
ax(1)=subplot(2,2,1);
loc=arrayfun(@(x) repmat(x.loc,length(x.sz),1),eeg2,'UniformOutput',0)';
loc=(vertcat(loc{:}));
for i=1:3
good_mid(i)=nnz(loc(~bad_nodes_mid)==i);
bad_mid(i)=nnz(loc(bad_nodes_mid)==i);
end


[~,chip]=chi2cont([good_mid',bad_mid']');


h=bar([(good_mid+bad_mid)./sum(good_mid+bad_mid);bad_mid./sum(bad_mid)],0.6,'stacked');
h(1).FaceColor=colors(1,:);
h(2).FaceColor=colors(2,:);
h(3).FaceColor=colors(3,:);

set(gca,'XTickLabel',{},'box','off','FontSize',10);
%pos=gca;
%pos.Position(3)=0.6;
l1=line([0,1],[0,1],'color',colors(1,:),'visible','off','LineWidth',5);
l2=line([0,1],[0,1],'color',colors(2,:),'visible','off','LineWidth',5);
l3=line([0,1],[0,1],'color',colors(3,:),'visible','off','LineWidth',5);
[leg,obj]=legend([l1,l2,l3],'Mesial','Lateral','Extra-temporal');
obj(4).XData=[0.2,0.3];
obj(6).XData=[0.2,0.3];
obj(8).XData=[0.2,0.3];
%obj(5).Children.Vertices([1,2,5],1)=0.2;
%obj(6).Children.Vertices([1,2,5],1)=0.2;
set(leg,'box','off','units','normalized','FontSize',8)
leg.Position(1:2)=[0.43,0.75];
xlim([0.5,2.5])
ylabel('Fraction')
text(-0.35,1.05,'a','Units','Normalized','FontSize',16,'FontWeight','bold')
ax(1).Position(3)=0.3;
text(1,-0.08,{'All';'Nodes'},'horizontalAlignment','center','FontSize',10)
text(2,-0.08,{'Extreme';'Nodes'},'horizontalAlignment','center','FontSize',10)
title('Mid-Seizure')

ax(2)=subplot(2,2,2);
odds=zeros(1,2);
se=odds;
err=zeros(2,2);

odds(1)=(bad_mid(1)/good_mid(1))/(bad_mid(3)/good_mid(3));
odds(2)=(bad_mid(2)/good_mid(2))/(bad_mid(3)/good_mid(3));

se(1)=sqrt(1/bad_mid(1)+1/good_mid(1)+1/bad_mid(3)+1/good_mid(3));
se(2)=sqrt(1/bad_mid(2)+1/good_mid(2)+1/bad_mid(3)+1/good_mid(3));

err(1:2,1)=exp(log(odds)-1.96*se);
err(1:2,2)=exp(log(odds)+1.96*se);

odplot1=ploterr(1,odds(1),[],{err(1,1),err(1,2)},'abshhy',0.1);
set(odplot1(1),'LineStyle','none');
set(odplot1(2),'color',colors(1,:),'LineWidth',2);
hold on
odplot2=ploterr(2,odds(2),[],{err(2,1),err(2,2)},'abshhy',0.1);
set(odplot2(1),'LineStyle','none');
set(odplot2(2),'color',colors(2,:),'LineWidth',2);
line(1,odds(1),'LineStyle','none','Marker','o','MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:))
line(2,odds(2),'LineStyle','none','Marker','o','MarkerFaceColor',colors(2,:),'MarkerEdgeColor',colors(2,:))
line([0,3],[1,1],'color','k','LineStyle','--')
ylim([0,4])
ylabel('Odds Ratio')
text(-0.1,1.05,'b','Units','Normalized','FontSize',16,'FontWeight','bold')
set(ax(2),'box','off','YTick',0:4,'XTick',1:2,'XTickLabel',{'Mesial','Lateral'},'YAxisLocation','right','color','none');
ax(2).Position([1,3])=[0.6,0.3];
sig1=text(1,3.2,'*','units','data','FontSize',20,'Fontweight','bold','horizontalAlignment','center');
sig2=text(2,3.2,'*','units','data','FontSize',20,'Fontweight','bold','horizontalAlignment','center');
sigboth=ploterr(1.5,3.7,0.5,[]);
set(sigboth(2),'color','k','LineWidth',2)
%sigboth(2).YData(4)=3.71;
%sigboth(2).YData(7)=3.71;
nonsig=text(1.52,3.9,'N.S.','units','data','FontSize',8,'Fontweight','bold','horizontalAlignment','center');
title('Mid-Seizure')

ax(1)=subplot(2,2,3);
loc=arrayfun(@(x) repmat(x.loc,length(x.sz),1),eeg2,'UniformOutput',0)';
loc=(vertcat(loc{:}));
for i=1:3
good_post(i)=nnz(loc(~bad_nodes_post)==i);
bad_post(i)=nnz(loc(bad_nodes_post)==i);
end


[~,chip]=chi2cont([good_post',bad_post']');


h=bar([(good_post+bad_post)./sum(good_post+bad_post);bad_post./sum(bad_post)],0.6,'stacked');
h(1).FaceColor=colors(1,:);
h(2).FaceColor=colors(2,:);
h(3).FaceColor=colors(3,:);

set(gca,'XTickLabel',{},'box','off','FontSize',10);
%pos=gca;
%pos.Position(3)=0.6;
xlim([0.5,2.5])
ylabel('Fraction')
text(-0.35,1.05,'c','Units','Normalized','FontSize',16,'FontWeight','bold')
ax(1).Position(3)=0.3;
text(1,-0.08,{'All';'Nodes'},'horizontalAlignment','center','FontSize',10)
text(2,-0.08,{'Extreme';'Nodes'},'horizontalAlignment','center','FontSize',10)
title('Post-Seizure')

ax(2)=subplot(2,2,4);
odds=zeros(1,2);
se=odds;
err=zeros(2,2);

odds(1)=(bad_post(1)/good_post(1))/(bad_post(3)/good_post(3));
odds(2)=(bad_post(2)/good_post(2))/(bad_post(3)/good_post(3));

se(1)=sqrt(1/bad_post(1)+1/good_post(1)+1/bad_post(3)+1/good_post(3));
se(2)=sqrt(1/bad_post(2)+1/good_post(2)+1/bad_post(3)+1/good_post(3));

err(1:2,1)=exp(log(odds)-1.96*se);
err(1:2,2)=exp(log(odds)+1.96*se);

odplot1=ploterr(1,odds(1),[],{err(1,1),err(1,2)},'abshhy',0.1);
set(odplot1(1),'LineStyle','none');
set(odplot1(2),'color',colors(1,:),'LineWidth',2);
hold on
odplot2=ploterr(2,odds(2),[],{err(2,1),err(2,2)},'abshhy',0.1);
set(odplot2(1),'LineStyle','none');
set(odplot2(2),'color',colors(2,:),'LineWidth',2);
line(1,odds(1),'LineStyle','none','Marker','o','MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:))
line(2,odds(2),'LineStyle','none','Marker','o','MarkerFaceColor',colors(2,:),'MarkerEdgeColor',colors(2,:))
line([0,3],[1,1],'color','k','LineStyle','--')
ylim([0,4])
ylabel('Odds Ratio')
text(-0.1,1.05,'d','Units','Normalized','FontSize',16,'FontWeight','bold')
set(ax(2),'box','off','YTick',0:4,'XTick',1:2,'XTickLabel',{'Mesial','Lateral'},'YAxisLocation','right','color','none');
ax(2).Position([1,3])=[0.6,0.3];
%sig1=text(1,3.2,'*','units','data','FontSize',20,'Fontweight','bold','horizontalAlignment','center');
sig2=text(2,3.2,'*','units','data','FontSize',20,'Fontweight','bold','horizontalAlignment','center');
sigboth=ploterr(1.5,3.7,0.5,[]);
set(sigboth(2),'color','k','LineWidth',2)
%sigboth(2).YData(4)=3.71;
%sigboth(2).YData(7)=3.71;
nonsig=text(1.52,3.9,'N.S.','units','data','FontSize',8,'Fontweight','bold','horizontalAlignment','center');
title('Post-Seizure')
uistack(leg,'top');
%%
obj(1).Color='k';
obj(2).Color='k';
obj(3).Color='k';
%%
colors=[0,118,192;163,2,52]/255;
colorshade=[186,207,236;228,184,180]/255;
subplot(2,3,1:2)
h1=bar([0:5;0:5]',[histcounts(pt_szo_overlap_mid(ptengel1a),0:6);histcounts(pt_szo_overlap_mid(~ptengel1a),0:6)]');
h1(1).FaceColor=colors(1,:);
h1(2).FaceColor=colors(2,:);
set(gca,'box','off','FontSize',10)
xlabel('Number of Overlapping Nodes')
ylabel('Number of Patients')
l1=line([0,1],[0,1],'color',colors(1,:),'visible','off','LineWidth',5);
l2=line([0,1],[0,1],'color',colors(2,:),'visible','off','LineWidth',5);
[leg,obj]=legend([l1,l2],'Engel 1a','Engel >1a');
obj(3).XData=[0.2,0.3];
obj(5).XData=[0.2,0.3];

set(leg,'box','off','units','normalized','FontSize',8)
text(-0.2,1.05,'a','Units','Normalized','FontSize',16,'FontWeight','bold')
title('Mid-Seizure')
subplot(2,3,3)
h2=bar(1,mean(pt_szo_overlap_mid(ptengel1a)),'FaceColor',colors(1,:));
hold on
h3=bar(2,mean(pt_szo_overlap_mid(~ptengel1a)),'FaceColor',colors(2,:));
e1=ploterr(1,mean(pt_szo_overlap_mid(ptengel1a)),[],std(pt_szo_overlap_mid(ptengel1a))/sqrt(nnz(ptengel1a)),'abshhy',0.1);
e2=ploterr(2,mean(pt_szo_overlap_mid(~ptengel1a)),[],std(pt_szo_overlap_mid(~ptengel1a))/sqrt(nnz(~ptengel1a)),'abshhy',0.1);

e1(1).Color=mean([0,0,0;colors(1,:)]);
e1(2).Color=mean([0,0,0;colors(1,:)]);
e1(1).LineWidth=2;
e1(2).LineWidth=2;

e2(1).Color=mean([0,0,0;colors(2,:)]);
e2(2).Color=mean([0,0,0;colors(2,:)]);
e2(1).LineWidth=2;
set(gca,'box','off','FontSize',10,'Xtick',[],'XTickLabel',{'Engel 1a','Engel >1a'},'YAxisLocation','right')
e2(2).LineWidth=2;

ylabel('Mean Overlap')
text(-0.1,1.05,'b','Units','Normalized','FontSize',16,'FontWeight','bold')
title('Mid-Seizure')

subplot(2,3,4:5)
h1=bar([0:5;0:5]',[histcounts(pt_szo_overlap_post(ptengel1a),0:6);histcounts(pt_szo_overlap_post(~ptengel1a),0:6)]');
h1(1).FaceColor=colors(1,:);
h1(2).FaceColor=colors(2,:);
set(gca,'box','off','FontSize',10)
xlabel('Number of Overlapping Nodes')
ylabel('Number of Patients')
l1=line([0,1],[0,1],'color',colors(1,:),'visible','off','LineWidth',5);
l2=line([0,1],[0,1],'color',colors(2,:),'visible','off','LineWidth',5);
text(-0.2,1.05,'c','Units','Normalized','FontSize',16,'FontWeight','bold')
title('Post-Seizure')
subplot(2,3,6)
h2=bar(1,mean(pt_szo_overlap_post(ptengel1a)),'FaceColor',colors(1,:));
hold on
h3=bar(2,mean(pt_szo_overlap_post(~ptengel1a)),'FaceColor',colors(2,:));
e1=ploterr(1,mean(pt_szo_overlap_post(ptengel1a)),[],std(pt_szo_overlap_post(ptengel1a))/sqrt(nnz(ptengel1a)),'abshhy',0.1);
e2=ploterr(2,mean(pt_szo_overlap_post(~ptengel1a)),[],std(pt_szo_overlap_post(~ptengel1a))/sqrt(nnz(~ptengel1a)),'abshhy',0.1);

e1(1).Color=mean([0,0,0;colors(1,:)]);
e1(2).Color=mean([0,0,0;colors(1,:)]);
e1(1).LineWidth=2;
e1(2).LineWidth=2;

e2(1).Color=mean([0,0,0;colors(2,:)]);
e2(2).Color=mean([0,0,0;colors(2,:)]);
e2(1).LineWidth=2;
set(gca,'box','off','FontSize',10,'Xtick',[],'XTickLabel',{'Engel 1a','Engel >1a'},'YAxisLocation','right')
e2(2).LineWidth=2;

ylabel('Mean Overlap')
text(-0.1,1.05,'d','Units','Normalized','FontSize',16,'FontWeight','bold')
title('Post-Seizure')

uistack(leg,'top');
%%
obj(1).Color='k';
obj(2).Color='k';

%%
p_opt=cell(36,1);
for m=1:36
m
p_opt_sz=zeros(length(eeg(m).sztype),5);
for k=1:length(eeg(m).sztype)
datacube=eeg(m).szraw(k).datacube;
for j=1:size(datacube,3)
[~,A]=arfit(datacube(:,:,j),1,20);
p_opt_sz(k,j)=size(A,2)/size(A,1);
end
end
p_opt{m}=p_opt_sz;
end

p_opt_bsl=cell(36,1);
for m=[1:32,34:36]
m
p_opt_bsl_sz=zeros(1,size(eeg(m).bslraw.datacube,3));
datacube=eeg(m).bslraw.datacube;
for j=1:size(datacube,3)
[~,A]=arfit(datacube(:,:,j),1,20);
p_opt_bsl_sz(j)=size(A,2)/size(A,1);
end
p_opt_bsl{m}=p_opt_bsl_sz;
end

%%

nnz([reshape((vertcat(p_opt{:})),1,[]),horzcat(p_opt_bsl{:})]==3)/length([reshape((vertcat(p_opt{:})),1,[]),horzcat(p_opt_bsl{:})]);

histogram([reshape((vertcat(p_opt{:})),1,[]),horzcat(p_opt_bsl{:})],'normalization','probability')
xlabel('Optimal Model Order')
ylabel('Fraction of Modeled Time Windows')
xlim([0,5])

set(gca,'FontSize',8,'box','off')