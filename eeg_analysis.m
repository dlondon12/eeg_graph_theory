%% standardizes all the seizure types
for i=1:length(eeg)
    for j=1:length(eeg(i).sz)
        if strcmp(eeg(i).sztype(j),'CP');
            eeg(i).sztype(j)={'cps'};
        elseif strcmp(eeg(i).sztype(j),'CPS');
            eeg(i).sztype(j)={'cps'};
        elseif strcmp(eeg(i).sztype(j),'focal');
            eeg(i).sztype(j)={'foc'};
        elseif strcmp(eeg(i).sztype(j),'sp');
            eeg(i).sztype(j)={'sps'};
        elseif strcmp(eeg(i).sztype(j),'subc')
            eeg(i).sztype(j)={'subcl'};
        elseif strcmp(eeg(i).sztype(j),'CPS2gen')
            eeg(i).sztype(j)={'cps2gen'};
        end
    end
end
sztype=unique(vertcat(eeg.sztype));
cmap=hsv(length(sztype));
%%
figure
for i=1:length(eeg)
    %figure
    %cmap=hsv(length(eeg(i).sz));
    for j=1:length(eeg(i).szraw)
        subplot(2,2,find(strcmp(sztype,eeg(i).sztype{j})))
        plot(2:6,eeg(i).szraw(j).glob_mat(:,2),'-','Color',cmap(find(strcmp(sztype,eeg(i).sztype{j})),:),'LineWidth',1);
        hold on
        %ylim([0 0.32])
    end
    xlabel('time');
    ylabel('sz');
    
    %zlabel('bsl');
    %legend(eeg(i).sztype)
    
    %cmap=hsv(length(eeg(i).bsl));
    %for j=1:length(eeg(i).bsl)
    
    %    plot(eeg(i).bsl(j).glob_mat(:,5),'Color',cmap(i,:))
    %    hold on
    
    %end
    
end
%% global metric in full network
engel=arrayfun(@(x) repmat((strcmp(x.engel,'1a')),length(x.sz),1),eeg,'UniformOutput',0)';
engel=vertcat(engel{:});
engel2=arrayfun(@(x) repmat((strcmp(x.engel,'1b')|(strcmp(x.engel,'1c'))|(strcmp(x.engel,'2a'))),length(x.sz),1),eeg,'UniformOutput',0)';
engel2=vertcat(engel2{:});
types=vertcat(eeg.sztype);
[sztypes,ia,ic]=unique(types);
figure
colors=get(gca,'colororder');
titles={'CPS','CPS with Gen.','Subclinical'};
metrics={'\sigma','\omega','Avg. Clust.','Path Length','Synchronizability'};
for m=1:5
    %subplot(2,3,m)
    szmet=[];
    for i=1:length(eeg)
        for j=1:length(eeg(i).szraw)
            szmet(size(szmet,1)+1,:)=eeg(i).szraw(j).glob_mat(:,m);
        end
    end
    %szmet=bsxfun(@times,szmet,1./szmet(:,1));
    
    szmet(isinf(szmet))=NaN;
    mszmet=zeros(length(sztypes),5);
    errszmet=mszmet;
    curcol=0;
    for i=[1,2,4]
        curcol=curcol+1;
        subplot(6,3,3*(m-1)+curcol)
        %boxplot(szmet(ic==i,:),'plotstyle','compact','positions',[(1+i/10):1:(5+i/10)],'colors',cmap(i,:))
        %hold on
        mszmet(i,:)=mean(szmet(ic==i & engel,:),'omitnan');
        errszmet(i,:)=std(szmet(ic==i & engel,:),'omitnan')/sqrt(nnz(ic==i & engel));
        
        errorbar(0.9:1:4.9,mszmet(i,:),...
            errszmet(i,:),errszmet(i,:),'LineWidth',2,'color',mean([0.75,0.75,0.75;colors(curcol,:)]),...
            'Marker','o','MarkerFaceColor',colors(curcol,:),'MarkerSize',3)
        
        hold on
        
        mszmet(i,:)=mean(szmet(ic==i & ~engel,:),'omitnan');
        errszmet(i,:)=std(szmet(ic==i & ~engel,:),'omitnan')/sqrt(nnz(ic==i & ~engel));
        
        errorbar(0.9:1:4.9,mszmet(i,:),...
            errszmet(i,:),errszmet(i,:),'LineWidth',2,'color',colors(curcol,:),...
            'Marker','o','MarkerFaceColor',colors(curcol,:),'MarkerSize',3,'LineStyle',':')
        
        
        %mszmet(i,:)=mean(szmet(ic==i & engel2,:),'omitnan');
        %errszmet(i,:)=std(szmet(ic==i & engel2,:),'omitnan')/sqrt(nnz(ic==i & engel2));
        
        %errorbar(1:5,mszmet(i,:),...
        %    errszmet(i,:),errszmet(i,:),'LineWidth',2,'color',colors(curcol,:),...
        %    'Marker','o','MarkerFaceColor',colors(curcol,:),'MarkerSize',3,'LineStyle','--')
        
        %mszmet(i,:)=mean(szmet(ic==i & ~(engel|engel2),:),'omitnan');
        %errszmet(i,:)=std(szmet(ic==i & ~(engel|engel2),:),'omitnan')/sqrt(nnz(ic==i & ~(engel|engel2)));
        %hold on
        
        %errorbar(1.1:1:5.1,mszmet(i,:),...
        %    errszmet(i,:),errszmet(i,:),'LineWidth',2,'color',mean([0,0,0;colors(curcol,:)]),...
        %    'Marker','o','MarkerFaceColor',colors(curcol,:),'MarkerSize',3,'LineStyle',':')
        set(gca,'XTick',1:5,'XTickLabel',{},'box','off','FontSize',10)
        xlim([0.7 5.3])
        if i==1
            ylab=ylabel(metrics(m));
            set(ylab,'units','normalized','position',[-0.3, 0.5, 0])
        else
            set(gca,'YTickLabel',{})
        end
        if m==1
            %ylim([1.5 2.5])
            title(titles(curcol),'FontSize',14)
        elseif m==2
            %ylim([-0.2 0.5])
        elseif m==3
            %ylim([1.5 3.2])
        elseif m==4
            %ylim([1 1.4])
        elseif m==5
            %ylim([0.03 0.15])
            set(gca,'box','off','XTick',1:5,'XTickLabel',{'Pre','Early','Mid','Late','Post'})
            
        end
    end
end
subplot(6,3,17)
line([0 1],[0 1],'color',[0.5 0.5 0.5],'LineWidth',2,'Visible','off')
hold on
%line([0 1],[0 1],'color',[0.25 0.25 0.25],'LineStyle','--','LineWidth',2,'Visible','off')
%hold on
line([0 1],[0 1],'color','k','LineStyle',':','LineWidth',2,'Visible','off')
set(gca,'Visible','off')
leg=legend('Engel 1a','Engel >1a');%, 'Engel 2b-4b'); %,...
%'FontSize',12,'position',[0.73 0.2 0.1375 0.0833]);
set(leg,'box','off','FontSize',12,'Location','South','Position',[0.45,0.1,0.15,0.1])
supt=suptitle('Across Seizures: Global Metrics in Full Network');
set(supt,'position',[0.5,-0.07,0],'FontSize',16,'FontWeight','bold')
%% global metrics in unresected network
types=vertcat(eeg.sztype);
[sztypes,ia,ic]=unique(types);
figure
colors=get(gca,'colororder');
titles={'CPS','CPS with Gen.','Subclinical'};
metrics={'SMI','Mod SMI','Avg. Clust.','Path Length','Synchronizability'};
for m=1:5
    %subplot(2,3,m)
    szmet_unres=[];
    for i=1:length(eeg)
        for j=1:length(eeg(i).sz)
            szmet_unres(size(szmet_unres,1)+1,:)=eeg(i).sz(j).glob_mat_unres(:,m);
        end
    end
    %szmet=bsxfun(@times,szmet,1./szmet(:,1));
    
    mszmet=zeros(length(sztypes),5);
    errszmet=mszmet;
    curcol=0;
    for i=[1,2,4]
        curcol=curcol+1;
        subplot(6,3,3*(m-1)+curcol)
        %boxplot(szmet(ic==i,:),'plotstyle','compact','positions',[(1+i/10):1:(5+i/10)],'colors',cmap(i,:))
        %hold on
        mszmet(i,:)=mean(szmet_unres(ic==i & engel,:));
        errszmet(i,:)=std(szmet_unres(ic==i & engel,:))/sqrt(nnz(ic==i & engel));
        
        errorbar([(0.8+curcol/10):1:(4.8+curcol/10)],mszmet(i,:),...
            errszmet(i,:),errszmet(i,:),'LineWidth',2,'color',colors(curcol,:),...
            'Marker','o','MarkerFaceColor',colors(curcol,:),'MarkerSize',3)
        
        mszmet(i,:)=mean(szmet_unres(ic==i & ~engel,:));
        errszmet(i,:)=std(szmet_unres(ic==i & ~engel,:))/sqrt(nnz(ic==i & ~engel));
        hold on
        
        errorbar([(0.85+curcol/10):1:(4.85+curcol/10)],mszmet(i,:),...
            errszmet(i,:),errszmet(i,:),'LineWidth',2,'color',colors(curcol,:),...
            'Marker','o','MarkerFaceColor',colors(curcol,:),'MarkerSize',3,'LineStyle','-.')
        set(gca,'XTick',1:5,'XTickLabel',{},'box','off','FontSize',12)
        xlim([0.7 5.3])
        if i==1
            ylab=ylabel(metrics(m));
            set(ylab,'units','normalized','position',[-0.3, 0.5, 0])
        else
            set(gca,'YTickLabel',{})
        end
        if m==1
            ylim([1.5 2.5])
            title(titles(curcol),'FontSize',14)
        elseif m==2
            ylim([-0.5 0])
        elseif m==3
            %ylim([0.3 0.6])
        elseif m==4
            ylim([1.8 3])
        elseif m==5
            ylim([0.0 0.15])
            set(gca,'box','off','XTick',1:5,'XTickLabel',{'Pre','Early','Mid','Late','Post'})
            %leg=legend(titles,'units','normalized','box','off',...
            %   'FontSize',12,'position',[0.73 0.2 0.1375 0.0833]);
        end
    end
end
subplot(6,3,17)
line([0 1],[0 1],'color','k','LineWidth',2,'Visible','off')
hold on
line([0 1],[0 1],'color','k','LineStyle','-.','LineWidth',2,'Visible','off')
set(gca,'Visible','off')
leg=legend('Engel 1','Engel>1');
set(leg,'box','off','FontSize',12,'Location','South','Position',[0.45,0.05,0.15,0.1])
supt=suptitle('Global Metrics in Unresected Network');
set(supt,'position',[0.5,-0.07,0])
%% averages with individual seizure lines
types=vertcat(eeg.sztype);
[sztypes,ia,sztypeind]=unique(types);
metrics={'SMI','Mod SMI','Avg. Clust.','Path Length','Synchronizability'};
titles={'CPS','CPS with Gen.','Subclinical'};
figure
colors=get(gca,'colororder');
for m=1:5
    szmet=[];
    for i=1:length(eeg)
        for j=1:length(eeg(i).sz)
            szmet(size(szmet,1)+1,:)=eeg(i).szraw(j).glob_mat(:,m);
        end
    end
    %szmet=bsxfun(@times,szmet,1./szmet(:,1));
    
    mszmet=zeros(length(sztypes),5);
    errszmet=mszmet;
    curcol=0;
    for i=[1,2,4]
        curcol=curcol+1;
        subplot(5,4,4*(m-1)+curcol)
        plot(1:5,szmet(ic==i,:)','color',mean([colors(curcol,:);1,1,1]),'LineWidth',1)
        %hold on
        %plot(1:5,szmet(ic==i & ~engel,:)','color',mean([colors(curcol,:);0,0,0]),'LineWidth',1)
        hold on
        %boxplot(szmet(ic==i,:),'plotstyle','compact','positions',[(1+i/10):1:(5+i/10)],'colors',cmap(i,:))
        %hold on
        mszmet(i,:)=mean(szmet(ic==i,:));
        errszmet(i,:)=std(szmet(ic==i,:))/sqrt(nnz(ic==i));
        
        errorbar([(1+i/20):1:(5+i/20)],mszmet(i,:),errszmet(i,:),errszmet(i,:),...
            'LineWidth',2,'color',colors(curcol,:),'Marker','o',...
            'MarkerFaceColor',colors(curcol,:),'MarkerSize',3)
        xlim([0.8 5.2])
        set(gca,'XTick',1:5,'XTickLabel',{},'box','off')
        
        if m==1
            title(titles(curcol))
            %ylim([1 3])
        elseif m==2
            %ylim([-1 1])
        elseif m==3
            %ylim([0 1])
        elseif m==4
            %ylim([1.5 3.5])
        elseif m==5
            %ylim([0 0.3])
            set(gca,'box','off','XTick',1:5,'XTickLabel',{'Pre','Early','Mid','Late','Post'})
        end
        if i==1
            set(gca,'Units','normalized')
            ylab(m)=ylabel(metrics(m));
        else
            set(gca,'YTickLabel',{})
        end
    end
    subplot(5,4,4*(m-1)+4)
    for i=[1,2,4]
        errorbar([(0.8+i/10):1:(4.8+i/10)],mszmet(i,:),errszmet(i,:),...
            errszmet(i,:),'LineWidth',2,'Marker','o','MarkerFaceColor',...
            colors(curcol,:),'MarkerSize',3)
        hold on
        set(gca,'XTick',1:5,'XTickLabel',{})
    end
    xlim([0.7 5.3])
    set(gca,'box','off')
    if m==5
        set(gca,'box','off','XTick',1:5,'XTickLabel',{'Pre','Early','Mid','Late','Post'})
    end
end
tightfig;
%% individual patient lines broken down by Engel class and seizure type
%% across patients global metrics in full network
type={'cps','cps2gen','subcl'};
titles={'CPS','CPS with Gen.','Subclinical'};
metrics={'\sigma','\omega','Average\newlineClust.','Path\newlineLength','S'};
p=zeros(5,5,3);
p2=p;
engel=arrayfun(@(x) strcmp(x.engel,'1a'),eeg,'UniformOutput',0)';
engel=vertcat(engel{:});
engel2=arrayfun(@(x) (strcmp(x.engel,'1b')|strcmp(x.engel,'1c')|strcmp(x.engel,'2a')),eeg,'UniformOutput',0)';
engel2=vertcat(engel2{:});
avgdistmat=[];
ptengelmat=[];
timepoint=[];
sztypemat=[];
ptmat=[];
figure
for k=1:3
    subclind=arrayfun(@(x) strcmp(x.sztype,type(k)),eeg,'UniformOutput',0);
    subclind=vertcat(subclind{:});
    ptsubclind=arrayfun(@(x,y) repmat(y,1,nnz(strcmp(x.sztype,type(k)))),eeg,1:length(eeg),'UniformOutput',0);
    ptsubclind=horzcat(ptsubclind{:})';
    ptengel=arrayfun(@(x) (strcmp(x.engel,'1a')),eeg);
    %colors=zeros(length(unique(ptsubclind)),3);
    %colors(:,3)=ptengel(unique(ptsubclind));
    %colors(:,1)=~ptengel(unique(ptsubclind));
    
    %p=zeros(5,5);
    %p2=p;
    avgdistmattype=[];
    ptengelmattype=[];
    timepoint=vertcat(timepoint,reshape(repmat(1:6,length(unique(ptsubclind)),1),[],1));
    sztypemat=vertcat(sztypemat,k*ones(6*length(unique(ptsubclind)),1));
    ptmat=vertcat(ptmat,repmat(unique(ptsubclind),6,1));
    for m=1:5
        %figure
        subplot(6,3,(m-1)*3+k)
        colors=get(gca,'colororder');
        szmet=[];
        for i=1:length(eeg)
            for j=1:length(eeg(i).sz)
                szmet(size(szmet,1)+1,:)=eeg(i).szraw(j).glob_mat(:,m);
            end
        end
        
        bslmet=[];
        for i=1:length(eeg)
            %glob_mat=cell2mat(eeg(i).bslraw.glob_mat(:,m));
            if size(eeg(i).bslraw.glob_mat{:})>1
                bslmet(size(bslmet,1)+1,:)=mean(eeg(i).bslraw.glob_mat{:}(:,m));
            else
                bslmet(size(bslmet,1)+1,:)=eeg(i).bslraw.glob_mat{:}(m);
            end
            
        end
        
        avgdist=NaN*ones(36,6);
        for i=1:5
            avgdist(:,i+1)=accumarray(ptsubclind,szmet(subclind,i),[36 1],@mean,NaN);
            %subplot(2,3,i)
            %plotSpread({avgdist(ptengel),avgdist(~ptengel)});
        end
        avgdist(:,1)=mean(bslmet,2);
        ptengel2=ptengel;
        ptengel2(setxor(1:36,unique(ptsubclind)))=[];
        avgdist(setxor(1:36,unique(ptsubclind)),:)=[];
        
        %ptengel2(isnan(avgdist(:,2)))=[];
        %avgdist(isnan(avgdist(:,2)),:)=[];
        avgdistmattype=horzcat(avgdistmattype,reshape(avgdist,[],1));
        ptengelmattype=repmat(ptengel2,1,6)';
        %ptengelmattype=horzcat(ptengelmattype,reshape(ptengel2,[],1));
        %plot(avgdist(ptengel2,:)','color',mean([1,1,1;colors(k,:)]),'LineWidth',1)
        %hold on
        %plot(avgdist(~ptengel2,:)','color',mean([0,0,0;colors(k,:)]),'LineWidth',1)
        %plot(mean(avgdist(ptengel2,:)),'color',colors(1,:),'LineWidth',3)
        %plot(mean(avgdist(~ptengel2,:)),'color',colors(2,:),'LineWidth',3)
        
        errorbar(0.95:1:5.95,mean(avgdist(ptengel2,:),'omitnan'),std(avgdist(ptengel2,:),'omitnan')/sqrt(nnz(ptengel2)),...
            std(avgdist(ptengel2,:),'omitnan')/sqrt(nnz(ptengel2)),'LineWidth',2,'color',colors(k,:),...
            'Marker','o','MarkerFaceColor',colors(k,:),'MarkerSize',3)
        hold on
        
        errorbar(1.05:1:6.05,mean(avgdist(~ptengel2,:),'omitnan'),std(avgdist(~ptengel2,:),'omitnan')/sqrt(nnz(~ptengel2)),...
            std(avgdist(~ptengel2,:),'omitnan')/sqrt(nnz(~ptengel2)),'LineWidth',2,'color',mean([0,0,0;colors(k,:)]),...
            'Marker','o','MarkerFaceColor',mean([0,0,0;colors(k,:)]),'MarkerSize',3,...
            'LineStyle',':')
        for i=1:6
            p(m,i,k)=poisstestu(avgdist(ptengel2,i),avgdist(~ptengel2,i),1000);
            p2(m,i,k)=ranksum(avgdist(ptengel2,i),avgdist(~ptengel2,i));
        end
        %(mean(avgdist(ptengel2,:))-mean(avgdist(~ptengel2,:)))./sqrt(var(avgdist(ptengel2,:))/nnz(ptengel2) + var(avgdist(~ptengel2,:))/nnz(~ptengel2))
        xlim([0.7 6.3])
        set(gca,'box','off','XTick',1:6,'XTickLabel',[],'FontSize',10)
        %h=breakxaxis([1.2 1.8],0.005);
        
        if k==1
            ylab=ylabel(metrics(m),'FontSize',14,'FontWeight','bold');
            set(ylab,'units','normalized','position',[-0.3, 0.5, 0])
        else
            %set(gca,'YTickLabel',{})
        end
        if m==1
            title(titles(k),'FontWeight','bold')
            %lim([1.5 2.5])
        elseif m==2
            %ylim([-0.5 0.5])
        elseif m==3
            ylim([1.5 3.5])
        elseif m==4
            %ylim([1 1.5])
        elseif m==5
            %ylim([5 40])
            set(gca,'XTick',1:6,'XTickLabel',{'Int','Pre','Early','Mid','Late','Post'},'XTickLabelRotation',90)
        end
        
        ymin=get(gca,'YLim');
        line([1.2 1.8],[ymin(1),ymin(1)],'LineStyle',':','color','w','LineWidth',2)
        
        %if pANCOVAN(1,m,k)<0.05/3 && pANCOVAN(1,m,k)>=0.01/3
        %    text(0.48,0.9,'*','Units','Normalized','FontSize',12,'FontWeight','bold')
        %elseif pANCOVAN(1,m,k)<0.01/3 && pANCOVAN(1,m,k)>=0.001/3
        %    text(0.46,0.9,'**','Units','Normalized','FontSize',12,'FontWeight','bold')
        %elseif pANCOVAN(1,m,k)<0.001/3
        %    text(0.44,0.9,'***','Units','Normalized','FontSize',12,'FontWeight','bold')
        %end
    end
    avgdistmat=vertcat(avgdistmat,avgdistmattype);
    ptengelmat=vertcat(ptengelmat,ptengelmattype);
end
subplot(6,3,17)
line([0 1],[0 1],'color',[0.5 0.5 0.5],'LineWidth',2,'Visible','off')
hold on
%line([0 1],[0 1],'color',[0.25 0.25 0.25],'LineStyle','--','LineWidth',2,'Visible','off')
%hold on
line([0 1],[0 1],'color','k','LineStyle',':','LineWidth',2,'Visible','off')
set(gca,'Visible','off')
leg=legend('Engel 1a','Engel >1a'); %,...
%'FontSize',12,'position',[0.73 0.2 0.1375 0.0833]);
set(leg,'box','off','FontSize',12,'Location','North','Position',[0.45,0.05,0.15,0.1])
supt=suptitle('Across Patients: Global Metrics in Full Network');
set(supt,'position',[0.5,-0.07,0],'FontSize',16,'FontWeight','bold')
%%
%prep=zeros(6,5,3,10);
%for a=1:10
    %a
h=zeros(6,5,3);
p=h;
clusts=cell(5,3);
for k=1:3
    k
    for m=1:5
        m
        clear x1
        clear x2
        for i=1:6
            x1(:,i)=avgdistmat(sztypemat==k & timepoint==i & ptengelmat==0,m);
            x2(:,i)=avgdistmat(sztypemat==k & timepoint==i & ptengelmat==1,m);
        end
        
        [h(:,m,k),p(:,m,k),clusts{m,k}]=clust_mass_1d_v2(x1,x2,0.2,0.05/15,1000);
    end
end
%prep(:,:,:,a)=p;
%end
%% MANOVAN
%[T,p,FANCOVAN,pANCOVAN,stats]=mancovan(avgdistmat,[sztypemat,ptengelmat],timepoint,{'verbose','group-covariate','group-group'});

[~,~,~,pANCOVAN(:,:,1)]=mancovan(avgdistmat(sztypemat==1,:),[ptengelmat(sztypemat==1),timepoint(sztypemat==1)],[],{'verbose'});
[~,~,~,pANCOVAN(:,:,2)]=mancovan(avgdistmat(sztypemat==2,:),[ptengelmat(sztypemat==2),timepoint(sztypemat==2)],[],{'verbose'});
[~,~,~,pANCOVAN(:,:,3)]=mancovan(avgdistmat(sztypemat==3,:),[ptengelmat(sztypemat==3),timepoint(sztypemat==3)],[],{'verbose'});
%%
p=zeros(3,5,6);
for m=1:5
for i=1:6
p(1,m,i)=poisstestu(squeeze(avgdiststruct(1,m,:,i)),squeeze(avgdiststruct(2,m,:,i)),10000);
p(2,m,i)=poisstestu(squeeze(avgdiststruct(1,m,:,i)),squeeze(avgdiststruct(3,m,:,i)),10000);
p(3,m,i)=poisstestu(squeeze(avgdiststruct(2,m,:,i)),squeeze(avgdiststruct(3,m,:,i)),10000);
end
end
%% across patients global metrics in full network
type={'cps','cps2gen','subcl'};
titles={'CPS','CPS with Gen.','Subclinical'};
metrics={'\sigma','\omega','Average Clust.','Path Length','S'};
%p=zeros(5,6,3);
%p2=p;

engel=arrayfun(@(x) strcmp(x.engel,'1a'),eeg,'UniformOutput',0)';
engel=vertcat(engel{:});
engel2=arrayfun(@(x) (strcmp(x.engel,'1b')|strcmp(x.engel,'1c')|strcmp(x.engel,'2a')),eeg,'UniformOutput',0)';
engel2=vertcat(engel2{:});
avgdistmat=[];
ptengelmat=[];
timepoint=[];
sztypemat=[];
ptmat=[];

avgdiststruct=zeros(3,5,36,6);

figure
for k=1:3
    subclind=arrayfun(@(x) strcmp(x.sztype,type(k)),eeg,'UniformOutput',0);
    subclind=vertcat(subclind{:});
    ptsubclind=arrayfun(@(x,y) repmat(y,1,nnz(strcmp(x.sztype,type(k)))),eeg,1:length(eeg),'UniformOutput',0);
    ptsubclind=horzcat(ptsubclind{:})';
    ptengel=arrayfun(@(x) (strcmp(x.engel,'1a')),eeg);
   

    avgdistmattype=[];
    ptengelmattype=[];
    timepoint=vertcat(timepoint,reshape(repmat(1:6,length(unique(ptsubclind)),1),[],1));
    sztypemat=vertcat(sztypemat,k*ones(6*length(unique(ptsubclind)),1));
    ptmat=vertcat(ptmat,repmat(unique(ptsubclind),6,1));
    for m=1:5
        h(m)=subplot(2,3,m);
        colors=get(gca,'colororder');
        szmet=[];
        for i=1:length(eeg)
            for j=1:length(eeg(i).sz)
                szmet(size(szmet,1)+1,:)=eeg(i).szraw(j).glob_mat(:,m);
            end
        end
        
        bslmet=[];
        for i=1:length(eeg)
            if size(eeg(i).bslraw.glob_mat{:})>1
                bslmet(size(bslmet,1)+1,:)=mean(eeg(i).bslraw.glob_mat{:}(:,m));
            else
                bslmet(size(bslmet,1)+1,:)=eeg(i).bslraw.glob_mat{:}(m);
            end
            
        end
        
        avgdist=NaN*ones(36,6);
        for i=1:5
            avgdist(:,i+1)=accumarray(ptsubclind,szmet(subclind,i),[36 1],@mean,NaN);
        end
        avgdist(:,1)=mean(bslmet,2);
        
        avgdiststruct(k,m,:,:)=avgdist;
        
        ptengel2=ptengel;
        ptengel2(setxor(1:36,unique(ptsubclind)))=[];
        avgdist(setxor(1:36,unique(ptsubclind)),:)=[];
        
        
        avgdistmattype=horzcat(avgdistmattype,reshape(avgdist,[],1));
        ptengelmattype=repmat(ptengel2,1,6)';
        
        errorbar((0.95:1:5.95)+0.15*k,mean(avgdist,'omitnan'),std(avgdist,'omitnan')./sqrt(sum(~isnan(avgdist))),...
            std(avgdist,'omitnan')./sqrt(sum(~isnan(avgdist))),'LineWidth',2,'color',colors(k,:),...
            'Marker','o','MarkerFaceColor',colors(k,:),'MarkerSize',3)
        hold on
        
        xlim([0.7 6.7])
        set(gca,'box','off','XTick',1:6,'XTickLabel',[],'FontSize',10)
      
        
        set(gca,'XTick',1:6,'XTickLabel',{'Int','Pre','Early','Mid','Late','Post'},'XTickLabelRotation',90)

        ymin=get(gca,'YLim');
        line([1.2 1.8],[ymin(1),ymin(1)],'LineStyle',':','color','w','LineWidth',2)
        
        %if pANCOVAN(1,m,k)<0.05/3 && pANCOVAN(1,m,k)>=0.01/3
        %    text(0.48,0.9,'*','Units','Normalized','FontSize',12,'FontWeight','bold')
        %elseif pANCOVAN(1,m,k)<0.01/3 && pANCOVAN(1,m,k)>=0.001/3
        %    text(0.46,0.9,'**','Units','Normalized','FontSize',12,'FontWeight','bold')
        %elseif pANCOVAN(1,m,k)<0.001/3
        %    text(0.44,0.9,'***','Units','Normalized','FontSize',12,'FontWeight','bold')
        %end
        
        color_comb=zeros(3,3);
        color_comb(1,:)=[1,0,1];
        color_comb(2,:)=[0,0,0];
        color_comb(3,:)=[1,0.5,0];
        
        
        
        if k==3
            if m<3
                ylab=ylabel(metrics(m),'FontSize',16,'FontWeight','bold');
            else
                ylab=ylabel(metrics(m),'FontSize',14,'FontWeight','bold');
            end
            set(ylab,'units','normalized','position',[-0.15, 0.5, 0])
        end
    end
    avgdistmat=vertcat(avgdistmat,avgdistmattype);
    ptengelmat=vertcat(ptengelmat,ptengelmattype);
    
end

pmarks={'o','s','^'};

for k=1:3
    for m=1:5
        ymin=get(h(m),'YLim');
        for i=1:6
            if p(k,m,i)<0.05/18
                plot(h(m),i+(k-1)*0.35,ymin(1)+1*(ymin(2)-ymin(1)),'Marker',pmarks{k},'MarkerFaceColor',color_comb(k,:),'MarkerEdgeColor','none');
            end
        end
    end
end

subplot(2,3,6)
line([0 1],[0 1],'color',colors(1,:),'LineWidth',2,'Visible','on')
line([0 1],[0 1],'color',colors(2,:),'LineWidth',2,'Visible','on')
line([0 1],[0 1],'color',colors(3,:),'LineWidth',2,'Visible','on')
line([0 1],[0 1],'color',color_comb(1,:),'LineWidth',2,'LineStyle','none','Marker',pmarks{1},'MarkerFaceColor',color_comb(1,:),'MarkerEdgeColor','none','Visible','on')
line([0 1],[0 1],'color',color_comb(2,:),'LineWidth',2,'LineStyle','none','Marker',pmarks{2},'MarkerFaceColor',color_comb(2,:),'MarkerEdgeColor','none','Visible','on')
line([0 1],[0 1],'color',color_comb(3,:),'LineWidth',2,'LineStyle','none','Marker',pmarks{3},'MarkerFaceColor',color_comb(3,:),'MarkerEdgeColor','none','Visible','on')
ylim([1000,1001])
%hold on
%line([0 1],[0 1],'color',[0.25 0.25 0.25],'LineStyle','--','LineWidth',2,'Visible','off')
%hold on
%line([0 1],[0 1],'color','k','LineStyle',':','LineWidth',2,'Visible','off')
set(gca,'Visible','off')
leg=legend('CPS','CPS w/gen','Subclinical','CPS vs CPS w/ gen','CPS vs Subclinical','CPS w/ gen vs Subclinical'); %,...
%'FontSize',12,'position',[0.73 0.2 0.1375 0.0833]);
set(leg,'box','off','FontSize',12,'Position',[0.7,0.1,0.25,0.3])
supt=suptitle('Metrics by seizure type and timepoint');
set(supt,'position',[0.5,-0.07,0],'FontSize',16,'FontWeight','bold')

%% all seizures bar graphs for all metrics
p=zeros(5,10);
p2=p;
for m=1:5
    szmet=[];
    for i=1:length(eeg)
        for j=1:length(eeg(i).sz)
            szmet(size(szmet,1)+1,:)=eeg(i).sz(j).glob_mat(:,m);
        end
    end
    szmet_unres=[];
    for i=1:length(eeg)
        for j=1:length(eeg(i).sz)
            szmet_unres(size(szmet_unres,1)+1,:)=eeg(i).sz(j).glob_mat_unres(:,m);
        end
    end
    
    szmet_rand=[];
    if m==5
        for i=1:length(eeg)
            for j=1:length(eeg(i).sz)
                szmet_rand(size(szmet_rand,1)+1,:)=mean(real(eeg(i).sz(j).glob_mat_rand),3); %wtf
            end
        end
    else
        szmet_rand=zeros(size(szmet_unres));
    end
    
    subplot(2,3,m)
    %boxplot([szmet(engel,:),szmet_unres(engel,:),szmet_rand(engel,:)],'positions',[1:5:21,2.5:5:22.5,4:5:24],'plotstyle','compact')
    %hold on;
    %boxplot([szmet(~engel,:),szmet_unres(~engel,:),szmet_rand(~engel,:)],'positions',[1.5:5:21.5,3:5:23,4.5:5:24.5],'plotstyle','compact','colors','r')
    %xlim([0 16])
    
    engm=mean(szmet(engel,:));
    engerr=std(szmet(engel,:))/sqrt(sum(engel));
    
    engm_unres=mean(szmet_unres(engel,:));
    engerr_unres=std(szmet_unres(engel,:))/sqrt(sum(engel));
    
    engm_rand=mean(szmet_rand(engel,:));
    engerr_rand=std(szmet_rand(engel,:))/sqrt(sum(engel));
    
    nengm=mean(szmet(~engel,:));
    nengerr=std(szmet(~engel,:))/sqrt(sum(~engel));
    
    nengm_unres=mean(szmet_unres(~engel,:));
    nengerr_unres=std(szmet_unres(~engel,:))/sqrt(sum(~engel));
    
    nengm_rand=mean(szmet_rand(~engel,:));
    nengerr_rand=std(szmet_rand(~engel,:))/sqrt(sum(~engel));
    
    errorbar([1:5:21,2.5:5:22.5,4:5:24],[engm,engm_unres,engm_rand],[engerr,engerr_unres,engerr_rand],'.b');
    hold on;
    engbar=bar([1:5:21,2.5:5:22.5,4:5:24],diag([engm,engm_unres,engm_rand]),0.3,'stacked','b');
    for j=1:5
        set(engbar(j+5),'FaceColor',[0.3 0.3 1])
        set(engbar(j+10),'FaceColor',[0.6 0.6 1])
    end
    errorbar([1.5:5:21.5,3:5:23,4.5:5:24.5],[nengm,nengm_unres,nengm_rand],[nengerr,nengerr_unres,nengerr_rand],'.r')
    nengbar=bar([1.5:5:21.5,3:5:23,4.5:5:24.5],diag([nengm,nengm_unres,nengm_rand]),0.3,'stacked','r');
    for j=1:5
        set(nengbar(j+5),'FaceColor',[1 0.3 0.3])
        set(nengbar(j+10),'FaceColor',[1 0.6 0.6])
    end
    if m==1
        title('Small World','FontSize',14)
    elseif m==2
        title('Small World','FontSize',14)
        ylim([-0.4 0.6])
    elseif m==3
        title('Average Clustering Coefficient','FontSize',14)
    elseif m==4
        ylim([0 3])
        title('Average Clustering Coefficient','FontSize',14)
    elseif m==5
        title('Synchronizability','FontSize',14)
        ylim([0 0.16])
        leg=legend([engbar(1),engbar(6),nengbar(1),nengbar(6),engbar(11),nengbar(11)],...
            'Engel 1','Engel 1 Unresected','Engel >1','Engel >1 Unresected','Engel 1 Random Resection','Engel >1 Random Resection');
        %legend 'boxoff'
    end
    
    ymax=get(gca,'YLim');
    ymin=ymax(1);
    ymax=ymax(2)-ymax(1);
    set(gca,'XTick',2:5:22,'XTickLabel',{'Pre';'Early';'Middle';'Late';'Post'},'box','off')
    line([1:5:21;1.5:5:21.5],(0.9*ymax+ymin)*ones(2,5),'color','k')
    line([2.5:5:22.5;3:5:23],(0.9*ymax+ymin)*ones(2,5),'color','k')
    
    
    for i=1:5
        p(m,i)=poisstestu(szmet(engel,i),szmet(~engel,i),10000);
        if p(m,i)<0.01 && p(m,i)>=0.001
            text(1.25+5*(i-1),0.93*ymax+ymin,'*','horizontalAlignment','center',...
                'FontSize',20,'FontWeight','bold')
        elseif p(m,i)<0.001
            text(1.25+5*(i-1),0.93*ymax+ymin,'**','horizontalAlignment','center',...
                'FontSize',20,'FontWeight','bold')
        end
    end
    for i=1:5
        p(m,i+5)=poisstestu(szmet_unres(engel,i),szmet_unres(~engel,i),10000);
        if p(m,i+5)<0.01 && p(m,i+5)>=0.001
            text(2.75+5*(i-1),0.93*ymax+ymin,'*','horizontalAlignment','center',...
                'FontSize',20,'FontWeight','bold')
        elseif p(m,i+5)<0.001
            text(2.75+5*(i-1),0.93*ymax+ymin,'**','horizontalAlignment','center',...
                'FontSize',20,'FontWeight','bold')
        end
    end
    
    for i=1:5
        p2(m,i)=poisstest2(szmet(engel,i),szmet_unres(engel,i),10000);
        p2(m,i+5)=poisstest2(szmet(~engel,i),szmet_unres(~engel,i),10000);
    end
    
    
end
set(leg,'position',[0.70 0.2 0.2 0.2],'units','normalized','box','off')
%% synchronizability analysis
synrand={};
for i=1:length(eeg)
    synrand=horzcat(synrand,(eeg(i).sz.glob_mat_rand));
end
for i=1:length(synrand)
    synrand{i}=squeeze(synrand{i});
end
synrand=cellfun(@(x) x(2,:),synrand,'UniformOutput',0);
synrand=cell2mat(synrand');
synrand=real(synrand);

figure
h=subplot(2,2,[1 2]);

beng=boxplot([szmet(engel,:),szmet_unres(engel,:),szmet_rand(engel,:)],...
    'positions',[1:5:21,2.5:5:22.5,4:5:24],'plotstyle','compact',...
    'colors',repmat([0 0 1;0 0.4 1;0 0.7 1],5,1),'datalim',[0 0.4]);
hold on;
bneng=boxplot([szmet(~engel,:),szmet_unres(~engel,:),szmet_rand(~engel,:)],...
    'positions',[1.5:5:21.5,3:5:23,4.5:5:24.5],'plotstyle','compact',...
    'colors',repmat([1 0 0;1 0 0.4;1 0 0.7],5,1),'datalim',[0 0.4]);

h=get(gca,'child');
hrem=findobj(get(h(1),'child'),'type','text');
delete(hrem);
hrem=findobj(get(h(2),'child'),'type','text');
delete(hrem);
box off

set(gca,'XTick',2.75:5:22.75,'XTickLabel',{'Pre';'Early';'Middle';'Late';'Post'},'FontSize',14)
ylabel('Synchronizability')

leg=legend([beng(2,1),beng(2,2),beng(2,3),bneng(2,1),bneng(2,2),bneng(2,3)],...
    'Engel 1','Engel 1 Unresected','Engel 1 Random Resection','Engel >1',...
    'Engel >1 Unresected','Engel >1 Random Resection','location','EastOutside');
set(leg,'box','off');


ymax=get(gca,'YLim');
ymin=ymax(1);
ymax=ymax(2)-ymax(1);
line([1:5:21;1.5:5:21.5],(0.68*ymax+ymin)*ones(2,5),'color','m','LineWidth',2,'Marker','+','MarkerSize',5)
line([2.5:5:22.5;3:5:23],(0.68*ymax+ymin)*ones(2,5),'color','m','LineWidth',2,'Marker','+','MarkerSize',5)

line([1:5:21;2.5:5:22.5],(0.78*ymax+ymin)*ones(2,5),'color','b','LineWidth',2,'Marker','+','MarkerSize',5)
line([1.5:5:21.5;3:5:23],(0.88*ymax+ymin)*ones(2,5),'color','r','LineWidth',2,'Marker','+','MarkerSize',5)

for i=1:5
    if p(m,i)<0.01 && p(m,i)>=0.001
        text(1.25+5*(i-1),0.7*ymax+ymin,'*','horizontalAlignment','center',...
            'FontSize',20,'FontWeight','bold','color','m')
    elseif p(m,i)<0.001
        text(1.25+5*(i-1),0.7*ymax+ymin,'**','horizontalAlignment','center',...
            'FontSize',20,'FontWeight','bold','color','m')
    end
end
for i=1:5
    if p(m,i+5)<0.01 && p(m,i+5)>=0.001
        text(2.75+5*(i-1),0.7*ymax+ymin,'*','horizontalAlignment','center',...
            'FontSize',20,'FontWeight','bold','color','m')
    elseif p(m,i+5)<0.001
        text(2.75+5*(i-1),0.7*ymax+ymin,'**','horizontalAlignment','center',...
            'FontSize',20,'FontWeight','bold','color','m')
    end
end

for i=1:5
    if p2(m,i)<0.01 && p(m,i)>=0.001
        text(1.75+5*(i-1),0.8*ymax+ymin,'*','horizontalAlignment','center',...
            'FontSize',20,'FontWeight','bold','color','b')
    elseif p(m,i)<0.001
        text(1.75+5*(i-1),0.8*ymax+ymin,'**','horizontalAlignment','center',...
            'FontSize',20,'FontWeight','bold','color','b')
    end
end
for i=1:5
    if p2(m,i+5)<0.01 && p(m,i+5)>=0.001
        text(2.25+5*(i-1),0.9*ymax+ymin,'*','horizontalAlignment','center',...
            'FontSize',20,'FontWeight','bold','color','r')
    elseif p(m,i+5)<0.001
        text(2.25+5*(i-1),0.9*ymax+ymin,'**','horizontalAlignment','center',...
            'FontSize',20,'FontWeight','bold','color','r')
    end
end


subplot 223
boxplot(synrand(engel,:)','positions',szmet_unres(engel,2),'plotstyle','compact','symbol','.w')
hold on;
boxplot(synrand(~engel,:)','positions',szmet_unres(~engel,2),'plotstyle','compact','colors','r','symbol','.w')
h=get(gca,'child');
hrem=findobj(get(h(1),'child'),'type','text');
delete(hrem);
hrem=findobj(get(h(2),'child'),'type','text');
delete(hrem);
%set(gca,'XTick',0:0.02:0.3,'XTickLabel',num2str((0:0.02:0.3)'))
set(gca,'XScale','log','YScale','log','FontSize',14)
xlim([0.001 0.3]);ylim([0.001 0.3]);
line([0.001 1],[0.001 1],'color','k');
set(gca,'XTick',get(gca,'YTick'),'XTickLabel',get(gca,'YTickLabel'))

xlabel('Unresected Synchronizability')
ylabel({'Random Resection';'Synchronizability'})
box off

subplot 224
scatter(szmet_unres(engel,2),min(synrand(engel,:)'),'b');
hold on;
scatter(szmet_unres(~engel,2),min(synrand(~engel,:)'),'r')
line([0.001 1],[0.001 1],'color','k')
xlim([0.001 0.3])
ylim([0.001 0.3])
set(gca,'XScale','log','YScale','log','FontSize',14)
xlabel('Unresected Synchronizability')
ylabel({'Minimum Random Resection';'Synchronizability'})
set(gca,'YAxisLocation','right')

btitle=uicontrol('style','text','units','normalized','string','Early Random Resection Analysis');
set(btitle,'position',[0.25,0.47,0.5,0.05],'FontWeight','bold','FontSize',14)
%%
figure;errorbar(1:5,engm,engerr,'b','LineWidth',3)
hold on; errorbar(1.05:1:5.05,engm_unres,engerr_unres,'c','LineWidth',3)
errorbar(1.1:1:5.1,nengm,nengerr,'r','LineWidth',3)
errorbar(1.15:1:5.15,nengm_unres,nengerr_unres,'m','LineWidth',3)
%%
figure
subplot 121
plot(szmet(engel,:)','b','LineWidth',1)
ylim([0.001 0.3])
set(gca,'YScale','log')
%hold on
subplot 122
plot(szmet(~engel,:)','r','LineWidth',1)
ylim([0.001 0.3])
set(gca,'YScale','log')
subplot 121
hold on
plot(szmet_unres(engel,:)','c','LineWidth',1)
ylim([0.001 0.3])
set(gca,'YScale','log')
subplot 122
hold on
plot(szmet_unres(~engel,:)','m','LineWidth',1)
ylim([0.001 0.3])
set(gca,'YScale','log')
%ylim([0.001 0.5])
%% boxplots
p=zeros(5,10);
p2=p;

ptengel=arrayfun(@(x) (strcmp(x.engel,'1a')),eeg);
inds=1:length(eeg);
ptind=arrayfun(@(x,y) repmat(y,1,length(x.sztype)),eeg,inds,'UniformOutput',0);
ptind=horzcat(ptind{:})';
figure
for m=1:5
    szmet=[];
    
    for i=1:length(eeg)
        for j=1:length(eeg(i).sz)
            szmet(size(szmet,1)+1,:)=eeg(i).sz(j).glob_mat(:,m);
        end
    end
    szmet_unres=[];
    for i=1:length(eeg)
        for j=1:length(eeg(i).sz)
            szmet_unres(size(szmet_unres,1)+1,:)=eeg(i).sz(j).glob_mat_unres(:,m);
        end
    end
    
    szmet_rand=[];
    if m==5
        for i=1:length(eeg)
            for j=1:length(eeg(i).sz)
                szmet_rand(size(szmet_rand,1)+1,:)=mean(real(eeg(i).sz(j).glob_mat_rand),3); %wtf
            end
        end
    else
        szmet_rand=zeros(size(szmet_unres));
    end
    
    szmet_frac=szmet_unres./szmet;
    szmet_min=szmet-szmet_unres;
    ptsz=zeros(length(eeg),5);
    ptsz_unres=ptsz;
    ptsz_frac=ptsz;
    ptsz_min=ptsz;
    ptsz_rand=ptsz;
    
    for i=1:5
        ptsz(:,i)=accumarray(ptind,szmet(:,i),[],@mean);
        ptsz_unres(:,i)=accumarray(ptind,szmet_unres(:,i),[],@mean);
        ptsz_frac(:,i)=accumarray(ptind,szmet_frac(:,i),[],@mean);
        ptsz_min(:,i)=accumarray(ptind,szmet_min(:,i),[],@mean);
        ptsz_rand(:,i)=accumarray(ptind,szmet_rand(:,i),[],@mean);
    end
    subplot(2,3,m)
    boxplot([ptsz(ptengel,:),ptsz_unres(ptengel,:),ptsz_rand(ptengel,:)],'positions',[1:5:21,2.5:5:22.5,4:5:24],'plotstyle','compact')
    hold on;
    boxplot([ptsz(~ptengel,:),ptsz_unres(~ptengel,:),ptsz_rand(~ptengel,:)],'positions',[1.5:5:21.5,3:5:23,4.5:5:24.5],'plotstyle','compact','colors','r')
    %xlim([0 16])
    
    engm=mean(ptsz(ptengel,:));
    engerr=std(ptsz(ptengel,:))/sqrt(sum(ptengel));
    
    engm_unres=mean(ptsz_unres(ptengel,:));
    engerr_unres=std(ptsz_unres(ptengel,:))/sqrt(sum(ptengel));
    
    engm_rand=mean(ptsz_rand(ptengel,:));
    engerr_rand=std(ptsz_rand(ptengel,:))/sqrt(sum(ptengel));
    
    nengm=mean(ptsz(~ptengel,:));
    nengerr=std(ptsz(~ptengel,:))/sqrt(sum(~ptengel));
    
    nengm_unres=mean(ptsz_unres(~ptengel,:));
    nengerr_unres=std(ptsz_unres(~ptengel,:))/sqrt(sum(~ptengel));
    
    nengm_rand=mean(ptsz_rand(~ptengel,:));
    nengerr_rand=std(ptsz_rand(~ptengel,:))/sqrt(sum(~ptengel));
    
    %errorbar([1:5:21,2.5:5:22.5,4:5:24],[engm,engm_unres,engm_rand],[engerr,engerr_unres,engerr_rand],'.b');
    %hold on;
    %engbar=bar([1:5:21,2.5:5:22.5,4:5:24],diag([engm,engm_unres,engm_rand]),0.3,'stacked','b');
    %for j=1:5
    %    set(engbar(j+5),'FaceColor',[0.3 0.3 1])
    %    set(engbar(j+10),'FaceColor',[0.6 0.6 1])
    %end
    %errorbar([1.5:5:21.5,3:5:23,4.5:5:24.5],[nengm,nengm_unres,nengm_rand],[nengerr,nengerr_unres,nengerr_rand],'.r')
    %nengbar=bar([1.5:5:21.5,3:5:23,4.5:5:24.5],diag([nengm,nengm_unres,nengm_rand]),0.3,'stacked','r');
    %for j=1:5
    %    set(nengbar(j+5),'FaceColor',[1 0.3 0.3])
    %    set(nengbar(j+10),'FaceColor',[1 0.6 0.6])
    %end
    if m==1
        title('Small World','FontSize',14)
    elseif m==2
        title('Small World','FontSize',14)
        ylim([-0.4 0.6])
    elseif m==3
        title('Average Clustering Coefficient','FontSize',14)
    elseif m==4
        ylim([0 3])
        title('Average Clustering Coefficient','FontSize',14)
    elseif m==5
        title('Synchronizability','FontSize',14)
        ylim([0 0.16])
        %leg=legend([engbar(1),engbar(6),nengbar(1),nengbar(6),engbar(11),nengbar(11)],...
        %    'Engel 1','Engel 1 Unresected','Engel >1','Engel >1 Unresected','Engel 1 Random Resection','Engel >1 Random Resection');
        %legend 'boxoff'
    end
    
    ymax=get(gca,'YLim');
    ymin=ymax(1);
    ymax=ymax(2)-ymax(1);
    set(gca,'XTick',2:5:22,'XTickLabel',{'Pre';'Early';'Middle';'Late';'Post'},'box','off')
    line([1:5:21;1.5:5:21.5],(0.9*ymax+ymin)*ones(2,5),'color','k')
    line([2.5:5:22.5;3:5:23],(0.9*ymax+ymin)*ones(2,5),'color','k')
    
    
    for i=1:5
        p(m,i)=poisstestu(ptsz(ptengel,i),ptsz(~ptengel,i),10000);
        if p(m,i)<0.01 && p(m,i)>=0.001
            text(1.25+5*(i-1),0.93*ymax+ymin,'*','horizontalAlignment','center',...
                'FontSize',20,'FontWeight','bold')
        elseif p(m,i)<0.001
            text(1.25+5*(i-1),0.93*ymax+ymin,'**','horizontalAlignment','center',...
                'FontSize',20,'FontWeight','bold')
        end
    end
    for i=1:5
        p(m,i+5)=poisstestu(ptsz_unres(ptengel,i),ptsz_unres(~ptengel,i),10000);
        if p(m,i+5)<0.01 && p(m,i+5)>=0.001
            text(2.75+5*(i-1),0.93*ymax+ymin,'*','horizontalAlignment','center',...
                'FontSize',20,'FontWeight','bold')
        elseif p(m,i+5)<0.001
            text(2.75+5*(i-1),0.93*ymax+ymin,'**','horizontalAlignment','center',...
                'FontSize',20,'FontWeight','bold')
        end
    end
    
    for i=1:5
        p2(m,i)=poisstest2(ptsz(ptengel,i),ptsz_unres(ptengel,i),10000);
        p2(m,i+5)=poisstest2(ptsz(~ptengel,i),ptsz_unres(~ptengel,i),10000);
    end
    
    
end
set(leg,'position',[0.70 0.2 0.2 0.2],'units','normalized','box','off')
%%
figure;errorbar(1:5,engm,engerr,'b','LineWidth',3)
hold on; errorbar(1.05:1:5.05,engm_unres,engerr_unres,'c','LineWidth',3)
errorbar(1.1:1:5.1,nengm,nengerr,'r','LineWidth',3)
errorbar(1.15:1:5.15,nengm_unres,nengerr_unres,'m','LineWidth',3)
%%
ptengel=arrayfun(@(x) (strcmp(x.engel,'1a')),eeg);
ptengel=ptengel+arrayfun(@(x) (2*strcmp(x.engel,'1b')),eeg);
ptengel=ptengel+arrayfun(@(x) (3*strcmp(x.engel,'2a')),eeg);
ptengel=ptengel+arrayfun(@(x) (4*strcmp(x.engel,'2b')),eeg);
ptengel=ptengel+arrayfun(@(x) (5*strcmp(x.engel,'3a')),eeg);
ptengel=ptengel+arrayfun(@(x) (6*strcmp(x.engel,'3b')),eeg);
ptengel=ptengel+arrayfun(@(x) (7*strcmp(x.engel,'4a')),eeg);
ptengel=ptengel+arrayfun(@(x) (8*strcmp(x.engel,'4b')),eeg);
%% ROC curve for synchronizability in early time point
[ptsz_unres_sort,sortind]=sort(ptsz_unres(:,2));
engelsort=ptengel(sortind);
for i=1:length(eeg)
    tpr(i)=sum(~engelsort(i:end))/nnz(~engelsort);
    fpr(i)=sum(engelsort(i:end))/nnz(engelsort);
end
figure
plot(fpr,tpr)
%% resected and non-resected averages
types=vertcat(eeg.sztype);
typeind=strcmp(types,'cps')+2*strcmp(types,'subcl');
typeind(typeind==0)=3;
engel=arrayfun(@(x) repmat((strcmp(x.engel,'1a')|strcmp(x.engel,'1b')),length(x.sz),1),eeg,'UniformOutput',0)';;
engel=vertcat(engel{:});

cmap=hsv(3);
cmap(3,:)=[1 1 1];
cmap=num2cell(cmap(typeind,:),2);
mtype='o+xsd';
%figure;
bay=[];
for m=1:5
    
    for i=1:length(eeg)
        
        res=logical(eeg(i).res);
        hub=zeros(length(eeg(i).leads),5,length(eeg(i).sz));
        for j=1:length(eeg(i).sz);
            hub(:,:,j)=squeeze(eeg(i).sz(j).metric(:,m,:));
        end
        
        hubr{i}=hub(res,:,:);
        hubu{i}=hub(~res,:,:);
        
        l=histc(hub(res,:,:),0:0.02:0.1)/sum(res);
        ev=histc(hub,0:0.02:0.1)/length(res);
        pr=sum(res)/length(res);
        
        bay=cat(3,bay,pr.*l./ev);
        
        %hubr{i}=reshape(hub(res,:,:),nnz(res),5*length(eeg(i).sz));
        %hubu{i}=reshape(hub(~res,:,:),nnz(~res),5*length(eeg(i).sz));
        
        xvalsu=reshape(bsxfun(@plus,[1:5]',0:0.2:(length(eeg(i).sz)-1)*0.2),1,5*length(eeg(i).sz));
        xvalsr=xvalsu+0.1;
        %if mod(i,4)==1;
        %    figure;
        %end
        %if mod(i,4)~=0
        %    subplot(1,4,mod(i,4))
        %else
        %    subplot(1,4,4)
        %end
        %boxplot(hubu{i},'plotStyle','compact','positions',xvalsu);
        %hold on
        %boxplot(hubr{i},'plotStyle','compact','positions',xvalsr,'colors','r');
        %line(repmat([1:5;1.2:5.2],1,length(eeg(i).sz)),[mhubu{i};mhubr{i}],'color','b','Marker','.')
        
    end
    
    mhubr=cellfun(@(x) reshape(mean(x),5,[]),hubr,'UniformOutput',0);
    mhubu=cellfun(@(x) reshape(mean(x),5,[]),hubu,'UniformOutput',0);
    
    mhubr=[mhubr{:}];
    mhubu=[mhubu{:}];
    
    for i=1:5
        %figure(i);
        subplot(2,3,m);
        a=line(repmat([i i+0.7]',1,131),[mhubu(i,:);mhubr(i,:)],'Color','b','Marker','.');
        hold on
        p=poisstest2(mhubu(i,:),mhubr(i,:),1000);
        text(0.19*i-0.1,0.9,num2str(p),'FontSize',12,'Units','Normalized')
        %for j=1:51
        %set(a,{'Color'},cmap);
        %end
        %scatter(mhubu(i,engel),mhubr(i,engel),'+','b');
        %scatter(mhubu(i,~engel),mhubr(i,~engel),'+','r');
        %maxm=max(get(gca,'XLim'),get(gca,'YLim'));
        %minm=min(get(gca,'XLim'),get(gca,'YLim'));
        %line([minm maxm],[minm maxm]);
    end
    %xlim([0 1])%;ylim([0 1])
    %xlim([0.7 6])
    %set(gca,'XTick',[1.4:5.4],'XTickLabel',{'pre';'early';'middle';'late';'post'},'FontSize',12)
    if m==1
        title('Degree','FontSize',14,'FontWeight','bold')
        ylabel('Degree','FontSize',14)
    elseif m==2
        title('Betweenness Centrality','FontSize',14,'FontWeight','bold')
        ylabel('Betweenness Centrality','FontSize',14)
    elseif m==3
        title('Clustering Coefficient','FontSize',14,'FontWeight','bold')
        ylabel('Clustering Coefficient','FontSize',14)
    elseif m==4
        title('Hubness','FontSize',14,'FontWeight','bold')
        ylabel('Hubness','FontSize',14)
    end
    hold on
end
h=suptitle('Unresected vs Resected Nodes Across All Seizures');
set(h,'FontWeight','bold','FontSize',15)
%% resected and not resected bayesian analysis
types={'cps','cps2gen','subcl'};
for k=1:4
    figure
    colors=get(gca,'colororder');
    for m=1:4
        hube=[];
        hubne=[];
        rese=[];
        resne=[];
        for i=1:length(eeg)
            %for i=[1:6,10:11,13:19]
            
            res=logical(eeg(i).res);
            %hub=zeros(length(eeg(i).leads),5,length(eeg(i).sz));
            for j=1:length(eeg(i).sz);
                %hub(:,:,j)=squeeze(eeg(i).sz(j).metric(:,m,:));
                if k~=1
                    if strcmp(eeg(i).sztype(j),types(k-1))
                        if engel(i)
                            rese=vertcat(rese,res);
                            hube=vertcat(hube,squeeze(eeg(i).sz(j).metric(:,m,:)));
                        else
                            resne=vertcat(resne,res);
                            hubne=vertcat(hubne,squeeze(eeg(i).sz(j).metric(:,m,:)));
                        end
                    end
                else
                    if engel(i)
                        rese=vertcat(rese,res);
                        hube=vertcat(hube,squeeze(eeg(i).sz(j).metric(:,m,:)));
                    else
                        resne=vertcat(resne,res);
                        hubne=vertcat(hubne,squeeze(eeg(i).sz(j).metric(:,m,:)));
                    end
                end
            end
        end
        
        
        rese=logical(rese);
        resne=logical(resne);
        if m==1
            centers=[6:2:14]';
            l=hist(hube(rese,:),centers);
            l=l/sum(rese);
        elseif m==2
            centers=[0:0.01:0.04]';
            l=hist(hube(rese,:),centers);
            l=l/sum(rese);
        elseif m==3
            centers=[0.5:0.1:0.9]';
            l=hist(hube(rese,:),centers);
            l=l/sum(rese);
        elseif m==4
            centers=[0:0.1:0.4]';
            l=hist(hube(rese,:),centers);
            l=l/sum(rese);
        elseif m==5
            centers=[0:0.07:0.28]';
            l=hist(hube(rese,:),centers);
            l=l/sum(rese);
        else
            %l=histc(hube(rese,:),0:0.2:1)/sum(rese);
            [l,centers]=hist(hube(rese,:),20);
            l=l/sum(rese);
        end
        %ev=histc(hube,0:0.2:1)/length(rese);
        ev=hist(hube,centers)/length(rese);
        pr=sum(rese)/length(rese);
        be=pr.*l./ev;
        
        %norme=be;
        %norme(isnan(norme))=0;
        %norme=sum(norme);
        %be=bsxfun(@times,be,1./(norme));
        
        %l=histc(hubne(resne,:),0:0.2:1)/sum(resne);
        l=hist(hubne(resne,:),centers)/sum(resne);
        %ev=histc(hubne,0:0.2:1)/length(resne);
        ev=hist(hubne,centers)/length(resne);
        pr=sum(resne)/length(resne);
        bne=pr.*l./ev;
        %normne=bne;
        %normne(isnan(normne))=0;
        %normne=sum(normne);
        %bne=bsxfun(@times,bne,1./(normne));
        
        e=hist(hube,centers);
        ne=hist(hubne,centers);
        
        re=be.*e;
        rne=bne.*ne;
        
        p=(re+rne)./(e+ne);
        re0=e.*p;
        rne0=ne.*p;
        
        observed=cat(3,re,e-re,rne,ne-rne);
        expected=cat(3,re0,e-re0,rne0,ne-rne0);
        
        chi2=sum(((observed-expected).^2)./expected,3);
        
        p=1-chi2cdf(chi2,1);
        
        belim=zeros(size(l,1),2,10);
        bnelim=belim;
        for i=1:5
            %[~,belim(:,:,i)]=binofit(be(:,i).*histc(hube(:,i),0:0.2:1),histc(hube(:,i),0:0.2:1));
            %[~,bnelim(:,:,i)]=binofit(bne(:,i).*histc(hubne(:,i),0:0.2:1),histc(hubne(:,i),0:0.2:1));
            [~,belim(:,:,i)]=binofit(round(re(:,i)),e(:,i));
            [~,bnelim(:,:,i)]=binofit(round(rne(:,i)),ne(:,i));
        end
        
        %be=reshape(be,6,1,5);
        %bne=reshape(bne,6,1,5);
        
        be=reshape(be,size(l,1),1,5);
        bne=reshape(bne,size(l,1),1,5);
        
        re=reshape(re,size(l,1),1,5);
        rne=reshape(re,size(l,1),1,5);
        
        rr=be./bne;
        ciu=rr.*reshape(exp(1.96*sqrt((1-be)./re + (1-bne)./rne)),size(l,1),1,5);
        cil=rr.*reshape(exp(-1.96*sqrt((1-be)./re + (1-bne)./rne)),size(l,1),1,5);
        
        %belim=abs(belim-[be,be]);
        %bnelim=abs(bnelim-[bne,bne]);
        
        
        %figure;
        for i=1:5
            subplot(5,4,4*(i-1)+m)
            %errorbar(0:0.2:1,be(:,1,i),belim(:,1,i),belim(:,2,i))
            
            patch([centers(~isnan(belim(:,1,i)))',fliplr(centers(~isnan(belim(:,2,i)))')],...
                [belim(~isnan(belim(:,1,i)),1,i)',fliplr(belim(~isnan(belim(:,2,i)),2,i)')],...
                colors(1,:),'FaceAlpha',0.5)
            hold on
            patch([centers(~isnan(bnelim(:,1,i)))',fliplr(centers(~isnan(bnelim(:,2,i)))')],...
                [bnelim(~isnan(bnelim(:,1,i)),1,i)',fliplr(bnelim(~isnan(bnelim(:,2,i)),2,i)')],...
                colors(2,:),'FaceAlpha',0.5)
            plot(centers,be(:,i),'LineWidth',3,'Marker','.','MarkerSize',20,'color',colors(1,:));
            plot(centers,bne(:,i),'LineWidth',3,'Marker','.','MarkerSize',20,'color',colors(2,:));
            %errorbar(0:0.2:1,bne(:,1,i),bnelim(:,1,i),bnelim(:,2,i),'r')
            ylim([0 0.6])
            %plot(centers,be(:,i)./bne(:,i),'g','LineWidth',3,'Marker','.','MarkerSize',25)
            %patch([centers(~(isnan(cil(:,1,i))|isinf(cil(:,1,i))))',fliplr(centers(~(isnan(ciu(:,1,i))|isinf(ciu(:,1,i))))')],[cil(~(isnan(cil(:,1,i))|isinf(cil(:,1,i))),1,i)',fliplr(ciu(~(isnan(ciu(:,1,i))|isinf(ciu(:,1,i))),1,i)')],[0 1 0])
            %alpha(0.5)
            set(gca,'box','off','XTick',centers)
            ymax=get(gca,'YLim');
            ymax=ymax(2);
            sig1=find(p(:,i)<0.01);
            sig2=find(p(:,i)<0.001);
            %sig3=find(p(:,i)<0.001);
            %sig2=setxor(sig2,sig3);
            sig1=setxor(sig1,sig2);
            for j=1:length(sig1)
                text(centers(sig1(j)),0.85*ymax,'*','FontSize',20,'horizontalAlignment','center','FontWeight','bold')
            end
            for j=1:length(sig2)
                text(centers(sig2(j)),0.85*ymax,{'**'},'FontSize',20,'horizontalAlignment','center','FontWeight','bold')
            end
            %for j=1:length(sig3)
            %    text(centers(sig3(j)),0.8*ymax,{'*';'*';'*'})
            %end
            if i==5
                if m==1
                    xlabel('Degree','FontSize',12)
                    set(gca','XTickLabel',{'6','','10','','14'})
                elseif m==2
                    xlabel({'Betweenness';'Centrality'},'FontSize',12)
                    set(gca','XTickLabel',{'0','','0.02','','0.04'})
                elseif m==3
                    xlabel({'Clustering';'Coefficient'},'FontSize',12)
                    set(gca','XTickLabel',{'0.5','','0.7','','0.9'})
                elseif m==4
                    xlabel('Hubness','FontSize',12)
                    set(gca','XTickLabel',{'0','','0.14','','0.28'})
                elseif m==5
                    xlabel('New Metric','FontSize',12)
                end
            else
                set(gca,'XTickLabel',{})
            end
            
            if m==1
                ylabel('Probability','FontSize',12)
            else
                set(gca,'YTickLabel',{})
            end
            if m==1
                xlim([2 16])
            elseif m==2
                xlim([-0.01 0.05])
            elseif m==3
                xlim([0.4 1])
            elseif m==4
                xlim([-0.1 0.5])
                if i==1
                    text(1.1,0.5,'Pre','Units','Normalized','Rotation',-90,'FontSize',16,'FontWeight','bold','horizontalAlignment','center')
                elseif i==2
                    text(1.1,0.5,'Early','Units','Normalized','Rotation',-90,'FontSize',16,'FontWeight','bold','horizontalAlignment','center')
                elseif i==3
                    text(1.1,0.5,'Middle','Units','Normalized','Rotation',-90,'FontSize',16,'FontWeight','bold','horizontalAlignment','center')
                elseif i==4
                    text(1.1,0.5,'Late','Units','Normalized','Rotation',-90,'FontSize',16,'FontWeight','bold','horizontalAlignment','center')
                elseif i==5
                    text(1.1,0.5,'Post','Units','Normalized','Rotation',-90,'FontSize',16,'FontWeight','bold','horizontalAlignment','center')
                end
            elseif m==5
                xlim([-0.07 0.35])
            end
            
        end
    end
    if k==1
        h=suptitle({'Resection Probability Given Metric';'All Seizures'});
    elseif k==2
        h=suptitle({'Resection Probability Given Metric';'CPS'});
    elseif k==3
        h=suptitle({'Resection Probability Given Metric';'CPS w/ Gen.'});
    else
        h=suptitle({'Resection Probability Given Metric';'Subclinical Seizures'});
    end
    set(h,'FontSize',18)
end

%%
figure;
for i=1:5
    subplot(2,3,i)
    %errorbar(0:0.2:1,be(:,1,i),belim(:,1,i),belim(:,2,i))
    plot(0:0.2:0.8,be(1:5,i),'LineWidth',3);
    patch([0:0.2:0.8,0.8:-0.2:0],[belim(1:5,1,i)',fliplr(belim(1:5,2,i)')],[0 0 1])
    hold on
    plot(0:0.2:0.8,bne(1:5,i),'r','LineWidth',3);
    patch([0:0.2:0.8,0.8:-0.2:0],[bnelim(1:5,1,i)',fliplr(bnelim(1:5,2,i)')],[1 0 0])
    alpha(0.5)
    %errorbar(0:0.2:1,bne(:,1,i),bnelim(:,1,i),bnelim(:,2,i),'r')
end
%% resected and not resected bayesian analysis for baseline
for m=1:4
    bay=[];
    for i=1:length(eeg)
        
        res=logical(eeg(i).res);
        hub=zeros(length(eeg(i).leads),1,2);
        for j=1:size(eeg(i).bsl.metric,3);
            hub(:,:,j)=squeeze(eeg(i).bsl.metric(:,m,j));
        end
        
        hubr{i}=hub(res,:,:);
        hubu{i}=hub(~res,:,:);
        
        l=histc(hub(res,:,:),0:0.02:0.1)/sum(res);
        ev=histc(hub,0:0.02:0.1)/length(res);
        pr=sum(res)/length(res);
        
        bay=cat(3,bay,pr.*l./ev);
        
        %hubr{i}=reshape(hub(res,:,:),nnz(res),5*length(eeg(i).sz));
        %hubu{i}=reshape(hub(~res,:,:),nnz(~res),5*length(eeg(i).sz));
        
        xvalsu=reshape(bsxfun(@plus,[1:5]',0:0.2:(length(eeg(i).sz)-1)*0.2),1,5*length(eeg(i).sz));
        xvalsr=xvalsu+0.1;
        %if mod(i,4)==1;
        %    figure;
        %end
        %if mod(i,4)~=0
        %    subplot(1,4,mod(i,4))
        %else
        %    subplot(1,4,4)
        %end
        %boxplot(hubu{i},'plotStyle','compact','positions',xvalsu);
        %hold on
        %boxplot(hubr{i},'plotStyle','compact','positions',xvalsr,'colors','r');
        %line(repmat([1:5;1.2:5.2],1,length(eeg(i).sz)),[mhubu{i};mhubr{i}],'color','b','Marker','.')
        
    end
    
    mhubr=cellfun(@(x) reshape(mean(x),2,[]),hubr,'UniformOutput',0);
    mhubu=cellfun(@(x) reshape(mean(x),2,[]),hubu,'UniformOutput',0);
    
    mhubr=[mhubr{:}];
    mhubu=[mhubu{:}];
    
    for i=1:2
        %figure(i);
        subplot(2,2,m);
        a=line(repmat([i i+0.7]',1,36),[mhubu(i,:);mhubr(i,:)],'Color','b','Marker','.');
        hold on
        p=poisstest2(mhubu(i,:),mhubr(i,:),1000);
        text(0.19*i-0.1,0.9,num2str(p),'FontSize',12,'Units','Normalized')
        %for j=1:51
        %set(a,{'Color'},cmap);
        %end
        %scatter(mhubu(i,engel),mhubr(i,engel),'+','b');
        %scatter(mhubu(i,~engel),mhubr(i,~engel),'+','r');
        %maxm=max(get(gca,'XLim'),get(gca,'YLim'));
        %minm=min(get(gca,'XLim'),get(gca,'YLim'));
        %line([minm maxm],[minm maxm]);
    end
    %xlim([0 1])%;ylim([0 1])
    %xlim([0.7 6])
    %set(gca,'XTick',[1.4:5.4],'XTickLabel',{'pre';'early';'middle';'late';'post'},'FontSize',12)
    if m==1
        title('Degree','FontSize',14,'FontWeight','bold')
        ylabel('Degree','FontSize',14)
    elseif m==2
        title('Betweenness Centrality','FontSize',14,'FontWeight','bold')
        ylabel('Betweenness Centrality','FontSize',14)
    elseif m==3
        title('Clustering Coefficient','FontSize',14,'FontWeight','bold')
        ylabel('Clustering Coefficient','FontSize',14)
    elseif m==4
        title('Hubness','FontSize',14,'FontWeight','bold')
        ylabel('Hubness','FontSize',14)
    end
    hold on
end
h=suptitle('Unresected vs Resected Nodes Across All Seizures');
set(h,'FontWeight','bold','FontSize',15)
%% seizure onset zones bayesian analysis
figure
types={'cps','cps2gen','subcl'};
colors=get(gca,'colororder');
colors(2:4,:)=colors(1:3,:);
colors(1,:)=[0.25,0.25,0.25];
for k=1:4
    for m=1:3
        hubcat=[];
        szocat=[];
        for i=1:length(eeg)
            %for i=[1:6,10:11,13:19]
            
            if k~=1
                typelog=strcmp(eeg(i).sztype,types(k-1));
                szo=logical(eeg(i).szo(:,typelog));
                szocat=vertcat(szocat,reshape(szo,[],1));
                %hub=zeros(length(eeg(i).leads),5,length(eeg(i).sz));
                for j=1:length(eeg(i).sz);
                    if strcmp(eeg(i).sztype(j),types(k-1))
                        hubcat=vertcat(hubcat,[mean(eeg(i).bsl.metric(:,m,:),3),squeeze(eeg(i).sz(j).metric(:,m,:))]);
                    end
                    
                end
            else
                szo=logical(eeg(i).szo);
                szocat=vertcat(szocat,reshape(szo,[],1));
                %hub=zeros(length(eeg(i).leads),5,length(eeg(i).sz));
                for j=1:length(eeg(i).sz)
                    hubcat=vertcat(hubcat,[mean(eeg(i).bsl.metric(:,m,:),3),squeeze(eeg(i).sz(j).metric(:,m,:))]);
                    
                end
            end
        end
        
        szocat=logical(szocat);
        
        %[l,centers]=hist(hubcat(szocat,:),10);
        %l=l/sum(szocat);
        %ev=histc(hube,0:0.2:1)/length(rese);
        
        if m==1
            centers=[6:2:14]';
            l=hist(hubcat(szocat,:),centers);
            l=l/sum(szocat);
        elseif m==2
            centers=[0:0.01:0.04]';
            l=hist(hubcat(szocat,:),centers);
            l=l/sum(szocat);
        elseif m==3
            centers=[0.5:0.1:0.9]';
            l=hist(hubcat(szocat,:),centers);
            l=l/sum(szocat);
        elseif m==4
            centers=[0:0.1:0.4]';
            l=hist(hubcat(szocat,:),centers);
            l=l/sum(szocat);
        elseif m==5
            centers=[0:0.07:0.28]';
            l=hist(hubcat(szocat,:),centers);
            l=l/sum(szocat);
        else
            %l=histc(hube(rese,:),0:0.2:1)/sum(rese);
            [l,centers]=hist(hubcat(szocat,:),20);
            l=l/sum(rese);
        end
        
        ev=hist(hubcat,centers)/length(szocat);
        pr=sum(szocat)/length(szocat);
        be=pr.*l./ev;
        
        e=hist(hubcat,centers);
        
        re=be.*e;
        %rne=bne.*ne;
        
        %p=(re+rne)./(e+ne);
        %re0=e.*p;
        %rne0=ne.*p;
        
        %observed=cat(3,re,e-re,rne,ne-rne);
        %expected=cat(3,re0,e-re0,rne0,ne-rne0);
        
        %chi2=sum(((observed-expected).^2)./expected,3);
        
        %p=1-chi2cdf(chi2,1);
        
        belim=zeros(size(l,1),2,6);
        bnelim=belim;
        for i=1:6
            %[~,belim(:,:,i)]=binofit(be(:,i).*histc(hube(:,i),0:0.2:1),histc(hube(:,i),0:0.2:1));
            %[~,bnelim(:,:,i)]=binofit(bne(:,i).*histc(hubne(:,i),0:0.2:1),histc(hubne(:,i),0:0.2:1));
            [~,belim(:,:,i)]=binofit(round(re(:,i)),e(:,i));
            %[~,bnelim(:,:,i)]=binofit(round(rne(:,i)),ne(:,i));
        end
        
        %be=reshape(be,6,1,5);
        %bne=reshape(bne,6,1,5);
        
        be=reshape(be,size(l,1),1,6);
        %bne=reshape(bne,size(l,1),1,5);
        
        re=reshape(re,size(l,1),1,6);
        %rne=reshape(re,size(l,1),1,5);
        
        %rr=be./bne;
        %ciu=rr.*reshape(exp(1.96*sqrt((1-be)./re + (1-bne)./rne)),size(l,1),1,5);
        %cil=rr.*reshape(exp(-1.96*sqrt((1-be)./re + (1-bne)./rne)),size(l,1),1,5);
        
        %belim=abs(belim-[be,be]);
        %bnelim=abs(bnelim-[bne,bne]);
        
        
        %figure;
        for i=1:6
            if k==1
                subplot(6,3,3*(i-1)+m)
            else
                subplot(7,3,3*(i-1)+m)
            end
            %errorbar(centers+0.1*centers(1)*(k-2),be(:,1,i),belim(:,1,i),belim(:,2,i),'LineWidth',2,'Marker','o','MarkerSize',3,'MarkerFaceColor',colors(k,:))
            if k==1
                plot(centers,be(:,i),'LineWidth',3,'Marker','.','MarkerSize',20,'color','k');
            else
                plot(centers,be(:,i),'LineWidth',3,'Marker','.','MarkerSize',20);
            end
            pat=patch([centers(~isnan(belim(:,1,i)))',fliplr(centers(~isnan(belim(:,2,i)))')],...
                [belim(~isnan(belim(:,1,i)),1,i)',fliplr(belim(~isnan(belim(:,2,i)),2,i)')],...
                colors(k,:),'EdgeColor','none');
            hAnnotation = get(pat,'Annotation');
            hLegendEntry = get(hAnnotation','LegendInformation');
            set(hLegendEntry,'IconDisplayStyle','off')
            hold on
            %plot(centers,bne(:,i),'r','LineWidth',3,'Marker','.','MarkerSize',25);
            %patch([centers(~isnan(bnelim(:,1,i)))',fliplr(centers(~isnan(bnelim(:,2,i)))')],[bnelim(~isnan(bnelim(:,1,i)),1,i)',fliplr(bnelim(~isnan(bnelim(:,2,i)),2,i)')],[1 0 0])
            %errorbar(0:0.2:1,bne(:,1,i),bnelim(:,1,i),bnelim(:,2,i),'r')
            ylim([0.005 0.4])
            %plot(centers,be(:,i)./bne(:,i),'g','LineWidth',3,'Marker','.','MarkerSize',25)
            %patch([centers(~(isnan(cil(:,1,i))|isinf(cil(:,1,i))))',fliplr(centers(~(isnan(ciu(:,1,i))|isinf(ciu(:,1,i))))')],[cil(~(isnan(cil(:,1,i))|isinf(cil(:,1,i))),1,i)',fliplr(ciu(~(isnan(ciu(:,1,i))|isinf(ciu(:,1,i))),1,i)')],[0 1 0])
            alpha(0.5)
            ymax=get(gca,'YLim');
            ymax=ymax(2);
            set(gca,'box','off','XTick',centers)
            set(gca,'YScale','log','YTick',[0.01,0.1])
            %sig1=find(p(:,i)<0.05);
            %sig2=find(p(:,i)<0.01);
            %sig3=find(p(:,i)<0.001);
            %sig2=setxor(sig2,sig3);
            %sig1=setxor(sig1,sig2);
            %for j=1:length(sig1)
            %    text(centers(sig1(j)),0.8*ymax,'*')
            %end
            %for j=1:length(sig2)
            %    text(centers(sig2(j)),0.8*ymax,{'*';'*'})
            %end
            %for j=1:length(sig3)
            %    text(centers(sig3(j)),0.8*ymax,{'*';'*';'*'})
            %end
            if i==6
                if m==1
                    xlabel('Degree','FontSize',12)
                    set(gca','XTickLabel',{'6','','10','','14'})
                elseif m==2
                    xlabel({'Betweenness';'Centrality'},'FontSize',12)
                    set(gca','XTickLabel',{'0','','0.02','','0.04'})
                elseif m==3
                    xlabel({'Clustering';'Coefficient'},'FontSize',12)
                    set(gca','XTickLabel',{'0.5','','0.7','','0.9'})
                elseif m==4
                    xlabel('Hubness','FontSize',12)
                    set(gca','XTickLabel',{'0','','0.14','','0.28'})
                end
            else
                set(gca,'XTickLabel',{})
            end
            
            if m==1
                ylabel('Probability','FontSize',12)
            else
                set(gca,'YTickLabel',{})
            end
            if m==1
                xlim([2 16])
            elseif m==2
                xlim([-0.01 0.05])
            elseif m==3
                xlim([0.4 1])
                %elseif m==4
                %    xlim([-0.1 0.5])
                if i==1
                    text(1.1,0.5,'Int','Units','Normalized','Rotation',-90,'FontSize',16,'FontWeight','bold','horizontalAlignment','center')
                elseif i==2
                    text(1.1,0.5,'Pre','Units','Normalized','Rotation',-90,'FontSize',16,'FontWeight','bold','horizontalAlignment','center')
                elseif i==3
                    text(1.1,0.5,'Early','Units','Normalized','Rotation',-90,'FontSize',16,'FontWeight','bold','horizontalAlignment','center')
                elseif i==4
                    text(1.1,0.5,'Middle','Units','Normalized','Rotation',-90,'FontSize',16,'FontWeight','bold','horizontalAlignment','center')
                elseif i==5
                    text(1.1,0.5,'Late','Units','Normalized','Rotation',-90,'FontSize',16,'FontWeight','bold','horizontalAlignment','center')
                elseif i==6
                    text(1.1,0.5,'Post','Units','Normalized','Rotation',-90,'FontSize',16,'FontWeight','bold','horizontalAlignment','center')
                end
            end
            %set(gca,'box','off')
        end
    end
    if k==1
        h=suptitle({'Seizure Onset Probability Given Metric';'All Seizures'});
        set(h,'FontSize',18)
        figure
    end
end
h=suptitle({'Seizure Onset Probability Given Metric';'By Seizure Type'});
set(h,'FontSize',18)
leg=legend('CPS','CPS w/ Gen','Subclinical');
set(leg,'box','off','FontSize',12,'Location','South','Position',[0.45,0.05,0.15,0.1])
%% highest nodes


for i=1:length(eeg)
    topleads=cell(5,4,5,length(eeg(i).sz));
    for j=1:length(eeg(i).sz)
        
        [sortmet,ind]=sort(eeg(i).sz(j).metric,'descend');
        topleads(:,:,:,j)=eeg(i).leads(ind(1:5,:,:));
        
    end
    eeg(i).topleads=topleads;
end

list=struct('topleads',[]);
k=1;
for i=1:length(eeg)
    for j=1:length(eeg(i).sz)
        list(k).topleads=eeg(i).topleads(:,1,1,j); %%%%change index here
        %second index is metric, third is time period
        list(k).pt=eeg(i).name;
        list(k).sznum=j;
        k=k+1;
    end
end

list2=cell(length(list),7);
for i=1:length(list)
    list2(i,1:5)=list(i).topleads;
    list2(i,6)={list(i).pt};
    list2(i,7)={list(i).sznum};
end


%% resected as function of outcome

ptengel=arrayfun(@(x) (strcmp(x.engel,'1a')),eeg);
ptengel=ptengel';

figure;
colors=get(gca,'colororder');

h(1)=subplot(3,4,1);
restot=arrayfun(@(x) nnz(x.res),eeg);
bar(1,mean(restot(ptengel)),'FaceColor',[0.5 0.5 0.5])
hold on
bar(2,mean(restot(~ptengel)),'FaceColor',[0.5 0.5 0.5])
hold on
errorbar(1:2,[mean(restot(ptengel)),mean(restot(~ptengel))],[std(restot(ptengel))...
    /sqrt(nnz(ptengel)),std(restot(~ptengel))/sqrt(nnz(~ptengel))],...
    'LineStyle','none','LineWidth',2,'color','k')

ylabel('Electrode Areas Resected')
%xlabel('Engel')
ylim([0 17])
xlim([0 3])
set(gca,'XTick',1:2,'XTickLabel',{'1a','>1a'},'box','off','FontSize',10)

ymax=get(gca,'YLim');
%line([1,2],[0.92*ymax(2),0.92*ymax(2)],'LineWidth',3,'Color','k')
%text(1.5,0.97*ymax(2),'N.S.','horizontalAlignment','center')

h(2)=subplot(3,4,2);
restotper=arrayfun(@(x) nnz(x.res)/length(x.res),eeg);
%colors=get(gca,'colororder');
bar(1,mean(restotper(ptengel)),'FaceColor',[0.5 0.5 0.5])
hold on
bar(2,mean(restotper(~ptengel)),'FaceColor',[0.5 0.5 0.5])
hold on
errorbar(1:2,[mean(restotper(ptengel)),mean(restotper(~ptengel))],[std(restotper(ptengel))...
    /sqrt(nnz(ptengel)),std(restotper(~ptengel))/sqrt(nnz(~ptengel))],...
    'LineStyle','none','LineWidth',2,'color','k')

ylabel({'Fraction of Electrode';'Areas Resected'})
%xlabel('Engel')
ylim([0 0.35])
xlim([0 3])
set(gca,'XTick',1:2,'XTickLabel',{'1a','>1a'},'box','off','FontSize',10,'YAxisLocation','right')

ymax=get(gca,'YLim');
%line([1,2],[0.92*ymax(2),0.92*ymax(2)],'LineWidth',3,'Color','k')
%text(1.5,0.97*ymax(2),'N.S.','horizontalAlignment','center')

h(3)=subplot(3,4,5:8);

ptszo=arrayfun(@(x) logical(sum(x.szo,2)),eeg,'UniformOutput',0);
res=arrayfun(@(x) x.res,eeg,'UniformOutput',0);
szotot=arrayfun(@(x) nnz(x.szo),eeg);

resszo=NaN*ones(1,length(eeg));
for i=1:length(eeg)
    resszo(i)=nnz(ptszo{i}&res{i});
    szoper(i)=resszo(i)/nnz(ptszo{i});
    szofrac(i)=nnz(ptszo{i})./length(ptszo{i});
end

%szofrac=arrayfun(@(x)  nnz(x.szo)/length())
%szoper=resszo./ptszo;

bar(1,mean(szoper(ptengel),'omitnan'),'FaceColor',[0.5 0.5 0.5])
hold on
bar(2,mean(szoper(~ptengel),'omitnan'),'FaceColor',[0.5 0.5 0.5])
hold on
errorbar(1:2,[mean(szoper(ptengel),'omitnan'),mean(szoper(~ptengel),'omitnan')],[std(szoper(ptengel),'omitnan')...
    /sqrt(nnz(ptengel)),std(szoper(~ptengel),'omitnan')/sqrt(nnz(~ptengel))],...
    'LineStyle','none','LineWidth',2,'color','k')

ylabel({'Fraction of SOZ Resected'})
%xlabel('Engel')
ylim([0 1])
%xlim([0 3])
%set(gca,'XTick',1:2,'XTickLabel',{'1a','>1a'},'box','off','FontSize',10)

ymax=get(gca,'YLim');
%line([1,2],[0.92*ymax(2),0.92*ymax(2)],'LineWidth',3,'Color','k')
%text(1.5,0.97*ymax(2),'N.S.','horizontalAlignment','center')


type={'cps','cps2gen','subcl'};
p=zeros(1,3);
szofractype=NaN*ones(3,length(eeg));
pttypevec=NaN*ones(3,length(eeg));
for k=1:3
    %subplot(3,4,k+5)
    pttype=arrayfun(@(x) strcmp(x.sztype,type(k)),eeg,'UniformOutput',0);
    for i=1:length(eeg)
        if nnz(pttype{i})~=0
            pttypevec(k,i)=1;
        end
    end
    resszotype=NaN*ones(1,length(eeg));
    ptszotype=cell(1,length(eeg));
    for i=1:length(eeg)
        ptszotype{i}=logical(sum(eeg(i).szo(:,pttype{i}),2));
        resszotype(i)=nnz(ptszotype{i}&res{i});
    end
    szotottype=cellfun(@nnz,ptszotype);
    szofractype(k,:)=cellfun(@(x) nnz(x)./length(x),ptszotype).*pttypevec(k,:);
    %szofractype(k,:)=szotottype./length(ptszotype);
    szopertype=resszotype./szotottype;
    
    bar(1+3*k,mean(szopertype(ptengel),'omitnan'),'FaceColor',colors(k,:));
    hold on
    bar(2+3*k,mean(szopertype(~ptengel),'omitnan'),'FaceColor',colors(k,:))
    hold on
    errorbar(1+3*k:2+3*k,[mean(szopertype(ptengel),'omitnan'),mean(szopertype(~ptengel),'omitnan')],...
        [std(szopertype(ptengel),'omitnan')/sqrt(nnz(ptengel)),std(szopertype(~ptengel),...
        'omitnan')/sqrt(nnz(~ptengel))],'LineStyle','none','LineWidth',2,'color','k')
    
    szopertypee=szopertype(ptengel);
    szopertypee=szopertypee(isnan(szopertypee)==0);
    szopertypeecell{k}=szopertypee;
    szopertypene=szopertype(~ptengel);
    szopertypene=szopertypene(isnan(szopertypene)==0);
    szopertypenecell{k}=szopertypene;
    p(k)=poisstestu(szopertypee,szopertypene,1000);
    
    %ylabel({'Fraction of Seizure';'Onset Zone Resected'})
    %xlabel('Engel')
    ylim([0 1])
    %xlim([0 3])
    %set(gca,'XTick',1:2,'XTickLabel',{'1a','>1a'},'box','off','FontSize',10)
    
    ymax=get(gca,'YLim');
    %line([1,2],[0.92*ymax(2),0.92*ymax(2)],'LineWidth',3,'Color','k')
    %text(1.5,0.97*ymax(2),'N.S.','horizontalAlignment','center')
end

set(gca,'XTick',[1,2,4,5,7,8,10,11],'XTickLabel',{'1a','>1a','1a','>1a','1a','>1a','1a','>1a'},'box','off','FontSize',10)

h(4)=subplot(3,4,4);
bar(1,'FaceColor',[0.5 0.5 0.5],'Visible','off');
hold on
for k=1:3
    bar(1,'FaceColor',colors(k,:),'Visible','off')
end
set(gca,'Visible','off')
[leg,obj]=legend('All','CPS','CPS w/ Gen.','Subclinical');
set(leg,'box','off','location','east','FontSize',6)

h(5)=subplot(3,4,9:12);

bar(1,mean(szofrac(ptengel),'omitnan'),'FaceColor',[0.5,0.5,0.5]);
hold on
bar(2,mean(szofrac(~ptengel),'omitnan'),'FaceColor',[0.5,0.5,0.5])
hold on
errorbar(1:2,[mean(szofrac(ptengel),'omitnan'),mean(szofrac(~ptengel),'omitnan')],...
    [std(szofrac(ptengel),'omitnan')/sqrt(nnz(ptengel)),std(szofrac(~ptengel),...
    'omitnan')/sqrt(nnz(~ptengel))],'LineStyle','none','LineWidth',2,'color','k')

%xlabel('Engel')
ylabel('Fraction of Areas in SOZ')
%ylim([0 0.3])
%xlim([0 3])
%set(gca,'XTick',1:2,'XTickLabel',{'1a','>1a'},'box','off','FontSize',10)

pfrac=zeros(1,3);
for k=1:3
    
    %subplot(3,4,9+k)
    bar(1+3*k,mean(szofractype(k,ptengel),'omitnan'),'FaceColor',colors(k,:));
    hold on
    bar(2+3*k,mean(szofractype(k,~ptengel),'omitnan'),'FaceColor',colors(k,:))
    hold on
    errorbar(1+3*k:2+3*k,[mean(szofractype(k,ptengel),'omitnan'),mean(szofractype(k,~ptengel),'omitnan')],...
        [std(szofractype(k,ptengel),'omitnan')/sqrt(nnz(~isnan(pttypevec(k,:)))),std(szofractype(k,~ptengel),...
        'omitnan')/sqrt(nnz(~isnan(pttypevec(k,:))))],'LineStyle','none','LineWidth',2,'color','k')
    
    
    %xlim([0 3])
    %set(gca,'XTick',1+3*k:2+3*k,'XTickLabel',{'1a','>1a'},'box','off','FontSize',10)
    
    szofractypee=szofractype(k,ptengel);
    szofractypee=szofractypee(isnan(szofractypee)==0);
    szofractypeecell{k}=szofractypee;
    szofractypene=szofractype(k,~ptengel);
    szofractypene=szofractypene(isnan(szofractypene)==0);
    szofractypenecell{k}=szofractypene;
    pfrac(k)=poisstestu(szofractypee,szofractypene,1000);
    
    
end

xlabel('Engel Classification')
ylim([0 0.43])
set(gca,'XTick',[1,2,4,5,7,8,10,11],'XTickLabel',{'1a','>1a','1a','>1a','1a','>1a','1a','>1a'},'box','off','FontSize',10)
line([7.5,10.5],[0.38,0.38],'LineWidth',2,'color','k')
line([6.5 8.5],[0.33,0.33],'LineWidth',1,'color','k')
line([6.5,6.5],[0.31,0.33],'LineWidth',1,'color','k')
line([8.5,8.5],[0.31,0.33],'LineWidth',1,'color','k')
line([9.5 11.5],[0.33,0.33],'LineWidth',1,'color','k')
line([9.5,9.5],[0.31,0.33],'LineWidth',1,'color','k')
line([11.5,11.5],[0.31,0.33],'LineWidth',1,'color','k')
line([7.5,7.5],[0.34,0.38],'LineWidth',2,'color','k')
line([10.5,10.5],[0.34,0.38],'LineWidth',2,'color','k')
text(9,0.39,'*','horizontalAlignment','center','FontSize',30,'FontWeight','bold')

engelcat=[repmat({'1a'},1,30),repmat({'>1a'},1,22)];
seizurecat=[repmat({'cps'},1,17),repmat({'cps2gen'},1,5),repmat({'subcl'},1,8),repmat({'cps'},1,11),repmat({'cps2gen'},1,3),repmat({'subcl'},1,8)];
grouping={engelcat,seizurecat};
%[~,~,stats]=anovan([horzcat(szofractypeecell{:}),horzcat(szofractypenecell{:})],grouping,'model','interaction','varnames',{'Engel','Type'})
%tightfig;
set(gcf,'position',[440.0000 378 600 500.8000])
h(3).Position(1)=0.15;
h(5).Position(1)=0.15;
h(1).Position(1)=0.15;
h(4).Position(1)=0.85;
h(1).Position(3)=h(3).Position(3)/4-0.04;
h(2).Position([1,3])=[h(3).Position(1)+h(3).Position(3)/4,h(3).Position(3)/4];
text(-0.9,1.07,'A','units','normalized','FontSize',16,'FontWeight','bold','Parent',h(1))
text(-0.02,1.07,'B','units','normalized','FontSize',16,'FontWeight','bold','Parent',h(2))
text(-0.18,1.07,'C','units','normalized','FontSize',16,'FontWeight','bold','Parent',h(3))
text(-0.18,1.07,'D','units','normalized','FontSize',16,'FontWeight','bold','Parent',h(5))
%sep=line([2.6 3],[0,0],'Parent',h(1),'LineWidth',3,'Color','w')
obj(5).Children.Vertices([1,2,5])=0.25;
obj(6).Children.Vertices([1,2,5])=0.25;
obj(7).Children.Vertices([1,2,5])=0.25;
obj(8).Children.Vertices([1,2,5])=0.25;
leg.Position(1)=0.63;
%% metric-defined SOZ
type={'cps','cps2gen','subcl'};
for k=3:3
    fracsoz=NaN.*zeros(20,length(eeg),17,5);
    fracsoz2=fracsoz;
    subclind=arrayfun(@(x) strcmp(x.sztype,type(k)),eeg,'UniformOutput',0);
    ptengel=arrayfun(@(x) (strcmp(x.engel,'1a')),eeg);
    metrics={'Degree','B/w-ness Cent.','Avg. Clust.','Hubness'};
    %subclind=vertcat(subclind{:});
    for i=1:length(eeg)
        if nnz(subclind{i})>0
            metricsmat=cat(4,eeg(i).szraw(subclind{i}).part); %concatenate metrics from different seizures
            metricsmat=max(metricsmat,[],4); %keep only max values for each lead
            
            for m=1:17
                
                for j=1:5
                    
                    [~,maxinds]=sort(metricsmat(:,m,j),'descend');
                    ressort=logical(eeg(i).res(maxinds));
                    fracsoz(:,i,m,j)=cumsum(ressort(1:20))./((1:20)');
                    fracsozcell{i}(:,m,j)=cumsum(ressort)./((1:length(ressort))');
                    
                    [~,maxinds]=sort(metricsmat(:,m,j),'ascend');
                    ressort=logical(eeg(i).res(maxinds));
                    fracsoz2(:,i,m,j)=cumsum(ressort(1:20))./((1:20)');
                    
                    
                end
                
            end
        end
    end
    
    figure
    p=zeros(20,17,5);
    tic
    for m=1:17
        m
        for j=1:5
            
            subplot(5,17,17*(j-1)+m)
            errorbar(1:20,mean(fracsoz(:,ptengel,m,j),2,'omitnan'),std(fracsoz(:,ptengel,m,j),[],2,'omitnan')...
                /sqrt(nnz(~isnan(fracsoz(1,ptengel,m,j)))),'LineWidth',2)
            hold on
            errorbar(1.2:1:20.2,mean(fracsoz(:,~ptengel,m,j),2,'omitnan'),std(fracsoz(:,~ptengel,m,j),[],2,'omitnan')...
                /sqrt(nnz(~isnan(fracsoz(1,~ptengel,m,j)))),'LineWidth',2)
            
            xlim([0 21])
            
            %plot(1:20,fracsoz(:,ptengel,m,j),'b');
            %hold on
            %plot(1:20,fracsoz(:,~ptengel,m,j),'r');
            
            %for i=1:20
            %    p(i,m,j)=poisstestu(fracsoz(i,ptengel,m,j),fracsoz(i,~ptengel,m,j),1000);
            %end
            for i=1:20
                fracsoz_nonan_e=fracsoz(i,ptengel,m,j);
                fracsoz_nonan_e(isnan(fracsoz_nonan_e))=[];
                
                fracsoz_nonan_ne=fracsoz(i,~ptengel,m,j);
                fracsoz_nonan_ne(isnan(fracsoz_nonan_ne))=[];
                
                p(i,m,j)=poisstestu(fracsoz_nonan_e,fracsoz_nonan_ne,1000);
            end
            
            set(gca,'box','off')
            if m==1
                ylim([0 0.8])
                ylabel({'Fraction';'Removed'})
            elseif m==2
                ylim([0 0.6])
            elseif m==3
                ylim([0 0.4])
            elseif m==4
                ylim([0 0.8])
            end
            ylim([0 0.8])
            
            if j==1
                %title(metrics{m});
            elseif j==5
                xlabel('SOZ Size')
            end
            
            if j~=5
                set(gca,'XTickLabel',{})
            end
            
            if m~=1
                set(gca,'YTickLabel',{})
            end
        end
        toc
    end
    leg=legend('1a','>1a');
    set(gcf,'position',[440,378,560,600])
    set(leg,'position',[0.47,0.01,0.1036,0.0583],'box','off')
end
%%
figure
for m=1:4
    for j=1:5
        
        subplot(5,4,4*(j-1)+m)
        errorbar(1:20,mean(fracsoz2(:,ptengel,m,j),2,'omitnan'),std(fracsoz2(:,ptengel,m,j),[],2,'omitnan')...
            /sqrt(nnz(ptengel)),'LineWidth',2)
        hold on
        errorbar(1.2:1:20.2,mean(fracsoz2(:,~ptengel,m,j),2,'omitnan'),std(fracsoz2(:,~ptengel,m,j),[],2,'omitnan')...
            /sqrt(nnz(~ptengel)),'LineWidth',2)
        
        xlim([0 21])
        %plot(1:20,fracsoz(:,ptengel,m,j),'b');
        %hold on
        %plot(1:20,fracsoz(:,~ptengel,m,j),'r');
    end
end
%%
for m=1:4
    for j=1:5
        subplot(5,4,4*(j-1)+m)
        hold on
        for i=1:length(eeg)
            h=plot(fracsozcell{i}(:,m,j));
            if ptengel(i)==1
                set(h,'color','b')
            else
                set(h,'color','r')
            end
        end
        
    end
end
%%
szconnmat=[];
for i=1:length(eeg)
    numleads=length(eeg(i).leads);
    for j=1:length(eeg(i).sztype)
        szconn=cellfun(@(x) length(x),eeg(i).sz(j).conn_cube)/numleads;
        szconnmat=vertcat(szconnmat,szconn);
    end
end
%% resected as function of outcome

ptengel=arrayfun(@(x) (strcmp(x.engel,'1a')),eeg);
ptengel=ptengel';

f=figure;
%colors=get(gca,'colororder');
colors=[0,118,192;163,2,52]/255;

h(1)=subplot(1,4,1);
restot=arrayfun(@(x) nnz(x.res),eeg);
bar(1,mean(restot(ptengel)),'FaceColor',colors(1,:),'EdgeColor','none')
hold on
bar(2,mean(restot(~ptengel)),'FaceColor',colors(2,:),'EdgeColor','none')
hold on
er1=ploterr([1,2],[mean(restot(ptengel)),mean(restot(~ptengel))],[],[std(restot(ptengel))...
    /sqrt(nnz(ptengel)),std(restot(~ptengel))/sqrt(nnz(~ptengel))],'abshhy',0.15);%,...
%'LineStyle','none','LineWidth',2,'color','k');
set(er1(1),'LineStyle','none')
set(er1(2),'LineWidth',1,'Color','k');

title({'Areas';'Resected'},'FontWeight','normal','FontSize',8)
text(-0.8,1.17,'a','FontSize',16,'FontWeight','bold','units','normalized')
ylabel('Number')
ylim([0 17])
xlim([0 3])
set(gca,'XTick',[1,2],'XTickLabel',{'\color{white}>\color{black}1a','>1a'},'box','off','FontSize',8,'XTickLabelRotation',90,'TickLabelInterpreter','tex')

ymax=get(gca,'YLim');
h(2)=subplot(1,4,2);
restotper=arrayfun(@(x) nnz(x.res)/length(x.res),eeg);
bar(1,mean(restotper(ptengel)),'FaceColor',colors(1,:),'EdgeColor','none')
hold on
bar(2,mean(restotper(~ptengel)),'FaceColor',colors(2,:),'EdgeColor','none')
hold on
er2=ploterr(1:2,[mean(restotper(ptengel)),mean(restotper(~ptengel))],[],[std(restotper(ptengel))...
    /sqrt(nnz(ptengel)),std(restotper(~ptengel))/sqrt(nnz(~ptengel))],'abshhy',0.15);%,...
%'LineStyle','none','LineWidth',2,'color','k')
set(er2(1),'LineStyle','none')
set(er2(2),'LineWidth',1,'Color','k');

title({'Areas';'Resected'},'FontWeight','normal','FontSize',8)
text(-0.35,1.17,'b','FontSize',16,'FontWeight','bold','units','normalized')

ylim([0 0.8])
xlim([0 3])
set(gca,'XTick',1:2,'XTickLabel',{'\color{white}>\color{black}1a','>1a'},'box','off','FontSize',8,'YAxisLocation','right','YTickLabel',{},'XTickLabelRotation',90)

ymax=get(gca,'YLim');
h(3)=subplot(1,4,3);

ptszo=arrayfun(@(x) logical(sum(x.szo,2)),eeg,'UniformOutput',0);
res=arrayfun(@(x) x.res,eeg,'UniformOutput',0);
szotot=arrayfun(@(x) nnz(x.szo),eeg);

resszo=NaN*ones(1,length(eeg));
for i=1:length(eeg)
    resszo(i)=nnz(ptszo{i}&res{i});
    szoper(i)=resszo(i)/nnz(ptszo{i});
    szofrac(i)=nnz(ptszo{i})./length(ptszo{i});
end


bar(1,mean(szoper(ptengel),'omitnan'),'FaceColor',colors(1,:),'EdgeColor','none')
hold on
bar(2,mean(szoper(~ptengel),'omitnan'),'FaceColor',colors(2,:),'EdgeColor','none')
hold on
er3=ploterr(1:2,[mean(szoper(ptengel),'omitnan'),mean(szoper(~ptengel),'omitnan')],[],[std(szoper(ptengel),'omitnan')...
    /sqrt(nnz(ptengel)),std(szoper(~ptengel),'omitnan')/sqrt(nnz(~ptengel))],'abshhy',0.15);%,...
%'LineStyle','none','LineWidth',2,'color','k')
set(er3(1),'LineStyle','none')
set(er3(2),'LineWidth',1,'Color','k');

title({'SOZ';'Resected'},'FontWeight','normal','FontSize',8)
text(-0.35,1.17,'c','FontSize',16,'FontWeight','bold','units','normalized')
ylim([0 0.8])
xlim([0 3])

set(gca,'XTick',1:2,'XTickLabel',{'\color{white}>\color{black}1a','>1a'},'box','off','FontSize',8,'YAxisLocation','right','YTickLabel',{},'XTickLabelRotation',90)


h(4)=subplot(1,4,4);

bar(1,mean(szofrac(ptengel),'omitnan'),'FaceColor',colors(1,:),'EdgeColor','none');
hold on
bar(2,mean(szofrac(~ptengel),'omitnan'),'FaceColor',colors(2,:),'EdgeColor','none')
hold on
er4=ploterr(1:2,[mean(szofrac(ptengel),'omitnan'),mean(szofrac(~ptengel),'omitnan')],[],...
    [std(szofrac(ptengel),'omitnan')/sqrt(nnz(ptengel)),std(szofrac(~ptengel),...
    'omitnan')/sqrt(nnz(~ptengel))],'abshhxy',0.15);%,'LineStyle','none','LineWidth',2,'color','k')
set(er4(1),'LineStyle','none')
set(er4(2),'LineWidth',1,'Color','k');

title({'Areas';'in SOZ'},'FontWeight','normal','FontSize',8)
text(-0.35,1.17,'d','FontSize',16,'FontWeight','bold','units','normalized')
ylabel('Fraction')
xlim([0 3])
ylim([0 0.8])
set(gca,'XTick',1:2,'XTickLabel',{'\color{white}>\color{black}1a','>1a'},'box','off','FontSize',8,'YAxisLocation','right','XTickLabelRotation',90)

margin=0.12;
ymargin=0.13;
space=0.07;
height=0.7;

width=(1-2*margin-3*space)/4;
h(1).Position([1,2,3,4])=[margin,ymargin,width,height];
h(2).Position([1,2,3,4])=[margin+width+space,ymargin,width,height];
h(3).Position([1,2,3,4])=[margin+2*(width+space),ymargin,width,height];
h(4).Position([1,2,3,4])=[margin+3*(width+space),ymargin,width,height];
