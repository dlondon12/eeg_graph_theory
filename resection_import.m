import=importdata('/Users/Dennis/Dropbox/Epilepsy Graph Theory/Demographic and Clinical/Electrodes Resected Edit_code3.xlsx','\t',1);
respts=fieldnames(import.textdata); %names of patients with resection data
eegpts={eeg.name}; %names of patients with eeg data

for i=1:length(respts)
%for i=1:1
    
    ptind=find(strcmpi(respts{i},eegpts)); %for each resected patient find the index in the eeg data
    if ~isempty(ptind)
        eegleads=eeg(ptind).leads;
        eegleads=strrep(eegleads,' ',''); %remove spaces from lead names
        nresleads=nnz(~isnan(import.data.(respts{i})(:,1))); %number of leads in resected file
        resleads=import.textdata.(respts{i})(2:(nresleads+1)); %resected leads
        resleads=strrep(resleads,' ',''); %remove spaces from resected leads
        [~,eegind,resind]=intersect(eegleads,resleads); %find leads common to eeg and resected file leads
        res=NaN.*ones(length(eegleads),1); %create vector of NaNs for resected
        res(eegind)=import.data.(respts{i})(resind,1); %populate vector in order of the eeg leads from the resected leads
        eeg(ptind).res=res; %import into eeg structure
        colhead=import.textdata.(respts{i})(1,:); %column headers
        szocol=find(~cellfun(@isempty,strfind(colhead,'SOZ'))); %find seizure origination columns
        clear szlabel
        for j=1:length(eeg(ptind).sztype)
            szlabel{j}=[eeg(ptind).sztype{j},eeg(ptind).szlabel{j}]; %restore seizure labels like cpsA
        end
        for j=1:length(szocol)
            szolabel=char(colhead(szocol(j))); %list of seizures onsets
            szind=find(strcmpi(szolabel(1:length(szolabel)-4),szlabel));
            if ~isempty(szind) && size(import.data.(respts{i}),2)>=length(szocol)
                szo=NaN.*ones(length(eegleads),1);
                szo(eegind)=import.data.(respts{i})(resind,szocol(j)-1);
                eeg(ptind).szo(:,szind)=szo;
            end
        end
    end
    
    
end



import=importdata('/Users/Dennis/Dropbox/Epilepsy Graph Theory/Demographic and Clinical/Table 1 data.xlsx','\t',1);
respts=import.textdata(1,2:end);
eegpts={eeg.name}; %names of patients with eeg data


for i=1:length(respts)
    
    ptind=find(strcmpi(respts{i},eegpts));
    if ~isempty(ptind)
        eeg(ptind).engel=import.textdata(15,i+1);
    end
    
end

import=importdata('/Users/Dennis/Dropbox/Epilepsy Graph Theory/Demographic and Clinical/Electrode Classification.xlsx','\t',1);
respts=fieldnames(import.textdata); %names of patients with resection data
eegpts={eeg.name}; %names of patients with eeg data

for i=1:length(respts)
%for i=1:1
    
    ptind=find(strcmpi(respts{i},eegpts)); %for each resected patient find the index in the eeg data
    if ~isempty(ptind)
        eegleads=eeg(ptind).leads;
        eegleads=strrep(eegleads,' ',''); %remove spaces from lead names
        nresleads=nnz(~isnan(import.data.(respts{i})(:,1))); %number of leads in resected file
        resleads=import.textdata.(respts{i})(2:(nresleads+1)); %resected leads
        resleads=strrep(resleads,' ',''); %remove spaces from resected leads
        [~,eegind,resind]=intersect(eegleads,resleads); %find leads common to eeg and resected file leads
        loc=NaN.*ones(length(eegleads),1); %create vector of NaNs for resected
        loc(eegind)=import.data.(respts{i})(resind,2); %populate vector in order of the eeg leads from the resected leads
        eeg(ptind).loc=loc; %import into eeg structure
    end
    
    
end
