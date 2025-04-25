function [rCond,subList,runCondAcqList,runCondStimList] = set2cond5(runSet,rCond,sesPhys,volAnat)
if ~exist('sesPhys','var'); sesPhys = []; end
if ~exist('volAnat','var'); volAnat = []; end




%%
runSetTmp = cat(2,runSet{:})';
runCondTmp = cat(2,rCond{:})';
runCondTmp = [runCondTmp{:}]';

runCondTmpTmpTmp = [];

for rs = 1:length(runSetTmp)
    if isempty(runSetTmp{rs}.fList); continue; end

        
    ind = ismember({runCondTmp.sub}',runSetTmp{rs}.sub) & ...
        ismember({runCondTmp.ses}',runSetTmp{rs}.ses) & ...
        ismember(strcat('acq-',{runCondTmp.acq}','_prsc-',{runCondTmp.prsc}'),runSetTmp{rs}.label) & ...
        ~cellfun('isempty',{runCondTmp.fList})'; 

    runCondTmpTmp = runCondTmp(ind);
    for rc = 1:length(runCondTmpTmp)
        runCondTmpTmp(rc).wd           = runSetTmp{rs}.finalFiles.wd;
        runCondTmpTmp(rc).info         = runSetTmp{rs}.finalFiles.info;
        runCondTmpTmp(rc).ppLabelList  = runSetTmp{rs}.finalFiles.ppLabelList;
        runCondTmpTmp(rc).dataType     = runSetTmp{rs}.finalFiles.dataType;
        
        [~,condBids] = fileparts(replace(runCondTmpTmp(rc).fList(:,1),'.nii.gz',''));
        [~,setBids] = fileparts(fileparts(runSetTmp{rs}.finalFiles.fPreprocList(:,1)));
        [a,b] = ismember(cellstr(setBids),cellstr(condBids));

        fieldList = {'fPreprocList' 'fTransList' 'fTransCatList' 'bidsList' 'nFrame' 'nFrameOrig' 'vSize' 'acqTime' 'fOrigList'};
        for i = 1:length(fieldList)
            tmp = runSetTmp{rs}.finalFiles.(fieldList{i})(a,:,:);
            runCondTmpTmp(rc).(fieldList{i}) = tmp(b(a),:,:);
        end

        
        % add preproc mask
        fMaskList = {};
        tmpField = {'wrMocoFiles' 'brMocoFiles' 'bsMocoFiles'};
        for i = 1:length(tmpField)
            if isfield(runSetTmp{rs},tmpField{i})
                if strcmp(tmpField{i},'bsMocoFiles')
                    curfMask = repmat(runSetTmp{rs}.(tmpField{i}).fMaskList,size(runSetTmp{rs}.(tmpField{i}).fList,1),1);
                else
                    curfMask = runSetTmp{rs}.(tmpField{i}).fMaskList(a,:,:);
                end
            else
                curfMask = repmat({''},size(runSetTmp{rs}.finalFiles.fPreprocList,1),1);
            end
            fMaskList = cat(3,fMaskList,curfMask(b(a),:,:));
        end
        runCondTmpTmp(rc).fPreprocMaskList = fMaskList;
    end

    runCondTmpTmpTmp = cat(1,runCondTmpTmpTmp,runCondTmpTmp);
end

%% Extract venc from filename (ideally from dcm)
for rc = 1:length(runCondTmpTmpTmp)
    venc1 = strsplit(runCondTmpTmpTmp(rc).fPreprocList{1},'_');
    venc1 = venc1(contains(venc1,'acq-pcVenc'));
    if isempty(venc1)
        runCondTmpTmpTmp(rc).vencAcq = 'none';
    else
        venc2 = cell(1,size(runCondTmpTmpTmp(rc).fPreprocList,2));
        for i = 1:size(runCondTmpTmpTmp(rc).fPreprocList,2)
            venc = strsplit(runCondTmpTmpTmp(rc).fPreprocList{1,i},'_');
            venc2(i) = venc(contains(venc,'rec-'));
        end
        venc1 = replace(venc1,'acq-','');
        venc2 = replace(venc2,'rec-','');
        runCondTmpTmpTmp(rc).vencAcq = char(venc1);
        % runCondTmpTmpTmp(rc).vencRec = venc2;
    end
end



%% Display summary
[{runCondTmpTmpTmp.sub}
    {runCondTmpTmpTmp.ses}
    {runCondTmpTmpTmp.acq}
    {runCondTmpTmpTmp.prsc}
    {runCondTmpTmpTmp.vencAcq}
    {runCondTmpTmpTmp.task}
    cellstr(num2str(cellfun('size',{runCondTmpTmpTmp.fList},1)'))']';




tmp = ...
    [{runCondTmpTmpTmp.sub}
    {runCondTmpTmpTmp.ses}
    {runCondTmpTmpTmp.acq}
    {runCondTmpTmpTmp.prsc}
    {runCondTmpTmpTmp.vencAcq}
    {runCondTmpTmpTmp.task}
    cellstr(num2str(cellfun('size',{runCondTmpTmpTmp.fList},1)'))']';
% cellstr(num2str(~cellfun('isempty',{runCondTmpTmpTmp.physSes})'))']';

tmp = ...
    [{runCondTmpTmpTmp.sub}
    {runCondTmpTmpTmp.ses}
    {runCondTmpTmpTmp.acq}
    {runCondTmpTmpTmp.prsc}
    {runCondTmpTmpTmp.vencAcq}
    {runCondTmpTmpTmp.task}
    cellstr(num2str(cellfun('size',{runCondTmpTmpTmp.fList},1)'))']';
% nPhysRuns']';

tmp(cellfun('isempty',tmp)) = {'-'};

disp(table(char(tmp(~ismember(tmp(:,7),'0'),1)),char(tmp(~ismember(tmp(:,7),'0'),2)),char(tmp(~ismember(tmp(:,7),'0'),3)),char(tmp(~ismember(tmp(:,7),'0'),4)),char(tmp(~ismember(tmp(:,7),'0'),5)),char(tmp(~ismember(tmp(:,7),'0'),6)),char(tmp(~ismember(tmp(:,7),'0'),7)),'VariableNames',{'sub' 'ses' 'mriCond' 'prsc' 'vencAcq' 'stimCond' 'mriRuns'}))
% disp(...
%     [{'sub' 'ses' 'mriCond' 'stimCond' 'mriRuns'}
%     tmp(~ismember(tmp(:,5),'0'),:)]...
%     )
% disp(...
%     [{'sub' 'ses' 'mriCond' 'stimCond' 'mriRuns' 'physRunsInSes'}
%     tmp(~ismember(tmp(:,5),'0'),:)]...
%     )





%% 

subList = unique({runCondTmpTmpTmp.sub})';
rCond = cell(size(subList));
for S = 1:length(subList)
    % if S~=2; continue; end
    % sort subjects
    ind = ismember({runCondTmpTmpTmp.sub},subList{S});
    tmp = runCondTmpTmpTmp(ind);

    % if any(~ismember({tmp.vencAcq},'none'))
    %     keyboard
    % end

    

    % [{tmp.sub}' {tmp.ses}' {tmp.labelAcq}' {tmp.label}' cellstr(num2str(cellfun('size',{tmp.fList},1)'))]
    acqAll     = {tmp.acq};
    prscAll    = {tmp.prsc};
    vencAcqAll = {tmp.vencAcq};
    runCondAcqListAll = strcat(...
    strcat('acq-'    ,acqAll)'    ,'_',...
    strcat('prsc-'   ,prscAll)'   ,'_',...
    strcat('vencAcq-',vencAcqAll)');
    [runCondAcqList,b] = unique(runCondAcqListAll);
    acq     = acqAll(b);
    prsc    = prscAll(b);
    vencAcq = vencAcqAll(b);
    
    for ac = 1:length(runCondAcqList)
        % if ac~=2; continue; end
        % sort acquisition conditions
        ind = ismember(runCondAcqListAll,runCondAcqList(ac));
        tmp2 = tmp(ind);
        [runCondStimList,~,runCondStimInd] = unique({tmp2.task});

  

        for rsc = 1:length(runCondStimList)
            % if strcmp('task_10sPrd1sDur',runCondStimList{rsc}) && strcmp('vfMRI',runCondAcqList{ac}) && S==1
            %     keyboard
            % end
            


            
            % sort stimulus conditions
            tmp3 = tmp2(rsc==runCondStimInd);
            indX = false(size(tmp3));
            for i = 1:length(tmp3)
                indX(i) = isempty(tmp3(i).fList);
            end
            tmp3(indX) = [];
            if isempty(tmp3); continue; end


            % Trim down bidsList elements that are not common across
            % sessions
            bidsList = cell(size(tmp3));
            for i = 1:length(tmp3)
                bidsList{i} = permute(tmp3(i).bidsList(1,1,1:end-1),[1 3 2]);
                for ii = 1:length(bidsList{i})
                    bidsList{i}{ii} = strsplit(bidsList{i}{ii},'-');
                    bidsList{i}{ii} = bidsList{i}{ii}{1};
                end
            end
            for i = 2:length(bidsList)
                bidsList{1} = intersect(bidsList{1},bidsList{i});
            end
            bidsList = bidsList{1};
            for i = 1:length(bidsList)
                bidsList{i} = [bidsList{i} '-'];
            end
            for i = 1:length(tmp3)
                ind = contains(tmp3(i).bidsList(1,1,:),bidsList);
                tmp3(i).bidsList = cat(3,tmp3(i).bidsList(:,:,ind),tmp3(i).bidsList(:,:,end));
            end
            bidsList = squeeze(tmp3(1).bidsList(1,1,1:end-1))';
            for i = 1:length(bidsList)
                bidsList{i} = strsplit(bidsList{i},'-');
                bidsList{i} = [bidsList{i}{1} '-'];
            end
            for i = 1:length(tmp3)
                bidsListX = squeeze(tmp3(i).bidsList(1,1,1:end-1))';
                for ii = 1:length(bidsListX)
                    bidsListX{ii} = strsplit(bidsListX{ii},'-');
                    bidsListX{ii} = [bidsListX{ii}{1} '-'];
                end
                [~,b] = ismember(bidsListX,bidsList);
                tmp3(i).bidsList = tmp3(i).bidsList(:,:,[b end]);
            end




            %% Copy all relevant fields
            tmp4 = tmp3(1);
            tmp4.ses          = [];
            tmp4.wd           = {};
            tmp4.bidsDir      = {};
            tmp4.bidsDerivDir = {};
            for i = 1:length(tmp3)
                tmp4.ses          = cat(1,tmp4.ses         ,repmat(        tmp3(i).ses          ,size(tmp3(i).fPreprocList,1),1));
                tmp4.wd           = cat(1,tmp4.wd          ,repmat(cellstr(tmp3(i).wd          ),size(tmp3(i).fPreprocList,1),1));
                % tmp4.bidsDir      = cat(1,tmp4.bidsDir     ,repmat(cellstr(tmp3(i).bidsDir)     ,size(tmp3(i).fPreprocList)));
                % tmp4.bidsDerivDir = cat(1,tmp4.bidsDerivDir,repmat(cellstr(tmp3(i).bidsDerivDir),size(tmp3(i).fPreprocList)));
            end
            % fieldList = {'fList' 'fOrigList' 'fPreprocList' 'fTransList' 'fTransCatList' 'bidsList' 'nFrame' 'vSize' 'acqTime' 'bhvr' 'nDummy' 'fPreprocUnderSesCatRunCatAvList' 'fPreprocUnderSesAvCatRunCatAvList' 'fPreprocUnderSesCatRunAvCatAvList' 'fPreprocUnderSesAvCatRunAvCatAvList'};
            fieldList = {'fList' 'fOrigList' 'fPreprocList' 'fPreprocMaskList' 'fTransList' 'fTransCatList' 'bidsList' 'tr' 'nFrame' 'nFrameOrig' 'vSize' 'date' 'acqTime' 'bhvr'};
            % fieldList = {'fList' 'fOrigList' 'fPreprocList' 'fTransList' 'fTransCatList'            'nFrame' 'vSize' 'date' 'acqTime' 'bhvr' 'nDummy'};
            for i = 1:length(fieldList)
                    tmp4(1).(fieldList{i}) = cat(1,tmp3(:).(fieldList{i}));
            end
            

            %% Compile
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %SPECIAL CASE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(vencAcq(ac),'pcVenc7z')
                vencAcq{ac} = 'pcVenc7ap';
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %SPECIAL CASE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            acqCond = strjoin(...
                [acq(ac)
                 prsc(ac)
                 vencAcq(ac)],'_');
            rCond{S}.(acqCond).(['task_' runCondStimList{rsc}]) = tmp4;
            % rCond{S}.(runCondAcqList{ac}).(['task_' runCondStimList{rsc}]) = tmp4;
        end
        % rCond{S}.(runCondAcqList{ac}).prcSmr = tmp3.QA;
    end
end




%% Match volAnat to mri subjects
tmp = [volAnat{:}];
for i = 1:length(tmp)
    S = ismember(subList,tmp{i}.sub);
    switch tmp{i}.label
        case 'fs'
            if isempty(tmp{i}.fsDir)
                continue
            end
            if ~isfield(rCond{S},'fs')
                rCond{S}.fs = tmp{i};
            else
                rCond{S}.fs(end+1,1) = tmp{i};
            end
        case 'avMap'
            if isempty(tmp{i}.fList)
                continue
            end
            if ~isfield(rCond{S},'avMap')
                rCond{S}.avMap = tmp{i};
            else
                rCond{S}.avMap(end+1,1) = tmp{i};
            end
    end
end



%% Match physio sessions to mri subjects
sesPhysSub = cell(size(sesPhys));
sesPhysSes = cell(size(sesPhys));
for s = 1:length(sesPhys)
    if isempty(sesPhys{s})
        sesPhysSub{s} = '';
        sesPhysSes{s} = '';
    else
        sesPhysSub{s} = unique({sesPhys{s}.mriRuns.sub}');
        sesPhysSes{s} = unique({sesPhys{s}.mriRuns.ses}');
        {sesPhys{s}.mriRuns.task}';
        {sesPhys{s}.mriRuns.acq}';
        if length(sesPhysSub{s})>1; dbstack; keyboard; error('physio sessions are confused'); end
        if length(sesPhysSes{s})>1; dbstack; keyboard; error('physio sessions are confused'); end
        sesPhysSub(s) = sesPhysSub{s};
        sesPhysSes(s) = sesPhysSes{s};
    end
end

for S = 1:length(rCond)
    Sind = ismember(sesPhysSub,subList{S});
    if ~any(Sind)
        rCond{S}.phs = [];
        continue
    end
    rCond{S}.phs = [sesPhys{Sind}];
end


%% Extract physio runs and put it in runCond
for S = 1:length(rCond)
    if isempty(rCond{S}.phs); continue; end
    % if length(rCond{S}.phs)>1; dbstack; error('more than one physio session, code that'); end
    
    physRun = [];
    for p = 1:length(rCond{S}.phs)
        physRun = cat(1,physRun,physSes2run(rCond{S}.phs(p)));
    end
    physRunMriFile = [physRun.mri]; physRunMriFile = {physRunMriFile.fspec}';
    
    runCondAcqList = fields(rCond{S});
    runCondAcqList(ismember(runCondAcqList,{'phs' 'fs' 'avMap'})) = [];
    for rc = 1:length(runCondAcqList)
        runCondStimList = fields(rCond{S}.(runCondAcqList{rc}));
        runCondStimList(~contains(runCondStimList,'task_')) = [];
        for sc = 1:length(runCondStimList)
            mriFile = rCond{S}.(runCondAcqList{rc}).(runCondStimList{sc}).fList(:,1);
            physInd = ismember(physRunMriFile,mriFile);
            mriInd = ismember(mriFile,physRunMriFile);
            rCond{S}.(runCondAcqList{rc}).(runCondStimList{sc}).phs         = repmat(physSes2run,size(mriInd));
            rCond{S}.(runCondAcqList{rc}).(runCondStimList{sc}).phs(mriInd) = physRun(physInd);
        end
    end
end



%% Summarize
subList2      = {};
sesList2      = {};
acqCondList2  = {};
prscCondList2 = {};
vencCondList2 = {};
stimCondList2 = {};
n             = {};
nPhys         = {};
% nPhysSes      = {};
for S = 1:length(rCond)
    acqTime       = {};
    runCondAcqList = fields(rCond{S});
    runCondAcqList(ismember(runCondAcqList,{'phs' 'fs' 'avMap'})) = [];
    for ac = 1:length(runCondAcqList)
        runCondStimList = fields(rCond{S}.(runCondAcqList{ac}));
        runCondStimList(~contains(runCondStimList,'task_')) = [];
        for rsc = 1:length(runCondStimList)
            subList2{end+1}      = rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).sub;
            sesList2{end+1}      = rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).ses';
            acqCondList2{end+1}  = rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).acq;
            prscCondList2{end+1} = rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).prsc;
            vencCondList2{end+1} = rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).vencAcq;
            stimCondList2{end+1} = rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).task;
            n{end+1}             = size(rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).fPreprocList,1);
            acqTime{end+1}       = rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).acqTime;
            nPhys{end+1}         = 0;
            if ~isempty(rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).phs)
                nPhys{end}       = rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).ses(~[rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).phs.isempty])';
                % nPhys{end}       = nnz(~[rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).phs.isempty]);
            end
        end
    end
end



[~,b] = sort(acqCondList2);
acqCondList2 = acqCondList2(b);
prscCondList2 = prscCondList2(b);
vencCondList2 = vencCondList2(b);
subList2 = subList2(b);
sesList2 = sesList2(b);
stimCondList2 = stimCondList2(b);
n = n(b);
nPhys = nPhys(b);

[~,b] = sort(subList2);
acqCondList2 = acqCondList2(b);
prscCondList2 = prscCondList2(b);
vencCondList2 = vencCondList2(b);
subList2 = subList2(b);
sesList2 = sesList2(b);
stimCondList2 = stimCondList2(b);
n = n(b);
nPhys = nPhys(b);

disp([{'mriCond' 'sub' 'mriSes' 'prsc' 'venc' 'stimCond' 'nRun' 'physSes'}
      {'-------' '---' '------' '----' '----' '--------' '----' '-------'}
[acqCondList2
subList2
sesList2
prscCondList2
vencCondList2
stimCondList2
n
nPhys]'])

runCondStimList = strcat('task_',unique(stimCondList2)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPECIAL CASE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vencCondList2(ismember(vencCondList2,'pcVenc7z')) = {'pcVenc7ap'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPECIAL CASE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runCondAcqList  = unique(strcat(acqCondList2','_',prscCondList2','_',vencCondList2'));
% runCondPrscList = unique(prscCondList2)';
% runCondVencList = unique(vencCondList2)';
