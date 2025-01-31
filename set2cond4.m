function [rCond,subList,runCondAcqList,runCondStimList] = set2cond4(runSet,rCond,sesPhys,volAnat)
force = 0;
if ~exist('sesPhys','var'); sesPhys = []; end
if ~exist('volAnat','var'); volAnat = []; end

% %% Insert phys data to mri (runCond) data
% for s = 1:length(runCond)
%     if ~isempty(sesPhys{s}) && ~isfield(sesPhys{s},'mriRuns'); dbstack; warning('can''t insert phys data into mri data'); end
%     for rc = 1:length(runCond{s})
%         % Skip if no phys data for that session
%         if isempty(sesPhys{s})
%             runCond{s}{rc}.physSes = [];
%             continue
%         end
% 
%         % find phys runs corresponding to current mri runs
%         ind = ...
%             ismember({sesPhys{s}.mriRuns.sub}'     ,runCond{s}{rc}.sub) & ...
%             ismember({sesPhys{s}.mriRuns.ses}'     ,runCond{s}{rc}.ses) & ...
%             ismember({sesPhys{s}.mriRuns.label}'   ,runCond{s}{rc}.label) & ...
%             ismember({sesPhys{s}.mriRuns.labelAcq}',runCond{s}{rc}.labelAcq);
% 
%         % match individual runs
%         if length(runCond{s}{rc}.fList) ~= length([sesPhys{s}.mriRuns(ind).fList]'); dbstack; error('number of phys runs doesn''t match number of mri runs...'); end
%         if isempty(runCond{s}{rc}.fList)
%             runCond{s}{rc}.physSes = [];
%             continue
%         end
%         [~,b] = ismember(runCond{s}{rc}.fList,[sesPhys{s}.mriRuns(ind).fList]');
%         if ~all(b'==1:nnz(ind)); dbstack; error('order of phys runs and mri runs doesn''t match. CODE THAT'); end
% 
%         % insert phys to runCond
%         runCond{s}{rc}.physSes = sesPhys{s};
%         % runCond{s}{rc}.phys = rmfield(sesPhys{s},{'data' 'mriRunFiles' 'mriRunTimes' 'mriRuns' 'runTimes' 'info'});
%         % runCond{s}{rc}.phys.mriRuns = sesPhys{s}.mriRuns(ind);
%         % runCond{s}{rc}.phys.mriRunTimes = sesPhys{s}.mriRunTimes(ind);
%         % runCond{s}{rc}.phys.mriRunFiles = sesPhys{s}.mriRunFiles(ind);
%         % runCond{s}{rc}.phys.runTimes = sesPhys{s}.runTimes(ind);
%         % runCond{s}{rc}.phys.info = rmfield(sesPhys{s}.info,'mri');
%         % runCond{s}{rc}.phys.info = sesPhys{s}.info.mri(ind);
%     end
% end


% 
% 
% if ~isempty(sesPhys)
%     % sesPhys is a vector of sessions while runCond is a vector of subject.
%     % First, insert phys sessions into corresponding subjects
%     runCondSub = cell(size(runCond));
%     for s = 1:length(runCond)
%         tmp = fields(runCond{s});
%         tmp2 = fields(runCond{s}.(tmp{1}));
%         runCondSub{s} = runCond{s}.(tmp{1}).(tmp2{1}).sub;
%     end
% 
%     sesPhysSub = cell(size(sesPhys));
%     for s = 1:length(sesPhys)
%         if ~isempty(sesPhys{s})
%             sub = unique({sesPhys{s}.mriRuns.sub});
%             if length(sub)~=1; dbstack; error('x'); end
%         else
%             sub = {''};
%         end
%         sesPhysSub(s) = sub ;
%     end
% 
%     for s = 1:length(runCond)
%         runCond{s}.sesPhys = [sesPhys{ismember(sesPhysSub,runCondSub{s})}]';
%     end
% 
%     % Then for each subject, insert physio data to the corresponding
%     % session and run
%     for s = 1:length(runCond)
%         mri = phys2mri(runCond{s});
%         sesPhys
% 
%     end
% 
% end








%%
runSetTmp = cat(2,runSet{:})';
runCondTmp = cat(2,rCond{:})';
runCondTmp = [runCondTmp{:}]';

runCondTmpTmpTmp = [];

for rs = 1:length(runSetTmp)
    if isempty(runSetTmp{rs}.fList); continue; end

    ind = ismember({runCondTmp.sub}',runSetTmp{rs}.sub) & ...
        ismember({runCondTmp.ses}',runSetTmp{rs}.ses) & ...
        ismember({runCondTmp.acq}',runSetTmp{rs}.label) & ...
        ~cellfun('isempty',{runCondTmp.fList})'; 

    runCondTmpTmp = runCondTmp(ind);
    for rc = 1:length(runCondTmpTmp)
        runCondTmpTmp(rc).wd           = runSetTmp{rs}.finalFiles.wd;
        runCondTmpTmp(rc).info         = runSetTmp{rs}.finalFiles.info;
        % runCondTmpTmp(rc).bidsDir      = runSetTmp{rs}.finalFiles.bidsDir;
        % runCondTmpTmp(rc).bidsDerivDir = runSetTmp{rs}.finalFiles.bidsDerivDir;
        runCondTmpTmp(rc).ppLabelList  = runSetTmp{rs}.finalFiles.ppLabelList;
        runCondTmpTmp(rc).dataType     = runSetTmp{rs}.finalFiles.dataType;
        % if isfield(runSetTmp{rs},'avMap')
        %     runCondTmpTmp(rc).avMap = runSetTmp{rs}.avMap;
        % else
        %     runCondTmpTmp(rc).avMap = [];
        % end

        try
            [~,condBids] = fileparts(replace(runCondTmpTmp(rc).fList(:,1),'.nii.gz',''));
        catch
            keyboard
        end
        [~,setBids] = fileparts(fileparts(runSetTmp{rs}.finalFiles.fPreprocList(:,1)));
        [a,b] = ismember(cellstr(setBids),cellstr(condBids));

        fieldList = {'fPreprocList' 'fTransList' 'fTransCatList' 'bidsList' 'nFrame' 'nFrameOrig' 'vSize' 'acqTime' 'fOrigList' 'nDummy'};
        for i = 1:length(fieldList)
            tmp = runSetTmp{rs}.finalFiles.(fieldList{i})(a,:,:);
            runCondTmpTmp(rc).(fieldList{i}) = tmp(b(a),:,:);
        end

        %QA and summaries
        % runCondTmpTmp(rc).QA     = runSetTmp{rs}.QA;
        % runCondTmpTmp(rc).QA.smr = runSetTmp{rs}.finalFiles.fPreprocSmr;


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
                % curfMask = {''};
                curfMask = repmat({''},size(runSetTmp{rs}.finalFiles.fPreprocList,1),1);
            end
            fMaskList = cat(3,fMaskList,curfMask(b(a),:,:));
        end
        runCondTmpTmp(rc).fPreprocMaskList = fMaskList;
    end

    runCondTmpTmpTmp = cat(1,runCondTmpTmpTmp,runCondTmpTmp);
end


%% Display summary
[{runCondTmpTmpTmp.sub}
    {runCondTmpTmpTmp.ses}
    {runCondTmpTmpTmp.acq}
    {runCondTmpTmpTmp.task}
    cellstr(num2str(cellfun('size',{runCondTmpTmpTmp.fList},1)'))']';


% nPhysRuns = cell(size(runCondTmpTmpTmp));
% for i = 1:length(runCondTmpTmpTmp)
%     if isempty(runCondTmpTmpTmp(i).physSes)
%         nPhysRuns{i} = num2str(0);
%     else
%         nPhysRuns{i} = num2str(length(runCondTmpTmpTmp(i).physSes.mriRuns));
%     end
% end


tmp = ...
    [{runCondTmpTmpTmp.sub}
    {runCondTmpTmpTmp.ses}
    {runCondTmpTmpTmp.acq}
    {runCondTmpTmpTmp.task}
    cellstr(num2str(cellfun('size',{runCondTmpTmpTmp.fList},1)'))']';
% cellstr(num2str(~cellfun('isempty',{runCondTmpTmpTmp.physSes})'))']';

tmp = ...
    [{runCondTmpTmpTmp.sub}
    {runCondTmpTmpTmp.ses}
    {runCondTmpTmpTmp.acq}
    {runCondTmpTmpTmp.task}
    cellstr(num2str(cellfun('size',{runCondTmpTmpTmp.fList},1)'))']';
% nPhysRuns']';

disp(table(tmp(~ismember(tmp(:,5),'0'),1),tmp(~ismember(tmp(:,5),'0'),2),tmp(~ismember(tmp(:,5),'0'),3),tmp(~ismember(tmp(:,5),'0'),4),tmp(~ismember(tmp(:,5),'0'),5),'VariableNames',{'sub' 'ses' 'mriCond' 'stimCond' 'mriRuns'}))
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

    % [{tmp.sub}' {tmp.ses}' {tmp.labelAcq}' {tmp.label}' cellstr(num2str(cellfun('size',{tmp.fList},1)'))]

    runCondAcqList = unique({tmp.acq});
    for ac = 1:length(runCondAcqList)
        % if ac~=2; continue; end
        % sort acquisition conditions
        ind = ismember({tmp.acq},runCondAcqList(ac));
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
            fieldList = {'fList' 'fOrigList' 'fPreprocList' 'fPreprocMaskList' 'fTransList' 'fTransCatList' 'bidsList' 'nFrame' 'nFrameOrig' 'vSize' 'date' 'acqTime' 'bhvr' 'nDummy'};
            % fieldList = {'fList' 'fOrigList' 'fPreprocList' 'fTransList' 'fTransCatList'            'nFrame' 'vSize' 'date' 'acqTime' 'bhvr' 'nDummy'};
            for i = 1:length(fieldList)
                % try
                    tmp4(1).(fieldList{i}) = cat(1,tmp3(:).(fieldList{i}));
                % catch
                %     sz = [0 max(cellfun('size',{tmp3(:).(fieldList{i})},2))];
                %     tmpX = cell(sz);
                %     for ii = 1:numel(tmp3)
                %         tmpXX = tmp3(ii).(fieldList{i});
                %         tmpXX = cat(2,tmpXX,repmat({''},[size(tmpXX,1) sz(2) - size(tmpXX,2)]));
                %         tmpX  = cat(1,tmpX,tmpXX);
                %     end
                %     tmp4.(fieldList{i}) = tmpX;
                % end
            end
            


            % % % % % %% Define underlays (catenate and average across sessions)
            % % % % % disp(['writing underlays for subj' num2str(S) '/' num2str(length(subList)) ', acqCond' num2str(ac) '/' num2str(length(runCondAcqList)) ', stimCond' num2str(rsc) '/' num2str(length(runCondStimList))])
            % % % % % 
            % % % % % %%% across runs, within stimulus conditions, within sessions
            % % % % % [bidsDerivDir,~,bidsDerivDirInd] = unique(tmp4.bidsDerivDir);
            % % % % % fSesRunCatAv   = cell(size(bidsDerivDir));
            % % % % % fSesRunAvCatAv = cell(size(bidsDerivDir));
            % % % % % for s = 1:length(bidsDerivDir)
            % % % % %     fList = dir(char(fullfile(bidsDerivDir{s},['*task-' runCondStimList{rsc} '*'],'av_preproc_volTs.nii.gz')));
            % % % % %     fList = fullfile({fList.folder},{fList.name})';
            % % % % %     mri   = cell(size(fList));
            % % % % %     for i = 1:length(fList)
            % % % % %         mri{i} = MRIread(fList{i});
            % % % % %     end
            % % % % %     mri = cat(1,mri{:});
            % % % % %     mriSesRunCatAv(s) = mri(1);
            % % % % %     mriSesRunCatAv(s).vol = cat(4,mri.vol); clear mri
            % % % % %     fSesRunCatAv{s} = fullfile(bidsDerivDir{s},['task-' runCondStimList{rsc} '_cat_av_preproc_volTs.nii.gz']);
            % % % % %     if force || ~exist(fSesRunCatAv{s},'file')
            % % % % %         MRIwrite(mriSesRunCatAv(s),fSesRunCatAv{s});
            % % % % %     end
            % % % % %     mriSesRunAvCatAv(s) = mriSesRunCatAv(s);
            % % % % %     mriSesRunAvCatAv(s).vol = mean(mriSesRunAvCatAv(s).vol,4);
            % % % % %     fSesRunAvCatAv{s} = fullfile(bidsDerivDir{s},['task-' runCondStimList{rsc} '_av_cat_av_preproc_volTs.nii.gz']);
            % % % % %     if force || ~exist(fSesRunAvCatAv{s},'file')
            % % % % %         MRIwrite(mriSesRunAvCatAv(s),fSesRunAvCatAv{s});
            % % % % %     end
            % % % % % end
            % % % % % 
            % % % % % %%% across sessions
            % % % % % % catenate sessions then average runs
            % % % % % [a,b,c] = fileparts(fSesRunCatAv); b = strcat('sesCat_',b);
            % % % % % fSesCatRunCatAv = cellstr(fullfile(a,strcat(b,c)));
            % % % % % mriSesCatRunCatAv   = mriSesRunCatAv;
            % % % % % [a,b,c] = fileparts(fSesRunCatAv); b = strcat('sesAvCat_',b);
            % % % % % fSesAvCatRunCatAv = cellstr(fullfile(a,strcat(b,c)));
            % % % % % mriSesAvCatRunCatAv = mriSesRunCatAv;
            % % % % % for i = 1:length(fSesRunCatAv)
            % % % % %     mriSesCatRunCatAv(i).vol = cat(4,mriSesRunCatAv.vol);
            % % % % %     if force || ~exist(fSesCatRunCatAv{i},'file')
            % % % % %         MRIwrite(mriSesCatRunCatAv(i),fSesCatRunCatAv{i});
            % % % % %     end
            % % % % %     mriSesAvCatRunCatAv(i).vol = mean(mriSesCatRunCatAv(i).vol,4);
            % % % % %     if force || ~exist(fSesAvCatRunCatAv{i},'file')
            % % % % %         MRIwrite(mriSesAvCatRunCatAv(i),fSesAvCatRunCatAv{i});
            % % % % %     end
            % % % % % end
            % % % % % tmp4.fPreprocUnderCondSpecSesCatRunCatAvList   = fSesCatRunCatAv(bidsDerivDirInd);
            % % % % % tmp4.fPreprocUnderCondSpecSesAvCatRunCatAvList = fSesAvCatRunCatAv(bidsDerivDirInd);
            % % % % % 
            % % % % % % average runs then catenate sessions
            % % % % % [a,b,c] = fileparts(fSesRunAvCatAv); b = strcat('sesCat_',b);
            % % % % % fSesCatRunAvCatAv = cellstr(fullfile(a,strcat(b,c)));
            % % % % % mriSesCatRunAvCatAv   = mriSesRunAvCatAv;
            % % % % % [a,b,c] = fileparts(fSesRunAvCatAv); b = strcat('sesAvCat_',b);
            % % % % % fSesAvCatRunAvCatAv = cellstr(fullfile(a,strcat(b,c)));
            % % % % % mriSesAvCatRunAvCatAv = mriSesRunAvCatAv;
            % % % % % for i = 1:length(fSesRunAvCatAv)
            % % % % %     mriSesCatRunAvCatAv(i).vol = cat(4,mriSesRunAvCatAv.vol);
            % % % % %     if force || ~exist(fSesCatRunAvCatAv{i},'file')
            % % % % %         MRIwrite(mriSesCatRunAvCatAv(i),fSesCatRunAvCatAv{i});
            % % % % %     end
            % % % % %     mriSesAvCatRunAvCatAv(i).vol = mean(mriSesCatRunAvCatAv(i).vol,4);
            % % % % %     if force || ~exist(fSesAvCatRunAvCatAv{i},'file')
            % % % % %         MRIwrite(mriSesAvCatRunAvCatAv(i),fSesAvCatRunAvCatAv{i});
            % % % % %     end
            % % % % % end
            % % % % % tmp4.fPreprocUnderCondSpecSesCatRunAvCatAvList   = fSesCatRunAvCatAv(bidsDerivDirInd);
            % % % % % tmp4.fPreprocUnderCondSpecSesAvCatRunAvCatAvList = fSesAvCatRunAvCatAv(bidsDerivDirInd);

            %% Compile
            rCond{S}.(runCondAcqList{ac}).(['task_' runCondStimList{rsc}]) = tmp4;
        end
        % rCond{S}.(runCondAcqList{ac}).prcSmr = tmp3.QA;
    end
end

% for S = 1:length(subList)
%     runCondAcqList = fields(runCond{S});
%     for ac = 1:length(runCondAcqList)
%         runCondStimList = fields(runCond{S}.(runCondAcqList{ac}));
%         for rsc = 1:length(runCondStimList)
%             length(runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).fList) ~= length(runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).fPreprocList)
%             % if length(runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).fList) ~= length(runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).fPreprocList)
%             %     keyboard
%             %     runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc})
%             % end
%         end
%     end
% end


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
    % for s = 1:length(rCond{S}.phs)
    %     rCond{S}.phs(s).sub = sesPhysSub{Sind};
    %     rCond{S}.phs(s).info.sub = sesPhysSub{Sind};
    %     rCond{S}.phs(s).ses = sesPhysSes{Sind};
    %     rCond{S}.phs(s).info.ses = sesPhysSes{Sind};
    % end
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


% %% Anonymize
% for S = 1:length(runCond)
%     runCondAcqList = fields(runCond{S});
%     for ac = 1:length(runCondAcqList)
%         runCondStimList = fields(runCond{S}.(runCondAcqList{ac}));
%         for rsc = 1:length(runCondStimList)
%             if ~any([runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).phys.isempty]); keyboard; end
%         end
%     end
% end

% 
% %% Phys
% if ~isempty(sesPhys)
%     % sesPhys is a vector of sessions while runCond is a vector of subject.
%     % First, insert phys sessions into corresponding subjects
%     runCondSub = cell(size(runCond));
%     for s = 1:length(runCond)
%         tmp = fields(runCond{s});
%         tmp2 = fields(runCond{s}.(tmp{1}));
%         runCondSub{s} = runCond{s}.(tmp{1}).(tmp2{1}).sub;
%     end
% 
%     sesPhysSub = cell(size(sesPhys));
%     for s = 1:length(sesPhys)
%         if ~isempty(sesPhys{s})
%             sub = unique({sesPhys{s}.mriRuns.sub});
%             if length(sub)~=1; dbstack; error('x'); end
%         else
%             sub = {''};
%         end
%         sesPhysSub(s) = sub ;
%     end
% 
%     for s = 1:length(runCond)
%         runCond{s}.sesPhys = [sesPhys{ismember(sesPhysSub,runCondSub{s})}]';
%     end
% 
%     % Then for each subject, insert physio data to the corresponding
%     % session and run
%     for s = 1:length(runCond)
%         mri = phys2mri(runCond{s});
%         sesPhys
% 
%     end
% 
% end


%% Summarize (and anonymize)
subList2      = {};
sesList2      = {};
acqCondList2  = {};
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
            stimCondList2{end+1} = rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).task;
            n{end+1}             = size(rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).fPreprocList,1);
            acqTime{end+1}       = rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).acqTime;
            nPhys{end+1}         = 0;
            if ~isempty(rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).phs)
                nPhys{end}       = rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).ses(~[rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).phs.isempty])';
                % nPhys{end}       = nnz(~[rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).phs.isempty]);
            end
            % if nPhys{end}~=0
            %     ses = [rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).phs.info];
            %     nPhysSes{end+1} = [ses.ses];
            % else
            %     nPhysSes{end+1} = '';
            % end
        end
    end
    % % %anonymize
    % % firstAcqTime = min(cat(1,acqTime{:}));
    % % runCondAcqList = fields(rCond{S});
    % % runCondAcqList(ismember(runCondAcqList,{'phs' 'fs' 'avMap'})) = [];
    % % for ac = 1:length(runCondAcqList)
    % %     runCondStimList = fields(rCond{S}.(runCondAcqList{ac}));
    % %     for rsc = 1:length(runCondStimList)
    % %         rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).acqTime = rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).acqTime - firstAcqTime;
    % %         if isfield(rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}),'date')
    % %             rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}) = rmfield(rCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}),'date');
    % %         end
    % %     end
    % % end
end


% %% seperate phys runs
% for s = 1:length(runCond)
%     mriCondList = fields(runCond{s});
%     for mc = 1:length(mriCondList)
%         stimCondList = fields(runCond{s}.(mriCondList{mc}));
%         runDur = zeros(0,2);
%         sc = 1;
%         try
%             nChan = size(runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physSes.titles,1);
%         catch
%             keyboard;
%         end
%         Fs     = zeros(nChan,0);
%         for sc = 1:length(stimCondList)
%             if isfield(runCond{s}.(mriCondList{mc}).(stimCondList{sc}),'physSes') && ~isempty(runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physSes)
%                 runDur = cat(1,runDur,runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physSes.mriRunTimes);
%                 try
%                     segmentInd = runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physSes.info.segmentInd;
%                     Fs     = cat(2,Fs,runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physSes.samplerate(:,segmentInd));
%                 catch
%                     keyboard
%                 end
%             end
%         end
%         if isempty(runDur); continue; end
%         % Grrr, not always the same number of data points...
%         runDur = max(diff(runDur,[],2));
%         Fs    = unique(Fs); if length(Fs)~=1; dbstack; error('X'); end
%         nRun  = ceil(runDur*Fs);
%         for sc = 1:length(stimCondList)
%             if ~isfield(runCond{s}.(mriCondList{mc}).(stimCondList{sc}),'physSes');   runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physSes = []; end
%             if ~isfield(runCond{s}.(mriCondList{mc}).(stimCondList{sc}),'physRuns'); runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physRuns = []; end
%             if isempty(runCond{s}.(mriCondList{mc}).(stimCondList{sc}).fList); continue; end
%             if isempty(runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physSes); continue; end
% 
%             curPhysSes = runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physSes;
%             tPhys = sum([hour(curPhysSes.blocktimes)*60*60 minute(curPhysSes.blocktimes)*60 second(curPhysSes.blocktimes)]);
% 
%             % find phys runs corresponding to current mri runs
%             ind = ...
%                 ismember({curPhysSes.mriRuns.sub}'      ,runCond{s}.(mriCondList{mc}).(stimCondList{sc}).sub) & ...
%                 ismember({curPhysSes.mriRuns.ses}'      ,runCond{s}.(mriCondList{mc}).(stimCondList{sc}).ses) & ...
%                 ismember({curPhysSes.mriRuns.label}'    ,runCond{s}.(mriCondList{mc}).(stimCondList{sc}).label) & ...
%                 ismember({curPhysSes.mriRuns.labelAcq}' ,runCond{s}.(mriCondList{mc}).(stimCondList{sc}).labelAcq);
% 
%             % match individual runs
%             if length(runCond{s}.(mriCondList{mc}).(stimCondList{sc}).fList) ~= length([curPhysSes.mriRuns(ind).fList]'); dbstack; error('number of phys runs doesn''t match number of mri runs...'); end
%             [~,b] = ismember(runCond{s}.(mriCondList{mc}).(stimCondList{sc}).fList,[curPhysSes.mriRuns(ind).fList]');
%             if ~all(b'==1:nnz(ind)); dbstack; error('order of phys runs and mri runs doesn''t match. CODE THAT'); end
% 
% 
%             % insert phys to runCond
%             runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physSes  = curPhysSes;
% 
%             curPhysSes.data = [];
%             curPhysSes.mriRuns(~ind)     = [];
%             curPhysSes.mriRunTimes(~ind,:) = [];
%             curPhysSes.mriRunFiles(~ind) = [];
%             % curPhysSes.runTimes(~ind)    = [];
%             curPhysSes.info.mri(~ind)    = [];
% 
%             runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physRuns = curPhysSes;
% 
%             % extract data from relevant runs only
%             data  = runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physSes.data;
% 
%             segmentInd = curPhysSes.info.segmentInd;
%             if segmentInd~=1; dbstack; error('physio data segment is not the first one, time will be messed up. Code that'); end
%             chanList = {'cardiac' 'resp' 'trigger'}; chanInd = ismember(chanList,curPhysSes.titles1); chanList = chanList(chanInd);
% 
%             runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physRuns.vec = zeros(nRun,length(chanList),1,length(curPhysSes.mriRuns));
%             runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physRuns.t   = zeros(nRun,1               ,1,length(curPhysSes.mriRuns));
%             runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physRuns.vecInfo = strjoin({'time/freq' 'vox/chan' 'taper/mode' 'run'},' x ');
%             runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physRuns.vecLabel = chanList;
% 
%             for chan = 1:length(chanList)
%                 chanInd = ismember(curPhysSes.titles1,chanList(chan));
%                 curData = data(curPhysSes.datastart(chanInd,segmentInd):curPhysSes.dataend(chanInd,segmentInd));
%                 t = (0:1/Fs:(size(curData,2)-1)/Fs) + tPhys;
% 
%                 for r = 1:length(curPhysSes.mriRuns)
%                     iS = find(t>curPhysSes.mriRunTimes(r,1),1,'first');
%                     iE = iS+nRun-1;
%                     runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physRuns.vec(:,chan,1,r) = curData(iS:iE);
%                     runCond{s}.(mriCondList{mc}).(stimCondList{sc}).physRuns.t(:,1,1,r)      = t(iS:iE);
%                 end
%             end
%         end
%     end
% end




[~,b] = sort(acqCondList2);
acqCondList2 = acqCondList2(b);
subList2 = subList2(b);
sesList2 = sesList2(b);
stimCondList2 = stimCondList2(b);
n = n(b);
nPhys = nPhys(b);

[~,b] = sort(subList2);
acqCondList2 = acqCondList2(b);
subList2 = subList2(b);
sesList2 = sesList2(b);
stimCondList2 = stimCondList2(b);
n = n(b);
nPhys = nPhys(b);

disp([{'mriCond' 'sub' 'mriSes' 'stimCond' 'nRun' 'physSes'}
      {'-------' '---' '------' '--------' '----' '-------'}
[acqCondList2
subList2
sesList2
stimCondList2
n
nPhys]'])

runCondStimList = strcat('task_',unique(stimCondList2)');
runCondAcqList  = unique(acqCondList2)';
