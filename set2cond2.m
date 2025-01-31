function [runCond,subList,runCondAcqList,runCondStimList] = set2cond2(runSet,runCond,sesPhys,volAnat)
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
runCondTmp = cat(2,runCond{:})';
runCondTmp = [runCondTmp{:}]';

runCondTmpTmpTmp = [];

for rs = 1:length(runSetTmp)
    if isempty(runSetTmp{rs}.fList); continue; end

    if isfield(runCondTmp,'labelAcq')
        ind = ismember({runCondTmp.sub}',runSetTmp{rs}.sub) & ...
            ismember({runCondTmp.ses}',runSetTmp{rs}.ses) & ...
            ismember({runCondTmp.labelAcq}',runSetTmp{rs}.label);
    else
        ind = ismember({runCondTmp.sub}',runSetTmp{rs}.sub) & ...
            ismember({runCondTmp.ses}',runSetTmp{rs}.ses) & ...
            ismember({runCondTmp.acq}',runSetTmp{rs}.label);
    end

    runCondTmpTmp = runCondTmp(ind);
    for rc = 1:length(runCondTmpTmp)
        runCondTmpTmp(rc).wd           = runSetTmp{rs}.finalFiles.wd;
        runCondTmpTmp(rc).bidsDir      = runSetTmp{rs}.finalFiles.bidsDir;
        runCondTmpTmp(rc).bidsDerivDir = runSetTmp{rs}.finalFiles.bidsDerivDir;
        runCondTmpTmp(rc).ppLabelList  = runSetTmp{rs}.finalFiles.ppLabelList;
        runCondTmpTmp(rc).dataType     = runSetTmp{rs}.finalFiles.dataType;
        % if isfield(runSetTmp{rs},'avMap')
        %     runCondTmpTmp(rc).avMap = runSetTmp{rs}.avMap;
        % else
        %     runCondTmpTmp(rc).avMap = [];
        % end

        [~,condBids] = fileparts(replace(runCondTmpTmp(rc).fList,'.nii.gz',''));
        [~,setBids] = fileparts(fileparts(runSetTmp{rs}.finalFiles.fPreprocList));
        [a,b] = ismember(cellstr(setBids),cellstr(condBids));

        fieldList = {'fPreprocList' 'fTransList' 'fTransCatList' 'bidsList' 'nFrame' 'vSize' 'acqTime' 'fOrigList' 'nDummy'};
        for i = 1:length(fieldList)
            tmp = runSetTmp{rs}.finalFiles.(fieldList{i})(a,:);
            runCondTmpTmp(rc).(fieldList{i}) = tmp(b(a),:);
        end
    end

    runCondTmpTmpTmp = cat(1,runCondTmpTmpTmp,runCondTmpTmp);
end


%% Display summary
if isfield(runCondTmp,'labelAcq') && isfield(runCondTmp,'label')
    [{runCondTmpTmpTmp.sub}
        {runCondTmpTmpTmp.ses}
        {runCondTmpTmpTmp.labelAcq}
        {runCondTmpTmpTmp.label}
        cellstr(num2str(cellfun('size',{runCondTmpTmpTmp.fList},1)'))']';
else
    [{runCondTmpTmpTmp.sub}
        {runCondTmpTmpTmp.ses}
        {runCondTmpTmpTmp.acq}
        {runCondTmpTmpTmp.stim}
        cellstr(num2str(cellfun('size',{runCondTmpTmpTmp.fList},1)'))']';
end


% nPhysRuns = cell(size(runCondTmpTmpTmp));
% for i = 1:length(runCondTmpTmpTmp)
%     if isempty(runCondTmpTmpTmp(i).physSes)
%         nPhysRuns{i} = num2str(0);
%     else
%         nPhysRuns{i} = num2str(length(runCondTmpTmpTmp(i).physSes.mriRuns));
%     end
% end

if isfield(runCondTmp,'labelAcq') && isfield(runCondTmp,'label')
    tmp = ...
        [{runCondTmpTmpTmp.sub}
        {runCondTmpTmpTmp.ses}
        {runCondTmpTmpTmp.labelAcq}
        {runCondTmpTmpTmp.label}
        cellstr(num2str(cellfun('size',{runCondTmpTmpTmp.fList},1)'))']';
    % cellstr(num2str(~cellfun('isempty',{runCondTmpTmpTmp.physSes})'))']';

    tmp = ...
        [{runCondTmpTmpTmp.sub}
        {runCondTmpTmpTmp.ses}
        {runCondTmpTmpTmp.labelAcq}
        {runCondTmpTmpTmp.label}
        cellstr(num2str(cellfun('size',{runCondTmpTmpTmp.fList},1)'))']';
    % nPhysRuns']';

    disp(...
        [{'sub' 'ses' 'mriCond' 'stimCond' 'mriRuns'}
        tmp(~ismember(tmp(:,5),'0'),:)]...
        )
else
    tmp = ...
        [{runCondTmpTmpTmp.sub}
        {runCondTmpTmpTmp.ses}
        {runCondTmpTmpTmp.acq}
        {runCondTmpTmpTmp.stim}
        cellstr(num2str(cellfun('size',{runCondTmpTmpTmp.fList},1)'))']';
    % cellstr(num2str(~cellfun('isempty',{runCondTmpTmpTmp.physSes})'))']';

    tmp = ...
        [{runCondTmpTmpTmp.sub}
        {runCondTmpTmpTmp.ses}
        {runCondTmpTmpTmp.acq}
        {runCondTmpTmpTmp.stim}
        cellstr(num2str(cellfun('size',{runCondTmpTmpTmp.fList},1)'))']';
    % nPhysRuns']';

    disp(...
        [{'sub' 'ses' 'mriCond' 'stimCond' 'mriRuns'}
        tmp(~ismember(tmp(:,5),'0'),:)]...
        )
end
% disp(...
%     [{'sub' 'ses' 'mriCond' 'stimCond' 'mriRuns' 'physRunsInSes'}
%     tmp(~ismember(tmp(:,5),'0'),:)]...
%     )





%%  

subList = unique({runCondTmpTmpTmp.sub})';
runCond = cell(size(subList));
for S = 1:length(subList)
    % if S~=2; continue; end
    % sort subjects
    ind = ismember({runCondTmpTmpTmp.sub},subList{S});
    tmp = runCondTmpTmpTmp(ind);

    % [{tmp.sub}' {tmp.ses}' {tmp.labelAcq}' {tmp.label}' cellstr(num2str(cellfun('size',{tmp.fList},1)'))]

    if isfield(runCondTmp,'labelAcq') && isfield(runCondTmp,'label')
        runCondAcqList = unique({tmp.labelAcq});
    else
        runCondAcqList = unique({tmp.acq});
    end
    for ac = 1:length(runCondAcqList)
        % if ac~=2; continue; end
        % sort acquisition conditions
        if isfield(runCondTmp,'labelAcq') && isfield(runCondTmp,'label')
            ind = ismember({tmp.labelAcq},runCondAcqList(ac));
            tmp2 = tmp(ind);
            [runCondStimList,~,runCondStimInd] = unique({tmp2.label});
        else
            ind = ismember({tmp.acq},runCondAcqList(ac));
            tmp2 = tmp(ind);
            [runCondStimList,~,runCondStimInd] = unique({tmp2.stim});
        end

        %% Define underlays (catenate and average across sessions)
        disp(['writing underlays for subj' num2str(S) '/' num2str(length(subList)) ', acqCond' num2str(ac) '/' num2str(length(runCondAcqList)) ', stimCondCombined'])
        %%% across sessions
        [bidsDerivDir,~,bidsDerivDirInd] = unique({tmp2.bidsDerivDir});
        fSesRunCatAv     = cell(size(bidsDerivDir));
        mriSesRunCatAv   = cell(size(bidsDerivDir));
        fSesRunAvCatAv   = cell(size(bidsDerivDir));
        mriSesRunAvCatAv = cell(size(bidsDerivDir));
        for s = 1:length(bidsDerivDir)
            fSesRunCatAv{s} = fullfile(bidsDerivDir{s},'cat_av_preproc_volTs.nii.gz');
            mriSesRunCatAv{s} = MRIread(fSesRunCatAv{s});
            fSesRunAvCatAv{s} = fullfile(bidsDerivDir{s},'av_cat_av_preproc_volTs.nii.gz');
            mriSesRunAvCatAv{s} = MRIread(fSesRunAvCatAv{s});
        end
        mriSesRunCatAv   = cat(1,mriSesRunCatAv{:});
        mriSesRunAvCatAv = cat(1,mriSesRunAvCatAv{:});

        % catenate sessions then average runs
        [a,b,c] = fileparts(fSesRunCatAv); b = strcat('sesCat_',b);
        fSesCatRunCatAv = cellstr(fullfile(a,strcat(b,c)));
        mriSesCatRunCatAv   = mriSesRunCatAv;
        [a,b,c] = fileparts(fSesRunCatAv); b = strcat('sesAvCat_',b);
        fSesAvCatRunCatAv = cellstr(fullfile(a,strcat(b,c)));
        mriSesAvCatRunCatAv = mriSesRunCatAv;
        for i = 1:length(fSesRunCatAv)
            mriSesCatRunCatAv(i).vol = cat(4,mriSesRunCatAv.vol);
            if force || ~exist(fSesCatRunCatAv{i},'file')
                MRIwrite(mriSesCatRunCatAv(i),fSesCatRunCatAv{i});
            end
            mriSesAvCatRunCatAv(i).vol = mean(mriSesCatRunCatAv(i).vol,4);
            if force || ~exist(fSesAvCatRunCatAv{i},'file')
                MRIwrite(mriSesAvCatRunCatAv(i),fSesAvCatRunCatAv{i});
            end
        end
        fSesCatRunCatAv   = fSesCatRunCatAv(bidsDerivDirInd)';
        fSesAvCatRunCatAv = fSesAvCatRunCatAv(bidsDerivDirInd)';
        for i = 1:length(tmp2)
            tmp2(i).fPreprocUnderSesCatRunCatAvList   = repmat(fSesCatRunCatAv(i),size(tmp2(i).fPreprocList));
            tmp2(i).fPreprocUnderSesAvCatRunCatAvList = repmat(fSesAvCatRunCatAv(i),size(tmp2(i).fPreprocList));
        end
        
        % average runs then catenate sessions
        [a,b,c] = fileparts(fSesRunAvCatAv); b = strcat('sesCat_',b);
        fSesCatRunAvCatAv = cellstr(fullfile(a,strcat(b,c)));
        mriSesCatRunAvCatAv   = mriSesRunAvCatAv;
        [a,b,c] = fileparts(fSesRunAvCatAv); b = strcat('sesAvCat_',b);
        fSesAvCatRunAvCatAv = cellstr(fullfile(a,strcat(b,c)));
        mriSesAvCatRunAvCatAv = mriSesRunAvCatAv;
        for i = 1:length(fSesRunAvCatAv)
            mriSesCatRunAvCatAv(i).vol = cat(4,mriSesRunAvCatAv.vol);
            if force || ~exist(fSesCatRunAvCatAv{i},'file')
                MRIwrite(mriSesCatRunAvCatAv(i),fSesCatRunAvCatAv{i});
            end
            mriSesAvCatRunAvCatAv(i).vol = mean(mriSesCatRunAvCatAv(i).vol,4);
            if force || ~exist(fSesAvCatRunAvCatAv{i},'file')
                MRIwrite(mriSesAvCatRunAvCatAv(i),fSesAvCatRunAvCatAv{i});
            end
        end
        fSesCatRunAvCatAv   = fSesCatRunAvCatAv(bidsDerivDirInd)';
        fSesAvCatRunAvCatAv = fSesAvCatRunAvCatAv(bidsDerivDirInd)';
        for i = 1:length(tmp2)
            tmp2(i).fPreprocUnderSesCatRunAvCatAvList   = repmat(fSesCatRunAvCatAv(i),size(tmp2(i).fPreprocList));
            tmp2(i).fPreprocUnderSesAvCatRunAvCatAvList = repmat(fSesAvCatRunAvCatAv(i),size(tmp2(i).fPreprocList));
        end


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

            %% Copy all relevant fields
            tmp4 = tmp3(1);
            tmp4.ses          = [];
            tmp4.wd           = {};
            tmp4.bidsDir      = {};
            tmp4.bidsDerivDir = {};
            for i = 1:length(tmp3)
                tmp4.ses          = cat(1,tmp4.ses         ,repmat(        tmp3(i).ses          ,size(tmp3(i).fPreprocList)));
                tmp4.wd           = cat(1,tmp4.wd          ,repmat(cellstr(tmp3(i).wd          ),size(tmp3(i).fPreprocList)));
                tmp4.bidsDir      = cat(1,tmp4.bidsDir     ,repmat(cellstr(tmp3(i).bidsDir)     ,size(tmp3(i).fPreprocList)));
                tmp4.bidsDerivDir = cat(1,tmp4.bidsDerivDir,repmat(cellstr(tmp3(i).bidsDerivDir),size(tmp3(i).fPreprocList)));
            end
            fieldList = {'fList' 'fOrigList' 'fPreprocList' 'fTransList' 'fTransCatList' 'bidsList' 'nFrame' 'vSize' 'acqTime' 'bhvr' 'nDummy' 'fPreprocUnderSesCatRunCatAvList' 'fPreprocUnderSesAvCatRunCatAvList' 'fPreprocUnderSesCatRunAvCatAvList' 'fPreprocUnderSesAvCatRunAvCatAvList'};
            for i = 1:length(fieldList)
                try
                    tmp4(1).(fieldList{i}) = cat(1,tmp3(:).(fieldList{i}));
                catch
                    sz = [0 max(cellfun('size',{tmp3(:).(fieldList{i})},2))];
                    tmpX = cell(sz);
                    for ii = 1:numel(tmp3)
                        tmpXX = tmp3(ii).(fieldList{i});
                        tmpXX = cat(2,tmpXX,repmat({''},[size(tmpXX,1) sz(2) - size(tmpXX,2)]));
                        tmpX  = cat(1,tmpX,tmpXX);
                    end
                    tmp4.(fieldList{i}) = tmpX;
                end
            end

            %% Define underlays (catenate and average across sessions)
            disp(['writing underlays for subj' num2str(S) '/' num2str(length(subList)) ', acqCond' num2str(ac) '/' num2str(length(runCondAcqList)) ', stimCond' num2str(rsc) '/' num2str(length(runCondStimList))])
            
            %%% across runs, within stimulus conditions, within sessions
            [bidsDerivDir,~,bidsDerivDirInd] = unique(tmp4.bidsDerivDir);
            fSesRunCatAv   = cell(size(bidsDerivDir));
            fSesRunAvCatAv = cell(size(bidsDerivDir));
            for s = 1:length(bidsDerivDir)
                fList = dir(char(fullfile(bidsDerivDir{s},['*task-' runCondStimList{rsc} '*'],'av_preproc_volTs.nii.gz')));
                fList = fullfile({fList.folder},{fList.name})';
                mri   = cell(size(fList));
                for i = 1:length(fList)
                    mri{i} = MRIread(fList{i});
                end
                mri = cat(1,mri{:});
                mriSesRunCatAv(s) = mri(1);
                mriSesRunCatAv(s).vol = cat(4,mri.vol); clear mri
                fSesRunCatAv{s} = fullfile(bidsDerivDir{s},['task-' runCondStimList{rsc} '_cat_av_preproc_volTs.nii.gz']);
                if force || ~exist(fSesRunCatAv{s},'file')
                    MRIwrite(mriSesRunCatAv(s),fSesRunCatAv{s});
                end
                mriSesRunAvCatAv(s) = mriSesRunCatAv(s);
                mriSesRunAvCatAv(s).vol = mean(mriSesRunAvCatAv(s).vol,4);
                fSesRunAvCatAv{s} = fullfile(bidsDerivDir{s},['task-' runCondStimList{rsc} '_av_cat_av_preproc_volTs.nii.gz']);
                if force || ~exist(fSesRunAvCatAv{s},'file')
                    MRIwrite(mriSesRunAvCatAv(s),fSesRunAvCatAv{s});
                end
            end

            %%% across sessions
            % catenate sessions then average runs
            [a,b,c] = fileparts(fSesRunCatAv); b = strcat('sesCat_',b);
            fSesCatRunCatAv = cellstr(fullfile(a,strcat(b,c)));
            mriSesCatRunCatAv   = mriSesRunCatAv;
            [a,b,c] = fileparts(fSesRunCatAv); b = strcat('sesAvCat_',b);
            fSesAvCatRunCatAv = cellstr(fullfile(a,strcat(b,c)));
            mriSesAvCatRunCatAv = mriSesRunCatAv;
            for i = 1:length(fSesRunCatAv)
                mriSesCatRunCatAv(i).vol = cat(4,mriSesRunCatAv.vol);
                if force || ~exist(fSesCatRunCatAv{i},'file')
                    MRIwrite(mriSesCatRunCatAv(i),fSesCatRunCatAv{i});
                end
                mriSesAvCatRunCatAv(i).vol = mean(mriSesCatRunCatAv(i).vol,4);
                if force || ~exist(fSesAvCatRunCatAv{i},'file')
                    MRIwrite(mriSesAvCatRunCatAv(i),fSesAvCatRunCatAv{i});
                end
            end
            tmp4.fPreprocUnderCondSpecSesCatRunCatAvList   = fSesCatRunCatAv(bidsDerivDirInd);
            tmp4.fPreprocUnderCondSpecSesAvCatRunCatAvList = fSesAvCatRunCatAv(bidsDerivDirInd);

            % average runs then catenate sessions
            [a,b,c] = fileparts(fSesRunAvCatAv); b = strcat('sesCat_',b);
            fSesCatRunAvCatAv = cellstr(fullfile(a,strcat(b,c)));
            mriSesCatRunAvCatAv   = mriSesRunAvCatAv;
            [a,b,c] = fileparts(fSesRunAvCatAv); b = strcat('sesAvCat_',b);
            fSesAvCatRunAvCatAv = cellstr(fullfile(a,strcat(b,c)));
            mriSesAvCatRunAvCatAv = mriSesRunAvCatAv;
            for i = 1:length(fSesRunAvCatAv)
                mriSesCatRunAvCatAv(i).vol = cat(4,mriSesRunAvCatAv.vol);
                if force || ~exist(fSesCatRunAvCatAv{i},'file')
                    MRIwrite(mriSesCatRunAvCatAv(i),fSesCatRunAvCatAv{i});
                end
                mriSesAvCatRunAvCatAv(i).vol = mean(mriSesCatRunAvCatAv(i).vol,4);
                if force || ~exist(fSesAvCatRunAvCatAv{i},'file')
                    MRIwrite(mriSesAvCatRunAvCatAv(i),fSesAvCatRunAvCatAv{i});
                end
            end
            tmp4.fPreprocUnderCondSpecSesCatRunAvCatAvList   = fSesCatRunAvCatAv(bidsDerivDirInd);
            tmp4.fPreprocUnderCondSpecSesAvCatRunAvCatAvList = fSesAvCatRunAvCatAv(bidsDerivDirInd);

            %% Compile
            runCond{S}.(runCondAcqList{ac}).(['task_' runCondStimList{rsc}]) = tmp4;
        end
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
            if ~isfield(runCond{S},'fs')
                runCond{S}.fs = tmp{i};
            else
                runCond{S}.fs(end+1,1) = tmp{i};
            end
        case 'avMap'
            if isempty(tmp{i}.fList)
                continue
            end
            if ~isfield(runCond{S},'avMap')
                runCond{S}.avMap = tmp{i};
            else
                runCond{S}.avMap(end+1,1) = tmp{i};
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
        if isfield(runCondTmp,'labelAcq') && isfield(runCondTmp,'label')
            {sesPhys{s}.mriRuns.label}';
            {sesPhys{s}.mriRuns.labelAcq}';
        else
            {sesPhys{s}.mriRuns.stim}';
            {sesPhys{s}.mriRuns.acq}';
        end
        if length(sesPhysSub{s})>1; dbstack; keyboard; error('physio sessions are confused'); end
        if length(sesPhysSes{s})>1; dbstack; keyboard; error('physio sessions are confused'); end
        sesPhysSub(s) = sesPhysSub{s};
        sesPhysSes(s) = sesPhysSes{s};
    end
end

for S = 1:length(runCond)
    Sind = ismember(sesPhysSub,subList{S});
    if ~any(Sind)
        runCond{S}.phys = [];
        continue
    end
    runCond{S}.phys = sesPhys{Sind};
    for s = 1:length(runCond{S}.phys)
        runCond{S}.phys(s).sub = sesPhysSub{Sind};
        runCond{S}.phys(s).info.sub = sesPhysSub{Sind};
        runCond{S}.phys(s).ses = sesPhysSes{Sind};
        runCond{S}.phys(s).info.ses = sesPhysSes{Sind};
    end
end


%% Extract physio runs and put it in runCond
for S = 1:length(runCond)
    if isempty(runCond{S}.phys); continue; end
    if length(runCond{S}.phys)>1; dbstack; error('more than one physio session, code that'); end
    
    physRun = physSes2run(runCond{S}.phys);
    physRunMriFile = [physRun.mri]; physRunMriFile = {physRunMriFile.fspec}';
    
    runCondAcqList = fields(runCond{S});
    runCondAcqList(ismember(runCondAcqList,{'phys' 'fs' 'avMap'})) = [];
    for rc = 1:length(runCondAcqList)
        runCondStimList = fields(runCond{S}.(runCondAcqList{rc}));
        for sc = 1:length(runCondStimList)
            mriFile = runCond{S}.(runCondAcqList{rc}).(runCondStimList{sc}).fList;
            physInd = ismember(physRunMriFile,mriFile);
            mriInd = ismember(mriFile,physRunMriFile);
            runCond{S}.(runCondAcqList{rc}).(runCondStimList{sc}).phys         = repmat(physSes2run,size(mriInd));
            runCond{S}.(runCondAcqList{rc}).(runCondStimList{sc}).phys(mriInd) = physRun(physInd);
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
nPhysSes      = {};
for S = 1:length(runCond)
    acqTime       = {};
    runCondAcqList = fields(runCond{S});
    runCondAcqList(ismember(runCondAcqList,{'phys' 'fs' 'avMap'})) = [];
    for ac = 1:length(runCondAcqList)
        runCondStimList = fields(runCond{S}.(runCondAcqList{ac}));
        for rsc = 1:length(runCondStimList)
            subList2{end+1}      = runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).sub;
            sesList2{end+1}      = runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).ses';
            if isfield(runCondTmp,'labelAcq') && isfield(runCondTmp,'label')
                acqCondList2{end+1}  = runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).labelAcq;
                stimCondList2{end+1} = runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).label;
            else
                acqCondList2{end+1}  = runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).acq;
                stimCondList2{end+1} = runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).stim;
            end
            n{end+1}             = length(runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).fPreprocList);
            acqTime{end+1}       = runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).acqTime;
            nPhys{end+1}         = 0;
            if isfield(runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}),'phys')
                nPhys{end}         = nnz(~[runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).phys.isempty]);
            end
            if nPhys{end}~=0
                ses = [runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).phys.info];
                nPhysSes{end+1} = [ses.ses];
            else
                nPhysSes{end+1} = '';
            end
        end
    end
    %anonymize
    firstAcqTime = min(cat(1,acqTime{:}));
    runCondAcqList = fields(runCond{S});
    runCondAcqList(ismember(runCondAcqList,{'phys' 'fs' 'avMap'})) = [];
    for ac = 1:length(runCondAcqList)
        runCondStimList = fields(runCond{S}.(runCondAcqList{ac}));
        for rsc = 1:length(runCondStimList)
            runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).acqTime = runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).acqTime - firstAcqTime;
            if isfield(runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}),'date')
                runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}) = rmfield(runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}),'date');
            end
        end
    end
end

runCond{S}.phys

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
nPhysSes = nPhysSes(b);

[~,b] = sort(subList2);
acqCondList2 = acqCondList2(b);
subList2 = subList2(b);
sesList2 = sesList2(b);
stimCondList2 = stimCondList2(b);
n = n(b);
nPhysSes = nPhysSes(b);

disp([{'mriCond' 'sub' 'mriSes' 'stimCond' 'nRun' 'physSes'}
      {'-------' '---' '------' '--------' '----' '-------'}
[acqCondList2
subList2
sesList2
stimCondList2
n
nPhysSes]'])

runCondStimList = strcat('task_',unique(stimCondList2)');
runCondAcqList  = unique(acqCondList2)';
