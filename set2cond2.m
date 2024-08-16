function [runCond,subList,runCondAcqList,runCondStimList] = set2cond2(runSet,runCond)


runSetTmp = cat(2,runSet{:})';
runCondTmp = cat(2,runCond{:})';
runCondTmp = [runCondTmp{:}]';

runCondTmpTmpTmp = [];

for rs = 1:length(runSetTmp)
    if isempty(runSetTmp{rs}.fList); continue; end

    ind = ismember({runCondTmp.sub}',runSetTmp{rs}.sub) & ...
        ismember({runCondTmp.ses}',runSetTmp{rs}.ses) & ...
        ismember({runCondTmp.labelAcq}',runSetTmp{rs}.label);

    runCondTmpTmp = runCondTmp(ind);
    for rc = 1:length(runCondTmpTmp)
        runCondTmpTmp(rc).wd           = runSetTmp{rs}.finalFiles.wd;
        runCondTmpTmp(rc).bidsDir      = runSetTmp{rs}.finalFiles.bidsDir;
        runCondTmpTmp(rc).bidsDerivDir = runSetTmp{rs}.finalFiles.bidsDerivDir;
        runCondTmpTmp(rc).ppLabelList  = runSetTmp{rs}.finalFiles.ppLabelList;
        runCondTmpTmp(rc).dataType     = runSetTmp{rs}.finalFiles.dataType;
        if isfield(runSetTmp{rs},'avMap')
            runCondTmpTmp(rc).avMap = runSetTmp{rs}.avMap;
        else
            runCondTmpTmp(rc).avMap = [];
        end

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

[{runCondTmpTmpTmp.sub}
{runCondTmpTmpTmp.ses}
{runCondTmpTmpTmp.labelAcq}
{runCondTmpTmpTmp.label}
cellstr(num2str(cellfun('size',{runCondTmpTmpTmp.fList},1)'))']'




subList = unique({runCondTmpTmpTmp.sub})';
runCond = cell(size(subList));
for S = 1:length(subList)
    % if S~=2; continue; end
    % sort subjects
    ind = ismember({runCondTmpTmpTmp.sub},subList{S});
    tmp = runCondTmpTmpTmp(ind);

    % [{tmp.sub}' {tmp.ses}' {tmp.labelAcq}' {tmp.label}' cellstr(num2str(cellfun('size',{tmp.fList},1)'))]

    runCondAcqList = unique({tmp.labelAcq});
    for ac = 1:length(runCondAcqList)
        % if ac~=2; continue; end
        % sort acquisition conditions
        ind = ismember({tmp.labelAcq},runCondAcqList(ac));
        tmp2 = tmp(ind);
        [runCondStimList,~,runCondStimInd] = unique({tmp2.label});

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
            MRIwrite(mriSesCatRunCatAv(i),fSesCatRunCatAv{i});
            mriSesAvCatRunCatAv(i).vol = mean(mriSesCatRunCatAv(i).vol,4);
            MRIwrite(mriSesAvCatRunCatAv(i),fSesAvCatRunCatAv{i});
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
            MRIwrite(mriSesCatRunAvCatAv(i),fSesCatRunAvCatAv{i});
            mriSesAvCatRunAvCatAv(i).vol = mean(mriSesCatRunAvCatAv(i).vol,4);
            MRIwrite(mriSesAvCatRunAvCatAv(i),fSesAvCatRunAvCatAv{i});
        end
        fSesCatRunAvCatAv   = fSesCatRunAvCatAv(bidsDerivDirInd)';
        fSesAvCatRunAvCatAv = fSesAvCatRunAvCatAv(bidsDerivDirInd)';
        for i = 1:length(tmp2)
            tmp2(i).fPreprocUnderSesCatRunAvCatAvList   = repmat(fSesCatRunAvCatAv(i),size(tmp2(i).fPreprocList));
            tmp2(i).fPreprocUnderSesAvCatRunAvCatAvList = repmat(fSesAvCatRunAvCatAv(i),size(tmp2(i).fPreprocList));
        end


        for rsc = 1:length(runCondStimList)
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
            fieldList = {'fPreprocList' 'fTransList' 'fTransCatList' 'bidsList' 'nFrame' 'vSize' 'acqTime' 'fOrigList' 'nDummy' 'fPreprocUnderSesCatRunCatAvList' 'fPreprocUnderSesAvCatRunCatAvList' 'fPreprocUnderSesCatRunAvCatAvList' 'fPreprocUnderSesAvCatRunAvCatAvList'};
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
                MRIwrite(mriSesRunCatAv(s),fSesRunCatAv{s});
                mriSesRunAvCatAv(s) = mriSesRunCatAv(s);
                mriSesRunAvCatAv(s).vol = mean(mriSesRunAvCatAv(s).vol,4);
                fSesRunAvCatAv{s} = fullfile(bidsDerivDir{s},['task-' runCondStimList{rsc} '_av_cat_av_preproc_volTs.nii.gz']);
                MRIwrite(mriSesRunAvCatAv(s),fSesRunAvCatAv{s});
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
                MRIwrite(mriSesCatRunCatAv(i),fSesCatRunCatAv{i});
                mriSesAvCatRunCatAv(i).vol = mean(mriSesCatRunCatAv(i).vol,4);
                MRIwrite(mriSesAvCatRunCatAv(i),fSesAvCatRunCatAv{i});
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
                MRIwrite(mriSesCatRunAvCatAv(i),fSesCatRunAvCatAv{i});
                mriSesAvCatRunAvCatAv(i).vol = mean(mriSesCatRunAvCatAv(i).vol,4);
                MRIwrite(mriSesAvCatRunAvCatAv(i),fSesAvCatRunAvCatAv{i});
            end
            tmp4.fPreprocUnderCondSpecSesCatRunAvCatAvList   = fSesCatRunAvCatAv(bidsDerivDirInd);
            tmp4.fPreprocUnderCondSpecSesAvCatRunAvCatAvList = fSesAvCatRunAvCatAv(bidsDerivDirInd);

            %% Compile
            runCond{S}.(runCondAcqList{ac}).(['task_' runCondStimList{rsc}]) = tmp4;
        end
    end
end


subList2      = {};
sesList2      = {};
acqCondList2  = {};
stimCondList2 = {};
n             = {};
for S = 1:length(runCond)
    runCondAcqList = fields(runCond{S});
    for ac = 1:length(runCondAcqList)
        runCondStimList = fields(runCond{S}.(runCondAcqList{ac}));
        for rsc = 1:length(runCondStimList)
            subList2{end+1}      = runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).sub;
            sesList2{end+1}      = runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).ses';
            acqCondList2{end+1}  = runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).labelAcq;
            stimCondList2{end+1} = runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).label;
            n{end+1}             = length(runCond{S}.(runCondAcqList{ac}).(runCondStimList{rsc}).fPreprocList);
        end
    end
end


[~,b] = sort(acqCondList2);
acqCondList2 = acqCondList2(b);
subList2 = subList2(b);
sesList2 = sesList2(b);
stimCondList2 = stimCondList2(b);
n = n(b);

[~,b] = sort(subList2);
acqCondList2 = acqCondList2(b);
subList2 = subList2(b);
sesList2 = sesList2(b);
stimCondList2 = stimCondList2(b);
n = n(b);

[acqCondList2
subList2
sesList2
stimCondList2
n]'

