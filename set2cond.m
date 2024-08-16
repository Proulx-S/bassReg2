function [runCond,subList,runCondAcqList,runCondStimList] = set2cond(runSet,runCond)


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
    % sort subjects
    ind = ismember({runCondTmpTmpTmp.sub},subList{S});
    tmp = runCondTmpTmpTmp(ind);

    runCondAcqList = unique({tmp.labelAcq});
    for ac = 1:length(runCondAcqList)
        % sort acquisition conditions
        ind = ismember({tmp.labelAcq},runCondAcqList(ac));
        tmp2 = tmp(ind);
        [runCondStimList,~,runCondStimInd] = unique({tmp2.label});
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
            % if length(tmp3)>1; keyboard; end
            % tmp3.bidsDerivDir
            % bidsDerivDir
            % ses
            fieldList = {'fPreprocList' 'fTransList' 'fTransCatList' 'bidsList' 'nFrame' 'vSize' 'acqTime' 'fOrigList' 'nDummy'};
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
            [bidsDerivDir,~,bidsDerivDirInd] = unique(tmp4.bidsDerivDir);
            disp(['writing underlays for subj' num2str(S) '/' num2str(length(subList)) ', cond' num2str(rsc) '/' num2str(length(runCondStimList))])
            % catenate sessions then average runs
            fSesRunCatAv = fullfile(bidsDerivDir,'cat_av_preproc_volTs.nii.gz');
            for i = 1:length(fSesRunCatAv)
                mriSesRunCatAv(i) = MRIread(fSesRunCatAv{i});
            end
            runCondStimList{rsc}
            fSesCatRunCatAv     = fullfile(bidsDerivDir,['task-' runCondStimList{rsc} '_sesCat_cat_av_preproc_volTs.nii.gz']);
            mriSesCatRunCatAv   = mriSesRunCatAv;
            fSesAvCatRunCatAv   = fullfile(bidsDerivDir,['task-' runCondStimList{rsc} '_sesAvCat_cat_av_preproc_volTs.nii.gz']);
            mriSesAvCatRunCatAv = mriSesRunCatAv;
            for i = 1:length(fSesRunCatAv)
                mriSesCatRunCatAv(i).vol = cat(4,mriSesRunCatAv.vol);
                MRIwrite(mriSesCatRunCatAv(i),fSesCatRunCatAv{i});
                mriSesAvCatRunCatAv(i).vol = mean(mriSesCatRunCatAv(i).vol,4);
                MRIwrite(mriSesAvCatRunCatAv(i),fSesAvCatRunCatAv{i});
            end
            tmp4.fPreprocUnderSesCatRunCatAvList   = fSesCatRunCatAv(bidsDerivDirInd);
            tmp4.fPreprocUnderSesAvCatRunCatAvList = fSesAvCatRunCatAv(bidsDerivDirInd);

            % average runs then catenate sessions
            fSesRunAvCatAv = fullfile(bidsDerivDir,'av_cat_av_preproc_volTs.nii.gz');
            for i = 1:length(fSesRunAvCatAv)
                mriSesRunAvCatAv(i) = MRIread(fSesRunAvCatAv{i});
            end
            fSesCatRunAvCatAv     = fullfile(bidsDerivDir,['task-' runCondStimList{rsc} '_sesCat_av_cat_av_preproc_volTs.nii.gz']);
            mriSesCatRunAvCatAv   = mriSesRunAvCatAv;
            fSesAvCatRunAvCatAv   = fullfile(bidsDerivDir,['task-' runCondStimList{rsc} '_sesAvCat_av_cat_av_preproc_volTs.nii.gz']);
            mriSesAvCatRunAvCatAv = mriSesRunAvCatAv;
            for i = 1:length(fSesRunAvCatAv)
                mriSesCatRunAvCatAv(i).vol = cat(4,mriSesRunAvCatAv.vol);
                MRIwrite(mriSesCatRunAvCatAv(i),fSesCatRunAvCatAv{i});
                mriSesAvCatRunAvCatAv(i).vol = mean(mriSesCatRunAvCatAv(i).vol,4);
                MRIwrite(mriSesAvCatRunAvCatAv(i),fSesAvCatRunAvCatAv{i});
            end
            tmp4.fPreprocUnderSesCatRunAvCatAvList   = fSesCatRunAvCatAv(bidsDerivDirInd);
            tmp4.fPreprocUnderSesAvCatRunAvCatAvList = fSesAvCatRunAvCatAv(bidsDerivDirInd);

            %% Compile
            runCond{S}.(runCondAcqList{ac}).(['task_' runCondStimList{rsc}]) = tmp4;
        end
    end

    runCond
end





