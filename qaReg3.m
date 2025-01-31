function cmd = qaReg3(xSet,ppLabel,acqLabel,xFlag)
global srcFs srcAfni
if ~exist('ppLabel','var');  ppLabel  = ''    ; end
if isempty(ppLabel);         ppLabel  = 'full'; end % 'full' 'init' 'final' 'bRun'
if ~exist('acqLabel','var'); acqLabel = ''    ; end
if ~exist('xFlag','var');       xFlag = []    ; end
if isempty(xFlag);              xFlag = 0     ; end

%% Refactor xSet
% should be a vector of subjects where each cell containing a vector of
% acquisition set
if isstruct(xSet); xSet = {xSet}; end
for s = 1:length(xSet); if isstruct(xSet{s}); xSet{s} = {xSet{s}}; end; end

subList   = {};
labelList = {};
sesList   = {};
dummy     = {};
nFrame    = {};
outDir    = {};
for s = 1:length(xSet)
    subList{end+1} = {};
    labelList{end+1} = {};
    for ss = 1:length(xSet{s})
        subList{end}{end+1} = xSet{s}{ss}.sub;
        labelList{end}{end+1} = xSet{s}{ss}.label;
    end
    if ~isempty(acqLabel)
        xSet{s} = xSet{s}{ismember(labelList{end},acqLabel)};
        labelList{end} = labelList{end}(ismember(labelList{end},acqLabel));
    end
    sesList{end+1,1} = xSet{s}.ses;
    dummy{end+1,1}   = xSet{s}.nDummy;
    nFrame{end+1,1}  = xSet{s}.initFiles.nFrame;
    outDir{end+1,1}    = xSet{s}.bidsDir;
end
subList   = unique([subList{:}]);   if length(subList)>1; dbstack; error('qa only one subject at a time'); end
labelList = unique([labelList{:}]); if length(labelList)>1 && isempty(acqLabel); dbstack; error('more than one acquisition types, please provide acqLabel'); end
dummy     = cat(1,dummy {:});
nFrame    = cat(1,nFrame{:});

outDir = strsplit(outDir{1},filesep); outDir = strjoin(outDir(1:find(ismember(outDir,['sub-' char(subList)]))),filesep);
outDir = fullfile(outDir,'derivatives'); if ~exist(outDir,'dir'); mkdir(outDir); end
outDir = fullfile(outDir,'qaReg');       if ~exist(outDir,'dir'); mkdir(outDir); end
outDir = fullfile(outDir,acqLabel);      if ~exist(outDir,'dir'); mkdir(outDir); end
outDir = fullfile(outDir,ppLabel);       if ~exist(outDir,'dir'); mkdir(outDir); end

% get all files
fBefore     = {};
fBidsBefore = {};
fAfter      = {};
fBidsAfter  = {};
sList   = {};
mList   = {};
tList   = {};
for s = 1:length(xSet)
    sList = cat(1,sList,repmat({xSet{s}.ses},size(xSet{s}.fList)));
    tList = cat(1,tList,xSet{s}.initFiles.acqTime);
    switch ppLabel
        case 'full'
            fBefore = cat(1,fBefore,xSet{s}.initFiles.fPlumbList        );
            bids = cell(size(xSet{s}.initFiles.bidsList,1),1);
            for f = 1:size(xSet{s}.initFiles.bidsList,1)
                bids{f,1} = strjoin(xSet{s}.initFiles.bidsList(f,:),'_');
            end
            fBidsBefore = cat(1,fBidsBefore,bids);
            
            fAfter  = cat(1,fAfter ,xSet{s}.finalFiles.fPreprocList);
            bids = cell(size(xSet{s}.finalFiles.bidsList,1),1);
            for f = 1:size(xSet{s}.finalFiles.bidsList,1)
                bids{f,1} = strjoin(xSet{s}.finalFiles.bidsList(f,:),'_');
            end
            fBidsAfter = cat(1,fBidsAfter,bids);
        case 'final'
            dbstack; error('code that');
            fAfter = cat(1,fAfter,xSet{s}.finalFiles.fPreprocList);
            bids = cell(size(xSet{s}.finalFiles.bidsList,1),1);
            for f = 1:size(xSet{s}.finalFiles.bidsList,1)
                bids{f,1} = strjoin(xSet{s}.finalFiles.bidsList(f,:),'_');
            end
            fBids = cat(1,fBids,bids);
        case 'bRun'
            dbstack; error('code that');
            fBase = xSet{s}.brMocoFiles.fBase;
            fBefore = cat(1,fBefore,xSet{s}.brMocoFiles.fEstimList);
            fAfter = cat(1,fAfter,xSet{s}.brMocoFiles.fMocoList);
            bids = cell(size(xSet{s}.brMocoFiles.bidsList,1),1);
            for f = 1:size(xSet{s}.brMocoFiles.bidsList,1)
                bids{f,1} = strjoin(xSet{s}.brMocoFiles.bidsList(f,:),'_');
            end
            fBids = cat(1,fBids,bids);
        case 'init'
            dbstack; error('code that');
            fAfter = cat(1,fAfter,xSet{s}.initFiles.fPlumbList);
            bids = cell(size(xSet{s}.finalFiles.bidsList,1),1);
            for f = 1:size(xSet{s}.finalFiles.bidsList,1)
                bids{f,1} = strjoin(xSet{s}.initFiles.bidsList(f,:),'_');
            end
            fBids = cat(1,fBids,bids);
        otherwise
            dbstack; error('x');
    end
    % get mask
    if isfield(xSet{s},'wrMocoFiles')
        mList = cat(1,mList,xSet{s}.wrMocoFiles.fMaskList);
    elseif isfield(xSet{s},'brMocoFiles')
        mList = cat(1,mList,xSet{s}.brMocoFiles.fMaskList);
    else
        dbstack; error('cannot find mask files')
    end
end
% [~,b] = sort(tList);
% fBids = fBids(b);
% fList = fList(b);
% sList = sList(b);
% mList = mList(b);
% tList = tList(b);

% %% Compute costs
% [costBefore,costLabelBefore] = getAllCost(fBefore);
% [costAfter ,costLabelAfter ] = getAllCost(fAfter );
% save(fullfile(outDir,'allcost'),'costBefore','costLabelBefore','costAfter','costLabelAfter')



%% Summarize data before
[im,imMask,cmd.allBefore,cmd.allFMLbefore,cmd.allMeansBefore] = summarizeData(fBefore,mList,dummy,outDir,fBidsBefore,'before');
correlateFrames(im,logical(1-imMask),subList,labelList,ppLabel,nFrame,fBidsBefore,fBefore,outDir)

%% Summarize after after
[im,imMask,cmd.allAfter,cmd.allFMLAfter,cmd.allMeansAfter] = summarizeData(fAfter,mList,dummy,outDir,fBidsAfter,'after');
correlateFrames(im,logical(1-imMask),subList,labelList,ppLabel,nFrame,fBidsAfter,fAfter,outDir)


function correlateFrames(im,imMask,subList,labelList,ppLabel,nFrame,fBids,f,outDir)
% correlate frames
for r = 1:length(im)
    im{r} = permute(im{r},[4 1 2 3]);
    im{r} = permute(im{r}(:,imMask(:,:,:,:,r)),[2 1]);
end
im = cat(2,im{:});
% im = permute(cat(4,im{:}),[4 1 2 3]);
% im = permute(im(:,:),[2 1]);
% im(any(isnan(im),2),:) = [];
rho = corr(im);

% sz = size(im);
% im = im(:,:,:,:);
% im = permute(im,[4 1 2 3]);
% im = permute(im(:,:),[2 1]); % vox x time
% im(any(isnan(im),2),:) = [];
% rho = corr(im);

figure('WindowStyle','docked');
imagesc(rho,[0 1]);
ax = gca;
ax.DataAspectRatio = [1 1 1];
ylabel(colorbar,'cross-frame correlation')
title(strjoin([subList labelList {ppLabel}],'; '))

dTick = round(cumsum(nFrame) - nFrame/2);
ax.XTick = dTick;
ax.YTick = dTick;
ax.TickLabelInterpreter = 'none';
ax.XTickLabels = fBids;
ax.YTickLabels = fBids;

% if length(rho)>length(f)
%     xline(dTick+0.5+nFrame/2)
%     yline(dTick+0.5+nFrame/2)
% end
drawnow

out = fullfile(outDir,[ppLabel '_xFrameCorrelation.fig']);
disp(out)
saveas(gcf,out);




return



cmd = {srcFs};
cmd{end+1} = 'freeview \';

for S = 1:length(xSet)

    switch ppLabel

        case 'final'
            xSet{S}.finalFiles

            fieldList = fields(xSet{S}.finalFiles);

            wd = cell(size(xSet{S}.finalFiles.fList));
            for r = 1:size(xSet{S}.finalFiles.fList,1)
                wd{r,1} = fullfile(xSet{S}.finalFiles.wd,strjoin(xSet{S}.finalFiles.bidsList(r,:),'_'));
            end
            %% All within-run averages, catenated across runs
            dir(xSet{S}.finalFiles.wd)




            fieldList(contains(fieldList,'rainMask')) = [];
            fieldList(contains(fieldList,'EchoCat')) = [];

            if any(contains(fieldList,'EchoRms'))
                fieldList(~contains(fieldList,'EchoRms')) = [];
            end
            if any(contains(fieldList,'fCatAv'))
                fieldList = fieldList(contains(fieldList,'fCatAv'));
            elseif any(contains(fieldList,'fAv'))
                fieldList = fieldList(contains(fieldList,'fAv'));
            elseif any(ismember(fieldList,'f')) || any(ismember(fieldList,'fEchoRms'))
                ismember(fieldList,'f') | ismember(fieldList,'fEchoRms')
                fieldList = fieldList(ismember(fieldList,'f') | ismember(fieldList,'fEchoRms'));
                % fieldList = fieldList(contains(fieldList,'fAv'));
            end

            if length(fieldList)~=1; dbstack; error('can''t find what to show'); end

            if xSet{S}.nVenc
                cmd{end+1} = [char(xSet{S}.finalFiles.(char(fieldList)){1}) ':resample=cubic:name=' replace(xSet{S}.label,'set-','') ' \'];
            else
                cmd{end+1} = [char(xSet{S}.finalFiles.(char(fieldList))) ':resample=cubic:name=' replace(xSet{S}.label,'set-','') ' \'];
            end



        case 'preproc'

            fieldList = fields(xSet{S}.preprocFiles);
            if any(contains(fieldList,'EchoRms'))
                fieldList(~contains(fieldList,'EchoRms')) = [];
            end
            if any(contains(fieldList,'fCorrectedCatAv'))
                fieldList = fieldList(contains(fieldList,'fCorrectedCatAv'));
            elseif any(contains(fieldList,'fCorrectedAv'))
                fieldList = fieldList(contains(fieldList,'fCorrectedAv'));
            elseif any(contains(fieldList,'fCorrected'))
                fieldList = fieldList(contains(fieldList,'fCorrected'));
            end

            if length(fieldList)>1; dbstack; error('can''t find what to show'); end

            if xSet{S}.nVenc
                cmd{end+1} = [xSet{S}.preprocFiles.(char(fieldList)){1} ':resample=cubic:name=' replace(xSet{S}.label,'set-','') ' \'];
            else
                cmd{end+1} = [char(xSet{S}.preprocFiles.(char(fieldList))) ':resample=cubic:name=' replace(xSet{S}.label,'set-','') ' \'];
            end


        otherwise
            dbstack; error('please specify valid fileLabel')
    end
end



switch ppLabel

    case 'final'
        cmd{end+1} = [xSet{1}.finalFiles.manBrainMaskInv ':name=mask &'];

    case 'preproc'
        cmd{end+1} = [xSet{1}.preprocFiles.manBrainMaskInv ':name=mask &'];

end

cmd = strjoin(cmd,newline);
disp('QA registration')
disp(cmd);
clipboard('copy',cmd);
disp('command copied to clipboard');


function [cost,costLabel] = getCost(fList,nFrame)
global srcAfni
% costList = cell(length(fList));
cost = cell(length(fList),length(fList));
for r1 = 1:length(fList)
    for r2 = 1:length(fList)
        tic
        if r2>r1; break; end
        disp(['cost run' num2str(r2) '->run' num2str(r1)])
        % costList{r1,r2} = cell(nFrame(r1),nFrame(r2));
        cost{r1,r2} = nan(nFrame(r1),nFrame(r2),14);
        for f1 = 1:nFrame(r1)
            disp(['cost run' num2str(r2) '->run' num2str(r1) 'frame' num2str(f1)])
            % for f2 = 1:nFrame(r2)
            cmd = {srcAfni};
            % costList{r1,r2}{f1,f2} = tempname;
            cmd{end+1} = '3dAllineate -overwrite \';
            cmd{end+1} = ['-base   ' fList{r1} '[' num2str(f1-1) '] \'];
            % cmd{end+1} = ['-source ' fList{r2} '[' num2str(f2-1) '] \'];
            % cmd{end+1} = ['-allcostX1D "IDENTITY" ' costList{r1,r2}{f1,f2} ' > /dev/null 2>&1'];
            cmd{end+1} = ['-source ' fList{r2} ' \'];
            % cmd{end+1} = ['-allcostX "IDENTITY" ' costList{r1,r2}{f1,f2}];
            cmd{end+1} = '-allcostX';
            [~,cmdout] = system(strjoin(cmd,newline));
            costTmp = cmdout;
            costTmp = strsplit(costTmp,'++ allcost output:'); costTmp(1) = [];
            for i = 1:length(costTmp)
                costTmp{i} = strsplit(costTmp{i},newline)';
                costTmp{i} = costTmp{i}(2:15);
                if all([r1 r2 f1]==1)
                    costLabel = cell(length(costTmp{i}),1);
                end
                for ii = 1:length(costTmp{i})
                    costTmp{i}{ii} = strsplit(costTmp{i}{ii},'=');
                    if all([r1 r2 f1]==1)
                        costLabel{ii} = char(replace(costTmp{i}{ii}(1),' ',''));
                    end
                    cost{r1,r2}(f1,i,ii) = str2double(costTmp{i}{ii}{2});
                end
            end
        end
        toc
    end
end

% cost = cell2mat(cost);
% for c = 1:size(cost,3)
%     cost(:,:,c) = triu(cost(:,:,c)) + triu(cost(:,:,c),1)';
% end







% function [cost,costLabel] = getCost2(fList,nFrame)
% global srcAfni
% % costList = cell(length(fList));
% cost = cell(length(fList),length(fList));
% for r1 = 1:length(fList)
%     for r2 = 1:length(fList)
%         tic
%         if r2>r1; break; end
%         disp(['cost run' num2str(r2) '->run' num2str(r1)])
%         % costList{r1,r2} = cell(nFrame(r1),nFrame(r2));
%         cost{r1,r2} = nan(nFrame(r1),nFrame(r2),14);
%         for f1 = 1:nFrame(r1)
%             disp(['cost run' num2str(r2) '->run' num2str(r1) 'frame' num2str(f1)])
%             % for f2 = 1:nFrame(r2)
%             cmd = {srcAfni};
%             % costList{r1,r2}{f1,f2} = tempname;
%             cmd{end+1} = '3dAllineate -overwrite \';
%             cmd{end+1} = ['-base   ' fList{r1} '[' num2str(f1-1) '] \'];
%             % cmd{end+1} = ['-source ' fList{r2} '[' num2str(f2-1) '] \'];
%             % cmd{end+1} = ['-allcostX1D "IDENTITY" ' costList{r1,r2}{f1,f2} ' > /dev/null 2>&1'];
%             cmd{end+1} = ['-source ' fList{r2} ' \'];
%             param1D = [tempname '.1D']; writematrix(repmat([0 0 0  0 0 0  1 1 1  0 0 0],[nFrame(r2) 1]),param1D,'FileType','text');
%             cost1D = [tempname '.1D'];
%             cmd{end+1} = ['-allcostX1D ' param1D ' ' cost1D ' > /dev/null 2>&1'];
%             % cmd{end+1} = ['-allcostX "IDENTITY" ' costList{r1,r2}{f1,f2}];
%             % cmd{end+1} = '-allcostX';
%             [~,cmdout] = system(strjoin(cmd,newline));
%             costTmp = cmdout;
%             costTmp = strsplit(costTmp,'++ allcost output:'); costTmp(1) = [];
%             for i = 1:length(costTmp)
%                 costTmp{i} = strsplit(costTmp{i},newline)';
%                 costTmp{i} = costTmp{i}(2:15);
%                 if all([r1 r2 f1]==1)
%                     costLabel = cell(length(costTmp{i}),1);
%                 end
%                 for ii = 1:length(costTmp{i})
%                     costTmp{i}{ii} = strsplit(costTmp{i}{ii},'=');
%                     if all([r1 r2 f1]==1)
%                         costLabel{ii} = char(replace(costTmp{i}{ii}(1),' ',''));
%                     end
%                     cost{r1,r2}(f1,i,ii) = str2double(costTmp{i}{ii}{2});
%                 end
%             end
%         end
%         toc
%     end
% end
% 
% % cost = cell2mat(cost);
% % for c = 1:size(cost,3)
% %     cost(:,:,c) = triu(cost(:,:,c)) + triu(cost(:,:,c),1)';
% % end







function [im,imMask,cmdAll,cmdFML,cmdMeans] = summarizeData(fList,mList,dummy,outDir,fBids,ppLabel)
global srcFs
xFlag = 0;

mriHd   = MRIread(fList{1},1);
mriMean = repmat(mriHd,size(fList));
fFML    = cell(size(fList));
nFrames = zeros(size(fList));
% im      = zeros([mriHd.volsize mriHd.nframes length(fList)],'single');
im      = cell(size(fList));
imMask  = zeros([mriHd.volsize 1             length(fList)],'single');
for f = 1:length(fList)
    disp(['reading and summarizing (' num2str(f) '/' num2str(length(fList)) ') ' fList{f} ])
    %read
    mri     = MRIread(fList{f});
    mriMask = MRIread(mList{f});


    % get all frames
    nFrames(f) = mri.nframes;
    tmp = permute(single(mri.vol),[4 1 2 3]);
    tmp(:,logical(mriMask.vol)) = nan;
    % im(:,:,:,:,f) = permute(tmp,[2 3 4 1]);
    im{f} = permute(tmp,[2 3 4 1]);

    imMask(:,:,:,:,f) = single(mriMask.vol);



    %get mean
    mriMean(f).vol = mean(mri.vol,4);


    %get first,middle and last frames
    if mri.nframes>1
        ind = [dummy(f) floor(mean([dummy(f) mri.nframes])) mri.nframes];
        mri.vol = mri.vol(:,:,:,ind);
        out = fullfile(outDir,'fml'); if ~exist(out,'dir'); mkdir(out); end
        fFML{f} = fullfile(out,[ppLabel '_' fBids{f} '.nii.gz']);

        tmp = mri; tmp.vol(isnan(tmp.vol)) = 0;
        disp(fFML{f})
        MRIwrite(tmp,fFML{f});
    else
        fFML{f} = '';
    end
end

% write catenated means
mriMean(1).vol = cat(4,mriMean.vol); mriMean(2:end) = [];
fMeans = fullfile(outDir,[ppLabel '_means.nii.gz']);
tmp = mriMean; tmp.vol(isnan(tmp.vol)) = 0;
disp(fMeans)
MRIwrite(tmp,fMeans); clear mriMean mri



% Produce visualization commands
% construct command

%means
cmd = {srcFs};
% switch acqLabel
%     case 'vfMRI'
        cmd{end+1} = 'fslview -m single \';
        cmd{end+1} = [fMeans ' -b 0,800 &'];
%     otherwise
%         cmd{end+1} = 'fslview \';
%         cmd{end+1} = [fMeans ' &'];
% end
cmdMeans = strjoin(cmd,newline);

%first, middle and last frames
cmdFML = cell(size(fFML));
for r = 1:length(fFML)
    cmd = {srcFs};
    % switch acqLabel
    %     case 'vfMRI'
            cmd{end+1} = 'fslview -m single \';
            cmd{end+1} = [strjoin(fFML(r),' -b 0,800 ') ' -b 0,800 &'];
    %     otherwise
    %         cmd{end+1} = 'fslview \';
    %         cmd{end+1} = [strjoin(fFML(r),' ') ' &'];
    % end
    cmdFML{r} = strjoin(cmd,newline);
end
cmd = {srcFs};
% switch acqLabel
%     case 'vfMRI'
        cmd{end+1} = 'fslview -m single \';
        cmd{end+1} = [strjoin(flip(fFML),' -b 0,800 ') ' -b 0,800 &'];
%     otherwise
%         cmd{end+1} = 'fslview \';
%         cmd{end+1} = [strjoin(flip(fFML),' ') ' &'];
% end
cmdFML{end+1} = strjoin(cmd,newline);

%all frames
cmdAll = cell(size(fList));
for r = 1:length(fList)
    cmd = {srcFs};
    % switch acqLabel
    %     case 'vfMRI'
            cmd{end+1} = 'fslview -m single \';
            cmd{end+1} = [strjoin(fList(r),' -b 0,800 ') ' -b 0,800 &'];
    %     otherwise
    %         cmd{end+1} = 'fslview \';
    %         cmd{end+1} = [strjoin(f(r),' ') ' &'];
    % end
    cmdAll{r} = strjoin(cmd,newline);
end

cmd = {srcFs};
% switch acqLabel
%     case 'vfMRI'
        cmd{end+1} = 'fslview -m single \';
        cmd{end+1} = [strjoin(flip(fList),' -b 0,800 ') ' -b 0,800 &'];
%     otherwise
%         cmd{end+1} = 'fslview \';
%         cmd{end+1} = [strjoin(flip(f),' ') ' &'];
% end
cmdAll{end+1} = strjoin(cmd,newline);

%execute commands
if xFlag
    system(cmdMeans);
    system(cmdFML{1});
    system(cmdAll{1});
end



