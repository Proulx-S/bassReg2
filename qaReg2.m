function [cmdMeans,cmdFML,cmdAll] = qaReg2(xSet,ppLabel,acqLabel,xFlag)
global srcFs srcAfni
if ~exist('ppLabel','var');  ppLabel  = ''     ; end
if isempty(ppLabel);         ppLabel  = 'final'; end % 'init' 'final'
if ~exist('acqLabel','var'); acqLabel = ''     ; end
if ~exist('xFlag','var');       xFlag = []     ; end
if isempty(xFlag);              xFlag = 0      ; end

%% Refactor xSet
% should be a vector of subjects where each cell containing a vector of
% acquisition set
if isstruct(xSet); xSet = {xSet}; end
for s = 1:length(xSet); if isstruct(xSet{s}); xSet{s} = {xSet{s}}; end; end

subList   = {};
labelList = {};
sesList   = {};
dummy     = {};
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
    outDir{end+1,1}    = xSet{s}.bidsDir;
end
subList   = unique([subList{:}]);   if length(subList)>1; dbstack; error('qa only one subject at a time'); end
labelList = unique([labelList{:}]); if length(labelList)>1 && isempty(acqLabel); dbstack; error('more than one acquisition types, please provide acqLabel'); end
dummy     = cat(1,dummy{:});

outDir = strsplit(outDir{1},filesep); outDir = strjoin(outDir(1:find(ismember(outDir,['sub-' char(subList)]))),filesep);
outDir = fullfile(outDir,'derivatives'); if ~exist(outDir,'dir'); mkdir(outDir); end
outDir = fullfile(outDir,'qaReg');       if ~exist(outDir,'dir'); mkdir(outDir); end
outDir = fullfile(outDir,acqLabel);      if ~exist(outDir,'dir'); mkdir(outDir); end
outDir = fullfile(outDir,ppLabel);       if ~exist(outDir,'dir'); mkdir(outDir); end

% get all files
fBids = {};
fEstimList = {};
fList = {};
sList = {};
mList = {};
tList = {};
for s = 1:length(xSet)
    sList = cat(1,sList,repmat({xSet{s}.ses},size(xSet{s}.fList)));
    tList = cat(1,tList,xSet{s}.initFiles.acqTime);
    switch ppLabel
        case 'final'
            fList = cat(1,fList,xSet{s}.finalFiles.fPreprocList);
            bids = cell(size(xSet{s}.finalFiles.bidsList,1),1);
            for f = 1:size(xSet{s}.finalFiles.bidsList,1)
                bids{f,1} = strjoin(xSet{s}.finalFiles.bidsList(f,:),'_');
            end
            fBids = cat(1,fBids,bids);
        case 'bRun'
            fBase = xSet{s}.brMocoFiles.fBase;
            fEstimList = cat(1,fEstimList,xSet{s}.brMocoFiles.fEstimList);
            fList = cat(1,fList,xSet{s}.brMocoFiles.fMocoList);
            bids = cell(size(xSet{s}.brMocoFiles.bidsList,1),1);
            for f = 1:size(xSet{s}.brMocoFiles.bidsList,1)
                bids{f,1} = strjoin(xSet{s}.brMocoFiles.bidsList(f,:),'_');
            end
            fBids = cat(1,fBids,bids);
        case 'init'
            fList = cat(1,fList,xSet{s}.initFiles.fList);
            bids = cell(size(xSet{s}.finalFiles.bidsList,1),1);
            for f = 1:size(xSet{s}.finalFiles.bidsList,1)
                bids{f,1} = strjoin(xSet{s}.initFiles.bidsList(f,:),'_');
            end
            fBids = cat(1,fBids,bids);
        otherwise
            dbstack; error('x');
    end
    % get mask
    mList = cat(1,mList,xSet{s}.wrMocoFiles.fMaskList);
end
% [~,b] = sort(tList);
% fBids = fBids(b);
% fList = fList(b);
% sList = sList(b);
% mList = mList(b);
% tList = tList(b);

%% Summarize data
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
    tmp = permute(single(mri.vol(:,:,:,dummy(f)+1)),[4 1 2 3]);
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
        MRIwrite(tmp,fFML{f});
    else
        fFML{f} = '';
    end
end

% write catenated means
mriMean(1).vol = cat(4,mriMean.vol); mriMean(2:end) = [];
fMeans = fullfile(outDir,[ppLabel '_means.nii.gz']);
tmp = mriMean; tmp.vol(isnan(tmp.vol)) = 0;
MRIwrite(tmp,fMeans); clear mriMean mri

% % get costs
% costList = cell(length(fList));
% cmd = {srcAfni};
% for r1 = 1:length(fList)
%     for r2 = 1:length(fList)
%         if r2<=r1; continue; end
%         costList{r1,r2} = tempname;
%         cmd{end+1} = '3dAllineate -overwrite \';
%         cmd{end+1} = ['-base   ' fList{r1} ' \'];
%         cmd{end+1} = ['-source ' fList{r2} ' \'];
%         % costList{end+1} = strsplit(fList{r},'['); costList{end} = replace(costList{end}{1},'.nii.gz','.cost');
%         cmd{end+1} = ['-allcostX1D "IDENTITY" ' costList{r1,r2}];
%     end
% end
% [status,cmdout] = system(strjoin(cmd,newline)); if status || contains(cmdout,'error','IgnoreCase',true); dbstack; error(cmdout); error('x'); end
% 
% r1 = 1; r2 = 2;
% fid = fopen(costList{r1,r2});
% tline = {fgetl(fid)}; tline = {fgetl(fid)};
% costLabel = replace(replace(replace(strsplit(char(tline),'___  ___'),'___',''),' ',''),'#','');
% fclose(fid);
% cost = cell(length(fList));
% for r1 = 1:length(fList)
%     for r2 = 1:length(fList)
%         if r2<=r1
%             cost{r1,r2} = zeros(1,1,length(costLabel));
%         else
%             fid = fopen(costList{r1,r2});
%             tline = {fgetl(fid)}; tline = {fgetl(fid)};
%             while ischar(tline{end})
%                 tline{end+1} = fgetl(fid);
%                 tline{end} = strsplit(tline{end},' '); tline{end}(1) = []; tline{end} = str2double(tline{end});
%             end
%             fclose(fid);
%             tline(1) = [];
%             cost{r1,r2} = permute(cat(1,tline{:}),[1 3 2]);
%         end
%     end
% end
% cost = cell2mat(cost);
% for c = 1:size(cost,3)
%     cost(:,:,c) = triu(cost(:,:,c)) + triu(cost(:,:,c),1)';
% end



% correlate frames
im = permute(cat(4,im{:}),[4 1 2 3]);
im = permute(im(:,:),[2 1]);
im(any(isnan(im),2),:) = [];
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

dTick = round(cumsum(nFrames) - nFrames/2);
ax.XTick = dTick;
ax.YTick = dTick;
ax.TickLabelInterpreter = 'none';
ax.XTickLabels = fBids;
ax.YTickLabels = fBids;

% df = size(rho,1)/length(fList);
% ax.TickLabelInterpreter = 'none';
% ax.XTick = df/2:df:size(rho,1);
% ax.YTick = df/2:df:size(rho,1);
% ax.XTickLabels = fBids;
% ax.YTickLabels = fBids;
if length(rho)>length(fList)
    xline(dTick+0.5+nFrames/2)
    yline(dTick+0.5+nFrames/2)
end
% xline(df+0.5:df:size(rho,1))
% yline(df+0.5:df:size(rho,1))
drawnow

out = fullfile(outDir,[ppLabel '_xFrameCorrelation.fig']);
disp(['saving QA to:' newline out])
saveas(gcf,out);





%% Produce visualization
% construct command

%means
cmd = {srcFs};
switch acqLabel
    case 'vfMRI'
        cmd{end+1} = 'fslview -m single \';
        cmd{end+1} = [fMeans ' -b 0,800 &'];
    otherwise
        cmd{end+1} = 'fslview \';
        cmd{end+1} = [fMeans ' &'];
end
cmdMeans = strjoin(cmd,newline);

%first, middle and last frames
cmdFML = cell(size(fFML));
for r = 1:length(fFML)
    cmd = {srcFs};
    switch acqLabel
        case 'vfMRI'
            cmd{end+1} = 'fslview -m single \';
            cmd{end+1} = [strjoin(fFML(r),' -b 0,800 ') ' -b 0,800 &'];
        otherwise
            cmd{end+1} = 'fslview \';
            cmd{end+1} = [strjoin(fFML(r),' ') ' &'];
    end
    cmdFML{r} = strjoin(cmd,newline);
end
cmd = {srcFs};
switch acqLabel
    case 'vfMRI'
        cmd{end+1} = 'fslview -m single \';
        cmd{end+1} = [strjoin(flip(fFML),' -b 0,800 ') ' -b 0,800 &'];
    otherwise
        cmd{end+1} = 'fslview \';
        cmd{end+1} = [strjoin(flip(fFML),' ') ' &'];
end
cmdFML{end+1} = strjoin(cmd,newline);

%all frames
cmdAll = cell(size(fList));
for r = 1:length(fList)
    cmd = {srcFs};
    switch acqLabel
        case 'vfMRI'
            cmd{end+1} = 'fslview -m single \';
            cmd{end+1} = [strjoin(fList(r),' -b 0,800 ') ' -b 0,800 &'];
        otherwise
            cmd{end+1} = 'fslview \';
            cmd{end+1} = [strjoin(fList(r),' ') ' &'];
    end
    cmdAll{r} = strjoin(cmd,newline);
end

cmd = {srcFs};
switch acqLabel
    case 'vfMRI'
        cmd{end+1} = 'fslview -m single \';
        cmd{end+1} = [strjoin(flip(fList),' -b 0,800 ') ' -b 0,800 &'];
    otherwise
        cmd{end+1} = 'fslview \';
        cmd{end+1} = [strjoin(flip(fList),' ') ' &'];
end
cmdAll{end+1} = strjoin(cmd,newline);

%execute commands
if xFlag
    system(cmdMeans);
    system(cmdFML{1});
    system(cmdAll{1});
end




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
