function runSet = initPreproc2(runSet,geomRef,param,skipMask,force,verbose)
global srcAfni

if ~isfield(runSet,'fOrigList'); runSet.fOrigList = []; end
if ~isfield(runSet,'label');         runSet.label = []; end
if ~isfield(runSet,'nDummy');        runSet.label = []; end
if isempty(runSet.fOrigList); runSet.fOrigList = runSet.fList; end
if isempty(runSet.label);         runSet.label = 'pp'; end
if isempty(runSet.nDummy);       runSet.nDummy = 0; end

if ~exist('force'   ,'var');        force    = []       ; end
if ~exist('skipMask','var');        skipMask = []       ; end
if ~isfield(runSet,'dataType'); dataType = {}     ; else; dataType = runSet.dataType; end
if ~isfield(param,'verbose');  verbose = []       ; else; verbose = param.verbose; end
if isempty(force);               force = 0        ; end
if isempty(skipMask);         skipMask = 0        ; end
if isempty(dataType);         dataType = {'volTs'}; elseif ischar(dataType); dataType = {dataType}; end
if isempty(verbose);           verbose = 1        ; end

if ~isfield(param,'duporigin'); duporigin = []    ; else; duporigin = param.duporigin; end
if isempty(duporigin);          duporigin = 0     ; end


% wd
runSet.wd   = fullfile(runSet.info.prcDir,['sub-' runSet.sub],['ses-' runSet.ses],['set-' runSet.label]);
% runSet.wd = fullfile(runSet.wd,['sub-' runSet.sub],['ses-' runSet.ses]);
% runSet.wd = fullfile(runSet.wd,['set-' runSet.label]);
if ~exist(runSet.wd,'dir'); mkdir(runSet.wd); end


% dataType: 'vol' 'volTs' 'singleEcho' 'multiEcho' 'pc'


% geomRef refers to the file to be used as the reference for the slice
% prescription in scanner space
% Hint: this should be a scan with a scanner-space slice presciption that
% matches that of the most important auxiliary scan, so they show nicely in
% freeview. This matters moslty when areg was used.


runSet.acqTime = runSet.date + getAcqTime(runSet.fOrigList(:,1));

[runSet.initFiles,fSort] = doIt(runSet,dataType,geomRef,param,skipMask,force,verbose);
runSet.fList     = runSet.fList(fSort,:);
runSet.date      = runSet.date(fSort);
runSet.nDummy    = runSet.nDummy(fSort);
runSet.fOrigList = runSet.fOrigList(fSort,:);
runSet.fMasks    = runSet.initFiles.fMasks;
% runSet = rmfield(runSet,'date');



function [runSet,fSort] = doIt(runSet,dataType,geomRef,param,skipMask,force,verbose)
global srcAfni

%% %%%%%%%%%%%%%%%
% Massage inputs %
%%%%%%%%%%%%%%% %%
%%% funcSet
%%%% sort by acquisition time
% runSet.acqTime = runSet.date + getAcqTime(runSet.fOrigList);
% runSet.acqTime = getAcqTime(runSet.fOrigList) + caldays(day(runSet.date)-1) + calmonths(month(runSet.date)-1) + calyears(year(runSet.date)-1);
% runSet = rmfield(runSet,'date');
[~,fSort] = sort(runSet.acqTime);
runSet.fOrigList = runSet.fOrigList(fSort,:);
runSet.fList     = runSet.fList(fSort,:);
runSet.acqTime   = runSet.acqTime(fSort);
runSet.nDummy    = runSet.nDummy(fSort);
%%%% bids
[~,bidsList,~] = fileparts(replace(runSet.fOrigList(:),'.nii.gz',''));
bidsList = cellstr(bidsList);
for R = 1:length(bidsList); bidsList{R} = strsplit(bidsList{R},'_'); end; bidsList = padcatcell(bidsList{:}); bidsList(cellfun('isempty',bidsList)) = {''};
runSet.bidsList = bidsList;
%%%% refractor
fOrigList = runSet.fOrigList;
if ~exist('param','var'); param = []; end; if ~isfield(param,'nDummy'); param.nDummy = []; end
if ~isempty(param.nDummy) && ~isempty(runSet.nDummy) && param.nDummy ~= runSet.nDummy; warning(['param.nDummy~=runSet.nDummy' newline 'using param.nDummy']); runSet.nDummy = param.nDummy; end
if isempty(param.nDummy); param.nDummy = runSet.nDummy; end; if isempty(runSet.nDummy); runSet.nDummy = param.nDummy; end

if ~isfield(param,'verbose') || isempty(param.verbose); param.verbose = 0; end


% bidsList = replace({xSet.fList}','.nii.gz',''); for R = 1:length(bidsList); bidsList{R} = strsplit(bidsList{R},'_'); end; bidsList = cat(1,bidsList{:})';
% dIn = xSet.files(1).folder;
% fOrig = fullfile(dIn,{xSet.files.name}');
% label = xSet.label;
% bidsList = replace({xSet.files.name}','.nii.gz',''); for R = 1:length(bidsList); bidsList{R} = strsplit(bidsList{R},'_'); end; bidsList = cat(1,bidsList{:})';
% %%% dOut
% dOut = fullfile(dOut,label); if ~exist(dOut,'dir'); mkdir(dOut); end
%%% geomRef
if ~exist('geomRef','var'); geomRef = []; end
if isempty(geomRef)
    if isfield(runSet,'fGeom') && ~isempty(runSet.fGeom)
        fRef = runSet.fGeom;
    else
        fRef = fOrigList{1};
    end
else
    if ~ischar(geomRef)
        fRef = fOrigList{geomRef};
    else
        fRef = geomRef;
    end
end
runSet.fGeom = fRef;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rewrite data at plumb without dummies %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%% Get oblique and plumb reference files
disp([' writing oblique and plumb references (for ' num2str(length(fOrigList)) ' files)'])
forceThis   = force;
verboseThis = verbose;
runSet = writeObliqueAndPlumbRef(runSet,forceThis,verboseThis);
%%% Get plumb timeseries files
if all(diff(runSet.nDummy)==0)
    disp([' rewriting data at plumb (for ' num2str(length(fOrigList)) ' files, removing ' num2str(runSet.nDummy(1)) ' dummies)'])
else
    disp([' rewriting data at plumb (for ' num2str(length(fOrigList)) ' files, removing ' num2str(runSet.nDummy') ' dummies)'])
end

forceThis   = force;
verboseThis = verbose;
runSet = writeObliqueAndPlumb(runSet,forceThis,verboseThis);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare data for motion estimation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
nRun = size(runSet.fPlumbList,1);
% if isfield(runSet,'nRun')
%     nRun = runSet.nRun;
% else
%     nRun = contains(squeeze(bidsList(1,:)),'run-');
%     if ~any(nRun); nRun = 1; else; nRun = length(unique(bidsList(nRun,:))); end
% end
% if isfield(runSet,'nEcho')
%     nEcho = runSet.nEcho;
% else
%     nRun = contains(squeeze(bidsList(:,1,1)),'run-');
%     if ~any(nRun); nRun = 1; else; nRun = length(unique(bidsList(nRun,:))); end
% end
% 
% % if isfield(xSet,'nVencDat')
% %     nVencDat = xSet.nVencDat;
% % else
%     nVencDat = length(fPlumb)/nEcho/nRun;
% % end
if isfield(runSet,'nEcho')
    nEcho = runSet.nEcho;
else
    nEcho = contains(bidsList,'echo-'); nEcho = unique(bidsList(nEcho,:));
    if isempty(nEcho)
        nEcho = 1;
    %     dataType{end+1} = 'singleEcho';
    else
        nEcho = length(nEcho);
        dbstack; error('double-check that')
        % dataType{end+1} = 'multiEcho';
    end
end
if nEcho==1
    dataType{end+1} = 'singleEcho';
else
    dataType{end+1} = 'multiEcho';
end
disp([dataType{end} ' data detected'])
runSet.dataType = dataType;

%% Simple single-echo functional timeseries
if all(ismember({'singleEcho' 'volTs'},dataType)) && ~any(ismember({'pc'},dataType))
    runSet.fEstimList = runSet.fPlumbList;
end

%% Multi-echo (cross-echo rms images for later motion/distortion estimation)
if all(ismember({'multiEcho' 'volTs'},dataType)) && ~any(ismember({'pc'},dataType))
    dbstack; error('code that')
end

%% Phase contrast data
if all(ismember({'pc'},dataType)) && ~any(ismember({''},dataType))
    dbstack; error('code that')
end




%% %%%%%%%%%%%%%%%
% Summarize data %
%%%%%%%%%%%%%%% %%
forceThis   = force;
verboseThis = verbose;
runSet.fPlumbSmr = summarizeVolTs2(runSet.fPlumbList,[],[],param.nDummy,[],forceThis,verboseThis);



%% %%%%%%%%%%%%%%%%%%
% Draw preproc mask %
%%%%%%%%%%%%%%%%%% %%

if ~skipMask
    forceThis = force;
    useSynth  = [];
    % if strcmp(runSet.label,'bold')
    %     useSynth = 1;
    % end

    [~,bidsSubFolder] = fileparts(fileparts(runSet.fOrigList{1}));
    fMaskInv = fullfile(runSet.info.bidsDir,'derivatives','manual',['sub-' runSet.info.sub],['ses-' runSet.info.ses],bidsSubFolder,['set-' runSet.label]);
    fMaskInv = fullfile(fMaskInv,['sub-' runSet.info.sub '_ses-' runSet.info.ses '_set-' runSet.label '_desc-forWRmoco_mask.nii.gz']);

    if exist(fMaskInv,'file') && ~forceThis
        %%% Get brain mask for preprocessing from bids derivatives if it exists
        runSet.fMasks.fUlay    = char(runSet.fPlumbSmr.sesCat.runAv.fList);
        runSet.fMasks.fMaskInv = replace(runSet.fMasks.fUlay,'_volTs.nii.gz','_volBrainMaskInv.nii.gz');
        copyfile(fMaskInv,runSet.fMasks.fMaskInv)
        runSet.fMasks.fMask    = replace(runSet.fMasks.fUlay,'_volTs.nii.gz','_volBrainMask.nii.gz'   );
        cmd = {srcAfni};
        cmd{end+1} = '3dcalc -overwrite \';
        cmd{end+1} = ['-prefix ' runSet.fMasks.fMask ' \'];
        cmd{end+1} = ['-a ' runSet.fMasks.fMaskInv ' \'];
        cmd{end+1} = '-expr ''-(a-1)''';
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end        

    else
        %%% Draw brain mask for preprocessing if it does not exist in bids derivatives
        runSet.fMasks = drawMask(char(runSet.fPlumbSmr.sesCat.runAv.fList),useSynth,forceThis);
        if ~exist(fileparts(fMaskInv),'dir'); mkdir(fileparts(fMaskInv)); end
        copyfile(runSet.fMasks.fMaskInv,fMaskInv)
    end




else
    runSet.fMasks = [];
end




    
return







%%% Detect multi-echo data if not specified in dataType


if nRun==1 && all(ismember({'singleEcho' 'pc' 'vol'},dataType))
    

    fOrigList = reshape(fOrigList,nRun,nEcho,nVencDat);
    fPlumb = reshape(fPlumb,nRun,nEcho,nVencDat);
    fEstim = fPlumb(:,:,contains(fPlumb(1,1,:),'proc-venc0_') & contains(fPlumb(1,1,:),'part-mag'));
    if nFrame>1
        fPlumbAv = reshape(fPlumbAv,nRun,nEcho,nVencDat);
        fEstimAv = fPlumbAv(:,:,contains(fPlumbAv(1,1,:),'proc-venc0_') & contains(fPlumbAv(1,1,:),'part-mag'));
        % fEstim = reshape(fPlumbAv,nRun,nEcho,nVencDat);   
    end
    % fEstim = fPlumb(contains(fPlumb,'proc-venc0_') & contains(fPlumb,'part-mag'));
    if size(fEstim,1)~=size(fOrigList,1); dbstack; error('double-check that'); end
    % bidsList;
    
    

elseif all(ismember({'singleEcho'},dataType))
    %%% Single echo data
    %%%% Write means (no rms) for later visualization
    if nRun>1
        cmd = {srcFs};
        if nFrame>1
            fIn = fPlumbAv;
        else
            fIn = fPlumb;
        end
        if any(contains(squeeze(bidsList(:,1,1)),'run-'))
            fOut = replace(fIn{1},bidsList{contains(squeeze(bidsList(:,1,1)),'run-'),1,1},'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        else
            fOut = fIn{1};
            if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        end
        fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
        if forceRewriteAtPlumb || ~exist(fOut,'file')
            cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
        end
        fPlumbCat = fOut;
        fIn = fOut;
        fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
        if forceRewriteAtPlumb || ~exist(fOut,'file')
            cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
        end
        fPlumbAvCat = fOut;

        disp(' averaging')
        if length(cmd)>1
            cmd = strjoin(cmd,newline); % disp(cmd)
            [status,cmdout] = system(cmd); if status; dbstack; error(cmdout); error('x'); end
            disp('  done')
        else
            disp('  already done, skipping')
        end
    end
    
    %%% Set fEstim
    fEstim = fPlumb;
    if nFrame>1
        fEstimAv = fPlumbAv;
    end
    if nRun>1
        fEstimCat = fPlumbCat;
        fEstimAvCat = fPlumbAvCat;
    end

elseif all(ismember({'multiEcho'},dataType))
    
    %%% Multi-echo data
    
    %%%% Reshape lists with echos in different columns
    echoList = bidsList(contains(bidsList(:,1,1),'echo-'),:,:);
    echoUniqueList = unique(echoList);
    
    bidsList = permute(reshape(bidsList,size(bidsList,1),length(echoUniqueList),length(echoList)/length(echoUniqueList)),[1 3 2]);
    fPlumb = reshape(fPlumb,length(echoUniqueList),length(echoList)/length(echoUniqueList))';
    if nFrame>1
        fPlumbAv = reshape(fPlumbAv,length(echoUniqueList),length(echoList)/length(echoUniqueList))';
    end
    fOrigList = reshape(fOrigList,length(echoUniqueList),length(echoList)/length(echoUniqueList))';

    disp('---------------------------------------------------')
    disp('---------------------------------------------------')
    warning(['Make sure echoes are properly sorted below' ...
        newline '(runs across rows and echos across columns)'])
    if size(bidsList,2)>1
        disp(fullfile(squeeze(bidsList(contains(bidsList(:,1,1),'run-'),:,:)),squeeze(bidsList(contains(bidsList(:,1,1),'echo-'),:,:))))
    else
        disp(fullfile(permute((bidsList(contains(bidsList(:,1,1),'echo-'),:,:)),[1 3 2])));
    end
    disp('---------------------------------------------------')
    disp('---------------------------------------------------')

    %%%% Write means (before rms) for later visualization
    if nRun>1
        cmd = {srcFs};
        if nFrame>1
            fPlumbCat = cell(1,size(fPlumbAv,2));
            fPlumbAvCat = cell(1,size(fPlumbAv,2));
        else
            fPlumbCat = cell(1,size(fPlumb,2));
            fPlumbAvCat = cell(1,size(fPlumb,2));
        end
        for E = 1:nEcho
            if nFrame>1
                fIn = fPlumbAv(:,E);
            else
                fIn = fPlumb(:,E);
            end
            fOut = replace(fIn{1},bidsList{contains(squeeze(bidsList(:,1,1)),'run-'),1,1},'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
            fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
            if forceRecomputeRMS || ~exist(fOut,'file')
                cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
            end
            fPlumbCat{E} = fOut;
            fIn = fOut;
            fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
            if forceRecomputeRMS || ~exist(fOut,'file')
                cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
            end
            fPlumbAvCat{E} = fOut;
        end
        disp(' averaging')
        if length(cmd)>1
            cmd = strjoin(cmd,newline); % disp(cmd)
            [status,cmdout] = system(cmd); if status; dbstack; error(cmdout); error('x'); end
            disp('  done')
        else
            disp('  already done, skipping')
        end
        disp(' non need for averaging (single-run)')
    end

    %%%% Compute and write rms
    disp('computing cross-echo rms')
    fEstim = cell(size(fPlumb,1),1);
    for R = 1:size(fPlumb,1)
        disp([' file' num2str(R) '/' num2str(size(fPlumb,1))])
        fOut = replace(fPlumb{R,1},'echo-1','echo-rms');
        if forceRecomputeRMS || ~exist(fOut,'file')
            [a,~] = fileparts(fOut); if ~exist(a,'dir'); mkdir(a); end

            mri = MRIread(fPlumb{R,1},1);
            mri.vol = zeros([mri.volsize mri.nframes]);
            for ii = 1:size(fPlumb,2)
                mriTmp = MRIread(fPlumb{R,ii});
                mri.vol = mri.vol + mriTmp.vol.^2;
            end
            mri.vol = sqrt(mri.vol./size(fPlumb,2));

            MRIwrite(mri,fOut);
            disp('  done')
        else
            disp('  already done, skipping')
        end
        fEstim{R} = fOut;
    end


    %%%% Write means (after rms) for later visualization
    cmd = {srcFs};
    if nFrame>1
        fEstimAv = cell(size(fEstim));
    end
    for R = 1:length(fEstim)
        fIn = fEstim{R};
        if nFrame>1
            fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
            if forceRecomputeRMS || ~exist(fOut,'file')
                cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
            end
            fEstimAv{R} = fOut;
        % else
        %     fEstimAv{R} = fIn;
        end
    end
    if nRun>1
        if nFrame>1
            fIn = fEstimAv;
        else
            fIn = fEstim;
        end
        fOut = replace(fIn{1},bidsList{contains(squeeze(bidsList(:,1,1)),'run-'),1,1},'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
        if forceRecomputeRMS || ~exist(fOut,'file')
            cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
        end
        fEstimCat = fOut;
        fIn = fOut;
        fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
        if forceRecomputeRMS || ~exist(fOut,'file')
            cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
        end
        fEstimAvCat = fOut;
    else
        fEstimCat = '';
        fEstimAvCat = '';
    end
    disp(' averaging')
    if length(cmd)>1
        cmd = strjoin(cmd,newline); % disp(cmd)
        [status,cmdout] = system(cmd); if status; dbstack; error(cmdout); error('x'); end
    else
        disp('  already done, skipping')
    end
end



%% Data on which to apply motion correction
fApply = fPlumb;



%% QA: Visualize motion
%%% between-run
if nRun>1 && ~isempty(fEstimCat)
    cmd = {srcFs};
    cmd{end+1} = ['fslview -m single ' fEstimCat ' &'];
    cmd = strjoin(cmd,newline); % disp(cmd)
    if param.verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    fFslviewBR = cmd;
% else
%     fFslviewBR = '';
end
%%% within-run
if nFrame>1
    cmd = {srcFs};
    cmd{end+1} = ['fslview -m single ' strjoin(fEstim,' ') ' &'];
    cmd = strjoin(cmd,newline); % disp(cmd)
    if param.verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    fFslviewWR = cmd;
    %%% within-run, first, middle and last frames
    fFslviewWRfstMdLst = qaFstMdLst(fEstim,forceRewriteAtPlumb || forceRecomputeRMS,param.verbose);
% else
%     fFslviewWR = '';
%     fFslviewWRfstMdLst = '';
end

if nRun>1
    runSet.initFiles.qaFiles.fFslviewBR = fFslviewBR; end
if nFrame>1
    runSet.initFiles.qaFiles.fFslviewWR = fFslviewWR;
    runSet.initFiles.qaFiles.fFslviewWRfstMdLst = fFslviewWRfstMdLst; end



%% Low SNR data (estimate motion on temporally smoothed data)


%%% Detect low SNR data it not specified in dataType
if any(ismember({'volTs'},dataType))
    if ~any(ismember({'lowSNR' 'highSNR'},dataType))
        dbstack; error(['automatic detection of lowSNR data not implemented' newline 'specify ''lowSNR'' or ''highSNR'' in dataType '])
    else
        disp([dataType{ismember(dataType,{'lowSNR' 'highSNR'})} ' data specified'])
    end
end


%%% Temporally smooth lowSNR data
if all(ismember({'volTs' 'lowSNR'},dataType))

    %%%% Visualize SNR before smoothing
    fEstimBeforeSm = fEstim{1};
    volTs = MRIread(fEstimBeforeSm);
    volAv = mean(volTs.vol,4);
    volEr = std(volTs.vol,[],4);
    if param.verbose
        hFig = figure('WindowStyle','docked');
    else
        hFig = figure('Visible','off','units','normalized','outerposition',[0 0 1 1]);
    end
    hTile = tiledlayout(2,3); hTile.TileSpacing = 'tight'; hTile.Padding = 'tight';
    ax = {};
    ax{end+1} = nexttile;
    imagesc(volAv(:,:,ceil(end/2)));
    ax{end}.DataAspectRatio = [1 1 1];
    ax{end}.XTick = []; ax{end}.YTick = [];
    ylabel(colorbar,'mean'); ax{end}.Colormap = gray;
    ylabel('without smoothing')
    ax{end+1} = nexttile;
    imagesc(volEr(:,:,ceil(end/2)));
    ax{end}.DataAspectRatio = [1 1 1];
    ax{end}.XTick = []; ax{end}.YTick = [];
    ylabel(colorbar,'std'); ax{end}.Colormap = gray;
    ax{end+1} = nexttile;
    imagesc(volAv(:,:,ceil(end/2))./volEr(:,:,ceil(end/2)));
    ax{end}.DataAspectRatio = [1 1 1];
    ax{end}.XTick = []; ax{end}.YTick = [];
    ylabel(colorbar,'tSNR'); ax{end}.Colormap = gray;
    if isempty(param) || ~isfield(param,'tSmWin_vol') || isempty(param.tSmWin_vol)
        dbstack; error('when dataType is lowSNR, must specify param.tSmWin_vol')
    end
    
    %%%% Smoothing
    disp('temporally smooth timeseries')
    cmd = {srcAfni};
    for R = 1:length(fEstim)
        fIn = fEstim{R};
        [d,fOut,~] = fileparts(replace(fIn,'.nii.gz',''));
        fOut = fullfile(d,[strjoin({['sm' num2str(param.tSmWin_vol)] fOut},'_') '.nii.gz']);
        cmd{end+1} = ['echo '' ''file' num2str(R) '/' num2str(length(fEstim))];
        if forceResmoothing || ~exist(fOut,'file')
            cmd{end+1} = '3dTsmooth -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = ['-hamming ' num2str(param.tSmWin_vol) ' \'];
            if verbose
                cmd{end+1} = fIn;
            else
                cmd{end+1} = [fIn ' > /dev/null 2>&1'];
            end
        else
            cmd{end+1} = ['echo ''  ''already smoothed, skipping'];
        end
        fEstim{R} = fOut;
    end
    cmd = strjoin(cmd,newline); % disp(cmd)
    [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end

    %%%% Visualize SNR after smoothing
    fEstimAfterSm = fEstim{1};
    volTs = MRIread(fEstimAfterSm);
    volAv = mean(volTs.vol,4);
    volEr = std(volTs.vol,[],4);
    ax{end+1} = nexttile;
    imagesc(volAv(:,:,ceil(end/2)));
    ax{end}.DataAspectRatio = [1 1 1];
    ax{end}.XTick = []; ax{end}.YTick = [];
    ylabel(colorbar,'mean'); ax{end}.Colormap = gray;
    ylabel(['with smoothing (' num2str(param.tSmWin_vol) 'volumes)'])
    ax{end+1} = nexttile;
    imagesc(volEr(:,:,ceil(end/2)));
    ax{end}.DataAspectRatio = [1 1 1];
    ax{end}.XTick = []; ax{end}.YTick = [];
    ylabel(colorbar,'std'); ax{end}.Colormap = gray;
    ax{end+1} = nexttile;
    imagesc(volAv(:,:,ceil(end/2))./volEr(:,:,ceil(end/2)));
    ax{end}.DataAspectRatio = [1 1 1];
    ax{end}.XTick = []; ax{end}.YTick = [];
    ylabel(colorbar,'tSNR'); ax{end}.Colormap = gray;
    for R = 1:3
        cLim = [min([ax{R}.CLim(1) ax{R+3}.CLim(1)]) max([ax{R}.CLim(2) ax{R+3}.CLim(2)])];
        ax{R}.CLim = cLim;
        ax{R+3}.CLim = cLim;
    end
    [~,fSort,~] = fileparts(fileparts(fEstimAfterSm));
    title(hTile,fSort,'interpreter','none')
    fFigSm = replace(fEstimAfterSm,'.nii.gz','.fig');
    set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    saveas(hFig,fFigSm)
    
    %%%% Visualize SNR before and after smoothing
    cmd = {srcFs};
    cmd{end+1} = ['fslview -m single ' fEstimBeforeSm ' ' fEstimAfterSm ' &'];
    cmd = strjoin(cmd,newline); % disp(cmd)
    fFslviewSm = cmd;
    if param.verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end

    %%%% Write means for later visualization
    cmd = {srcFs};
    fEstimAv = cell(size(fEstim));
    for R = 1:length(fEstim)
        fIn = fEstim{R};
        fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
        if forceResmoothing || ~exist(fOut,'file')
            cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
        end
        fEstimAv{R} = fOut;
    end
    fIn = fEstimAv;
    fOut = replace(fIn{1},bidsList{contains(squeeze(bidsList(:,1,1)),'run-'),1,1},'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
    fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
    if forceResmoothing || ~exist(fOut,'file')
        cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
    end
    fIn = fOut;
    fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
    if forceResmoothing || ~exist(fOut,'file')
        cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
    end
    disp(' averaging')
    if length(cmd)>1
        cmd = strjoin(cmd,newline); % disp(cmd)
        [status,cmdout] = system(cmd); if status; dbstack; error(cmdout); error('x'); end
    else
        disp('  already done, skipping')
    end

    %% QA: Visualize motion (with smoothing)
    %%% between-run
    cmd = {srcFs};
    cmd{end+1} = ['fslview -m single ' fEstimCat ' &'];
    cmd = strjoin(cmd,newline); % disp(cmd)
    if param.verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    fFslviewBRsm = cmd;
    %%% within-run
    cmd = {srcFs};
    cmd{end+1} = ['fslview -m single ' strjoin(fEstim,' ') ' &'];
    cmd = strjoin(cmd,newline); % disp(cmd)
    if param.verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    fFslviewWRsm = cmd;
    %%% within-run, first, middle and last frames
    fFslviewWRsmFstMdLst = qaFstMdLst(fEstim,forceRewriteAtPlumb||forceRecomputeRMS,param.verbose);

    runSet.initFiles.qaFilesSm.fFslviewBR = fFslviewBRsm;
    runSet.initFiles.qaFilesSm.fFslviewWR = fFslviewWRsm;
    runSet.initFiles.qaFilesSm.fFslviewWRfstMdLst = fFslviewWRsmFstMdLst;


end


%% Sumarize outputs
if ~exist('nRun','var') || isempty(nRun)
    nRun = size(fOrigList,1); end
if ~exist('nEcho','var') || isempty(nEcho)
    nEcho = size(fOrigList,2); end
% if any(ismember(dataType,'pc')) && (~exist('nVencDat','var') || isempty(nVencDat))
%     nVencDat = size(fOrig,3); end

runSet.dataType = dataType;
if ~isfield(runSet,'nRun');  runSet.nRun = nRun; end
if ~isfield(runSet,'nEcho'); runSet.nEcho = nEcho; end
if ~isfield(runSet,'nVenc'); runSet.nRun = nVenc; end
if ~isfield(runSet,'nVencDat'); runSet.nRun = nVencDat; end
if ~isfield(runSet,'nFrame'); runSet.nFrame = nFrame; end
runSet.initFiles.nRun   = runSet.nRun;
runSet.initFiles.nEcho  = runSet.nEcho;
runSet.initFiles.nVenc  = runSet.nVenc;
runSet.initFiles.nVencDat  = runSet.nVencDat;
runSet.initFiles.nFrame = runSet.nFrame;
% if exist('nVencDat','var')
%     xSet.initFiles.nVencDat = nVencDat; end
runSet.initFiles.fOrig = fOrigList;

runSet.initFiles.fPlumbRef = fPlumbRef;
runSet.initFiles.fPlumb = fPlumb;
% if exist('nVencDat','var')
if nFrame>1
    runSet.initFiles.fPlumbAv = fPlumbAv; end
if nRun>1
    runSet.initFiles.fPlumbCatAv = fPlumbCat;
    runSet.initFiles.fPlumbAvCatAv = fPlumbAvCat; end

runSet.initFiles.fEstim = fEstim;
if nFrame>1
    runSet.initFiles.fEstimAv = fEstimAv; end
if nRun>1
    if iscell(fEstimCat)
        runSet.initFiles.fEstimCatAv = fEstimCat;
    else
        runSet.initFiles.fEstimCatAv = {fEstimCat};
    end
    if iscell(fEstimAvCat)
        runSet.initFiles.fEstimAvCatAv = fEstimAvCat;
    else
        runSet.initFiles.fEstimAvCatAv = {fEstimAvCat};
    end
end

runSet.initFiles.fApply = fApply;

runSet.initFiles.fObliqueRef = fObliqueRef;
runSet.initFiles.fPlumbRef = fPlumbRef;

%%% QA
if exist('fFigSm','var')
    runSet.initFiles.qaFilesSm.fFigSm = fFigSm;
end
if exist('fFslviewSm','var')
    runSet.initFiles.qaFilesSm.fFslviewSm = fFslviewSm;
end
