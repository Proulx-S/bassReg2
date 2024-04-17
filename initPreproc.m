function xSet = initPreproc(xSet,dOut,geomRef,dataType,param,force)
if ~exist('force','var');    force    = []; end
if ~exist('dataType','var'); dataType = {}; end % 'vol' or 'volTs'

% override dataType
if isfield(xSet,'dataType'); dataType = xSet.dataType; end

if isempty(force);    force    = 0        ; end
if isempty(dataType); dataType = {'volTs'}; elseif ischar(dataType); dataType = {dataType}; end





forceRewriteAtPlumb = force;
forceRecomputeRMS = force;
forceResmoothing = force;
global srcFs srcAfni

% geomRef refers to the file to be used as the reference for the slice
% prescription in scanner space
% Hint: this should be a scan with a scanner-space slice presciption that
% matches that of the most important auxiliary scan, so they show nicely in
% freeview. This matters moslty when areg was used.




%% Massage inputs
%%% funcSet
%%%% sort by acquisition time
[~,b] = sort(getAcqTime(fullfile({xSet.files.folder},{xSet.files.name})));
xSet.files = xSet.files(b);
%%%% refactor
dIn = xSet.files(1).folder;
fOrig = fullfile(dIn,{xSet.files.name}');
label = xSet.label;
bidsList = replace({xSet.files.name}','.nii.gz',''); for i = 1:length(bidsList); bidsList{i} = strsplit(bidsList{i},'_'); end; bidsList = cat(1,bidsList{:})';
%%% dOut
dOut = fullfile(dOut,label); if ~exist(dOut,'dir'); mkdir(dOut); end
%%% geomRef
if ~exist('geomRef','var'); geomRef = []; end
if isempty(geomRef)
    if isfield(xSet,'fGeom') && ~isempty(xSet.fGeom)
        fRef = xSet.fGeom;
    else
        fRef = fOrig{1};
    end
else
    if ~ischar(geomRef)
        fRef = fOrig{geomRef};
    else
        fRef = geomRef;
    end
end
xSet.fGeom = fRef;
% param
if ~exist('param','var'); param = []; end
if ~isfield(param,'nDummy') || isempty(param.nDummy); param.nDummy = 0; end
if ~isfield(param,'verbose') || isempty(param.verbose); param.verbose = 1; end


%% Rewrite data at plumb without dummies
disp('Rewriting data at plumb (consider using 3drefit -deobllique "https://discuss.afni.nimh.nih.gov/t/3drefit-deoblique/3303/2")')

%%% Get oblique and plumb files
disp([' writing oblique and plumb references (for ' num2str(length(fOrig)) ' files)'])
cmd = {srcAfni};
%%%% reference
fIn = fRef;
nframes = MRIread(fIn,1); nframes = nframes.nframes; if nframes>8; nframes = 8; end
fPlumbRef = fullfile(dOut,'setPlumb_volRef.nii.gz');
if forceRewriteAtPlumb || ~exist(fPlumbRef,'file')
    cmd{end+1} = '3dNwarpApply -overwrite \';
    cmd{end+1} = ['-nwarp ''IDENT(' fIn ')'' \'];
    cmd{end+1} = ['-prefix ' fPlumbRef ' \'];
    cmd{end+1} = ['-source ' fIn '[0..' num2str(nframes-1) ']'];
end
fObliqueRef = fullfile(dOut,'setOblique_volRef.nii.gz');
if forceRewriteAtPlumb || ~exist(fObliqueRef,'file')
    cmd{end+1} = '3dcalc -overwrite \';
    cmd{end+1} = ['-a ' fIn '[0..' num2str(nframes-1) '] \'];
    cmd{end+1} = ['-expr a \'];
    cmd{end+1} = ['-prefix ' fObliqueRef];
end
%%%% individual files
for i = 1:length(fOrig)
    fIn = fOrig{i};
    nframes = MRIread(fIn,1); nframes = nframes.nframes; if nframes>8; nframes = 8; end
    fOut = replace(fIn,dIn,dOut); fOut = replace(fOut,'.nii.gz',''); if ~exist(fOut,'dir'); mkdir(fOut); end;
    fOut = fullfile(fOut,'runPlumb_volRef.nii.gz');
    if forceRewriteAtPlumb || ~exist(fOut,'file')
        cmd{end+1} = '3dNwarpApply -overwrite \';
        cmd{end+1} = ['-nwarp ''IDENT(' fIn ')'' \'];
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = ['-source ' fIn '[0..' num2str(nframes-1) ']'];
    end
    fOut = replace(fOut,'runPlumb_volRef.nii.gz','runOblique_volRef.nii.gz');
    if forceRewriteAtPlumb || ~exist(fOut,'file')
        cmd{end+1} = '3dcalc -overwrite \';
        cmd{end+1} = ['-a ' fIn '[0..' num2str(nframes-1) '] \'];
        cmd{end+1} = ['-expr a \'];
        cmd{end+1} = ['-prefix ' fOut];
    end
end
%%%% launch command
if length(cmd)>1
    cmd = strjoin(cmd,newline); % disp(cmd)
    [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
end

%%% Rewrite at plumb (and means for later visualization)
disp(' writing data (and means for later visualization)')
mri_setPlumbRef = MRIread(fPlumbRef,1);
fPlumb = cell(size(fOrig));
fPlumbAv = cell(size(fOrig));
for i = 1:numel(fOrig)
    disp(['  file' num2str(i) '/' num2str(numel(fOrig))])
    fIn = fOrig{i};
    fOut = replace(replace(fOrig{i},'.nii.gz',''),dIn,dOut); if ~exist(fOut,'dir'); mkdir(fOut); end
    if any(ismember(dataType,'volTs'))
        fOut_plumb = fullfile(fOut,'setPlumb_volTs.nii.gz');
        fOut_av_plumb = fullfile(fOut,'av_setPlumb_volTs.nii.gz');
        fOut_av_oblique = fullfile(fOut,'av_runOblique_volTs.nii.gz'); clear fOut
    elseif any(ismember(dataType,'vol'))
        fOut_plumb = fullfile(fOut,'setPlumb_vol.nii.gz');
        fOut_av_plumb = fullfile(fOut,'av_setPlumb_vol.nii.gz');
        fOut_av_oblique = fullfile(fOut,'av_runOblique_vol.nii.gz'); clear fOut
    else
        dbstack; error(['don''t know which suffix to give' newline 'please set ''vol'' or ''volTs'' as dataType'])
    end
    if forceRewriteAtPlumb || ...
            any([~exist(fOut_plumb,'file') ~exist(fOut_av_plumb,'file') ~exist(fOut_av_oblique,'file')])
        mri_runOblique = MRIread(fIn);
        mri_runOblique.vol = mri_runOblique.vol(:,:,:,param.nDummy+1:end);
        mri_setPlumbRef.vol = mri_runOblique.vol;
        MRIwrite(mri_setPlumbRef,fOut_plumb);
        mri_setPlumbRef.vol = mean(mri_setPlumbRef.vol,4);
        MRIwrite(mri_setPlumbRef,fOut_av_plumb);
        mri_runOblique.vol = mean(mri_runOblique.vol,4);
        MRIwrite(mri_runOblique,fOut_av_oblique);
        disp('   done')
    else
        disp('   already done, skipping')
    end
    fPlumb{i} = fOut_plumb;
    fPlumbAv{i} = fOut_av_plumb;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare data for motion estimation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Multi-echo (cross-echo rms images for later motion/distortion estimation)

%%% Detect multi-echo data if not specified in dataType
tmp = bidsList(contains(bidsList(:,1,1),'echo-'),:,:);
if ~isempty(tmp) && length(unique(tmp))>1
    % if nnz(contains(fPlumb,'echo-'))==length(fPlumb)
    dataType{end+1} = 'multiEcho';
    nEcho = length(tmp);
else
    dataType{end+1} = 'singleEcho';
    nEcho = 1;
end
disp([dataType{end} ' data detected'])
% if ~any(ismember({'multiEcho' 'singleEcho'},dataType))
%     tmp = bidsList(contains(bidsList(:,1,1),'echo-'),:,:);
%     if ~isempty(tmp) && length(unique(tmp))>1
%     % if nnz(contains(fPlumb,'echo-'))==length(fPlumb)
%         dataType{end+1} = 'multiEcho';
%     else
%         dataType{end+1} = 'singleEcho';
%     end
%     disp([dataType{end} ' data detected'])
% else
%     disp([dataType{ismember(dataType,{'multiEcho' 'singleEcho'})} ' data specified'])
% end

if ~any(ismember({'multiEcho'},dataType))
    %%% Single echo data
    %%%% Write means (no rms) for later visualization
    cmd = {srcFs};
    fIn = fPlumbAv;
    fOut = replace(fIn{1},bidsList{contains(squeeze(bidsList(:,1,1)),'run-'),1,1},'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
    fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
    if forceRewriteAtPlumb || ~exist(fOut,'file')
        cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
    end
    fPlumbCatAv = fOut;
    fIn = fOut;
    fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
    if forceRewriteAtPlumb || ~exist(fOut,'file')
        cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
    end
    fPlumbAvCatAv = fOut;
    
    disp(' averaging')
    if length(cmd)>1
        cmd = strjoin(cmd,newline); % disp(cmd)
        [status,cmdout] = system(cmd); if status; dbstack; error(cmdout); error('x'); end
        disp('  done')
    else
        disp('  already done, skipping')
    end
    %%% Set fEstim
    fEstim = fPlumb;
    fEstimAv = fPlumbAv;
    fEstimCatAv = fPlumbCatAv;
    fEstimAvCatAv = fPlumbAvCatAv;

else
    
    %%% Multi-echo data
    
    %%%% Reshape lists with echos in different columns
    echoList = bidsList(contains(bidsList(:,1,1),'echo-'),:,:);
    echoUniqueList = unique(echoList);
    
    bidsList = permute(reshape(bidsList,size(bidsList,1),length(echoUniqueList),length(echoList)/length(echoUniqueList)),[1 3 2]);
    fPlumb = reshape(fPlumb,length(echoUniqueList),length(echoList)/length(echoUniqueList))';
    fPlumbAv = reshape(fPlumbAv,length(echoUniqueList),length(echoList)/length(echoUniqueList))';
    fOrig = reshape(fOrig,length(echoUniqueList),length(echoList)/length(echoUniqueList))';

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
    cmd = {srcFs};
    fPlumbCatAv = cell(1,size(fPlumbAv,2));
    fPlumbAvCatAv = cell(1,size(fPlumbAv,2));
    if any(contains(squeeze(bidsList(:,1,1)),'run-'))
        for e = 1:size(fPlumbAv,2)
            fIn = fPlumbAv(:,e);
            fOut = replace(fIn{1},bidsList{contains(squeeze(bidsList(:,1,1)),'run-'),1,1},'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
            fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
            if forceRecomputeRMS || ~exist(fOut,'file')
                cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
            end
            fPlumbCatAv{e} = fOut;
            fIn = fOut;
            fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
            if forceRecomputeRMS || ~exist(fOut,'file')
                cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
            end
            fPlumbAvCatAv{e} = fOut;
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
    for i = 1:size(fPlumb,1)
        disp([' file' num2str(i) '/' num2str(size(fPlumb,1))])
        fOut = replace(fPlumb{i,1},'echo-1','echo-rms');
        if forceRecomputeRMS || ~exist(fOut,'file')
            [a,~] = fileparts(fOut); if ~exist(a,'dir'); mkdir(a); end

            mri = MRIread(fPlumb{i,1},1);
            mri.vol = zeros([mri.volsize mri.nframes]);
            for ii = 1:size(fPlumb,2)
                mriTmp = MRIread(fPlumb{i,ii});
                mri.vol = mri.vol + mriTmp.vol.^2;
            end
            mri.vol = sqrt(mri.vol./size(fPlumb,2));

            MRIwrite(mri,fOut);
            disp('  done')
        else
            disp('  already done, skipping')
        end
        fEstim{i} = fOut;
    end


    %%%% Write means (after rms) for later visualization
    cmd = {srcFs};
    fEstimAv = cell(size(fEstim));
    for i = 1:length(fEstim)
        fIn = fEstim{i};
        if nframes>1
            fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
            if forceRecomputeRMS || ~exist(fOut,'file')
                cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
            end
            fEstimAv{i} = fOut;
        else
            fEstimAv{i} = fIn;
        end
    end
    if any(contains(squeeze(bidsList(:,1,1)),'run-'))
        fIn = fEstimAv;
        fOut = replace(fIn{1},bidsList{contains(squeeze(bidsList(:,1,1)),'run-'),1,1},'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
        if forceRecomputeRMS || ~exist(fOut,'file')
            cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
        end
        fEstimCatAv = fOut;
        fIn = fOut;
        fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
        if forceRecomputeRMS || ~exist(fOut,'file')
            cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
        end
        fEstimAvCatAv = fOut;
    else
        fEstimCatAv = '';
        fEstimAvCatAv = '';
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
if ~isempty(fEstimCatAv)
    cmd = {srcFs};
    cmd{end+1} = ['fslview -m single ' fEstimCatAv ' &'];
    cmd = strjoin(cmd,newline); % disp(cmd)
    if param.verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    fFslviewBR = cmd;
else
    fFslviewBR = '';
end
%%% within-run
if nframes>1
    cmd = {srcFs};
    cmd{end+1} = ['fslview -m single ' strjoin(fEstim,' ') ' &'];
    cmd = strjoin(cmd,newline); % disp(cmd)
    if param.verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    fFslviewWR = cmd;
    %%% within-run, first, middle and last frames
    fFslviewWRfstMdLst = qaFstMdLst(fEstim,forceRewriteAtPlumb || forceRecomputeRMS,param.verbose);
else
    fFslviewWR = '';
    fFslviewWRfstMdLst = '';
end

xSet.initFiles.qaFiles.fFslviewBR = fFslviewBR;
xSet.initFiles.qaFiles.fFslviewWR = fFslviewWR;
xSet.initFiles.qaFiles.fFslviewWRfstMdLst = fFslviewWRfstMdLst;



%% Low SNR data (estimate motion on temporally smoothed data)


%%% Detect low SNR data it not specified in dataType
if ~any(ismember({'lowSNR' 'highSNR'},dataType))
    dbstack; error(['automatic detection of lowSNR data not implemented' newline 'specify ''lowSNR'' or ''highSNR'' in dataType '])
else
    disp([dataType{ismember(dataType,{'lowSNR' 'highSNR'})} ' data specified'])
end


%%% Temporally smooth lowSNR data
if any(ismember({'lowSNR'},dataType))

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
    for i = 1:length(fEstim)
        fIn = fEstim{i};
        [d,fOut,~] = fileparts(replace(fIn,'.nii.gz',''));
        fOut = fullfile(d,[strjoin({['sm' num2str(param.tSmWin_vol)] fOut},'_') '.nii.gz']);
        cmd{end+1} = ['echo '' ''file' num2str(i) '/' num2str(length(fEstim))];
        if forceResmoothing || ~exist(fOut,'file')
            cmd{end+1} = '3dTsmooth -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = ['-hamming ' num2str(param.tSmWin_vol) ' \'];
            cmd{end+1} = fIn;
        else
            cmd{end+1} = ['echo ''  ''already smoothed, skipping'];
        end
        fEstim{i} = fOut;
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
    for i = 1:3
        cLim = [min([ax{i}.CLim(1) ax{i+3}.CLim(1)]) max([ax{i}.CLim(2) ax{i+3}.CLim(2)])];
        ax{i}.CLim = cLim;
        ax{i+3}.CLim = cLim;
    end
    [~,b,~] = fileparts(fileparts(fEstimAfterSm));
    title(hTile,b,'interpreter','none')
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
    for i = 1:length(fEstim)
        fIn = fEstim{i};
        fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
        if forceResmoothing || ~exist(fOut,'file')
            cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
        end
        fEstimAv{i} = fOut;
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
    cmd{end+1} = ['fslview -m single ' fEstimCatAv ' &'];
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

    xSet.initFiles.qaFilesSm.fFslviewBR = fFslviewBRsm;
    xSet.initFiles.qaFilesSm.fFslviewWR = fFslviewWRsm;
    xSet.initFiles.qaFilesSm.fFslviewWRfstMdLst = fFslviewWRsmFstMdLst;


end


%% Sumarize outputs
xSet.nEcho = nEcho;
xSet.initFiles.fOrig = fOrig;

xSet.initFiles.fPlumbRef = fPlumbRef;
xSet.initFiles.fPlumb = fPlumb;
xSet.initFiles.fPlumbAv = fPlumbAv;
xSet.initFiles.fPlumbCatAv = fPlumbCatAv;
xSet.initFiles.fPlumbAvCatAv = fPlumbAvCatAv;

xSet.initFiles.fEstim = fEstim;
xSet.initFiles.fEstimAv = fEstimAv;
if iscell(fEstimCatAv)
    xSet.initFiles.fEstimCatAv = fEstimCatAv;
else
    xSet.initFiles.fEstimCatAv = {fEstimCatAv};
end
if iscell(fEstimAvCatAv)
    xSet.initFiles.fEstimAvCatAv = fEstimAvCatAv;
else
    xSet.initFiles.fEstimAvCatAv = {fEstimAvCatAv};
end

xSet.initFiles.fApply = fApply;

xSet.initFiles.fObliqueRef = fObliqueRef;
xSet.initFiles.fPlumbRef = fPlumbRef;

%%% QA
if exist('fFigSm','var')
    xSet.initFiles.qaFilesSm.fFigSm = fFigSm;
end
if exist('fFslviewSm','var')
    xSet.initFiles.qaFilesSm.fFslviewSm = fFslviewSm;
end
