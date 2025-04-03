function [outFile,hFig] = xCorrQA(fList,mList,nDummy,ttStr,outDir,force,verbose)

if ~exist('force','var');        force = []                ; end
if ~exist('verbose','var');    verbose = []                ; end
if ~exist('nDummy','var');      nDummy = []                ; end
if ~exist('ttStr','var');        ttStr = []                ; end
if ~exist('outDir','var')       outDir = []                ; end
if isempty(nDummy);             nDummy = zeros(size(fList)); end
if isempty(force);               force = 0                 ; end
if isempty(verbose);           verbose = 0                 ; end
if isempty(ttStr);               ttStr = 'QA'              ; end
if isempty(outDir);             outDir = ''                ; end
%% 
outFile = fullfile(outDir,[ttStr '.fig']);
if ~force && exist(outFile,'file')
    if verbose
        hFig = open(outFile);
    else
        hFig = [];
    end
    return
end


%% Read data and masks
mriHd    = MRIread(fList{1},1);
im       = false([mriHd.volsize 0              ]);
frameI   = false([1 1 1         0              ]);
fileI    = false([1 1 1         0              ]);
fileTR   = false([1 1 1         0              ]);
nDummy2  = false([1 1 1         0              ]);
imMask   = false([mriHd.volsize 1 size(mList,1)]);
nFrames  = zeros(size(fList,1),1);
bidsName = cell(size(fList,1),1);
for f = 1:size(fList,1)
    disp(['reading (' num2str(f) '/' num2str(size(fList,1)) ') ' fList{f,1} ])
    
    mri = MRIread(fList{f,1});
    im    = cat(4,im     ,single(mri.vol      ));

    fileI   = cat(4,fileI  ,repmat(          f   ,[1 1 1 mri.nframes]));
    fileTR  = cat(4,fileTR ,repmat(mri.tr/1000   ,[1 1 1 mri.nframes]));
    nDummy2 = cat(4,nDummy2,repmat(nDummy(f)     ,[1 1 1 mri.nframes]));
    frameI = cat(4,frameI,permute(1:mri.nframes,[1 3 4 2]));
    nFrames(f) = mri.nframes;
    
    mri = MRIread(mList{f,1});
    imMask(:,:,:,:,f) = logical(mri.vol);

    bidsName{f} = strsplit(fList{f,1},filesep);
    while ~contains(bidsName{f}{end},'sub-') & ~contains(bidsName{f}{end},'ses-')
        bidsName{f}(end) = [];
    end
    bidsName{f} = bidsName{f}{end};
end


%% Uniformize crop masks
[~,b] = fileparts(replace(mList{f,1},'.nii.gz',''));
if contains(b,'Inv')
    imMask = all(~imMask,5);
else
    imMask = all( imMask,5);
end
if size(imMask,3)>1; imMask(:,:,[1 end],:,:) = false; end


%% Correlate frames
im = permute(im,[4 1 2 3]);
rho = corr(permute(im(:,imMask),[2 1]));


%% Plot
if verbose
    hFig = figure('WindowStyle','docked');
else
    hFig = figure('WindowStyle','docked','Visible','off');
end
imagesc(rho,[0 1]);
ax = gca;
ax.DataAspectRatio = [1 1 1];
ylabel(colorbar,'cross-frame Pearson''s correlation')

dTick = round(cumsum(nFrames) - nFrames/2);
ax.XTick = dTick;
ax.YTick = dTick;
ax.TickLabelInterpreter = 'none';
ax.XTickLabels = [];
ax.YTickLabels = bidsName;

if length(nFrames)>1
    sLine = cumsum(nFrames)'; sLine = [0 sLine(1:end-1)]+0.5;
    xline(sLine,'r','LineWidth',eps)
    yline(sLine,'r','LineWidth',eps)
end
sub = strsplit(outDir,filesep);
sub = sub{contains(sub,'sub-')};
[~,b] = fileparts(replace(outDir,'_set',''));
title([sub '; set ' b newline ttStr],'Interpreter','none')
drawnow


%% Store metadata for later use
hFig.UserData.fileNames   = fList(:,1);
hFig.UserData.nDummy      = nDummy;
hFig.UserData.fileInd     = fileI;
hFig.UserData.frameNumber = frameI;
hFig.UserData.frameTime   = (frameI+nDummy2-1).*fileTR;

%% Save figure
if ~isempty(outDir)
    if ~verbose
        set(hFig, 'CreateFcn', @(obj, ~) set(obj, 'Visible', 'on'));
    end
    if ~exist(fileparts(outFile),'dir'); mkdir(fileparts(outFile)); end
    saveas(hFig,outFile);
end

