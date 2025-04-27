function [fFig,hFig] = xCorrQA(fList,mList,nDummy,ttStr,outDir,force,verbose)
global src
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



    %% Massage fList
    cmd = {src.afni};

    % run-by-run timeseries averages
    fAvList = cell(size(fList));
    for R = 1:size(fList,1)
        fAvList{R} = strsplit(fList{R,1,1},filesep); fAvList{R}{end} = ['av_' fAvList{R}{end}]; fAvList{R} = strjoin(fAvList{R},filesep);
        if force || ~exist(fAvList{R},'file')
            cmd{end+1} =  '3dTstat -overwrite -mean \';
            cmd{end+1} = ['-prefix ' fAvList{R} ' \'];
            cmd{end+1} =  fList{R};
        end
    end

    % cross-run catenation
    fName = strsplit(fAvList{1},filesep); fName = fName{end};
    venc = strsplit(fAvList{1},'_');
    switch nnz(contains(venc,'acq-pcVenc'))
        case 1
            venc = venc{contains(venc,'acq-pcVenc')};
            fName = replace(fName,'_volTs.nii.gz',['_' venc '_volTs.nii.gz']);
        case 0
            venc = '';
        otherwise
            error('Multiple venc strings found in fAvList');
    end
    fCatAv = unique(fileparts(fileparts(fAvList)));
    if length(fCatAv)>1
        fCatAv = strsplit(fCatAv{1},filesep);
        fCatAv{contains(fCatAv,'ses-')} = 'ses-cat';
        fCatAv = strjoin(fCatAv,filesep);
    else
        fCatAv = char(fCatAv);
    end
    if force || ~exist(fCatAv,'file')
        cmd{end+1} = '3dTcat -overwrite \';
        cmd{end+1} = ['-prefix ' fullfile(fCatAv,fName) ' \'];
        cmd{end+1} = strjoin(fAvList, ' ');
    end

    % cross-run arerage
    fCatAv = fullfile(fCatAv,['cat_' fName]);
    fAvCatAv = strsplit(fCatAv,filesep); fAvCatAv{end} = ['av_' fAvCatAv{end}]; fAvCatAv = strjoin(fAvCatAv,filesep);
    if force || ~exist(fAvCatAv,'file')
        cmd{end+1} = '3dTstat -overwrite -mean \';
        cmd{end+1} = ['-prefix ' fAvCatAv ' \'];
        cmd{end+1} = fCatAv;
    end
    
    if length(cmd)>1
        if verbose
            [status, result] = system(strjoin(cmd, newline), '-echo');
        else
            [status, result] = system(strjoin(cmd, newline));
        end
    end






    %%
    if isempty(outDir); outDir = fileparts(fAvCatAv); else; fFig = fullfile(outDir,[ttStr '.fig']); end
    if isempty(ttStr); fFig = fullfile(outDir,['xCorrQA.fig']); else; fFig = fullfile(outDir,[ttStr '.fig']); end
    if ~isempty(venc)
        fFig = strsplit(fFig,filesep); fFig{end} = [venc '_' fFig{end}]; fFig = strjoin(fFig,filesep);
    end
    
    
    if ~force && exist(fFig,'file')
        if verbose
            hFig = open(fFig);
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
        [~,b] = fileparts(replace(mList{f,1},'.nii.gz',''));
        if contains(b,'Inv')
            imMask(:,:,:,:,f) = ~logical(mri.vol);
        else
            imMask(:,:,:,:,f) =  logical(mri.vol);
        end

        bidsName{f} = strsplit(fList{f,1},filesep);
        while ~contains(bidsName{f}{end},'sub-') & ~contains(bidsName{f}{end},'ses-')
            bidsName{f}(end) = [];
        end
        bidsName{f} = bidsName{f}{end};
    end


    %% Uniformize crop masks
    imMask = all(imMask,5);
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
    axis image
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
    bids = strsplit(outDir,filesep);
    sub = bids{contains(bids,'sub-')};
    acq = bids{contains(bids,'acq-')};
    ttStr2 = [sub '; set: ' acq];
    if ~isempty(venc)
        ttStr2 = [ttStr2 '; ' venc];
    end
    title(ttStr2,'Interpreter','none')
    drawnow


    %% Store metadata for later use
    hFig.UserData.fileNames   = fList(:,1);
    hFig.UserData.fAvList     = fAvList;
    hFig.UserData.fCatAv      = fCatAv;
    hFig.UserData.fAvCatAv    = fAvCatAv;
    hFig.UserData.nDummy      = nDummy;
    hFig.UserData.fileInd     = fileI;
    hFig.UserData.frameNumber = frameI;
    hFig.UserData.frameTime   = (frameI+nDummy2-1).*fileTR;
    hFig.UserData.venc        = venc;

    %% Save figure
    if ~verbose
        set(hFig, 'CreateFcn', @(obj, ~) set(obj, 'Visible', 'on'));
    end
    if ~exist(fileparts(fFig),'dir'); mkdir(fileparts(fFig)); end
    saveas(hFig,fFig);


