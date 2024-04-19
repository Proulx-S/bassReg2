function [funcSet,param] = setPreproc(funcSet,param,funDsgn,force)
if ~exist('force','var');                         force = []        ; end
if ~isfield(funcSet,'inplaneFlag'); funcSet.inplaneFlag = []        ; end
if ~isfield(funcSet,'dataType');       funcSet.dataType = {}        ; end
if isempty(force);                                force = 0         ; end
if isempty(funcSet.inplaneFlag);    funcSet.inplaneFlag = 0         ; end
if isempty(funcSet.dataType);          funcSet.dataType = {'lowSNR'}; end

wDirFunc = param.wDirFunc;
subInd = param.subInd;
setInd = param.setInd;

%% Initiate data
sm = dir(fullfile(wDirFunc,funcSet.label,'*','*sm*.nii.gz')); sm = {sm.name}'; for i = 1:length(sm); sm{i} = strsplit(sm{i},'_'); sm(i) = sm{i}(contains(sm{i},'sm')); sm{i} = str2double(sm{i}(3:end)); end
% disp(unique([sm{:}])')

% switch subInd
%     case {1 2}
%         funDsgn.trMri = 3;
%     case 3
%         funDsgn.trMri = 2.5;
% end
% funDsgn.trMri;

mriInfo = MRIread(fullfile(funcSet.files(1).folder,funcSet.files(1).name),1);
funcSet.nFrame = mriInfo.nframes;
if funcSet.nFrame>1
    if ~any(ismember(funcSet.dataType,'volTs'))
        funcSet.dataType{end+1} = 'volTs';
    end
    funDsgn.trMri = mriInfo.tr/1000;
    funDsgn.durSeq;
    if any(ismember(funcSet.dataType,'lowSNR'))
        % sm = (18+39)/funDsgn.trMri;
        sm = (18)/funDsgn.trMri;
        sm = ceil((sm-1)/2)*2+1; sm(sm<0) = 1; % ensure sm is an odd number
    else
        sm = 0;
        if ~any(ismember(funcSet.dataType,'highSNR'))
            funcSet.dataType{end+1} = 'highSNR';
        end
    end
    param.tSmWin_vol = sm; % smoohting window size in number of volumes, must be odd
    % param.tSmWin_vol = ceil((18+39)/2.5); % smoohting window size in number of volumes, must be odd
    if isfield(funcSet,'nDummy') && ~isempty(funcSet.nDummy)
        param.nDummy = funcSet.nDummy;
    else
        mri = MRIread(fullfile(funcSet.files(1).folder,funcSet.files(1).name));
        ts = permute(mri.vol,[4 1 2 3]); ts = mean(ts(:,:),2);
        t = mri.tr/1000; t = 0:t:t*(mri.nframes-1);
        figure('WindowStyle','docked');
        plot(t,ts,'-o'); xlabel('sec'); ylabel('means signal'); grid on; title(funcSet.files(1).name,'Interpreter','none')
        dbstack; error(['please specify funcSet.nDummy for ' funcSet.label])
    end
    param.verbose = 0;
else
    if ~any(ismember(funcSet.dataType,'vol'))
        funcSet.dataType{end+1} = {'vol'};
    end
end

forceThis = force>1;%0;
funcSet = initPreproc(funcSet,wDirFunc,[],param,forceThis);


%% Draw brain mask
funcSet.initFiles = makeMask(funcSet.initFiles,funcSet.bidsDerivDir);

% %% Visualize data
% system(funcSet.initFiles.qaFiles.fFslviewBR);
% system(funcSet.initFiles.qaFiles.fFslviewWR);

% system(changeCLim(funcSet.initFiles.qaFiles.fFslviewWRfstMdLst,[100 800]));
% system(funcSet.initFiles.qaFilesSm.fFslviewBR);
% system(funcSet.initFiles.qaFilesSm.fFslviewWR);
% system(funcSet.initFiles.qaFilesSm.fFslviewWRfstMdLst);
% open(funcSet.initFiles.qaFilesSm.fFigSm)
% system(funcSet.initFiles.qaFilesSm.fFslviewSm);



%% Within-run motion correction
if funcSet.nFrame>1
    %%% Generate base for motion estimation
    param.baseType = 'first';%'first' (default), 'mid', 'last', 'av', 'mcAv';
    force;
    funcSet.wrMocoFiles = genBaseWR(funcSet.initFiles,param,force);
    %%% Estimate motion
    force;
    if any(ismember(funcSet.dataType,'highSNR'))
        param.spSmFac = 0;
    elseif any(ismember(funcSet.dataType,'lowSNR'))
        param.spSmFac = 2; % spatial smoothing for estimation of motion, in number of voxels
    end
    if any(ismember(funcSet.dataType,'inplane'))
        param.explicitlyFixThroughPlane = 1;
    else
        param.explicitlyFixThroughPlane = 0;
    end

    funcSet.wrMocoFiles = estimMotionWR(funcSet.wrMocoFiles,param,force);

    %%% Visualize
    % system(changeCLim(funcSet.initFiles.qaFiles.fFslviewWR,[0 600]));
    % system(changeCLim(funcSet.initFiles.qaFilesSm.fFslviewWR,[0 600]));
    % system(changeCLim(funcSet.wrMocoFiles.qaFiles.fFslviewWR,[0 600]));
    % system(changeCLim(funcSet.initFiles.qaFilesSm.fFslviewWRfstMdLst,[0 600]));
    % system(changeCLim(funcSet.wrMocoFiles.qaFiles.fFslviewWRfstMdLst,[0 600]));
end

%%% Between-run motion correction
if funcSet.nRun>1
    force;
    source = funcSet.wrMocoFiles.fMocoAv;
    base = funcSet.wrMocoFiles.fMocoAv{1}; %!!! should allow to input the index to use as the base !!!%
    mask = funcSet.wrMocoFiles.manBrainMaskInv;
    param.spSmFac = 0;
    funcSet.brMocoFiles = estimMotionBR(source,base,mask,param,force);
end


%%% Apply motion correction
force;
param.maskFile = funcSet.initFiles.manBrainMaskInv;
imFiles = funcSet.initFiles;
if funcSet.nFrame>1 || funcSet.nRun>1
    motFiles = funcSet.wrMocoFiles;
    if isfield(funcSet,'brMocoFiles')
        % catenate successive transformations in the 3rd dimension
        motFiles.fMocoMat = cat(10,motFiles.fMocoMat,funcSet.brMocoFiles.fMocoMat);
    end
else
    motFiles = [];
end
funcSet.preprocFiles = applyMotion(imFiles,motFiles,param,force);
% else
%     disp('applying motion correction: nothing to apply (probably nFrame=nRun=1)')
%     for i = 1:numel(funcSet.initFiles.fOrig)
%     end
% 
% 
%     funcSet.preprocFiles.fOrig = funcSet.initFiles.fOrig;
%     funcSet.preprocFiles.nRun = funcSet.initFiles.nRun;
%     funcSet.preprocFiles.nFrame = funcSet.initFiles.nFrame;
%     funcSet.preprocFiles.nEcho = funcSet.initFiles.nEcho;
% end
% funcSet.preprocFiles = addMaskToCmd(funcSet.preprocFiles,maskFile);

% system(changeCLim(funcSet.initFiles.qaFilesSm.fFslviewWRfstMdLst,[0 600]));
% system(changeCLim(funcSet.preprocFiles.qaFilesSm.fFslviewWRfstMdLst,[0 600]));

% system(changeCLim(funcSet.preprocFiles.qaFilesSm.fFslviewWR,[0 600]));

% system(changeCLim(funcSet.initFiles.qaFiles.fFslviewBR,[0 600]));
% system(changeCLim(funcSet.preprocFiles.qaFiles.fFslviewBR,[0 600]));





