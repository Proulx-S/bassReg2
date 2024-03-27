function [funcSet,param] = setPreproc(funcSet,param,funDsgn,force)
if ~exist('force','var'); force = []; end
if isempty(force); force = 0; end

wDirFunc = param.wDirFunc;
subInd = param.subInd;

%% Initiate data
sm = dir(fullfile(wDirFunc,funcSet.label,'*','*sm*.nii.gz')); sm = {sm.name}'; for i = 1:length(sm); sm{i} = strsplit(sm{i},'_'); sm(i) = sm{i}(contains(sm{i},'sm')); sm{i} = str2double(sm{i}(3:end)); end;
disp(unique([sm{:}])')

switch subInd
    case {1 2}
        funDsgn.trMri = 3;
    case 3
        funDsgn.trMri = 2.5;
end
funDsgn.trMri;
funDsgn.durSeq;
% sm = (18+39)/funDsgn.trMri;
sm = (18)/funDsgn.trMri;
sm = round((sm-1)/2)*2+1; sm(sm<0) = 1; % ensure sm is an odd number
param.tSmWin_vol = sm; % smoohting window size in number of volumes, must be odd
% param.tSmWin_vol = ceil((18+39)/2.5); % smoohting window size in number of volumes, must be odd
param.nDummy = 1;
param.verbose = 0;
dataType = 'lowSNR';
force;
funcSet = initPreproc(funcSet,wDirFunc,length(funcSet.files),dataType,param,force);

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
%%% Generate base for motion estimation
param.baseType = 'first';%'first' (default), 'mid', 'last', 'av', 'mcAv';
force;
funcSet.wrMocoFiles = genBaseWR(funcSet.initFiles,param,force);
%%% Estimate motion
force;;
param.spSmFac = 2; % spatial smoothing for estimation of motion, in number of voxels
funcSet.wrMocoFiles = estimMotionWR(funcSet.wrMocoFiles,param,force); %!!!! add base to qa cmd!!!!

%%% Visualize
% system(changeCLim(funcSet.initFiles.qaFilesSm.fFslviewWR,[0 600]));
% system(changeCLim(funcSet.wrMocoFiles.qaFiles.fFslviewWR,[0 600]));
% system(changeCLim(funcSet.initFiles.qaFilesSm.fFslviewWRfstMdLst,[0 600]));
% system(changeCLim(funcSet.wrMocoFiles.qaFiles.fFslviewWRfstMdLst,[0 600]));

%%% Between-run motion correction
if size(funcSet.initFiles.fPlumbList,1)>1
    force;
    source = funcSet.wrMocoFiles.fMocoAvList;
    base = funcSet.wrMocoFiles.fMocoAvList{1};
    mask = funcSet.wrMocoFiles.manBrainMaskInv;
    param.spSmFac = 0;
    funcSet.brMocoFiles = estimMotionBR(source,base,mask,param,force); %!!!! add base to qa cmd!!!!
end


%%% Apply motion correction
force;
param.maskFile = funcSet.initFiles.manBrainMaskInv;
imFiles = funcSet.initFiles.fPlumbList;
motFiles = funcSet.wrMocoFiles.fMocoMatList;
if isfield(funcSet,'brMocoFiles')
    motFiles = cat(2,motFiles,funcSet.brMocoFiles.fMocoMatList);
end
funcSet.preprocFiles = applyMotion(imFiles,motFiles,param,force);
% funcSet.preprocFiles = addMaskToCmd(funcSet.preprocFiles,maskFile);

% system(changeCLim(funcSet.initFiles.qaFilesSm.fFslviewWRfstMdLst,[0 600]));
% system(changeCLim(funcSet.preprocFiles.qaFilesSm.fFslviewWRfstMdLst,[0 600]));

% system(changeCLim(funcSet.preprocFiles.qaFilesSm.fFslviewWR,[0 600]));

% system(changeCLim(funcSet.initFiles.qaFiles.fFslviewBR,[0 600]));
% system(changeCLim(funcSet.preprocFiles.qaFiles.fFslviewBR,[0 600]));

%%% Put back to oblique (same for all run in set)
force;
funcSet.finalFiles = finalizePreproc(funcSet.preprocFiles,funcSet.initFiles,force);




