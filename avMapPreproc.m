function avMapSet = avMapPreproc(avMapSet,regBaseSet,force,verbose)
global srcFs srcAfni bidsDir

if ~isfield(avMapSet,'S');  S       = []; else; S = avMapSet.S; end
if ~exist('force','var');   force   = []; end
if ~exist('verbose','var'); verbose = []; end
% if ~isfield(avMapSet,'fRegBasePlumb'); fRegBasePlumb = []                    ; end
% if ~isfield(avMapSet,'fMaskPlumb');    fMaskPlumb    = []                    ; end

if isempty(S);              S       = 1 ; end
if isempty(force);          force   = 0 ; end
if isempty(verbose);        verbose = 1 ; end
% if isfield(avMapSet,'fRegBasePlumb');  fRegBasePlumb = avMapSet.fRegBasePlumb; end %must already be at plumb
% if isfield(avMapSet,'fMaskPlumb');     fMaskPlumb    = avMapSet.fMaskPlumb   ; end %must already be at plumb

wDirAnat    = avMapSet.wDirAnat;
bidsPrpList = avMapSet.bidsPrpList;

%% Initiate preproc
for i = 1:length(bidsPrpList)
    av = dir(fullfile(bidsDir,'anat',bidsPrpList{i}));
    % av = fullfile({av(:).folder},{av(:).name})';
    if ~isempty(av); break; end
end
avMapSet.files = av;
avMapSet.label
avMapSet.dataType = {'vol' 'highSNR'};

force;
avMapSet = initPreproc(avMapSet,wDirAnat,[],[],[],force);


%% Register
%%% estimate registration from avMap to functional
force;
param.verbose = verbose;
param.cost = 'nmi';
avMapSet.bsMocoFiles = estimMotionBS(avMapSet,regBaseSet{S},[],avMapSet.regTo_funcSet.runInd,param,force);

% avMapSet.fRegBasePlumb = funcSet{S}.preprocFiles.fCorrectedAvCatAvEchoRms; % pick deoblic file
% avMapSet.fMaskPlumb = funcSet{S}.initFiles.manBrainMaskInv; % pick deoblic file


%%%% add other sets to visualization command
indList = 1:length(regBaseSet); indList(indList==S) = [];
cmd = avMapSet.bsMocoFiles.qaFiles.fFslviewBS; cmd = strsplit(cmd,newline);
cmd1 = cmd(1:find(contains(cmd,'freeview')));
cmdX = {};
cmd2 = cmd(find(contains(cmd,'freeview'))+1:end);
for ind = indList
    cmdX{end+1} = [autoSelect(regBaseSet{ind}) ':name=' regBaseSet{ind}.label ' \'];
end
avMapSet.bsMocoFiles.qaFiles.fFslviewBS = strjoin([cmd1 cmdX cmd2],newline);

% clipboard('copy',avMapSet.bsMocoFiles.qaFiles.fFslviewBS)


%%% apply registration
imFiles = avMapSet.initFiles;
motFiles = avMapSet.bsMocoFiles;
avMapSet.preprocFiles = applyMotion(imFiles,motFiles,[],force);

% clipboard('copy',avMapSet.preprocFiles.qaFiles.fFslviewBS)


%% Finalize (put back to oblique)
avMapSet.finalFiles = finalizePreproc(avMapSet.preprocFiles,avMapSet.initFiles,force);


%% S0 and T2*
TE = fullfile({avMapSet.files.folder},{avMapSet.files.name})';
for i = 1:length(TE)
    [~,TE{i}] = system(['jq ''.EchoTime'' ' replace(TE{i},'.nii.gz','.json')]); TE{i} = str2num(replace(replace(TE{i},'[0;39m',''),['[0m' newline],''));
end
TE = cell2mat(TE);

nneg = 0;
fInList = avMapSet.finalFiles.fAvEchoCat;
avMapSet.finalFiles.fAvEchoS0 = cell(size(fInList,1),1,1);
avMapSet.finalFiles.fAvEchoT2s = cell(size(fInList,1),1,1);
for i = 1:size(fInList,1)
    fIn = fInList{i};
    avS0 = replace(fIn,'_vol.nii.gz','_S0.nii.gz');
    avT2s = replace(fIn,'_vol.nii.gz','_T2s.nii.gz');
    if force || ~exist(avS0,'file') || ~exist(avT2s,'file')
        av = MRIread(fIn);
        [t2star,S0] = calc_t2s_vol(av.vol(:,:,:,1:end), TE*1000, nneg, verbose); % last echo is actually the rms average
        av.vol = S0; if force || ~exist(avS0,'file'); MRIwrite(av,avS0); end
        av.vol = t2star; if force || ~exist(avT2s,'file'); MRIwrite(av,avT2s); end
    end
    avMapSet.finalFiles.fAvEchoS0{i} = avS0;
    avMapSet.finalFiles.fAvEchoT2s{i} = avT2s;
end

    


