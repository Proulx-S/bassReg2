function [xSet,pSetInd] = computeSatin(xSet,type)

%% Get set labels
labelList = cell(size(xSet))';
for S = 1:length(xSet)
    labelList{S} = xSet{S}.label;
end

%% Identify appropriate sets
switch type
    case 'abs'
        pSetInd = contains(labelList,'f0b0');
        mSetInd = contains(labelList,'f1b1');
    case 'dir'
        pSetInd = contains(labelList,'f1b0');
        mSetInd = contains(labelList,'f0b1');
end
if all([any(pSetInd) any(mSetInd)])
    pSet = xSet{pSetInd};
    mSet = xSet{mSetInd};
else
    return
end

if pSet.nEcho>1    
    %% %%%%%%%%%%%%%%%
    % For every echo %
    %%%%%%%%%%%%%%% %%
    
    %%% Extract approproate files
    fPList = pSet.finalFiles;
    candidate = {'fAvCatAvEchoCat' 'fAvEchoCat' 'fEchoCat'};
    candidate = candidate(ismember(candidate,fields(fPList))); candidate = candidate{1};
    fPList = fPList.(candidate);

    fMList = mSet.finalFiles;
    candidate = {'fAvCatAvEchoCat' 'fAvEchoCat' 'fEchoCat'};
    candidate = candidate(ismember(candidate,fields(fMList))); candidate = candidate{1};
    fMList = fMList.(candidate);
    label = type; label(1) = upper(label(1)); label(2:end) = lower(label(2:end));
    label = ['satIndex' label];

    %%% Compute satin
    xSet{pSetInd}.(label).fEchoCat = doIt(fPList,fMList,type);

    %%% Add some outputs
    xSet{pSetInd}.(label).manBrainMask    = xSet{pSetInd}.finalFiles.manBrainMask;
    xSet{pSetInd}.(label).manBrainMaskInv = xSet{pSetInd}.finalFiles.manBrainMaskInv;


    %% %%%%%%%%%%%%%%%%%%%
    % For cross-echo RMS %
    %%%%%%%%%%%%%%%%%%% %%
    
    %%% Extract approproate files
    fPList = pSet.finalFiles;
    candidate = {'fAvCatAvEchoRms' 'fAvEchoRms' 'fEchoRms'};
    candidate = candidate(ismember(candidate,fields(fPList))); candidate = candidate{1};
    fPList = fPList.(candidate);

    fMList = mSet.finalFiles;
    candidate = {'fAvCatAvEchoRms' 'fAvEchoRms' 'fEchoRms'};
    candidate = candidate(ismember(candidate,fields(fMList))); candidate = candidate{1};
    fMList = fMList.(candidate);
    label = type; label(1) = upper(label(1)); label(2:end) = lower(label(2:end));
    label = ['satIndex' label];

    %%% Compute satin
    xSet{pSetInd}.(label).fEchoRms = doIt(fPList,fMList,type);

    %%% Add some outputs
    xSet{pSetInd}.(label).manBrainMask    = xSet{pSetInd}.finalFiles.manBrainMask;
    xSet{pSetInd}.(label).manBrainMaskInv = xSet{pSetInd}.finalFiles.manBrainMaskInv;


    %% %%%%%%%
    % For S0 %
    %%%%%%% %%
    %%% Extract approproate files
    fPList = pSet.finalFiles;
    candidate = {'fAvCatAvS0' 'fAvS0' 'fS0'};
    candidate = candidate(ismember(candidate,fields(fPList))); candidate = candidate{1};
    fPList = fPList.(candidate);

    fMList = mSet.finalFiles;
    candidate = {'fAvCatAvS0' 'fAvS0' 'fS0'};
    candidate = candidate(ismember(candidate,fields(fMList))); candidate = candidate{1};
    fMList = fMList.(candidate);
    
    label = type; label(1) = upper(label(1)); label(2:end) = lower(label(2:end));
    label = ['satIndex' label];

    %%% Compute satin
    xSet{pSetInd}.(label).fS0 = doIt(fPList,fMList,type);

    %%% Add some outputs
    xSet{pSetInd}.(label).manBrainMask    = xSet{pSetInd}.finalFiles.manBrainMask;
    xSet{pSetInd}.(label).manBrainMaskInv = xSet{pSetInd}.finalFiles.manBrainMaskInv;




else
    %% %%%%%%%%%%%%%%%%%%%%%
    % For single-echo data %
    %%%%%%%%%%%%%%%%%%%%% %%
    dbstack; error('compute that')
end


function fOut = doIt(fPList,fMList,type)

if ~all(size(fPList)==size(fMList))
    dbstack; error('something wrong');
end


pList = repmat(MRIread(fPList{1}),size(fPList));
mList = repmat(MRIread(fMList{1}),size(fMList));
fOut  = cell(size(fPList));


for i = 1:numel(fPList)
    %% Read data
    fP = fPList{i};
    fM = fMList{i};
    pMri = MRIread(fP);
    p = pMri.vol;
    mMri = MRIread(fM);
    m = mMri.vol;


    switch type

        case 'abs'
            %% Compute satIndex
            satIndex = (p - m) ./ p;
            satIndex(isnan(satIndex)|satIndex==inf|satIndex==-inf) = 0;

        case 'dir'        
            %% Compute satIndex
            satIndex = (p - m) ./ max(cat(10,p,m),[],10);
            satIndex(isnan(satIndex)|satIndex==inf|satIndex==-inf) = 0;
    end
    
    %% Set outputs
    label = type; label(1) = upper(label(1)); label(2:end) = lower(label(2:end));
    fOut{i} = replace(replace(replace(fP,'_volTs.nii.gz',['_satIndex' label '.nii.gz']),'_vol.nii.gz',['_satIndex' label '.nii.gz']),'_S0.nii.gz',['_S0satIndex' label '.nii.gz']);
    if strcmp(fP,fOut{i}); dbstack; error('could not determine output filename from input filename'); end

    %% Write to file
    satIndexMri = pMri;
    satIndexMri.vol = satIndex;
    MRIwrite(satIndexMri,fOut{i});
end


