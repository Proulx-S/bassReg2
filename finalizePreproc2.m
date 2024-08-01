function finalPreprocFiles = finalizePreproc2(initFiles,preprocFiles,force,verbose)
global srcAfni
if ~exist('force','var');     force = []; end
if ~exist('verbose','var'); verbose = []; end
if isempty(force);            force = 0; end
if isempty(verbose);        verbose = 0; end

% if preprocFiles.updatedFlag;  force = 1; end % override force=0 when preproc was updated
bidsDerivDir = fullfile(fileparts(fileparts(initFiles.fOrigList{1})),'derivatives');

%% Combine transformation
finalPreprocFiles           = rmfield(initFiles,'fEstimList');
finalPreprocFiles.fTransList    = cell(size(finalPreprocFiles.fOrigList,1),length(preprocFiles));
finalPreprocFiles.fTransCatList = cell(size(finalPreprocFiles.fOrigList,1),1);
finalPreprocFiles.ppLabelList   = cell(1,length(preprocFiles));
for i = 1:length(preprocFiles)
    finalPreprocFiles.ppLabelList{1,i} = preprocFiles{i}.ppLabel;
    switch preprocFiles{i}.ppLabel
        case {'withinRunMoco' 'betweenRunMoco'}
            finalPreprocFiles.fTransList(:,i) = replace(preprocFiles{i}.fMocoList,'.nii.gz','.aff12.1D');
        otherwise
            dbstack; error('code that')
    end
end
cmd = {srcAfni};
for r = 1:size(finalPreprocFiles.fTransList,1)
    [~,b,~] = fileparts(fileparts(finalPreprocFiles.fTransList{r,1}));
    finalPreprocFiles.fTransCatList{r} = fullfile(bidsDerivDir,b);
    if ~exist(finalPreprocFiles.fTransCatList{r},'dir'); mkdir(finalPreprocFiles.fTransCatList{r}); end
    finalPreprocFiles.fTransCatList{r} = fullfile(finalPreprocFiles.fTransCatList{r},'transCat.aff12.1D');

    cmd{end+1} = ['rm -f ' finalPreprocFiles.fTransCatList{r}];
    cmd{end+1} = ['head -1 ' strjoin(finalPreprocFiles.fTransList(r,1),' ') ' > ' finalPreprocFiles.fTransCatList{r}];
    cmd{end+1} = ['cat_matvec ' strjoin(flip(finalPreprocFiles.fTransList(r,:)),' ') ' >> ' finalPreprocFiles.fTransCatList{r}];
end



%% Apply transformations
finalPreprocFiles.fPreprocList = cell(size(finalPreprocFiles.fTransCatList));
fPreprocUpdated                = false(size(finalPreprocFiles.fTransCatList));
for r = 1:size(finalPreprocFiles.fTransList,1)
    fIn    = finalPreprocFiles.fPlumbList{r};
    fTrans = finalPreprocFiles.fTransCatList{r};
    fOut   = fullfile(fileparts(finalPreprocFiles.fTransCatList{r}),'preproc_volTs.nii.gz');

    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dAllineate -overwrite -nocmass \';
        cmd{end+1} = ['-source ' fIn ' \'];
        cmd{end+1} = ['-1Dmatrix_apply ' fTrans ' \'];
        cmd{end+1} = ['-prefix ' fOut];
        fPreprocUpdated(r) = true;
    end

    finalPreprocFiles.fPreprocList{r} = fOut;
end


%% Run command
if verbose
    [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
else
    [status,cmdout] = system(strjoin(cmd,newline)); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
end


%% Rewrite at setOblique
disp(' rewrting preproc data at setOblique')
mriOblique = MRIread(initFiles.fGeom);
for r = 1:size(finalPreprocFiles.fPreprocList,1)
    disp(['file ' num2str(r) '/' num2str(size(finalPreprocFiles.fPreprocList,1))])
    if fPreprocUpdated(r)
        mri = MRIread(finalPreprocFiles.fPreprocList{r,1});
        mriOblique.vol = mri.vol;
        MRIwrite(mriOblique,finalPreprocFiles.fPreprocList{r,1});
        disp(' done')
    else
        disp(' already done, skipping')
    end
end


%% Summarize
forceThis   = force;
verboseThis = verbose;
summarizeVolTs(finalPreprocFiles.fPreprocList,[],initFiles.nFrame,[],initFiles.dataType,forceThis,verboseThis)
% fOrigList = replace(finalPreprocFiles.fPreprocList,'preproc_volTs.nii.gz','orig_volTs.nii.gz')
summarizeVolTs(initFiles.fOrigList,replace(finalPreprocFiles.fPreprocList,'preproc_volTs.nii.gz','orig_volTs.nii.gz'),initFiles.nFrame,initFiles.nDummy,initFiles.dataType,forceThis,verboseThis)



return

if strcmp(files.fTrans{R,1},'''IDENTITY''')
    fMotCat = '''IDENTITY''';
else
    fMotCat = replace(fOut,'.nii.gz','.aff12.1D');
    if size(files.fTrans,10)>1
        cmd{end+1} = ['rm -f ' fMotCat];
        cmd{end+1} = ['head -1 ' strjoin(files.fTrans(R,1,1,1,1,1,1,1,1,1),' ') ' > ' fMotCat];
        cmd{end+1} = ['cat_matvec ' strjoin(flip(files.fTrans(R,1,1,1,1,1,1,1,1,:)),' ') ' >> ' fMotCat];
    else
        fMotCat = files.fTrans{R,1,1,1,1,1,1,1,1,1};
    end
end

cmd{end+1} = '3dAllineate -overwrite \';
cmd{end+1} = ['-source ' fIn ' \'];
cmd{end+1} = ['-1Dmatrix_apply ' fMotCat ' \'];
cmd{end+1} = ['-prefix ' fOut];






%% Define all filenames
fieldList1 = {'fCorrected' 'fCorrectedAv'};
fieldList2 = {'f'          'fAv'};

fieldList2 = fieldList2(ismember(fieldList1,fields(preprocFiles)));
fieldList1 = fieldList1(ismember(fieldList1,fields(preprocFiles)));
for i = 1:length(fieldList1)
    finalFiles.(fieldList2{i}) = cell(size(preprocFiles.fOrig));
    for I = 1:numel(preprocFiles.fOrig)
        fIn = preprocFiles.(fieldList1{i}){I};
        fOut = replace(fIn,'setPlumb_',''); fOut = strsplit(fOut,filesep); fOut = strjoin([fOut(1:end-1) {'preproc'} fOut(end)],filesep); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        finalFiles.(fieldList2{i}){I} = fOut;

        fInListA{I}{i} = fIn;
        fOutListA{I}{i} = fOut;
    end
end




fieldList1 = {'fCorrectedEchoRms' 'fCorrectedAvEchoRms' 'fCorrectedEchoCat' 'fCorrectedAvEchoCat' 'fCorrectedEchoRms' 'fCorrectedAvEchoRms'};
fieldList2 = {'fEchoRms'          'fAvEchoRms'          'fEchoCat'          'fAvEchoCat'          'fEchoRms'          'fAvEchoRms'};
clear fInListA2 fOutListA2
for i = 1:length(fieldList1)
    if ~isfield(preprocFiles,fieldList1{i}); continue; end
    finalFiles.(fieldList2{i}) = cell(size(preprocFiles.(fieldList1{i})));
    for I = 1:size(preprocFiles.(fieldList1{i}),1)
        fIn = preprocFiles.(fieldList1{i}){I};
        fOut = replace(fIn,'setPlumb_',''); fOut = strsplit(fOut,filesep); fOut = strjoin([fOut(1:end-1) {'preproc'} fOut(end)],filesep); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        finalFiles.(fieldList2{i}){I} = fOut;

        fInListA2{I}{i} = fIn;
        fOutListA2{I}{i} = fOut;
    end
end
if exist('fInListA2','var') && exist('fOutListA2','var')
    for I = 1:length(fInListA2)
        ind = cellfun('isempty',fInListA2{I});
        fInListA2{I}(ind) = [];
        fOutListA2{I}(ind) = [];
    end

    fInListA = cat(2,fInListA,fInListA2);
    fOutListA = cat(2,fOutListA,fOutListA2);
end


fieldList1 = {'fCorrectedCatAv' 'fCorrectedAvCatAv'};
fieldList2 = {'fCatAv'          'fAvCatAv'};
fInListB = {};
fOutListB = {};
for i = 1:length(fieldList1)
    if ~isfield(preprocFiles,fieldList1{i}); continue; end
    for E = 1:size(preprocFiles.(fieldList1{i}),2)
        fIn = preprocFiles.(fieldList1{i}){E};
        fOut = replace(fIn,'setPlumb_',''); fOut = strsplit(fOut,filesep); fOut = strjoin([fOut(1:end-1) {'preproc'} fOut(end)],filesep); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        finalFiles.(fieldList2{i}){1,E} = fOut;

        fInListB{end+1} = fIn;
        fOutListB{end+1} = fOut;
    end
end


fieldList1 = {'fCorrectedAvCatAvEchoCat' 'fCorrectedCatAvEchoRms' 'fCorrectedAvCatAvEchoRms'};
fieldList2 = {'fAvCatAvEchoCat'          'fCatAvEchoRms'          'fAvCatAvEchoRms'         };
% fieldList1 = {'fCorrectedCatAvEchoRms' 'fCorrectedAvCatAvEchoRms' 'fCorrectedCatAvEchoCat' 'fCorrectedAvCatAvEchoCat'};
% fieldList2 = {'fCatAvEchoRms' 'fAvCatAvEchoRms' 'fCatAvEchoCat' 'fAvCatAvEchoCat'};
fInListB2 = {};
fOutListB2 = {};
for i = 1:length(fieldList1)
    if ~isfield(preprocFiles,fieldList1{i}); continue; end
        fIn = char(preprocFiles.(fieldList1{i}));
        fOut = replace(fIn,'setPlumb_',''); fOut = strsplit(fOut,filesep); fOut = strjoin([fOut(1:end-1) {'preproc'} fOut(end)],filesep); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        finalFiles.(fieldList2{i}) = {fOut};

        fInListB2{end+1} = fIn;
        fOutListB2{end+1} = fOut;
end

fInListB = cat(2,fInListB,fInListB2);
fOutListB = cat(2,fOutListB,fOutListB2);



fInListB{end+1} = preprocFiles.manBrainMask;
fOutListB{end+1} = replace(fInListB{end},'.nii.gz','Oblique.nii.gz');
finalFiles.manBrainMask = fOutListB{end};

fInListB{end+1} = preprocFiles.manBrainMaskInv;
fOutListB{end+1} = replace(fInListB{end},'.nii.gz','Oblique.nii.gz');
finalFiles.manBrainMaskInv = fOutListB{end};

%% Rewrite all to oblique
obliqueRef = MRIread(initFiles.fObliqueRef,1);
disp('writing back to set oblique (for visualization with other images)')
for I = 1:length(fInListA)
    disp([' file ' num2str(I) '/' num2str(length(fInListA))])
    for i = 1:length(fInListA{I})
        if force || ~exist(fOutListA{I}{i},'file')
            mriTmp     = MRIread(fInListA{I}{i});
            mri        = obliqueRef;
            mri.vol    = mriTmp.vol;

            mri.xsize  = mriTmp.xsize;
            mri.ysize  = mriTmp.ysize;
            mri.zsize  = mriTmp.zsize;
            mri.volres = mriTmp.volres;
            
            mri.height  = mriTmp.height;
            mri.width  = mriTmp.width;
            mri.depth  = mriTmp.depth;
            mri.volsize = mriTmp.volsize;
            
            mri.nframes = mriTmp.nframes;
            mri.nvoxels = mriTmp.nvoxels;

            mri.tr     = mriTmp.tr;
            mri.te     = mriTmp.te;
            mri.ti     = mriTmp.ti;
            MRIwrite(mri,fOutListA{I}{i});
        end
        disp(['  in : ' fInListA{I}{i}])
        disp(['  out: ' fOutListA{I}{i}])
    end
    disp(' done')
end
disp(' runCat')
for i = 1:length(fInListB)
    if force || ~exist(fOutListB{i},'file')
        mriTmp = MRIread(fInListB{i});
        obliqueRef.vol = mriTmp.vol;
        MRIwrite(obliqueRef,fOutListB{i});
    end
    disp(['  in : ' fInListB{i}])
    disp(['  out: ' fOutListB{i}])
end
disp(' done')

