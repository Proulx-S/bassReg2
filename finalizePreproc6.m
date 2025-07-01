function finalPreprocFiles = finalizePreproc6(initFiles,preprocFiles,force,verbose)
global src
if ~exist('force','var');     force = []; end
if ~exist('verbose','var'); verbose = []; end
if isempty(force);            force = 0; end
if isempty(verbose);        verbose = 0; end

[~,bidsDerivDir] = fileparts(fileparts(initFiles.fOrigList{1}));
bidsDerivDir = fullfile(initFiles.info.prcDir,'bids','derivatives',['sub-' initFiles.sub],['ses-' initFiles.ses],bidsDerivDir,initFiles.label);


%% Combine transformations
finalPreprocFiles           = rmfield(initFiles,'fEstimList');
finalPreprocFiles.fTransList    = cell(size(finalPreprocFiles.fPlumbList,1),1,size(preprocFiles,3));
finalPreprocFiles.fTransCatList = cell(size(finalPreprocFiles.fPlumbList,1),1,1                   );
finalPreprocFiles.ppLabelList   = cell(1                                   ,1,length(preprocFiles));
for i = 1:length(preprocFiles)
    if isempty(preprocFiles{i}); continue; end
    finalPreprocFiles.ppLabelList{1,1,i} = preprocFiles{1,1,i}.ppLabel;
    switch preprocFiles{i}.ppLabel
        case {'withinRunMoco' 'betweenRunMoco'}
            finalPreprocFiles.fTransList(:,:,i) = replace(preprocFiles{i}.fMocoList,'.nii.gz','.aff12.1D');
        case 'betweenSesMoco'
            if length(preprocFiles{i}.fMocoList)~=1; dbstack; error('figure that out'); end
            finalPreprocFiles.fTransList(:,i) = repmat(replace(preprocFiles{i}.fMocoList,'.nii.gz','.aff12.1D'),size(finalPreprocFiles.fTransList(:,i)));
        otherwise
            dbstack; error('code that')
    end
end
[finalPreprocFiles.fTransList{:,:,all(cellfun('isempty',finalPreprocFiles.fTransList),1)}] = deal('');
[finalPreprocFiles.ppLabelList{:,:,all(cellfun('isempty',finalPreprocFiles.ppLabelList),1)}] = deal('');

cmd = {src.afni};
for r = 1:size(finalPreprocFiles.fTransList,1)
    [~,b,~] = fileparts(fileparts(finalPreprocFiles.fTransList{r,1,1}));
    finalPreprocFiles.fTransCatList{r} = fullfile(bidsDerivDir,b);
    if ~exist(finalPreprocFiles.fTransCatList{r},'dir'); mkdir(finalPreprocFiles.fTransCatList{r}); end
    finalPreprocFiles.fTransCatList{r} = fullfile(finalPreprocFiles.fTransCatList{r},'transCat.aff12.1D');

    if force || ~exist(finalPreprocFiles.fTransCatList{r},'file')
        cmd{end+1} = ['rm -f ' finalPreprocFiles.fTransCatList{r}];
        cmd{end+1} = ['head -1 ' strjoin(finalPreprocFiles.fTransList(r,1,1),' ') ' > ' finalPreprocFiles.fTransCatList{r}];
        cmd{end+1} = ['cat_matvec ' strjoin(flip(finalPreprocFiles.fTransList(r,1,:)),' ') ' >> ' finalPreprocFiles.fTransCatList{r}];
    end
end



%% Apply transformations
finalPreprocFiles.fPreprocList = cell(size(finalPreprocFiles.fPlumbList));
fPreprocUpdated                = false(size(finalPreprocFiles.fPlumbList));
if size(finalPreprocFiles.fTransCatList,2)==1
    finalPreprocFiles.fTransCatList = repmat(finalPreprocFiles.fTransCatList,[1 size(finalPreprocFiles.fPreprocList,2)]);
end

for r = 1:numel(finalPreprocFiles.fPlumbList)
    fIn    = finalPreprocFiles.fPlumbList{r};
    fTrans = finalPreprocFiles.fTransCatList{r};
    fOut   = fullfile(fileparts(fIn),'preproc_volTs.nii.gz');

    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dAllineate -overwrite -nocmass -final wsinc5 \';
        cmd{end+1} = ['-source ' fIn ' \'];
        cmd{end+1} = ['-master ' fIn ' \'];
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
    [status,cmdout] = system(strjoin(cmd,newline)); if status || contains(cmdout,'error','IgnoreCase',true); dbstack; error(cmdout); error('x'); end
end


%% Rewrite at setOblique
%%% Pick the correct geometry
fGeomWR = [];
fGeomBR = [];
fGeomBS = [];
for pp = 1:length(preprocFiles)
    if ~isempty(preprocFiles{pp}) && strcmp(preprocFiles{pp}.ppLabel,'withinRunMoco')
        fGeomWR = preprocFiles{pp}.fGeom;
    end
    if ~isempty(preprocFiles{pp}) && strcmp(preprocFiles{pp}.ppLabel,'betweenRunMoco')
        if isfield(preprocFiles{pp},'fGeomSes1') && ~isempty(preprocFiles{pp}.fGeomSes1)
            fGeomBR = preprocFiles{pp}.fGeomSes1;
        else
            fGeomBR = preprocFiles{pp}.fGeom;
        end
    end
    if ~isempty(preprocFiles{pp}) && strcmp(preprocFiles{pp}.ppLabel,'betweenSesMoco')
        if isfield(preprocFiles{pp},'fGeomSes1') && ~isempty(preprocFiles{pp}.fGeomSes1)
            fGeomBS = preprocFiles{pp}.fGeomSes1;
        else
            fGeomBS = preprocFiles{pp}.fGeom;
        end
    end
end
fGeom = {initFiles.fGeom fGeomWR fGeomBR fGeomBS}; fGeom = fGeom(~cellfun('isempty',fGeom));
finalPreprocFiles.fGeom = fGeom{end}; clear fGeom fGeomWR fGeomBR fGeomBS

%%% Rewrite
disp(' rewriting preproc data at setOblique')
if any(fPreprocUpdated)
    mriOblique = MRIread(finalPreprocFiles.fGeom,1);
end
for r = 1:numel(finalPreprocFiles.fPreprocList)
    disp(['file ' num2str(r) '/' num2str(numel(finalPreprocFiles.fPreprocList))])
    if fPreprocUpdated(r)
        mri = MRIread(finalPreprocFiles.fPreprocList{r});
        mriOblique.vol = mri.vol;
        MRIwrite(mriOblique,finalPreprocFiles.fPreprocList{r});
        disp(' done')
    else
        disp(' already done, skipping')
    end
end


%% Summarize
forceThis   = force;
verboseThis = verbose;
finalPreprocFiles.fPreprocSmr = summarizeVolTs4(finalPreprocFiles.fPreprocList,0,finalPreprocFiles.dataType,forceThis,verboseThis);



%% Transform multivariate data
if ismember('PC',finalPreprocFiles.dataType)
    fReal = finalPreprocFiles.fPreprocList(:,contains(finalPreprocFiles.fPreprocList(1,:),'part-real'));
    fImag = finalPreprocFiles.fPreprocList(:,contains(finalPreprocFiles.fPreprocList(1,:),'part-imag'));
    [fMag,fPhs] = cmplx2plr(fReal,fImag,force);
    finalPreprocFiles.fPreprocList = cat(2,finalPreprocFiles.fPreprocList,fMag,fPhs);
end





