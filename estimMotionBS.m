function files = estimMotionBS(funcSetSource,funcSetBase,sourceRunInd,baseRunInd,param,force)
global srcAfni srcFs
if ~exist('sourceRunInd','var'); sourceRunInd  = []; end
if ~exist('baseRunInd','var');   baseRunInd    = []; end
if ~exist('param','var');        param         = []; end
if ~isfield(param,'verbose');    param.verbose = []; end
if ~isfield(param,'cost');       param.cost    = []; end
if ~exist('force','var');        force         = []; end

if isempty(sourceRunInd);  sourceRunInd = 1       ; end
if isempty(baseRunInd);    baseRunInd   = 1       ; end
if isempty(param.verbose); verbose      = 0       ; param.verbose = verbose; else; verbose  = param.verbose; end
if isempty(param.cost);    cost         = 'lpa+ZZ'; param.cost    = cost   ; else; cost     = param.cost   ; end
if isempty(force);         force        = 0       ; end
explicitlyFixThroughPlane = 0;
spSmFac = 0;


%% Regenerate preprocFiles to undo any bsMoco that may have been applied already
% Usefull when debugging.
% NOTE that this is done for funcSetSource, not funcSetBase.


robustFlag = 0; % robust detection (will sometime force unnecessary computaitons)
for i = 1:length(funcSetSource.preprocFiles.fTrans)
    robustFlag = robustFlag + isempty(dir(fullfile(fileparts(funcSetSource.preprocFiles.fTrans{i}),'mcBS_*.aff12.1D')));
end

if robustFlag || ...
        ( isfield(funcSetSource,'preprocFiles') && isfield(funcSetSource.preprocFiles,'bsMocoFlag') && funcSetSource.preprocFiles.bsMocoFlag ) % less robust detection
    warning(['Between-set registration may have already been applied' newline '-->Regenerating preprocFiles for this new between-set registration'])
    % dbstack; warning(['between-set registration already applied' newline 'must regenerate preprocFiles to do a new between-set registration'])
    % keyboard
    % % error(['between-set registration already applied' newline 'must regenerate preprocFiles to do a new between-set registration']);

    motFiles = {};
    candidateTrans = {'wrMocoFiles' 'brMocoFiles'};
    candidateTrans = candidateTrans(ismember(candidateTrans,fields(funcSetSource)));
    for i = 1:length(candidateTrans); motFiles{end+1} = funcSetSource.(candidateTrans{i}); end
    imFiles = funcSetSource.initFiles;

    funcSetSource.preprocFiles = applyMotion(imFiles,motFiles,[],1,verbose);
    funcSetSource.preprocFiles.bsMocoFlag = 0;
end


%% Define files and param
[fBase,fBaseMask,fBaseMaskInv] = autoSelect(funcSetBase,baseRunInd);
[fSource,fSourceMask,~]        = autoSelect(funcSetSource,sourceRunInd);
vsize = MRIread(fBase,1); vsize = mean([vsize.xsize vsize.ysize]);
fOut = strsplit(fSource,filesep); fOut{end} = ['mcBS_' fOut{end}]; fOut = strjoin(fOut,filesep);
fOutParam = replace(fOut,'.nii.gz','');


%% Generate qa files before
disp('Estimating motion between sets')
cmd = {srcAfni};
mriSource = MRIread(fSource,1);
mriBase =  MRIread(fBase,1);
sameGridFlag = all([mriSource.volsize mriSource.xsize mriSource.ysize mriSource.zsize]==[mriBase.volsize mriBase.xsize mriBase.ysize mriBase.zsize]);
if sameGridFlag
    cmd{end+1} = '3dTcat -overwrite \';
    fCatBefore = strsplit(fSource,filesep); fCatBefore{end} = ['bsCatBefore_' fCatBefore{end}]; fCatBefore = strjoin(fCatBefore,filesep);
    cmd{end+1} = ['-prefix ' fCatBefore ' \'];
    cmd{end+1} = [fSource ' ' fBase];
end


%% Moco
cmd{end+1} = '3dAllineate -overwrite \';
cmd{end+1} = ['-base ' fBase ' \'];
cmd{end+1} = ['-source ' fSource ' \'];
cmd{end+1} = ['-prefix ' fOut ' \'];
cmd{end+1} = ['-1Dparam_save ' fOutParam ' \'];
cmd{end+1} = ['-1Dmatrix_save ' fOutParam ' \'];

afni3dAlineateArg = {['-cost ' cost] '-interp quintic' '-final wsinc5'};
% % afni3dAlineateArg = {'-cost ls' '-interp quintic' '-final wsinc5'};
% % afni3dAlineateArg = {'-cost ls+' '-interp quintic' '-final wsinc5'};
% % afni3dAlineateArg = {'-cost ls+ZZ' '-interp quintic' '-final wsinc5'};
% %%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%
% %%% so far the best %%%
% % afni3dAlineateArg = {'-cost lpa' '-interp quintic' '-final wsinc5'};
% % afni3dAlineateArg = {'-cost lpa+' '-interp quintic' '-final wsinc5'};
% afni3dAlineateArg = {'-cost lpa+ZZ' '-interp quintic' '-final wsinc5'};
% %%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%
% % afni3dAlineateArg = {'-cost lpc' '-interp quintic' '-final wsinc5'};
% % afni3dAlineateArg = {'-cost lpc+' '-interp quintic' '-final wsinc5'};
% % afni3dAlineateArg = {'-cost lpc+ZZ' '-interp quintic' '-final wsinc5'};
% % afni3dAlineateArg = {'-cost mi' '-interp quintic' '-final wsinc5'};
% % afni3dAlineateArg = {'-cost mi+ZZ' '-interp quintic' '-final wsinc5'};
% % afni3dAlineateArg = {'-cost mi+' '-interp quintic' '-final wsinc5'};
% % afni3dAlineateArg = {'-cost nmi' '-interp quintic' '-final wsinc5'};
% % afni3dAlineateArg = {'-cost nmi+ZZ' '-interp quintic' '-final wsinc5'};
% % afni3dAlineateArg = {'-cost nmi+' '-interp quintic' '-final wsinc5'};
cmd{end+1} = [strjoin(afni3dAlineateArg,' ') ' \'];
if ~isempty(fBaseMaskInv)
    disp([' using exclusion mask on base: ' fBaseMaskInv])
    cmd{end+1} = ['-emask ' fBaseMaskInv ' \'];
else
    disp(' not using exclusion mask on base')
end
if ~isempty(fSourceMask)
    disp([' using inclusion mask on source: ' fSourceMask])
    cmd{end+1} = ['-source_mask ' fSourceMask ' \'];
else
    disp(' not using inclusion mask on source')
end
if explicitlyFixThroughPlane
    cmd{end+1} = '-parfix 2 0 -parfix 4 0 -parfix 5 0 \';
    disp(' explicitly enforcing no through-plane motion')
end
cmd{end+1} = '-nopad -conv 0 -nmatch 100% -onepass -nocmass \';
% cmd{end+1} = '-nopad -conv 0 -nmatch 100% -nocmass \';
% cmd{end+1} = '-maxrot 2 -maxshf 1 \';
if spSmFac>0
    cmd{end+1} = ['-fineblur ' num2str(vsize*spSmFac) ' \'];
end
cmd{end+1} = '-warp shift_rotate'; % cmd{end+1} = ['-warp shift_rotate -parfix 2 0 -parfix 4 0 -parfix 5 0'];
if verbose; disp(strjoin(cmd,newline)); end

%% Generate qa files after
if sameGridFlag
    cmd{end+1} = '3dTcat -overwrite \';
    fCatAfter = strsplit(fSource,filesep); fCatAfter{end} = ['bsCatAfter_' fCatAfter{end}]; fCatAfter = strjoin(fCatAfter,filesep);
    cmd{end+1} = ['-prefix ' fCatAfter ' \'];
    cmd{end+1} = [fOut ' ' fBase];
end

%%% with the spatial smoothing applied internaly within 3dAllineate
if sameGridFlag
    if spSmFac
        cmd{end+1} = '3dmerge -overwrite \';
        cmd{end+1} = ['-1blur_fwhm ' num2str(vsize*spSmFac) ' \'];
        cmd{end+1} = '-doall \';
        fCatAfterBlur = strsplit(fCatAfter,filesep); fCatAfterBlur{end} = ['blur_' fCatAfterBlur{end}]; fCatAfterBlur = strjoin(fCatAfterBlur,filesep);
        cmd{end+1} = ['-prefix ' fCatAfterBlur ' \'];
        cmd{end+1} = fCatAfter;

        cmd{end+1} = '3dmerge -overwrite \';
        cmd{end+1} = ['-1blur_fwhm ' num2str(vsize*spSmFac) ' \'];
        cmd{end+1} = '-doall \';
        fCatBeforeBlur = strsplit(fCatBefore,filesep); fCatBeforeBlur{end} = ['blur_' fCatBeforeBlur{end}]; fCatBeforeBlur = strjoin(fCatBeforeBlur,filesep);
        cmd{end+1} = ['-prefix ' fCatBeforeBlur ' \'];
        cmd{end+1} = fCatBefore;
    end
end

%% Run command
if force || ~exist([fOutParam '.aff12.1D'],'file')
    cmd = strjoin(cmd,newline); % disp(cmd)
    if verbose
        [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    else
        [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    end
    disp([' source       : ' fSource])
    disp([' base         : ' fBase])
    disp([' mask (base)  : ' fBaseMaskInv])
    disp([' mask (source): ' fSourceMask])
    disp([' outputs      : ' fOut])
    disp(['                ' fOutParam '.aff12.1D'])
    disp('done')
else
    disp([' source       : ' fSource])
    disp([' base         : ' fBase])
    disp([' mask (base)  : ' fBaseMaskInv])
    disp([' mask (source): ' fSourceMask])
    disp([' outputs      : ' fOut])
    disp(['                ' fOutParam '.aff12.1D'])
    disp(' already done, skipping')
end

%% Outputs
files.fMoco = {fSource};
files.fMocoParam = {[fOutParam '.param.1D']};
files.fMocoMat = {[fOutParam '.aff12.1D']};
files.fBase = {fBase};
files.manBrainMaskInv = fBaseMaskInv;
files.manBrainMask = fBaseMask;
files.manBrainMask_source = fSourceMask;
if sameGridFlag
    files.fCatBefore = fCatBefore;
    files.fCatAfter = fCatAfter;
end
if spSmFac && sameGridFlag
    files.fCatBeforeBlur = fCatBeforeBlur;
    files.fCatAfterBlur = fCatAfterBlur;
end

if sameGridFlag
    cmd = {srcFs};
    cmd{end+1} = 'fslview -m single \';
    if spSmFac
        cmd{end+1} = [fCatBeforeBlur ' \'];
        cmd{end+1} = [fCatAfterBlur ' \'];
    end
    cmd{end+1} = [fCatBefore ' \'];
    cmd{end+1} = [fCatAfter ' \'];
    cmd{end+1} = [fBaseMaskInv ' &'];
else
    cmd = {srcFs};
    cmd{end+1} = 'freeview \';
    cmd{end+1} = [fBase ':name=' funcSetBase.label ' \'];
    cmd{end+1} = [fSource ':name=' funcSetSource.label '_before:visible=0 \'];
    cmd{end+1} = [fOut ':name=' funcSetSource.label '_after \'];
    cmd{end+1} = [fBaseMaskInv ' &'];
end

files.qaFiles.fFslviewBS = strjoin(cmd,newline);
if verbose
    if any(contains(cmd,'freeview'))
        clipboard('copy',files.qaFiles.fFslviewBS)
        disp('QA command copied to clipboard')
    else
        system(files.qaFiles.fFslviewBS);
    end
end



