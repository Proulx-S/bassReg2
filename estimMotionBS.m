function files = estimMotionBS(funcSetSource,funcSetBase,sourceRunInd,baseRunInd,verbose,force)
global srcAfni srcFs
if ~exist('sourceRunInd','var'); sourceRunInd = []; end
if ~exist('baseRunInd','var'); baseRunInd = []; end
if ~exist('verbose','var'); verbose = []; end
if ~exist('force','var'); force = []; end
if isempty(sourceRunInd); sourceRunInd = 1; end
if isempty(baseRunInd); baseRunInd = 1; end
if isempty(verbose); verbose = 0; end
if isempty(force); force = 0; end
explicitlyFixThroughPlane = 0;
spSmFac = 0;

%% Define files and param
if isfield(funcSetBase.preprocFiles,'fCorrectedAvEchoRmsList')
    fBase = funcSetBase.preprocFiles.fCorrectedAvEchoRmsList{baseRunInd};
else
    fBase = funcSetBase.preprocFiles.fCorrectedAvList{baseRunInd};
end
if isfield(funcSetSource.preprocFiles,'fCorrectedAvEchoRmsList')
    fSource = funcSetSource.preprocFiles.fCorrectedAvEchoRmsList{sourceRunInd};
else
    fSource = funcSetSource.preprocFiles.fCorrectedAvList{sourceRunInd};
end
fMask = funcSetBase.initFiles.manBrainMaskInv;
vsize = MRIread(fBase,1); vsize = mean([vsize.xsize vsize.ysize]);

fOut = strsplit(fSource,filesep); fOut{end} = ['mcBS_' fOut{end}]; fOut = strjoin(fOut,filesep);
fOutParam = replace(fOut,'.nii.gz','');



disp('Estimating motion between sets')
cmd = {srcAfni};
%% Generate qa files before
cmd{end+1} = '3dTcat -overwrite \';
fCatBefore = strsplit(fSource,filesep); fCatBefore{end} = ['bsCatBefore_' fCatBefore{end}]; fCatBefore = strjoin(fCatBefore,filesep);
cmd{end+1} = ['-prefix ' fCatBefore ' \'];
cmd{end+1} = [fSource ' ' fBase];

%% Moco
cmd{end+1} = '3dAllineate -overwrite \';
cmd{end+1} = ['-base ' fBase ' \'];
cmd{end+1} = ['-source ' fSource ' \'];
cmd{end+1} = ['-prefix ' fOut ' \'];
cmd{end+1} = ['-1Dparam_save ' fOutParam ' \'];
cmd{end+1} = ['-1Dmatrix_save ' fOutParam ' \'];
% afni3dAlineateArg = {'-cost ls' '-interp quintic' '-final wsinc5'};
% afni3dAlineateArg = {'-cost ls+' '-interp quintic' '-final wsinc5'};
% afni3dAlineateArg = {'-cost ls+ZZ' '-interp quintic' '-final wsinc5'};
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%% so far the best %%%
% afni3dAlineateArg = {'-cost lpa' '-interp quintic' '-final wsinc5'};
% afni3dAlineateArg = {'-cost lpa+' '-interp quintic' '-final wsinc5'};
afni3dAlineateArg = {'-cost lpa+ZZ' '-interp quintic' '-final wsinc5'};
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
% afni3dAlineateArg = {'-cost lpc' '-interp quintic' '-final wsinc5'};
% afni3dAlineateArg = {'-cost lpc+' '-interp quintic' '-final wsinc5'};
% afni3dAlineateArg = {'-cost lpc+ZZ' '-interp quintic' '-final wsinc5'};
% afni3dAlineateArg = {'-cost mi' '-interp quintic' '-final wsinc5'};
% afni3dAlineateArg = {'-cost mi+ZZ' '-interp quintic' '-final wsinc5'};
% afni3dAlineateArg = {'-cost mi+' '-interp quintic' '-final wsinc5'};
% afni3dAlineateArg = {'-cost nmi' '-interp quintic' '-final wsinc5'};
% afni3dAlineateArg = {'-cost nmi+ZZ' '-interp quintic' '-final wsinc5'};
% afni3dAlineateArg = {'-cost nmi+' '-interp quintic' '-final wsinc5'};
cmd{end+1} = [strjoin(afni3dAlineateArg,' ') ' \'];
if ~isempty(fMask)
    disp([' using mask: ' fMask])
    cmd{end+1} = ['-emask ' fMask ' \'];
else
    disp(' not using mask')
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
cmd{end+1} = '3dTcat -overwrite \';
fCatAfter = strsplit(fSource,filesep); fCatAfter{end} = ['bsCatAfter_' fCatAfter{end}]; fCatAfter = strjoin(fCatAfter,filesep);
cmd{end+1} = ['-prefix ' fCatAfter ' \'];
cmd{end+1} = [fOut ' ' fBase];

%%% with the spatial smoothing applied internaly within 3dAllineate
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

%% Run command
if force || ~exist([fOutParam '.aff12.1D'],'file')
    cmd = strjoin(cmd,newline); % disp(cmd)
    if verbose
        [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    else
        [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    end
    disp([' source : ' fSource])
    disp([' base   : ' fBase])
    disp([' mask   : ' fMask])
    disp([' outputs: ' fOut])
    disp(['           ' fOutParam '.aff12.1D'])
    disp('done')
else
    disp([' source : ' fSource])
    disp([' base   : ' fBase])
    disp([' mask   : ' fMask])
    disp([' outputs: ' fOut])
    disp(['          ' fOutParam '.aff12.1D'])
    disp(' already done, skipping')
end

%% Outputs
files.fMoco = fSource;
files.fMocoParam = [fOutParam '.param.1D'];
files.fMocoMat = [fOutParam '.aff12.1D'];
files.fBase = fBase;
files.manBrainMaskInv = fMask;
files.fCatBefore = fCatBefore;
files.fCatAfter = fCatAfter;
if spSmFac
    files.fCatBeforeBlur = fCatBeforeBlur;
    files.fCatAfterBlur = fCatAfterBlur;
end


cmd = {srcFs};
cmd{end+1} = 'fslview -m single \';
if spSmFac
    cmd{end+1} = [fCatBeforeBlur ' \'];
    cmd{end+1} = [fCatAfterBlur ' \'];
end
cmd{end+1} = [fCatBefore ' \'];
cmd{end+1} = [fCatAfter ' \'];
cmd{end+1} = [fMask ' &'];


files.qaFiles.fFslviewBS = strjoin(cmd,newline);
if verbose
    system(files.qaFiles.fFslviewBS);
end
