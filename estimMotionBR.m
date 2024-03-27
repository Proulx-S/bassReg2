function files = estimMotionBR(fSource,fBase,fMask,param,force,verbose)
global srcAfni srcFs

if exist('param','var') && ~isempty(param) && isfield(param,'explicitlyFixThroughPlane'); explicitlyFixThroughPlane = param.explicitlyFixThroughPlane; else; explicitlyFixThroughPlane = []; end
if exist('param','var') && ~isempty(param) && isfield(param,'afni3dAlineateArg'); afni3dAlineateArg = param.afni3dAlineateArg; else; afni3dAlineateArg = {}; end
if exist('param','var') && ~isempty(param) && isfield(param,'spSmFac'); spSmFac = param.spSmFac; else; spSmFac = []; end

if ~exist('fMask','var'); fMask = []; end
if ~exist('force','var'); force = []; end
if ~exist('verbose','var'); verbose = []; end
if ~exist('spSmFac','var'); spSmFac = 0; end
if isempty(explicitlyFixThroughPlane); explicitlyFixThroughPlane = 0; end
if isempty(afni3dAlineateArg); afni3dAlineateArg = {'-cost ls' '-interp quintic' '-final wsinc5'}; end
if isempty(force); force = 0; end
if isempty(verbose); verbose = 0; end

disp(['estimating between-run motion (moco to ' param.baseType ')'])
files.fMocoList = cell(size(fSource));
files.fMocoAvList = cell(size(fSource));
files.fMocoParamList = cell(size(fSource));
files.fMocoMatList = cell(size(fSource));
files.fBaseList = cell(size(fSource));
files.manBrainMaskInv = fMask;
for I = 1:numel(fSource)
    disp([' run' num2str(I) '/' num2str(length(fSource))])
    %%% set filename
    fIn = fSource{I};
    vsize = MRIread(fIn,1); vsize = mean([vsize.xsize vsize.ysize]);
    fOut = strsplit(fIn,filesep); fOut{end} = ['mcBR_' fOut{end}]; fOut = strjoin(fOut,filesep);
    fOutWeights = strsplit(fIn,filesep); fOutWeights{end} = ['mcBR_' fOutWeights{end}]; fOutWeights{end} = strsplit(fOutWeights{end},'_'); fOutWeights{end}{end} = 'weights.nii.gz'; fOutWeights{end} = strjoin(fOutWeights{end},'_'); fOutWeights = strjoin(fOutWeights,filesep);
    fOutParam = replace(fOut,'.nii.gz','');
    % fOutAv = strsplit(fOut,filesep); fOutAv{end} = ['av_' fOutAv{end}]; fOutAv = strjoin(fOutAv,filesep);
    if force || ~exist(fOut,'file')
        cmd = {srcAfni};
        %%% moco
        cmd{end+1} = '3dAllineate -overwrite \';
        cmd{end+1} = ['-base ' fBase ' \'];
        cmd{end+1} = ['-source ' fIn ' \'];
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = ['-1Dparam_save ' fOutParam ' \'];
        cmd{end+1} = ['-1Dmatrix_save ' fOutParam ' \'];
        % afni3dAlineateArg = {'-cost ls' '-interp quintic' '-final wsinc5'};
        afni3dAlineateArg = {'-cost lpa+ZZ' '-interp quintic' '-final wsinc5'};
        cmd{end+1} = [strjoin(afni3dAlineateArg,' ') ' \'];
        if ~isempty(fMask)
            disp(['  using mask: ' fMask])
            cmd{end+1} = ['-emask ' fMask ' \'];
        else
            disp('  not using mask')
        end
        if explicitlyFixThroughPlane
            cmd{end+1} = '-parfix 2 0 -parfix 4 0 -parfix 5 0 \';
            disp('  explicitly enforcing no through-plane motion')
        end
        cmd{end+1} = '-nopad -conv 0 -nmatch 100% -onepass -nocmass \'; 
        % cmd{end+1} = '-maxrot 1 -maxshf 0.5 \'; 
        if spSmFac>0
            cmd{end+1} = ['-fineblur ' num2str(vsize*spSmFac) ' \'];
        end
        cmd{end+1} = ['-wtprefix ' fOutWeights ' \'];
        cmd{end+1} = '-warp shift_rotate'; % cmd{end+1} = ['-warp shift_rotate -parfix 2 0 -parfix 4 0 -parfix 5 0'];
        if verbose; disp(strjoin(cmd,newline)); end

        %%% detect smoothing
        sm = strsplit(fIn,filesep); sm = strsplit(sm{end},'_'); ind = ~cellfun('isempty',regexp(sm,'^sm\d+$')); if any(ind); sm = sm{ind}; else sm = 'sm1'; end; sm = str2num(sm(3:end));
        n = MRIread(fIn,1); n = n.nframes - 1;
        nLim = [0 n] + [1 -1].*((sm+1)/2-1);
        
        %%% execute command
        cmd = strjoin(cmd,newline); % disp(cmd)
        if verbose
            [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
        else
            [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
        end

        disp(['  outputs: ' fOutParam])
        disp('  done')
    else
        disp(['  outputs: ' fOutParam])
        disp('  already done, skipping')
    end

    %%% output files
    files.fMocoList{I} = fOut;
    files.fMocoParamList{I} = [fOutParam '.param.1D'];
    files.fMocoMatList{I} = [fOutParam '.aff12.1D'];
    files.fBaseList{I} = fBase;
end

%%% write means
cmd = {srcFs};
fIn = files.fMocoList;
fOut = replace(fIn{1},char(regexp(fIn{1},'run-\d+','match')),'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
% fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
if force || ~exist(fOut,'file')
    cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
end
files.fMocoCat = fOut;

fIn = fOut;
fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
if force || ~exist(fOut,'file')
    cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
end
files.fMocoAvCat = fOut;

disp('averaging')
if length(cmd)>1
    cmd = strjoin(cmd,newline); % disp(cmd)
    if verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    else
        [status,cmdout] = system(cmd); if status; dbstack; error(cmdout); error('x'); end
    end
    disp(' done')
else
    disp(' already done, skipping')
end

%% Outputs
%%% QA

%%%% between-run motion
cmd = {srcFs};
cmd{end+1} = ['fslview -m single ' files.fMocoCat ' \'];
if isfield(files,'manBrainMaskInv') && ~isempty(files.manBrainMaskInv)
    cmd{end+1} = [files.manBrainMaskInv ' &'];
end
cmd = strjoin(cmd,newline); % disp(cmd)
if verbose
    [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
end
files.qaFiles.fFslviewBR = cmd;

