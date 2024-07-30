function files = estimMotionWR(files,param,force,verbose)
global srcAfni srcFs

if exist('param','var') && ~isempty(param) && isfield(param,'explicitlyFixThroughPlane'); explicitlyFixThroughPlane = param.explicitlyFixThroughPlane; else; explicitlyFixThroughPlane = []; end
if exist('param','var') && ~isempty(param) && isfield(param,'afni3dAlineateArg'); afni3dAlineateArg = param.afni3dAlineateArg; else; afni3dAlineateArg = {}; end
if exist('param','var') && ~isempty(param) && isfield(param,'spSmFac'); spSmFac = param.spSmFac; else; spSmFac = []; end

if ~exist('force','var'); force = []; end
if ~exist('verbose','var'); verbose = []; end
if ~exist('spSmFac','var'); spSmFac = 0; end
if isempty(explicitlyFixThroughPlane); explicitlyFixThroughPlane = 0; end
if isempty(afni3dAlineateArg); afni3dAlineateArg = {'-cost ls' '-interp quintic' '-final wsinc5'}; end
if isempty(force); force = 0; end
if isempty(verbose); verbose = 0; end

nRun = size(files.fEstim,1);

disp(['estimating within-run motion (second-pass moco to ' files.param.baseType '; accounting for smoothing)'])
fMask = files.manBrainMaskInv;
files.fMoco = cell(size(files.fEstim));
files.fMocoAv = cell(size(files.fEstim));
files.fMocoParam = cell(size(files.fEstim));
files.fMocoMat = cell(size(files.fEstim));
% files.fBase = cell(size(files.fEstim));
for I = 1:length(files.fEstim)
    disp([' run' num2str(I) '/' num2str(length(files.fEstim))])
    %%% set filename
    fIn = files.fEstim{I};
    vsize = MRIread(fIn,1); vsize = mean([vsize.xsize vsize.ysize]);
    fBase = files.fEstimBase{I};
    fOut = strsplit(fIn,filesep); fOut{end} = ['mcWR_' fOut{end}]; fOut = strjoin(fOut,filesep);
    % fOutWeights = strsplit(fIn,filesep); fOutWeights{end} = ['mcWR_' fOutWeights{end}]; fOutWeights{end} = strsplit(fOutWeights{end},'_'); fOutWeights{end}{end} = 'weights.nii.gz'; fOutWeights{end} = strjoin(fOutWeights{end},'_'); fOutWeights = strjoin(fOutWeights,filesep);
    % fOutPear = strsplit(fIn,filesep); fOutPear{end} = ['mcWR_' fOutPear{end}]; fOutPear{end} = strsplit(fOutPear{end},'_'); fOutPear{end}{end} = 'pearCor.nii.gz'; fOutPear{end} = strjoin(fOutPear{end},'_'); fOutPear = strjoin(fOutPear,filesep);
    fOutParam = replace(fOut,'.nii.gz','');
    fOutAv = strsplit(fOut,filesep); fOutAv{end} = ['av_' fOutAv{end}]; fOutAv = strjoin(fOutAv,filesep);
    if force || ~exist(fOut,'file')
        cmd = {srcAfni};
        %%% moco
        cmd{end+1} = '3dAllineate -overwrite \';
        cmd{end+1} = ['-base ' fBase ' \'];
        % cmd{end+1} = ['-source ' fIn '[11..86] \'];
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
        % cmd{end+1} = ['-wtprefix ' fOutWeights ' \'];
        % cmd{end+1} = ['-PearSave ' tempname ' \']; % not working (because single slice ??)
        cmd{end+1} = '-warp shift_rotate'; % cmd{end+1} = ['-warp shift_rotate -parfix 2 0 -parfix 4 0 -parfix 5 0'];
        % disp(strjoin(cmd,newline))

        %%% detect smoothing
        sm = strsplit(fIn,filesep); sm = strsplit(sm{end},'_'); ind = ~cellfun('isempty',regexp(sm,'^sm\d+$')); if any(ind); sm = sm{ind}; else sm = 'sm1'; end; sm = str2num(sm(3:end));
        n = MRIread(fIn,1); n = n.nframes - 1;
        nLim = [0 n] + [1 -1].*((sm+1)/2-1);
        
        %%% average
        cmd{end+1} = '3dTstat -overwrite \';
        cmd{end+1} = ['-prefix ' fOutAv ' \'];
        cmd{end+1} = '-mean \';
        cmd{end+1} = [fOut '[' num2str(nLim(1)) '..' num2str(nLim(2)) ']'];
        
        %%% execute command
        cmd = strjoin(cmd,newline); % disp(cmd)
        if verbose
            [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
        else
            [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
        end

        %%% adjust motion estimates for smoothing effects
        editMocoParam([fOutParam '.param.1D'],sm)
        editMocoParam([fOutParam '.aff12.1D'],sm)

        
        disp(['  outputs: ' fOutParam])
        disp('  done')
    else
        disp(['  outputs: ' fOutParam])
        disp('  already done, skipping')
    end

    %%% output files
    files.fMoco{I} = fOut;
    files.fMocoAv{I} = fOutAv;
    files.fMocoParam{I} = [fOutParam '.param.1D'];
    files.fMocoMat{I} = [fOutParam '.aff12.1D'];
end

%%% write means
if nRun>1
    cmd = {srcFs};
    fIn = files.fMocoAv;
    fOut = replace(fIn{1},char(regexp(fIn{1},'run-\d+','match')),'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
    fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
    if force || ~exist(fOut,'file')
        cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
    end
    files.fMocoCatAv = fOut;

    fIn = fOut;
    fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
    if force || ~exist(fOut,'file')
        cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
    end
    files.fMocoAvCatAv = fOut;
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
end

%% Outputs
%%% QA

if nRun>1
    %%%% between-run motion
    cmd = {srcFs};
    cmd{end+1} = ['fslview -m single ' files.fMocoCatAv ' \'];
    if isfield(files,'manBrainMaskInv') && ~isempty(files.manBrainMaskInv)
        cmd{end+1} = [files.manBrainMaskInv ' &'];
    end
    cmd = strjoin(cmd,newline); % disp(cmd)
    if verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    files.qaFiles.fFslviewBR = cmd;
end
%%%% within-run motion
%%%%% rewrite files with [frame] because not working with fslview grrrrr
cmd = {srcAfni};
for i = 1:length(files.fEstimBase)
    if strcmp(files.fEstimBase{i}(end),']')
        fIn = files.fEstimBase{i};
        fOut = fIn; [~,fOut] = regexp(fOut,'\[.*\]','match','split'); fOut = strsplit(fOut{1},filesep); fOut{end} = ['mcRef-first_' fOut{end}]; fOut = strjoin(fOut,filesep);

        cmd{end+1} = '3dcalc -overwrite \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = ['-a ' fIn ' \'];
        cmd{end+1} = '-expr ''a''';
        
        files.fEstimBase{i} = fOut;
    end
end
if length(cmd)>1
    cmd = strjoin(cmd,newline); % disp(cmd)
    [status,cmdout] = system(cmd); if status; dbstack; error(cmdout); error('x'); end
end
%%%%% generate the command including the bases
cmd = {srcFs};
cmd{end+1} = ['fslview -m single ' strjoin([files.fEstimBase files.fMoco],' ') ' \'];
if isfield(files,'manBrainMaskInv') && ~isempty(files.manBrainMaskInv)
    cmd{end+1} = [files.manBrainMaskInv ' &'];
end
cmd = strjoin(cmd,newline); % disp(cmd)
if verbose
    [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
end
files.qaFiles.fFslviewWR = cmd;

%%%% motion between first, mid and last frames (acounting for smoothing) of each run
files.qaFiles.fFslviewWRfstMdLst = qaFstMdLst(files.fMoco,force,verbose);
files = addMaskToCmd(files);

files.param = param;



