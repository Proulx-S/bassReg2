function files = applyMotionWR(fMri,fMot,param,force,verbose)
global srcAfni srcFs

if ~exist('param','var'); param = []; end; if isempty(param); param.fMask = []; param.tSmWin_vol = []; end
if ~exist('force','var'); force = []; end
if ~exist('verbose','var'); verbose = []; end
if isempty(param.tSmWin_vol); param.tSmWin_vol = 0; end
if isempty(force); force = 0; end
if isempty(verbose); verbose = 0; end

files.fOrig = fMri;
files.fCorrectedList = cell(size(fMri));
files.fCorrectedAvList = cell(size(fMri));

disp('applying within-run motion correction')
for I = 1:length(files.fOrig)
    disp([' run' num2str(I) '/' num2str(length(files.fOrig))])
    %%% set filename
    fIn = files.fOrig{I};
    fOut = strsplit(fIn,filesep); fOut{end} = ['preproc_' fOut{end}]; fOut = strjoin(fOut,filesep);
    fOutAv = strsplit(fOut,filesep); fOutAv{end} = ['av_' fOutAv{end}]; fOutAv = strjoin(fOutAv,filesep);
    if force || ~exist(fOut,'file')
        %%% detect smoothing
        sm = strsplit(fIn,filesep); sm = strsplit(sm{end},'_'); ind = ~cellfun('isempty',regexp(sm,'^sm\d+$')); if any(ind); sm = sm{ind}; else sm = 'sm1'; end; sm = str2num(sm(3:end));
        n = MRIread(fIn,1); n = n.nframes - 1;
        nLim = [0 n] + [1 -1].*((sm+1)/2-1);

        cmd = {srcAfni};
        % %%% 3dNwarpApply
        % cmd{end+1} = '3dNwarpApply -overwrite \';
        % cmd{end+1} = ['-source ' fIn ' \'];
        % cmd{end+1} = ['-nwarp ''IDENT(' fIn ') ' fMot{I} ''' \'];
        % cmd{end+1} = ['-prefix ' fOut];
        %%% 3dAllineate
        cmd{end+1} = '3dAllineate -overwrite \';
        cmd{end+1} = ['-source ' fIn ' \'];
        cmd{end+1} = ['-1Dmatrix_apply ' fMot{I} ' \'];
        cmd{end+1} = ['-prefix ' fOut];

        %%% average
        cmd{end+1} = '3dTstat -overwrite \';
        cmd{end+1} = ['-prefix ' fOutAv ' \'];
        cmd{end+1} = '-mean \';
        cmd{end+1} = fOut;

        %%% output files
        files.fCorrectedList{I} = fOut;
        files.fCorrectedAvList{I} = fOutAv;

        %%% smooth for visualization only
        if param.tSmWin_vol>1
            fIn = fOut;
            fOut = strsplit(fIn,filesep); fOut{end} = ['sm' num2str(param.tSmWin_vol) '_' fOut{end}]; fOut = strjoin(fOut,filesep);
            cmd{end+1} = ['echo '' ''smoothing for visualization only'];
            if force || ~exist(fOut,'file')
                cmd{end+1} = '3dTsmooth -overwrite \';
                cmd{end+1} = ['-prefix ' fOut ' \'];
                cmd{end+1} = ['-hamming ' num2str(param.tSmWin_vol) ' \'];
                cmd{end+1} = fIn;
            else
                cmd{end+1} = ['echo ''  ''already smoothed, skipping'];
            end
            fCorrectedSm{I} = fOut;
        end

        %%% run command
        cmd = strjoin(cmd,newline); % disp(cmd)
        if verbose
            [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
        else
            [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
        end
        disp('  done')
    else
        disp('  already done, skipping')
    end

end

%%% write means
cmd = {srcFs};
fIn = files.fCorrectedAvList;
fOut = replace(fIn{1},char(regexp(fIn{1},'run-\d+','match')),'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
if force || ~exist(fOut,'file')
    cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
end
files.fCorrectedCatAv = fOut;

fIn = fOut;
fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
if force || ~exist(fOut,'file')
    cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
end
files.fCorrectedAvCatAv = fOut;

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
if isfield(param,'maskFile') && ~isempty(param.maskFile) && exist(param.maskFile,'file')
    maskFile = param.maskFile;
elseif isfield(files,'manBrainMaskInv') && ~isempty(files.manBrainMaskInv)
    maskFile = files.manBrainMaskInv;
else
    maskFile = [];
end

%%%% between-run motion
cmd = {srcFs};
if ~isempty(maskFile)
    cmd{end+1} = ['fslview -m single ' files.fCorrectedCatAv ' \'];
    cmd{end+1} = [maskFile ' &'];
else
    cmd{end+1} = ['fslview -m single ' files.fCorrectedCatAv ' &'];
end
cmd = strjoin(cmd,newline); % disp(cmd)
if verbose
    [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
end
files.qaFiles.fFslviewBR = cmd;

%%%% within-run motion
%%%%% generate the command including the bases
cmd = {srcFs};
if ~isempty(maskFile)
    cmd{end+1} = ['fslview -m single ' strjoin(files.fCorrectedList,' ') ' \'];
    cmd{end+1} = [maskFile ' &'];
else
    cmd{end+1} = ['fslview -m single ' strjoin(files.fCorrectedList,' ') ' &'];
end
cmd = strjoin(cmd,newline); % disp(cmd)
if verbose
    [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
end
files.qaFiles.fFslviewWR = cmd;
%%%%% with smoothing
if sm>1
    cmd = {srcFs};
    if ~isempty(maskFile)
        cmd{end+1} = ['fslview -m single ' strjoin(fCorrectedSm,' ') ' \'];
        cmd{end+1} = [maskFile ' &'];
    else
        cmd{end+1} = ['fslview -m single ' strjoin(fCorrectedSm,' ') ' &'];
    end
    cmd = strjoin(cmd,newline); % disp(cmd)
    if verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    files.qaFiles.fFslviewWRsm = cmd;
end

%%%% motion between first, mid and last frames (acounting for smoothing) of each run
files.qaFiles.fFslviewWRfstMdLst = qaFstMdLst(files.fCorrectedList,force,verbose);
%%%%% with smoothing
if sm>1
    files.qaFiles.fFslviewWRsmFstMdLst = qaFstMdLst(fCorrectedSm,force,verbose);
    if ~isempty(maskFile)
        files = addMaskToCmd(files,maskFile);
    end
end




