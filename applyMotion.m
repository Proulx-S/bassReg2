function files = applyMotion(fMri,fMot,param,force,verbose)
global srcAfni srcFs

if ~exist('param','var'); param = []; end; if isempty(param); param.fMask = []; param.tSmWin_vol = []; end
if ~exist('force','var'); force = []; end
if ~exist('verbose','var'); verbose = []; end
if ~isfield(param,'tSmWin_vol'); param.tSmWin_vol = []; end
if isempty(param.tSmWin_vol); param.tSmWin_vol = 0; end
if isempty(force); force = 0; end
if isempty(verbose); verbose = 0; end

files.fOrig = fMri;
files.fCorrectedList = cell(size(fMri));
files.fCorrectedAvList = cell(size(fMri));
files.fTransList = fMot;
files.fTransCatList = cell(size(fMri));

disp('applying motion correction')
for I = 1:size(files.fOrig,1)
    disp([' run' num2str(I) '/' num2str(size(files.fOrig,1))])

    for E = 1:size(files.fOrig,2)
        disp(['  echo' num2str(E) '/' num2str(size(files.fOrig,2))])

        %%% set filename
        fIn = files.fOrig{I,E};
        fOut = strsplit(fIn,filesep); fOut{end} = ['preproc_' fOut{end}]; fOut = strjoin(fOut,filesep);
        fOutAv = strsplit(fOut,filesep); fOutAv{end} = ['av_' fOutAv{end}]; fOutAv = strjoin(fOutAv,filesep);
        fMotCat = replace(fOut,'.nii.gz','.aff12.1D');
        cmd = {srcAfni};
        % %%% 3dNwarpApply
        % cmd{end+1} = '3dNwarpApply -overwrite \';
        % cmd{end+1} = ['-source ' fIn ' \'];
        % cmd{end+1} = ['-nwarp ''IDENT(' fIn ') ' fMot{I,E} ''' \'];
        % cmd{end+1} = ['-prefix ' fOut];
        %%% 3dAllineate
        %%%% cat transformations
        cmd{end+1} = ['rm -f ' fMotCat];
        cmd{end+1} = ['head -1 ' strjoin(fMot(I,1),' ') ' > ' fMotCat];
        cmd{end+1} = ['cat_matvec ' strjoin(flip(fMot(I,:)),' ') ' >> ' fMotCat];
        
        cmd{end+1} = '3dAllineate -overwrite \';
        cmd{end+1} = ['-source ' fIn ' \'];
        cmd{end+1} = ['-1Dmatrix_apply ' fMotCat ' \'];
        cmd{end+1} = ['-prefix ' fOut];

        %%% average
        cmd{end+1} = '3dTstat -overwrite \';
        cmd{end+1} = ['-prefix ' fOutAv ' \'];
        cmd{end+1} = '-mean \';
        cmd{end+1} = fOut;

        %%% output files
        files.fCorrectedList{I,E} = fOut;
        files.fCorrectedAvList{I,E} = fOutAv;
        files.fTransCatList{I,E} = fMotCat;

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


            fCorrectedSm{I,E} = fOut;
        end

        %%% run command
        if force ||...
        (exist('fCorrectedSm','var') && ~exist(fCorrectedSm{I,E},'file')) ||...
        (isfield(files,'fCorrectedList') && ~exist(files.fCorrectedList{I,E},'file')) ||...
        (isfield(files,'fCorrectedAvList') && ~exist(files.fCorrectedAvList{I,E},'file')) ||...
        (isfield(files,'fTransCatList') && ~exist(files.fTransCatList{I,E},'file'))
        % if force || (~exist(fCorrectedSm{I,E},'file') || ~exist(files.fCorrectedList{I,E},'file') || ~exist(files.fCorrectedAvList{I,E},'file') || ~exist(files.fTransCatList{I,E},'file'))
            cmd = strjoin(cmd,newline); % disp(cmd)
            if verbose
                [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
            else
                [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
            end
            disp(['   outputs: ' files.fCorrectedList{I,E}])
            disp('   done')
        else
            disp(['   outputs: ' files.fCorrectedList{I,E}])
            disp('   already done, skipping')
        end

    end

    %%% Cross-echo processing
    cmd = {srcAfni};
    %%%% rms
    if size(files.fCorrectedList,2)>1
        fOut = files.fCorrectedList{I,1}; fOut = strsplit(fOut,'_'); fOut{contains(fOut,'echo-')} = 'echo-rms'; fOut = strjoin(fOut,'_'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        if force || ~exist(fOut,'file')
            fIn = files.fCorrectedList(I,:);
            tmp = ('a':'z'); tmp = tmp(1:length(fIn));
            expr = 'sqrt('; for i = 1:length(fIn); expr = [expr tmp(i) '*' tmp(i) '+']; end; expr(end) = ')';
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            for i = 1:length(fIn)
                cmd{end+1} = ['-' tmp(i) ' ' fIn{i} ' \'];
            end
            cmd{end+1} = ['-expr ''' expr ''''];
        end
        files.fCorrectedEchoRmsList{I,1} = fOut;

        if param.tSmWin_vol
            fOut = fCorrectedSm{I,1}; fOut = strsplit(fOut,'_'); fOut{contains(fOut,'echo-')} = 'echo-rms'; fOut = strjoin(fOut,'_'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
            if force || ~exist(fOut,'file')
                fIn = fCorrectedSm(I,:);
                tmp = ('a':'z'); tmp = tmp(1:length(fIn));
                expr = 'sqrt('; for i = 1:length(fIn); expr = [expr tmp(i) '*' tmp(i) '+']; end; expr(end) = ')';
                cmd{end+1} = '3dcalc -overwrite \';
                cmd{end+1} = ['-prefix ' fOut ' \'];
                for i = 1:length(fIn)
                    cmd{end+1} = ['-' tmp(i) ' ' fIn{i} ' \'];
                end
                cmd{end+1} = ['-expr ''' expr ''''];
            end
            fCorrectedSmRms{I,1} = fOut;
        end

        fOut = files.fCorrectedAvList{I,1}; fOut = strsplit(fOut,'_'); fOut{contains(fOut,'echo-')} = 'echo-rms'; fOut = strjoin(fOut,'_'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        if force || ~exist(fOut,'file')
            fIn = files.fCorrectedAvList(I,:);
            tmp = ('a':'z'); tmp = tmp(1:length(fIn));
            expr = 'sqrt('; for i = 1:length(fIn); expr = [expr tmp(i) '*' tmp(i) '+']; end; expr(end) = ')';
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            for i = 1:length(fIn)
                cmd{end+1} = ['-' tmp(i) ' ' fIn{i} ' \'];
            end
            cmd{end+1} = ['-expr ''' expr ''''];
        end
        files.fCorrectedAvEchoRmsList{I,1} = fOut;

        %%%%cat
        fOut = files.fCorrectedAvList{I,1}; fOut = strsplit(fOut,'_'); fOut{contains(fOut,'echo-')} = 'echo-cat'; fOut = strjoin(fOut,'_'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        if force || ~exist(fOut,'file')
            fIn = files.fCorrectedAvList(I,:);
            cmd{end+1} = '3dTcat -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = strjoin(fIn,' ');
        end
        files.fCorrectedAvEchoCatList{I,1} = fOut;
        
        %%%% run command
        disp('  cross-echo rms and cat')
        if length(cmd)>1
            cmd = strjoin(cmd,newline); % disp(cmd)
            if verbose
                [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
            else
                [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
            end
            disp('   done')
        else
            disp('   already done, skipping')
        end
    end
end

%%% write means
cmd = {srcFs};
fInList = files.fCorrectedAvList;
for E = 1:size(fInList,2)
    fOut = replace(fInList{1,E},char(regexp(fInList{1,E},'run-\d+','match')),'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
    fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
    if force || ~exist(fOut,'file')
        cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fInList(:,E),' ')];
    end
    files.fCorrectedCatAv{1,E} = fOut;

    fIn = fOut;
    fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
    if force || ~exist(fOut,'file')
        cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
    end
    files.fCorrectedAvCatAv{1,E} = fOut;
end

%%%% multi-echo
fInList = files.fCorrectedAvEchoRmsList;
fOut = replace(fInList{1,1},char(regexp(fInList{1,1},'run-\d+','match')),'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
if force || ~exist(fOut,'file')
    cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fInList(:,1),' ')];
end
files.fCorrectedCatAvEchoRmsList = fOut;

fIn = fOut;
fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
if force || ~exist(fOut,'file')
    cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
end
files.fCorrectedAvCatAvEchoRmsList = fOut;


fInList = files.fCorrectedAvCatAv;
fOut = replace(fInList{1,1},char(regexp(fInList{1,1},'echo-\d+','match')),'echo-cat'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
if force || ~exist(fOut,'file')
    cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fInList,' ')];
end
files.fCorrectedAvCatAvEchoCat = fOut;

%%% run command
disp('cross-run catenating and averaging')
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
if iscell(files.fCorrectedCatAv)
    cmd{end+1} = ['fslview -m single ' strjoin(files.fCorrectedCatAv,' ') ' &'];
else
    cmd{end+1} = ['fslview -m single ' files.fCorrectedCatAv ' &'];
end
if ~isempty(maskFile)
    cmd{end} = replace(cmd{end},' &',' \');
    cmd{end+1} = [maskFile ' &'];
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
    cmd{end+1} = ['fslview -m single ' strjoin(files.fCorrectedList(:,1),' ') ' \'];
    cmd{end+1} = [maskFile ' &'];
else
    cmd{end+1} = ['fslview -m single ' strjoin(files.fCorrectedList(:,1),' ') ' &'];
end
cmd = strjoin(cmd,newline); % disp(cmd)
if verbose
    [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
end
files.qaFiles.fFslviewWR = cmd;
%%%%% with smoothing
if param.tSmWin_vol>1
    cmd = {srcFs};
    if ~isempty(maskFile)
        cmd{end+1} = ['fslview -m single ' strjoin(fCorrectedSm(:,1),' ') ' \'];
        cmd{end+1} = [maskFile ' &'];
    else
        cmd{end+1} = ['fslview -m single ' strjoin(fCorrectedSm(:,1),' ') ' &'];
    end
    cmd = strjoin(cmd,newline); % disp(cmd)
    if verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    files.qaFilesSm.fFslviewWR = cmd;
end

%%%% motion between first, mid and last frames (acounting for smoothing) of each run
files.qaFiles.fFslviewWRfstMdLst = qaFstMdLst(files.fCorrectedList(:,1),force,verbose);
%%%%% with smoothing
if param.tSmWin_vol>1
    files.qaFilesSm.fFslviewWRfstMdLst = qaFstMdLst(fCorrectedSm(:,1),force,verbose);
    if ~isempty(maskFile)
        files = addMaskToCmd(files,maskFile);
    end
end

%%%% between-set motion
if isfield(param,'bsMoco')
    fBase = param.bsMoco.fBase;
    if isfield(files,'fCorrectedAvEchoRmsList')
        fSource = files.fCorrectedAvEchoRmsList{param.bsMoco.iSource};
    else
        fSource = files.fCorrectedAvList{param.bsMoco.iSource,1};
    end
    %%%%% generate qa file
    cmd = {srcAfni};
    cmd{end+1} = '3dTcat -overwrite \';
    fCatAfter = tempname; mkdir(fCatAfter); fCatAfter = fullfile(fCatAfter,'bsCatAfter.nii.gz');
    cmd{end+1} = ['-prefix ' fCatAfter ' \'];
    cmd{end+1} = [fSource ' ' fBase];
    %%%%% generate qa command
    cmd{end+1} = srcFs;
    cmd{end+1} = 'fslview -m single \';
    cmd{end+1} = [fCatAfter ' \'];
    cmd{end+1} = [param.maskFile ' &'];

    cmd = strjoin(cmd,newline); % disp(cmd)
    if verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    files.qaFiles.fFslviewBS = cmd;
end



