function files = applyMotion(fMri,fMot,param,force,verbose)
global srcAfni srcFs

% if ~isstruct(fMri); tmp = fMri; fMri = []; fMri.fPlumb = tmp; end
% if ~isstruct(fMot); tmp = fMot; fMot = []; fMot.fMocoMat = tmp; end

if ~exist('param','var'); param = []; end; if isempty(param); param.fMask = []; param.tSmWin_vol = []; end
if ~exist('force','var'); force = []; end
if ~exist('verbose','var'); verbose = []; end
if ~isfield(param,'tSmWin_vol'); param.tSmWin_vol = []; end
if isempty(param.tSmWin_vol); param.tSmWin_vol = 0; end
if isempty(force); force = 0; end
if isempty(verbose); verbose = 0; end


files.fOrig = fMri.fPlumb;
[nRun,nEcho] = size(files.fOrig); nFrame = MRIread(files.fOrig{1,1}); nFrame = nFrame.nframes;
files.fCorrected = cell(size(fMri.fPlumb));
files.fCorrectedAv = cell(size(fMri.fPlumb));
if ~isempty(fMot)
    if iscell(fMot)
        sz = [1 1]; for i = 1:length(fMot); sz = max(sz,size(fMot{i}.fMocoMat)); end
        files.fTrans      = cell([sz length(fMot)]);
        files.fTransLabel = cell([1 1 length(fMot)]);
        for i = 1:length(fMot)
            files.fTrans(:,:,i) = fMot{i}.fMocoMat;
            [~,files.fTransLabel{1,1,i}] = fileparts(fMot{i}.fMocoMat{1});
            files.fTransLabel{1,1,i} = strsplit(files.fTransLabel{1,1,i},'_'); files.fTransLabel{1,1,i} = files.fTransLabel{1,1,i}{1};
        end
    else
        files.fTrans      = fMot.fMocoMat;
        files.fTransLabel = cell([1 1 size(fMot.fMocoMat,3)]);
        for i = 1:size(fMot.fMocoMat,3)
            [~,files.fTransLabel{1,1,i}] = fileparts(fMot.fMocoMat{1,1,i});
            files.fTransLabel{1,1,i} = strsplit(files.fTransLabel{1,1,i},'_'); files.fTransLabel{1,1,i} = files.fTransLabel{1,1,i}{1};
        end
    end
else
    sz = size(fMri.fPlumb); sz(2) = 1;
    files.fTrans = repmat({'''IDENTITY'''},sz);
    sz(1) = 1;
    files.fTransLabel = repmat({'''IDENTITY'''},sz);
end
files.fTransCat = cell(size(fMri.fPlumb));

disp('applying motion correction')
for I = 1:size(files.fOrig,1)
    disp([' run' num2str(I) '/' num2str(size(files.fOrig,1))])

    for E = 1:size(files.fOrig,2)
        disp(['  echo' num2str(E) '/' num2str(size(files.fOrig,2))])

        %%% set filename
        fIn = files.fOrig{I,E};
        fOut = strsplit(fIn,filesep); fOut{end} = ['preproc_' fOut{end}]; fOut = strjoin(fOut,filesep);
        fOutAv = strsplit(fOut,filesep); fOutAv{end} = ['av_' fOutAv{end}]; fOutAv = strjoin(fOutAv,filesep);
        cmd = {srcAfni};
        % %%% 3dNwarpApply
        % cmd{end+1} = '3dNwarpApply -overwrite \';
        % cmd{end+1} = ['-source ' fIn ' \'];
        % cmd{end+1} = ['-nwarp ''IDENT(' fIn ') ' fMot{I,E} ''' \'];
        % cmd{end+1} = ['-prefix ' fOut];
        %%% 3dAllineate
        %%%% cat transformations
        if strcmp(files.fTrans{I,1},'''IDENTITY''')
            fMotCat = '''IDENTITY''';
        else
            fMotCat = replace(fOut,'.nii.gz','.aff12.1D');
            if size(files.fTrans,3)>1
                cmd{end+1} = ['rm -f ' fMotCat];
                cmd{end+1} = ['head -1 ' strjoin(files.fTrans(I,1,1),' ') ' > ' fMotCat];
                cmd{end+1} = ['cat_matvec ' strjoin(flip(files.fTrans(I,1,:)),' ') ' >> ' fMotCat];
            else
                fMotCat = files.fTrans{I,1,1};
            end
        end

        cmd{end+1} = '3dAllineate -overwrite \';
        cmd{end+1} = ['-source ' fIn ' \'];
        cmd{end+1} = ['-1Dmatrix_apply ' fMotCat ' \'];
        cmd{end+1} = ['-prefix ' fOut];


        %%% average
        if nFrame>1
            cmd{end+1} = '3dTstat -overwrite \';
            cmd{end+1} = ['-prefix ' fOutAv ' \'];
            cmd{end+1} = '-mean \';
            cmd{end+1} = fOut;
        else
            fOutAv = fOut;
        end

        %%% output files
        files.fCorrected{I,E} = fOut;
        files.fCorrectedAv{I,E} = fOutAv;
        files.fTransCat{I,E} = fMotCat;

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
        (isfield(files,'fCorrected') && ~exist(files.fCorrected{I,E},'file')) ||...
        (isfield(files,'fCorrectedAv') && ~exist(files.fCorrectedAv{I,E},'file')) ||...
        (isfield(files,'fTransCat') && ~exist(files.fTransCat{I,E},'file'))
            updatedFlag = 1;
        % if force || (~exist(fCorrectedSm{I,E},'file') || ~exist(files.fCorrected{I,E},'file') || ~exist(files.fCorrectedAv{I,E},'file') || ~exist(files.fTransCat{I,E},'file'))
            cmd = strjoin(cmd,newline); % disp(cmd)
            if verbose
                [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
            else
                [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
            end
            
            disp(['   outputs: ' files.fCorrected{I,E}])
            if ischar(fMotCat) && strcmp(fMotCat,'''IDENTITY''')
                disp('   done (IDENTITY)')
            else
                disp('   done')
            end
        else
            updatedFlag = 0;
            disp(['   outputs: ' files.fCorrected{I,E}])
            disp('   already done, skipping')
        end

    end

    %%% Cross-echo processing
    cmd = {srcAfni};
    %%%% rms
    if nEcho>1
        fOut = files.fCorrected{I,1}; fOut = strsplit(fOut,'_'); fOut{contains(fOut,'echo-')} = 'echo-rms'; fOut = strjoin(fOut,'_'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        if force || ~exist(fOut,'file')
            fIn = files.fCorrected(I,:);
            tmp = ('a':'z'); tmp = tmp(1:length(fIn));
            expr = 'sqrt('; for i = 1:length(fIn); expr = [expr tmp(i) '*' tmp(i) '+']; end; expr(end) = ')';
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            for i = 1:length(fIn)
                cmd{end+1} = ['-' tmp(i) ' ' fIn{i} ' \'];
            end
            cmd{end+1} = ['-expr ''' expr ''''];
        end
        files.fCorrectedEchoRms{I,1} = fOut;

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

        %%% rms the average only when relevant
        if nFrame>1
            fOut = files.fCorrectedAv{I,1}; fOut = strsplit(fOut,'_'); fOut{contains(fOut,'echo-')} = 'echo-rms'; fOut = strjoin(fOut,'_'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
            if force || ~exist(fOut,'file')
                fIn = files.fCorrectedAv(I,:);
                tmp = ('a':'z'); tmp = tmp(1:length(fIn));
                expr = 'sqrt('; for i = 1:length(fIn); expr = [expr tmp(i) '*' tmp(i) '+']; end; expr(end) = ')';
                cmd{end+1} = '3dcalc -overwrite \';
                cmd{end+1} = ['-prefix ' fOut ' \'];
                for i = 1:length(fIn)
                    cmd{end+1} = ['-' tmp(i) ' ' fIn{i} ' \'];
                end
                cmd{end+1} = ['-expr ''' expr ''''];
            end
            files.fCorrectedAvEchoRms{I,1} = fOut;
        end

        %%%%cat
        fOut = files.fCorrectedAv{I,1}; fOut = strsplit(fOut,'_'); fOut{contains(fOut,'echo-')} = 'echo-cat'; fOut = strjoin(fOut,'_'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        if force || ~exist(fOut,'file')
            fIn = files.fCorrectedAv(I,:);
            cmd{end+1} = '3dTcat -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = strjoin(fIn,' ');
        end
        files.fCorrectedAvEchoCat{I,1} = fOut;
        
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
if nRun>1
    cmd = {srcFs};
    fInList = files.fCorrectedAv;
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
    fInList = files.fCorrectedAvEchoRms;
    fOut = replace(fInList{1,1},char(regexp(fInList{1,1},'run-\d+','match')),'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
    fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
    if force || ~exist(fOut,'file')
        cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fInList(:,1),' ')];
    end
    if iscell(fOut)
        files.fCorrectedCatAvEchoRms = fOut;
    else
        files.fCorrectedCatAvEchoRms = {fOut};
    end

    fIn = fOut;
    fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
    if force || ~exist(fOut,'file')
        cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
    end
    if iscell(fOut)
        files.fCorrectedAvCatAvEchoRms = fOut;
    else
        files.fCorrectedAvCatAvEchoRms = {fOut};
    end


    fInList = files.fCorrectedAvCatAv;
    fOut = replace(fInList{1,1},char(regexp(fInList{1,1},'echo-\d+','match')),'echo-cat'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
    if force || ~exist(fOut,'file')
        cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fInList,' ')];
    end
    if iscell(fOut)
        files.fCorrectedAvCatAvEchoCat = fOut;
    else
        files.fCorrectedAvCatAvEchoCat = {fOut};
    end

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
end

%% Outputs
if isfield(fMri,'manBrainMask')
    files.manBrainMask = fMri.manBrainMask;
elseif iscell(fMot)
    files.manBrainMask = fMot{1}.manBrainMask;
elseif isstruct(fMot)
    files.manBrainMask = fMot.manBrainMask;
end
if isfield(fMri,'manBrainMaskInv')
    files.manBrainMaskInv = fMri.manBrainMaskInv;
elseif iscell(fMot)
    files.manBrainMaskInv = fMot{1}.manBrainMaskInv;
elseif isstruct(fMot)
    files.manBrainMaskInv = fMot.manBrainMaskInv;
end
files.updatedFlag         = updatedFlag;


%% QA
if isfield(param,'maskFile') && ~isempty(param.maskFile) && exist(param.maskFile,'file')
    maskFile = param.maskFile;
elseif isfield(files,'manBrainMaskInv') && ~isempty(files.manBrainMaskInv)
    maskFile = files.manBrainMaskInv;
else
    maskFile = [];
end

%%% between-run motion
if nRun>1
    cmd = {srcFs};
    cmd{end+1} = ['fslview -m single ' strjoin(files.fCorrectedCatAv,' ') ' &'];
    if ~isempty(maskFile)
        cmd{end} = replace(cmd{end},' &',' \');
        cmd{end+1} = [maskFile ' &'];
    end
    cmd = strjoin(cmd,newline); % disp(cmd)
    if verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    files.qaFiles.fFslviewBR = cmd;
end

%%% within-run motion
%%%% generate the command including the bases
if nFrame>1
    cmd = {srcFs};
    if ~isempty(maskFile)
        cmd{end+1} = ['fslview -m single ' strjoin(files.fCorrected(:,1),' ') ' \'];
        cmd{end+1} = [maskFile ' &'];
    else
        cmd{end+1} = ['fslview -m single ' strjoin(files.fCorrected(:,1),' ') ' &'];
    end
    cmd = strjoin(cmd,newline); % disp(cmd)
    if verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    files.qaFiles.fFslviewWR = cmd;
    %%%% with smoothing
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

    %%% motion between first, mid and last frames (acounting for smoothing) of each run
    files.qaFiles.fFslviewWRfstMdLst = qaFstMdLst(files.fCorrected(:,1),force,verbose);
    %%%% with smoothing
    if param.tSmWin_vol>1
        files.qaFilesSm.fFslviewWRfstMdLst = qaFstMdLst(fCorrectedSm(:,1),force,verbose);
        if ~isempty(maskFile)
            files = addMaskToCmd(files,maskFile);
        end
    end
end

%%% between-set motion
bsInd = find(contains(files.fTransLabel,'mcBS'));
if ~isempty(bsInd)
% bsFlag = files.fTrans; if iscell(bsFlag); bsFlag = bsFlag{1}; end; bsFlag = strsplit(bsFlag,filesep); bsFlag = bsFlag{end};
% bsFlag = contains(bsFlag,'mcBS_');
% if bsFlag
    fBase = fMot{bsInd}.fBase{1};
    fSource = fMot{bsInd}.fMoco{1};
    if nEcho>1
        fOut = {'fCorrectedAvEchoRms' 'fCorrectedEchoRms'};
    else
        fOut = {'fCorrectedAv' 'fCorrected'};
    end
    fOut = fOut(ismember(fOut,fields(files))); fOut = char(files.(fOut{1}));
    fMask = fMot{bsInd}.manBrainMaskInv;
    
    mriBase =  MRIread(fBase,1); mriSource = MRIread(fSource,1);
    sameGridFlag = all([mriSource.volsize mriSource.xsize mriSource.ysize mriSource.zsize]==[mriBase.volsize mriBase.xsize mriBase.ysize mriBase.zsize]);

    % fBase = param.bsMoco.fBase;
    % if isfield(files,'fCorrectedAvEchoRms')
    %     fSource = files.fCorrectedAvEchoRms{param.bsMoco.iSource};
    % else
    %     fSource = files.fCorrectedAv{param.bsMoco.iSource,1};
    % end

    %%%% generate qa file
    if sameGridFlag
        %%%%% before
        %too complicated
        %%%%% after
        cmd = {srcAfni};
        cmd{end+1} = '3dTcat -overwrite \';
        fCatAfter = tempname; mkdir(fCatAfter); fCatAfter = fullfile(fCatAfter,'bsCatAfter.nii.gz');
        cmd{end+1} = ['-prefix ' fCatAfter ' \'];
        cmd{end+1} = [fSource ' ' fBase];
        %%%% generate qa command
        cmd{end+1} = srcFs;
        cmd{end+1} = 'fslview -m single \';
        cmd{end+1} = [fCatAfter ' \'];
        cmd{end+1} = [param.maskFile ' &'];
    else
        cmd = {srcFs};
        cmd{end+1} = 'freeview \';
        cmd{end+1} = [fBase ' \'];
        cmd{end+1} = [fSource ':visible=0 \'];
        cmd{end+1} = [fOut ' \'];
        cmd{end+1} = [fMask ' &'];
    end

    cmd = strjoin(cmd,newline); % disp(cmd)
    if verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    files.qaFiles.fFslviewBS = cmd;
end


