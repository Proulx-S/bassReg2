function files = funcAnaWithRespMode(fFunc,param,force,verbose)
global srcAfni srcFs
if ~exist('force'  ,'var');   force = []; end
if ~exist('verbose','var'); verbose = []; end
if isempty(force)         ;   force = 0 ; end
if isempty(verbose)       ; verbose = 0 ; end

sz = size(fFunc.f,[1 2]);
nRun = sz(1); nEcho = sz(2);


fMask = fFunc.manBrainMask;
% if contains(fMask,'Inv.nii.gz')
%     fMask = replace(fMask,'Inv.nii.gz','.nii.gz');
%     if exist(fMask,'file')
%         warning(['inverse of brain mask is provided?' newline 'inverting it back']);
%     else
%         dbstack; error('code that');
%     end
% end

%% Run afni's 3dDeconvolve

%%% Per echo analysis
disp(['Functional analysis (custom response model): ' param.label])
[fRun,fSes,param] = runAfni(fFunc.f,param,fMask,force,verbose); % analysis performed on each echoe within that function
if isempty(fSes); fSes = fRun; end

%%% Cross-echo RMS analysis
if nEcho>1
    disp(['Functional analysis (cross-echo rms): ' param.label])
    [fRun_echoRms,fSes_echoRms] = runAfni(fFunc.fEchoRms,param,fMask,force,verbose);
    if isempty(fSes_echoRms); fSes_echoRms = fRun_echoRms; end
    % fRun = cat(2,fRun,fRun_echoRms);
    % fSes = cat(2,fSes,fSes_echoRms);
end


%% Simplify stat outputs
label = strsplit(fFunc.f{1},filesep); label = replace(label{contains(label,'set-')},'set-','');
disp('Simplifying stat outputs')
forceThis = force;
cmd = {srcAfni};

%%% Individual runs
for R = 1:nRun

    %%%% Individual echo
    for E = 1:nEcho
        fIn = fRun(R,E).fStat;

        %%%%% t-stat
        fOut = replace(fIn,'_stats.nii.gz','_model-custom_respT.nii.gz');
        fRun(R,E).fRespT = fOut;
        if forceThis || ~exist(fOut,'file')
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fIn '[' label '#0_Tstat] 2> /dev/null'];
        end

        %%%%% amplitude
        fOut = replace(fIn,'_stats.nii.gz','_model-custom_respAmp.nii.gz');
        fRun(R,E).fRespAmp = fOut;
        if forceThis || ~exist(fOut,'file')
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fIn '[' label '#0_Coef] 2> /dev/null'];
        end

        %%%%% baseline
        fOut = replace(fIn,'_stats.nii.gz','_model-custom_base.nii.gz');
        fRun(R,E).fBase = fOut;
        if forceThis || ~exist(fOut,'file')
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fIn '[Run#1Pol#0_Coef] 2> /dev/null'];
        end
    end

    %%%% cross-echo rms
    if nEcho>1
        fIn = fRun_echoRms(R,1).fStat;

        %%%%% t-stat
        fOut = replace(fIn,'_stats.nii.gz','_model-custom_respT.nii.gz');
        fRun_echoRms(R,1).fRespT = fOut;
        if forceThis || ~exist(fOut,'file')
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fIn '[' label '#0_Tstat] 2> /dev/null'];
        end

        %%%%% amplitude
        fOut = replace(fIn,'_stats.nii.gz','_model-custom_respAmp.nii.gz');
        fRun_echoRms(R,1).fRespAmp = fOut;
        if forceThis || ~exist(fOut,'file')
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fIn '[' label '#0_Coef] 2> /dev/null'];
        end

        %%%%% baseline
        fOut = replace(fIn,'_stats.nii.gz','_model-custom_base.nii.gz');
        fRun_echoRms(R,1).fBase = fOut;
        if forceThis || ~exist(fOut,'file')
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fIn '[Run#1Pol#0_Coef] 2> /dev/null'];
        end
    end
end

%%% Whole session
if nRun>1
    R = 1;
    %%%% Individual echo
    for E = 1:nEcho
        fIn = fSes(R,E).fStat;

        %%%%% t-stat
        fOut = replace(fIn,'_stats.nii.gz','_model-custom_respT.nii.gz');
        fSes(R,E).fRespT = fOut;
        if forceThis || ~exist(fOut,'file')
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fIn '[' label '#0_Tstat] 2> /dev/null'];
        end

        %%%%% amplitude
        fOut = replace(fIn,'_stats.nii.gz','_model-custom_respAmp.nii.gz');
        fSes(R,E).fRespAmp = fOut;
        if forceThis || ~exist(fOut,'file')
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fIn '[' label '#0_Coef] 2> /dev/null'];
        end

        %%%%% baseline
        fOut = replace(fIn,'_stats.nii.gz','_model-custom_base.nii.gz');
        fSes(R,E).fBase = fOut;
        if forceThis || ~exist(fOut,'file')
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            buck = num2str(1:nRun,'Run#%iPol#0_Coef,'); buck(end) = [];
            cmd{end+1} = [fIn '[' buck '] 2> /dev/null'];
            cmd{end+1} = '3dTstat -overwrite -mean \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fOut ' 2> /dev/null'];
        end
    end

    %%%% cross-echo rms
    if nEcho>1
        E = 1;
        fIn = fSes_echoRms(R,E).fStat;

        %%%%% t-stat
        fOut = replace(fIn,'_stats.nii.gz','_model-custom_respT.nii.gz');
        fSes_echoRms(R,E).fRespT = fOut;
        if forceThis || ~exist(fOut,'file')
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fIn '[' label '#0_Tstat] 2> /dev/null' ];
        end

        %%%%% amplitude
        fOut = replace(fIn,'_stats.nii.gz','_model-custom_respAmp.nii.gz');
        fSes_echoRms(R,E).fRespAmp = fOut;
        if forceThis || ~exist(fOut,'file')
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fIn '[' label '#0_Coef] 2> /dev/null'];
        end

        %%%%% baseline
        fOut = replace(fIn,'_stats.nii.gz','_model-custom_base.nii.gz');
        fSes_echoRms(R,E).fBase = fOut;
        if forceThis || ~exist(fOut,'file')
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            buck = num2str(1:nRun,'Run#%iPol#0_Coef,'); buck(end) = [];
            cmd{end+1} = [fIn '[' buck '] 2> /dev/null'];
            cmd{end+1} = '3dTstat -overwrite -mean \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fOut ' 2> /dev/null'];
        end
    end
end

%%% Run command
if length(cmd)>1
    if verbose
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
    else
        [status,cmdout] = system(strjoin(cmd,newline)); if status; dbstack; error(cmdout); error('x'); end
        % [status,cmdout] = system(strjoin(cmd(1:4),newline)); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    end
    disp(' done')
else
    disp(' already done, skipping')
end





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add baseline to estimated responses (for relaxation fitting) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
disp('Adding baseline to response amplitude')
forceThis = force;
% useFittedBaseline = 1;
cmd = {srcAfni};

%%% Each run
for R = 1:nRun
    for E = 1:nEcho+1
        %%%% Define file names
        if E>nEcho
            if nEcho==1; continue; end
            E = 1;
            fAmp = fRun_echoRms(R,E).fRespAmp;
            fBase = fRun_echoRms(R,E).fBase;
            fAmpOnBase = replace(fAmp,'_respAmp.nii.gz','_respAmpOnBase.nii.gz');
            fRun_echoRms(R,E).fRespAmpOnBase = fAmpOnBase;
            fStat = fRun_echoRms(R,E).fStat;
        else
            fAmp = fRun(R,E).fRespAmp;
            fBase = fRun(R,E).fBase;
            fAmpOnBase = replace(fAmp,'_respAmp.nii.gz','_respAmpOnBase.nii.gz');
            fRun(R,E).fRespAmpOnBase = fAmpOnBase;
            fStat = fRun(R,E).fStat;
        end

        %%% Average fitted baselines across runs
        if forceThis || ~exist(fBase,'file')
            cmd{end+1} = '3dTstat -overwrite \';
            cmd{end+1} = '-mean \';
            cmd{end+1} = ['-prefix ' fBase ' \'];
            buck = 'Run#1Pol#0_Coef';
            cmd{end+1} = [fStat '[' buck '] 2> /dev/null'];
        end
        
        %%% Add baseline to response
        if forceThis || ~exist(fAmpOnBase,'file')
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-prefix ' fAmpOnBase ' \'];
            cmd{end+1} = ['-a ' fBase ' \'];
            cmd{end+1} = ['-b ' fAmp ' \'];
            cmd{end+1} = '-expr ''a+b'' 2> /dev/null';
        end
    end
end

%%% Whole session
if nRun>1
    R = 1;
    for E = 1:nEcho+1
        %%% Define file names
        if E>nEcho
            E = 1;
            if nEcho==1; continue; end
            fAmp = fSes_echoRms(R,E).fRespAmp;
            fBase = fSes_echoRms(R,E).fBase;
            fAmpOnBase = replace(fAmp,'_respAmp.nii.gz','_respAmpOnBase.nii.gz');
            fSes_echoRms(R,E).fRespAmpOnBase = fAmpOnBase;
            fStat = fSes_echoRms(R,E).fStat;
        else
            fAmp = fSes(R,E).fRespAmp;
            fBase = fSes(R,E).fBase;
            fAmpOnBase = replace(fAmp,'_respAmp.nii.gz','_respAmpOnBase.nii.gz');
            fSes(R,E).fRespAmpOnBase = fAmpOnBase;
            fStat = fSes(R,E).fStat;
        end
        %%% Average fitted baselines across runs
        if forceThis || ~exist(fBase,'file')
            cmd{end+1} = '3dTstat -overwrite \';
            cmd{end+1} = '-mean \';
            cmd{end+1} = ['-prefix ' fBase ' \'];
            buck = num2str(1:size(fRun,1),'Run#%iPol#0_Coef,'); buck(end) = [];
            cmd{end+1} = [fStat '[' buck '] 2> /dev/null'];
        end

        %%% Add baseline to response
        if forceThis || ~exist(fAmpOnBase,'file')
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-prefix ' fAmpOnBase ' \'];
            cmd{end+1} = ['-a ' fBase ' \'];
            cmd{end+1} = ['-b ' fAmp ' \'];
            cmd{end+1} = '-expr ''a+b'' 2> /dev/null';
        end
    end
end

%%% run command
if length(cmd)>1
    if verbose
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
    else
        [status,cmdout] = system(strjoin(cmd,newline)); if status; dbstack; error(cmdout); error('x'); end
        % [status,cmdout] = system(strjoin(cmd(1:4),newline)); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    end
    disp(' done')
else
    disp(' already done, skipping')
end


%% Reconcile output with convention
forceThis = force;
disp('Catenating functional analyses')
cmd = {srcAfni};
fieldListIn =  {'fRespT' 'fRespAmp' 'fBase' 'fRespAmpOnBase' 'fStat' 'fMat'   };
fieldListOut = {'respT'  'respAmp'  'base'  'respAmpOnBase'  'stat'  'dsgnMat'};
for i = 1:length(fieldListIn)
    fieldIn = fieldListIn{i};
    fieldOut = fieldListOut{i};
    if isfield(fRun,fieldIn)
        files.(fieldOut).f = reshape({fRun.(fieldIn)},size(fRun));
    end
    if nEcho>1 && isfield(fRun_echoRms,fieldIn)
        files.(fieldOut).fEchoRms = reshape({fRun_echoRms.(fieldIn)},size(fRun_echoRms));
    end
    if nRun>1 && isfield(fSes,fieldIn)
        files.(fieldOut).fSes = reshape({fSes.(fieldIn)},size(fSes));
    end
    if nRun>1 && nEcho>1 && isfield(fSes_echoRms,fieldIn)
        files.(fieldOut).fSesEchoRms = reshape({fSes_echoRms.(fieldIn)},size(fSes_echoRms));
    end

    %%% catenate echo
    if nEcho>1
        if any(ismember(fieldOut,{'respT' 'respAmp'}))
            if isfield(fRun,fieldIn)
                fIn = reshape({fRun.(fieldIn)},size(fRun));
                files.(fieldOut).fEchoCat = cell([size(fIn,1) 1]);
                for R = 1:size(fIn,1)
                    fOut = fIn{R,1}; fOut = strsplit(fOut,'_'); fOut{contains(fOut,'echo-')} = 'echo-cat'; fOut = {strjoin(fOut,'_')};
                    files.(fieldOut).fEchoCat(R,1) = fOut;
                    if forceThis || ~exist(char(fOut),'file')
                        if any(ismember(fieldOut,{'respF'}))
                            cmd{end+1} = '3dbucket -overwrite \';
                        else
                            cmd{end+1} = '3dTcat -overwrite \';
                        end
                        cmd{end+1} = ['-prefix ' char(fOut) ' \'];
                        cmd{end+1} = [strjoin(fIn(R,:)) ' 2> /dev/null'];
                    end
                end
            end
            if isfield(fSes,fieldIn)
                fIn = {fSes.(fieldIn)};
                fOut = fIn{1}; fOut = strsplit(fOut,'_'); fOut{contains(fOut,'echo-')} = 'echo-cat'; fOut = {strjoin(fOut,'_')};
                files.(fieldOut).fSesEchoCat = fOut;
                if forceThis || ~exist(char(fOut),'file')
                    if any(ismember(fieldOut,{'respF'}))
                        cmd{end+1} = '3dbucket -overwrite \';
                    else
                        cmd{end+1} = '3dTcat -overwrite \';
                    end
                    cmd{end+1} = ['-prefix ' char(fOut) ' \'];
                    cmd{end+1} = [strjoin(fIn) ' 2> /dev/null'];
                end
            end
        end
    end

    %%% catenate runs
    if nRun>1
        if any(ismember(fieldOut,{'respT' 'respAmp'}))
            if isfield(fRun,fieldIn)
                fIn = reshape({fRun.(fieldIn)},size(fRun));
                files.(fieldOut).fCat = cell([1 size(fIn,2)]);
                for E = 1:size(fIn,2)
                    fOut = fIn{1,E}; fOut = strsplit(fOut,'_'); fOut{contains(fOut,'run-')} = 'run-cat'; fOut = {strjoin(fOut,'_')};
                    fOut = strsplit(char(fOut),filesep); fOut{end} = ['run-cat_' fOut{end}]; fOut = {strjoin(fOut,filesep)};
                    files.(fieldOut).fCat(1,E) = fOut;
                    if forceThis || ~exist(char(fOut),'file')
                        if any(ismember(fieldOut,{'respF'}))
                            cmd{end+1} = '3dbucket -overwrite \';
                        else
                            cmd{end+1} = '3dTcat -overwrite \';
                        end
                        cmd{end+1} = ['-prefix ' char(fOut) ' \'];
                        cmd{end+1} = [strjoin(fIn(:,E)) ' 2> /dev/null'];
                    end
                end
                if nEcho>1
                    fIn = reshape({fRun_echoRms.(fieldIn)},size(fRun_echoRms));
                    fOut = fIn{1}; fOut = strsplit(fOut,'_'); fOut{contains(fOut,'run-')} = 'run-cat'; fOut = {strjoin(fOut,'_')};
                    fOut = strsplit(char(fOut),filesep); fOut{end} = ['run-cat_' fOut{end}]; fOut = {strjoin(fOut,filesep)};
                    files.(fieldOut).fCatEchoRms = fOut;
                    if forceThis || ~exist(char(fOut),'file')
                        if any(ismember(fieldOut,{'respF'}))
                            cmd{end+1} = '3dbucket -overwrite \';
                        else
                            cmd{end+1} = '3dTcat -overwrite \';
                        end
                        cmd{end+1} = ['-prefix ' char(fOut) ' \'];
                        cmd{end+1} = [strjoin(fIn(:)) ' 2> /dev/null'];
                    end
                end
            end
        end
    end
end
% files.cmd = strjoin(cmd,newline);

if length(cmd)>1
    if verbose
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
    else
        [status,cmdout] = system(strjoin(cmd,newline)); if status; dbstack; error(cmdout); error('x'); end
        % [status,cmdout] = system(strjoin(cmd(1:4),newline)); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    end
    disp(' done')
else
    disp(' already done, skipping')
end






function [fRun,fSes,param] = runAfni(fList,param,fMask,force,verbose)
global srcAfni srcFs
if ~exist('fMask','var'); fMask = []; end
if ~exist('force','var'); force = []; end
if ~exist('verbose','var'); verbose = []; end
if isempty(force); force = 0; end
if isempty(verbose); verbose = 0; end


%% %%%%%%%%%%%%%%%%%%
% Functional design %
%%%%%%%%%%%%%%%%%% %%
k = param.funDsgn.k;
trStim = param.funDsgn.trStim;
durSeq =  param.funDsgn.durSeq;
condSeq = param.funDsgn.condSeq;
startSeq = cumsum(durSeq)-durSeq;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed-effect model on individual runs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
cmd = {srcAfni};
for E = 1:size(fList,2)
    for I = 1:size(fList,1)

        %% Define files
        fIn = fList{I,E};
        fOut = fIn;
        fStat = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_model-custom_stats.nii.gz');
        fMat = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_model-custom_stats.xmat.1D');
        fStim = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_model-custom_startTime.1D');
        fResp = [];
        fFit = [];
        fResid = [];

        
        fRun(I,E).fStat = fStat;
        fRun(I,E).fMat = fMat;


        %% Contruct afni command
        cmdTmp = {};
        if size(fList,2)>1
            cmdTmp{end+1} = ['echo ''  ''run' num2str(I) '/' num2str(size(fList,1)) ' -- echo' num2str(E) '/' num2str(size(fList,2))];
        else
            cmdTmp{end+1} = ['echo ''  ''run' num2str(I) '/' num2str(size(fList,1))];
        end
        if ~exist(fStat,'file') || force
            cmdTmp = [cmdTmp afniCmd(fIn,fMask,param,fStim,startSeq,condSeq,trStim,fResp,fFit,fResid,fMat,fStat,verbose)];

            cmdTmp{end+1} = ['echo ''   ''' fStat];
            cmdTmp{end+1} = ['echo ''   ''' fMat];
            cmdTmp{end+1} = 'echo ''   ''done';
        else
            cmdTmp{end+1} = ['echo ''   ''' fStat];
            cmdTmp{end+1} = ['echo ''   ''' fMat];
            cmdTmp{end+1} = 'echo ''   ''already done, skipping';
        end
        fRun(I,E).cmd =  strjoin(cmdTmp,newline);

        cmd = [cmd cmdTmp];
        % [status,cmdout] = system(strjoin(cmd,newline),'-echo');
    end
end


%% Run afni command
cmd = strjoin(cmd,newline); % disp(strjoin(cmd,newline))
[status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end

% Extract and plot design matrix
if verbose>1
    mri = MRIread(fIn,1); trMri = mri.tr/1000; nFrame = mri.nframes; clear mri
    
    dt = trStim;
    cmdX = {srcAfni};
    cmdX{end+1} = ['1dcat ' fMat];
    [status,cmdout] = system(strjoin(cmdX,newline));
    param.funDsgn.mat = str2num(cmdout);
    val = param.funDsgn.mat(:,end);
    t = 0:trMri:trMri*(size(param.funDsgn.mat,1)-1);
    figure('WindowStyle','docked');
    plot(t,val); hold on
    grid on
    grid minor

    startTimes = param.funDsgn.startSeq(logical(param.funDsgn.condSeq));
    stopTimes = param.funDsgn.startSeq(~logical(param.funDsgn.condSeq));
    xline(startTimes,'r')
    xline(stopTimes,'g')

    for i = 1:length(startTimes)
        plot(...
            param.funDsgn.hrf(:,1)+startTimes(i),...
            param.funDsgn.hrf(:,2),'m')
    end

end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed-effect model on concatenated runs (separate baselines) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
cmd = {srcAfni};
if size(fList,1)>1
    for E = 1:size(fList,2)

        %% Define files
        fIn = fList(:,E);
        tmp = strsplit(fIn{1},'_'); tmp = char(tmp(contains(tmp,'run-')));
        fOut = replace(fIn{1},tmp,'run-cat'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end

        fStat = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_model-custom_stats.nii.gz');
        fMat = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_model-custom_stats.xmat.1D');
        fStim = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_model-custom_startTime.1D');
        fResp = [];
        fFit = [];
        fResid = [];

        % fSes(1,E).fResp = fResp;
        % fSes(1,E).fFit = fFit;
        % fSes(1,E).fResid = fResid;
        fSes(1,E).fStat = fStat;
        fSes(1,E).fMat = fMat;


        %% Contruct afni command
        cmdTmp = {};
        if size(fList,2)>1
            cmdTmp{end+1} = ['echo ''  ''runCat/' num2str(size(fList,1)) ' -- echo' num2str(E) '/' num2str(size(fList,2))];
        else
            cmdTmp{end+1} = ['echo ''  ''runCat/' num2str(size(fList,1))];
        end
        if ~exist(fStat,'file') || force
            cmdTmp = [cmdTmp afniCmd(fIn,fMask,param,fStim,startSeq,condSeq,trStim,fResp,fFit,fResid,fMat,fStat,verbose)];

            % cmdTmp{end+1} = ['echo ''   ''' fResp];
            cmdTmp{end+1} = ['echo ''   ''' fStat];
            cmdTmp{end+1} = ['echo ''   ''' fMat];
            cmdTmp{end+1} = 'echo ''   ''done';
        else
            % cmdTmp{end+1} = ['echo ''   ''' fResp];
            cmdTmp{end+1} = ['echo ''   ''' fStat];
            cmdTmp{end+1} = ['echo ''   ''' fMat];
            cmdTmp{end+1} = 'echo ''   ''already done, skipping';
        end
        fSes(1,E).cmd =  strjoin(cmdTmp,newline);

        cmd = [cmd cmdTmp];
        % [status,cmdout] = system(strjoin(cmdTmp,newline),'-echo');
    end

    %% Run afni command
    cmd = strjoin(cmd,newline); % disp(strjoin(cmd,newline))
    if verbose
        [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    else
        [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    end


    % Extract and plot design matrix
    if verbose>1
        mri = MRIread(fIn{1},1); trMri = mri.tr/1000; nFrame = mri.nframes; clear mri

        dt = trStim;
        cmdX = {srcAfni};
        cmdX{end+1} = ['1dcat ' fMat];
        [status,cmdout] = system(strjoin(cmdX,newline));
        param.funDsgn.mat = str2num(cmdout);
        val = param.funDsgn.mat(:,end);
        t = 0:trMri:trMri*(size(param.funDsgn.mat,1)-1);
        figure('WindowStyle','docked');
        plot(t,val); hold on
        grid on
        grid minor

        startTimes = param.funDsgn.startSeq(logical(param.funDsgn.condSeq));
        stopTimes = param.funDsgn.startSeq(~logical(param.funDsgn.condSeq));
        xline(startTimes,'r')
        xline(stopTimes,'g')

        for i = 1:length(startTimes)
            plot(...
                param.funDsgn.hrf(:,1)+startTimes(i),...
                param.funDsgn.hrf(:,2),'m')
        end


        dt = trStim;
        cmdX = {srcAfni};
        cmdX{end+1} = ['1dcat ' fMat];
        [status,cmdout] = system(strjoin(cmdX,newline));
        param.funDsgn.mat = str2num(cmdout);
        figure('WindowStyle','docked');
        imagesc(param.funDsgn.mat); colormap gray
        ax = gca;
        ax.XTick = ( 1:size(param.funDsgn.mat(:,4:end),2) )+size(fIn,1)*3;
        ax.XTickLabel = num2str(((ax.XTick-size(fIn,1)*3)*dt)','%f');
        % clim([-0.001 0.001])
    end

else
    fSes = [];
end




function cmd = afniCmd(fIn,fMask,param,fStim,startSeq,condSeq,trStim,fResp,fFit,fResid,fMat,fStat,verbose)
cmd = {'3dDeconvolve -overwrite \'};
% cmdTmp{end+1} = ['-force_TR ' num2str(trStim) ' \'];
if iscell(fIn)
    cmd{end+1} = ['-input ' strjoin(fIn,' ') ' \'];
else
    cmd{end+1} = ['-input ' fIn ' \'];
end
if ~isempty(fMask)
    cmd{end+1} = ['-mask ' fMask ' \'];
end
cmd{end+1} = '-polort A \';
% if ~exist('trMri','var') || isempty(trMri) || ~exist('nFrame','var') || isempty(nFrame)
if iscell(fIn)
    mri = MRIread(fIn{1},1);
else
    mri = MRIread(fIn,1);
end
trMri = mri.tr/1000;
cmd{end+1} = ['-stim_times_subtract ' num2str(trMri*param.nDummy,'%f') ' \'];
cmd{end+1} = '-num_stimts 1 \';
cmd{end+1} = ['-stim_label 1 ' replace(param.label,'set-','') ' \'];

% write design to file
fido = fopen(fStim, 'w');
for i = 1:length(fIn)
    fprintf(fido,'%.3f ',startSeq(condSeq==1));
    fprintf(fido,'\n');
end
fclose(fido);

% Define response model
k = 1;
cmd{end+1} = ['-stim_times ' num2str(k) ' ' fStim ' \'];

t = param.funDsgn.hrf(:,1);
hrf = param.funDsgn.hrf(:,2);
dt = mode(diff(t));
b = t(1,1);
c = t(end,1) + dt*2;
expr = {};
for i = 1:length(hrf)
    expr{end+1} = [num2str(hrf(i),'%f') '*TENT(t-' num2str(t(i)+2*dt,'%f') ')'];
end
expr = ['''EXPR(' num2str(b) ',' num2str(c) ') ' strjoin(expr,'+') ''''];

cmd{end+1} = [expr ' \'];
cmd{end+1} = '-tout \';
cmd{end+1} = '-bout \';
cmd{end+1} = ['-x1D ' fMat ' \'];
if verbose>1
    cmd{end+1} = ['-bucket ' fStat];
else
    cmd{end+1} = ['-bucket ' fStat ' 2>/dev/null'];
end


