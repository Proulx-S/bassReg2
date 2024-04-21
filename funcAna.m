function [fRun,fSes,fSes_echoCat,param] = funcAna(fFunc,param,fAnat,fFmap,force)
global srcAfni srcFs
if ~exist('force','var'); force = []; end
if ~exist('fAnat','var'); fAnat = []; end
if ~exist('fFmap','var'); fFmap = []; end
if ~isfield(param,'verbose'); param.verbose = []; end
if isempty(force); force = 0; end
if isempty(param.verbose); param.verbose = 0; end; verbose = param.verbose;

meFlag = isfield(fFunc,'fEchoRms') && ~isempty(fFunc.fEchoRms);

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
disp(['Functional analysis: ' param.label])
[fRun,fSes,param] = runAfni(fFunc.f,param,fMask,force,verbose); % analysis performed on each echoe within that function
if isempty(fSes); fSes = fRun; end
if meFlag % relaunch the function to perform analysis on cross-echo rms
    disp(['Functional analysis (cross-echo rms): ' param.label])
    [fRun_echoRms,fSes_echoRms] = runAfni(fFunc.fEchoRms,param,fMask,force,verbose);
    if isempty(fSes_echoRms); fSes_echoRms = fRun_echoRms; end
    % fRun = cat(2,fRun,fRun_echoRms);
    % fSes = cat(2,fSes,fSes_echoRms);
end


%% Add baseline to estimated responses (for later visualization as movies)
forceThis = 0;
% useFittedBaseline = 1;
cmd = {srcAfni};
for E = 1:size(fSes,2)+1
    %%% Define file names
    if E>size(fSes,2)
        fResp = fSes_echoRms.fResp;
        fBase = replace(fResp,'_resp.nii.gz','_base.nii.gz');
        fRespOnBase = replace(fResp,'_resp.nii.gz','_respOnBase.nii.gz');
        fSes_echoRms.fBase = fBase;
        fSes_echoRms.fRespOnBase = fRespOnBase;
        % if useFittedBaseline
            % Use fitted baseline
            fStat = fSes_echoRms.fStat;
        % else
        %     % Use plain average as baseline
        %     copyfile(fFunc.fAvCatAvEchoRms,fBase);
        % end
    else
        fResp = fSes(1,E).fResp;
        fBase = replace(fResp,'_resp.nii.gz','_base.nii.gz');
        fRespOnBase = replace(fResp,'_resp.nii.gz','_respOnBase.nii.gz');
        fSes(1,E).fBase = fBase;
        fSes(1,E).fRespOnBase = fRespOnBase;
        % if useFittedBaseline
            % Use fitted baseline
            fStat = fSes(1,E).fStat;
        % else
        %     % Use plain average as baseline
        %     copyfile(fFunc.fAvCatAv{E},fBase);
        % end
    end
    
    % if useFittedBaseline
        %%% Average fitted baselines across runs
        cmd{end+1} = '3dTstat -overwrite \';
        cmd{end+1} = '-mean \';
        cmd{end+1} = ['-prefix ' fBase ' \'];
        buck = num2str(1:size(fRun,1),'Run#%iPol#0_Coef,'); buck(end) = [];
        cmd{end+1} = [fStat '[' buck ']'];
    % end

    %%% Add baseline to response
    cmd{end+1} = '3dcalc -overwrite \';
    cmd{end+1} = ['-prefix ' fRespOnBase ' \'];
    cmd{end+1} = ['-a ' fBase ' \'];
    cmd{end+1} = ['-b ' fResp ' \'];
    cmd{end+1} = '-expr ''a+b''';
end

%%% Run command
disp('Adding estimated baseline and responses for visualization as movies')
cmd = strjoin(cmd,newline); % disp(cmd)
if forceThis || ~exist(fRespOnBase,'file') || ~exist(fBase,'file')
    if verbose
        [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    else
        [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    end
    disp(' done')
else
    disp(' already done, skipping')
end

%%% Make the movies
forceThis = 0;
disp('Making movies')
nLoop = 4;
mask = MRIread(fFunc.manBrainMask); mask = logical(mask.vol);
for E = 1:size(fSes,2)+1
    disp([' ' num2str(E) '/' num2str(size(fSes,2)+1)])
    if E>size(fSes,2)
        fIn = fSes_echoRms.fRespOnBase;
        fOut = replace(fIn,'.nii.gz','');
        fSes_echoRms.fRespOnBaseMovie = [fOut '.avi'];
        fSes_echoRms.fRespOnBaseMovieHighBit = [fOut '.mj2'];
    else
        fIn = fSes(E).fRespOnBase;
        fOut = replace(fIn,'.nii.gz','');
        fSes(E).fRespOnBaseMovie = [fOut '.avi'];
        fSes(E).fRespOnBaseMovieHighBit = [fOut '.mj2'];
    end
    if forceThis || ~exist([fOut '.avi'],'file')
        vOut = VideoWriter(fOut,'Uncompressed AVI');
    end
    if forceThis || ~exist([fOut '.mj2'],'file')
        vOutHighBit = VideoWriter(fOut,'Archival');
    end

    if forceThis || ~exist([fOut '.avi'],'file') || ~exist([fOut '.mj2'],'file')
        resp = MRIread(fIn);
        % Scale
        mask = repmat(mask,[1 1 resp.depth resp.nframes]);
        resp.vol(mask) = resp.vol(mask) - min(resp.vol(mask));
        resp.vol(mask) = resp.vol(mask) ./ max(resp.vol(mask));
        mask = mask(:,:,1,1);
        % Crop
        resp.vol(all(~mask,2),:,:,:) = [];
        resp.vol(:,all(~mask,1),:,:) = [];
        % Upsample
        resp.vol = imresize(resp.vol,3,'nearest');
        % Set frame rate to 1cycle/sec
        vOutHighBit.FrameRate = resp.nframes;
        vOut.FrameRate = resp.nframes;

        % Write
        open(vOut)
        open(vOutHighBit)
        for L = 1:nLoop
            for f = 1:resp.nframes
                writeVideo(vOut,resp.vol(:,:,1,f));
                writeVideo(vOutHighBit,uint16(resp.vol(:,:,1,f)*(2^16-1)));
            end
        end
        close(vOut)
        close(vOutHighBit)

        disp('  done')
    else
        disp('  already done,skipping')
    end
end





%% Generate visualization commands
if ~isempty(fAnat) && ~isempty(fAnat{1})
    fT1w = fAnat{contains(fAnat,'_T1w.nii.gz')};

    ind = contains(fAnat,'echo-rms_satIndex.nii.gz');
    if any(ind); fSatinIndex_echoRms = fAnat{ind}; else; fSatinIndex_echoRms = []; end
    ind = contains(fAnat,'echo-rms_satIndexNumOnBack.nii.gz');
    if any(ind); fSatinIndexPlus_echoRms = fAnat{ind}; else; fSatinIndexPlus_echoRms = []; end
    ind = contains(fAnat,'echo-rms_satIndexNumAbs.nii.gz');
    if any(ind); fSatinIndexAbs_echoRms = fAnat{ind}; else; fSatinIndexAbs_echoRms = []; end
    ind = contains(fAnat,'echo-cat_satIndex.nii.gz');
    if any(ind); fSatinIndex_echoCat = fAnat{ind}; else; fSatinIndex_echoCat = []; end
    ind = contains(fAnat,'echo-cat_satIndexNumOnBack.nii.gz');
    if any(ind); fSatinIndexPlus_echoCat = fAnat{ind}; else; fSatinIndexPlus_echoCat = []; end
    ind = contains(fAnat,'echo-cat_satIndexNumAbs.nii.gz');
    if any(ind); fSatinIndexAbs_echoCat = fAnat{ind}; else; fSatinIndexAbs_echoCat = []; end
else
    fT1w = [];
    fSatinIndex_echoRms = [];
    fSatinIndexPlus_echoRms = [];
    fSatinIndexAbs_echoRms = [];
    fSatinIndex_echoCat = [];
    fSatinIndexPlus_echoCat = [];
    fSatinIndexAbs_echoCat = [];
end
if ~isempty(fFmap) && ~isempty(fFmap{1})
    ind = contains(fFmap,'rec-FA_TB1SRGE.nii.gz');
    if any(ind); fB1 = fFmap{ind}; else; fB1 = []; end
    % fB1 = fFmap{contains(fFmap,'rec-FA_TB1SRGE.nii.gz')};
else
    fB1 = [];
end

%%% individual runs
for I = 1:size(fRun,1)
    E = 1; % multi-echo visualization not implemented
    fRun(I,E).fUnder = fFunc.fAv{I,E};

    cmd = {srcAfni};
    if isfield(param,'layout') && ~isempty(param.layout) && exist(param.layout,'file')
        cmd{end+1} = ['afni -layout ' param.layout ' \'];
    else
        cmd{end+1} = 'afni \';
    end
    cmd{end+1} = ['-tbar run' num2str(I) 'echo ' num2str(E) ' \'];
    cmd{end+1} = [fRun(I,E).fUnder ' \'];
    cmd{end+1} = [fRun(I,E).fStat ' \'];
    cmd{end+1} = [fRun(I,E).fResp ' &'];
    cmd = strjoin(cmd,newline); % disp(cmd)
    fRun(I,E).cmdVisAfni = cmd;

    if verbose; disp(cmd); end
    if verbose>1; [status,cmdout] = system(cmd); end

    % https://afni.nimh.nih.gov/pub/dist/doc/program_help/README.driver.html
    % https://afni.nimh.nih.gov/pub/dist/doc/program_help/plugout_drive.html
    % afni -yesplugouts
    % plugout_drive  -com 'SWITCH_SESSION A.afni'                       \
    % -com 'OPEN_WINDOW A.axialimage geom=600x600+416+44 \
    % ifrac=0.8 opacity=9'                         \
    % -com 'OPEN_WINDOW A.sagittalimage geom=+45+430     \
    % ifrac=0.8 opacity=9'                         \
    % -com 'SWITCH_UNDERLAY anat'                        \
    % -com 'SWITCH_OVERLAY strip'                        \
    % -com 'SEE_OVERLAY +'                               \
    % -com 'SET_DICOM_XYZ 7 12 2'                        \
    % -com 'OPEN_WINDOW A.axialimage keypress=v'         \
    % -quit
    %
    % SET_SUBBRICKS



    % fT1w = dir(fullfile(bidsDir,'anat','*proc-RMS_T1w.nii.gz')); fT1w = fullfile(fT1w.folder,fT1w.name);
    % fB1 = dir(fullfile(bidsDir,'fmap','*rec-FA_TB1SRGE.nii.gz')); fB1 = fullfile(fB1.folder,fB1.name);
    cmd = {srcFs};
    cmd{end+1} = 'freeview \';
    cmd{end+1} = [fRun(I).fUnder ' \'];
    if ~isempty(fT1w)
        cmd{end+1} = [fT1w ':resample=cubic \'];
    end
    if ~isempty(fB1)
        cmd{end+1} = [fB1 ':resample=cubic:visible=0 \'];
    end
    cmd{end+1} = [fRun(I).fResp ':resample=cubic:visible=0 \'];
    if ~isempty(fSatinIndexPlus_echoRms)
        cmd{end+1} = [fSatinIndexPlus_echoRms ':resample=cubic:visible=0 \'];
    end
    thresh = abs(norminv(0.95));
    if ~isempty(fSatinIndex_echoRms)
        cmd{end+1} = [fSatinIndex_echoRms ':colormap=heat:resample=cubic:visible=0 \'];
    end
    cmd{end+1} = [fRun(I).fStat ':colormap=heat:resample=cubic:visible=0 &'];
    cmd = strjoin(cmd,newline);

    fRun(I).cmdVisFs = cmd;
end


%%% full session
if isempty(fSes)
    fSes = fRun;
end
if meFlag
    if isempty(fSes_echoRms)
        fSes_echoRms = fRun_echoRms;
    end
    %%%% Multi echo
    %%%%% catenate echo
    fSes_echoCat = fSes(:,1);
    fSes_echoCat.fResp = {fSes.fResp};
    fSes_echoCat.fBase = {fSes.fBase};
    fSes_echoCat.fRespOnBase = {fSes.fRespOnBase};
    fSes_echoCat.fFit = {fSes.fFit};
    fSes_echoCat.fResid = {fSes.fResid};

    fOut = strsplit(fSes(1).fStat,'_'); fOut{contains(fOut,'echo-')} = 'echo-cat';
    fOut = replace(strjoin(fOut,'_'),'_stats.nii.gz','_fullF.nii.gz'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
    
    cmd = {srcAfni};
    cmd{end+1} = '3dbucket -overwrite \';
    cmd{end+1} = ['-prefix ' fOut ' \'];
    cmd{end+1} = [strjoin({fSes.fStat},'[0] ') '[0]'];
    % cmd = {srcAfni};
    % cmd{end+1} = '3dTcat -overwrite \';
    % cmd{end+1} = ['-prefix ' fOut ' \'];
    % cmd{end+1} = [strjoin({fSes.fStat},'[0] ') '[0]'];
    cmd = strjoin(cmd,newline);
    [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    fSes_echoCat.cmd = cmd;
    fSes_echoCat.fStat = fOut;
    
    cadidateField = {'fAvCatAvEchoCat' 'fAvEchoCat' 'fAvCatAv' 'fAv'};
    cadidateField = cadidateField(ismember(cadidateField,fields(fFunc)));
    fSes_echoCat.fUnder = fFunc.(cadidateField{1});

    %%%%% catenate cross-echo rms
    fSes_echoCat.rms.fResp = fSes_echoRms.fResp;
    fSes_echoCat.rms.fBase = fSes_echoRms.fBase;
    fSes_echoCat.rms.fRespOnBase = fSes_echoRms.fRespOnBase;
    fSes_echoCat.rms.fRespOnBaseMovie = fSes_echoRms.fRespOnBaseMovie;
    fSes_echoCat.rms.fRespOnBaseMovieHighBit = fSes_echoRms.fRespOnBaseMovieHighBit;
    fSes_echoCat.rms.fFit = fSes_echoRms.fFit;
    fSes_echoCat.rms.fResid = fSes_echoRms.fResid;
    fSes_echoCat.rms.fStat = fSes_echoRms.fStat;

    cadidateField = {'fAvCatAvEchoRms' 'fAvEchoRms' 'fAvCatAv' 'fAv'};
    cadidateField = cadidateField(ismember(cadidateField,fields(fFunc)));
    fSes_echoCat.rms.fUnder = fFunc.(cadidateField{1});





    % %%%%% visulaize with afni
    % cmd = {srcAfni};
    % if isfield(param,'layout') && ~isempty(param.layout) && exist(param.layout,'file')
    %     cmd{end+1} = ['afni -layout ' param.layout ' \'];
    % else
    %     cmd{end+1} = 'afni \';
    % end
    % cmd{end+1} = ['-tbar ' num2str(size(fRun,1)) 'runs \'];
    % cmd{end+1} = [fSes_echoCat.fUnder ' \'];
    % cmd{end+1} = [fSes_echoCat.rms.fUnder ' \'];
    % cmd{end+1} = [fSes_echoCat.fStat ' \'];
    % cmd{end+1} = [fSes_echoCat.rms.fStat ' \'];
    % cmd{end+1} = [strjoin([fSes_echoCat.fResp fSes_echoCat.rms.fResp],' ') ' &'];
    % cmd = strjoin(cmd,newline); % disp(cmd)
    % fSes_echoCat.cmdVisAfni = cmd;
    % 
    % if verbose; disp(cmd); end 
    % if verbose>1; [status,cmdout] = system(cmd); end

    %%%%% visulaize with freeview
    cmd = {srcFs};
    cmd{end+1} = 'freeview \';
    cmd{end+1} = '-timecourse \';
    % cmd{end+1} = '-subtitle \';
    cmd{end+1} = [char(fSes_echoCat.fUnder) ':name=echo-cat_funcAv:visible=0 \'];
    cmd{end+1} = [char(fSes_echoCat.rms.fUnder) ':name=echo-rms_funcAv:visible=1 \'];
    if ~isempty(fT1w)
        cmd{end+1} = [fT1w ':resample=cubic:name=T1w:visible=0 \'];
    end
    if ~isempty(fB1)
        cmd{end+1} = [fB1 ':colormap=turbo:resample=cubic:name=B1map:colorscale=50,130:visible=0 \'];
    end
    for E = 1:length(fSes_echoCat.fResp)
        tmp = strsplit(fSes_echoCat.fResp{E},'_');
        cmd{end+1} = [fSes_echoCat.fResp{E} ':name=' tmp{contains(tmp,'echo-')} '_resp:visible=0 \'];
        cmd{end+1} = [fSes_echoCat.fBase{E} ':name=' tmp{contains(tmp,'echo-')} '_base:visible=0 \'];
        cmd{end+1} = [fSes_echoCat.fRespOnBase{E} ':name=' tmp{contains(tmp,'echo-')} '_respOnBase:visible=0 \'];
    end
    cmd{end+1} = [fSes_echoCat.rms.fResp ':name=echo-rms_resp:visible=0 \'];
    cmd{end+1} = [fSes_echoCat.rms.fBase ':name=echo-rms_base:visible=0 \'];
    cmd{end+1} = [fSes_echoCat.rms.fRespOnBase ':name=echo-rms_respOnBase:visible=0 \'];
    if ~isempty(fSatinIndexPlus_echoRms)
        cmd{end+1} = [fSatinIndexPlus_echoRms ':resample=cubic:name=echo-rms_satIndexNumPlusBack:visible=0 \'];
    end
    if ~isempty(fSatinIndexPlus_echoCat)
        cmd{end+1} = [fSatinIndexPlus_echoCat ':resample=cubic:name=echo-cat_satIndexNumPlusBack:visible=0 \'];
    end
    if ~isempty(fSatinIndexAbs_echoRms)
        cmd{end+1} = [fSatinIndexAbs_echoRms ':resample=cubic:name=echo-rms_satIndexNumAbs:visible=0 \'];
    end
    if ~isempty(fSatinIndexAbs_echoCat)
        cmd{end+1} = [fSatinIndexAbs_echoCat ':resample=cubic:name=echo-cat_satIndexNumAbs:visible=0 \'];
    end
    if ~isempty(fSatinIndex_echoRms)
        cmd{end+1} = [fSatinIndex_echoRms ':colormap=heat:resample=cubic:name=echo-rms_satIndex:visible=0 \'];
    end
    if ~isempty(fSatinIndex_echoCat)
        cmd{end+1} = [fSatinIndex_echoCat ':colormap=heat:resample=cubic:name=echo-cat_satIndex:visible=0 \'];
    end

    for E = 1:size(fSes,2)
        cmdTmp = {srcAfni};
        cmdTmp{end+1} = ['fdrval -qinput ' fSes(E).fStat ' 0 0.05'];
        [~,tmp] = system(strjoin(cmdTmp,newline));
        tmp = strsplit(tmp,newline);
        thresh(E) = str2num(tmp{end-1});

        cmdTmp = {srcAfni};
        cmdTmp{end+1} = ['3dBrickStat ' fSes(E).fStat];
        [~,tmp] = system(strjoin(cmdTmp,newline));
        tmp = strsplit(tmp,newline);        
        maxF(E) = str2num(tmp{end-1});
        if maxF(E)>20; maxF(E) = 20; end
        cmd{end+1} = [fSes(E).fStat ':colormap=heat:heatscale=' num2str(thresh(E)) ',' num2str(maxF(E)) ':name=echo-' num2str(E) '_fullF:visible=0 \'];
    end
    thresh = mean(thresh);
    maxF = min(maxF); if maxF>20; maxF = 20; end
    cmd{end+1} = [fSes_echoCat.fStat ':colormap=heat:heatscale=' num2str(thresh) ',' num2str(maxF) ':name=echo-cat_fullF:visible=0 \'];

    cmdTmp = {srcAfni};
    cmdTmp{end+1} = ['fdrval -qinput ' fSes_echoCat.rms.fStat ' 0 0.05'];
    [~,tmp] = system(strjoin(cmdTmp,newline));
    tmp = strsplit(tmp,newline);
    thresh = str2num(tmp{end-1});

    cmdTmp = {srcAfni};
    cmdTmp{end+1} = ['3dBrickStat ' fSes_echoCat.rms.fStat];
    [~,tmp] = system(strjoin(cmdTmp,newline));
    tmp = strsplit(tmp,newline);
    maxF = str2num(tmp{end-1});
    maxF = min(maxF); if maxF>20; maxF = 20; end
    cmd{end+1} = [fSes_echoCat.rms.fStat ':colormap=heat:heatscale=' num2str(mean(thresh)) ',' num2str(maxF) ':name=echo-rms_fullF &'];
    
    cmd = strjoin(cmd,newline);
    fSes_echoCat.cmdVisFs = cmd;
else
    %%%% single-echo
    fSes_echoCat = [];
    E = 1;

    %%%%% visualize with afni
    fSes(1,E).fUnder = fFunc.fAvCatAv{1,E};
    cmd = {srcAfni};
    if isfield(param,'layout') && ~isempty(param.layout) && exist(param.layout,'file')
        cmd{end+1} = ['afni -layout ' param.layout ' \'];
    else
        cmd{end+1} = 'afni \';
    end
    cmd{end+1} = ['-tbar ' num2str(size(fRun,1)) 'runs \'];
    cmd{end+1} = [fSes(1,E).fUnder ' \'];
    cmd{end+1} = [fSes(1,E).fStat ' \'];
    cmd{end+1} = [fSes(1,E).fResp ' &'];
    cmd = strjoin(cmd,newline); %disp(cmd)
    fSes(1,E).cmdVisAfni = cmd;

    if verbose; disp(cmd); end
    if verbose>1; [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end; end

    %%%%% visualize with freeview
    cmd = {srcFs};
    cmd{end+1} = 'freeview \';
    cmd{end+1} = [fSes(1,E).fUnder ' \'];
    if ~isempty(fT1w)
        cmd{end+1} = [fT1w ':resample=cubic \'];
    end
    if ~isempty(fB1)
        cmd{end+1} = [fB1 ':resample=cubic:visible=0 \'];
    end
    cmd{end+1} = [fSes(1,E).fResp ':resample=cubic:visible=0 \'];
    if ~isempty(fSatinIndexPlus_echoRms)
        cmd{end+1} = [fSatinIndexPlus_echoRms ':resample=cubic:visible=0 \'];
    end
    thresh = abs(norminv(0.95));
    if ~isempty(fSatinIndex_echoRms)
        cmd{end+1} = [fSatinIndex_echoRms ':colormap=heat:resample=cubic:visible=0 \'];
    end
    cmd{end+1} = [fSes(1,E).fStat ':colormap=heat:resample=cubic &'];
    cmd = strjoin(cmd,newline);

    fSes(1,E).cmdVisFs = cmd;
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
        fStat = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_stats.nii.gz');
        fMat = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_stats.xmat.1D');
        fStim = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_startTime.1D');
        fResp = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_resp.nii.gz');
        fFit = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_fit.nii.gz');
        fResid = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_resid.nii.gz');

        fRun(I,E).fResp = fResp;
        fRun(I,E).fFit = fFit;
        fRun(I,E).fResid = fResid;
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

            cmdTmp{end+1} = ['echo ''   ''' fResp];
            cmdTmp{end+1} = ['echo ''   ''' fStat];
            cmdTmp{end+1} = ['echo ''   ''' fMat];
            cmdTmp{end+1} = 'echo ''   ''done';
        else
            cmdTmp{end+1} = ['echo ''   ''' fResp];
            cmdTmp{end+1} = ['echo ''   ''' fStat];
            cmdTmp{end+1} = ['echo ''   ''' fMat];
            cmdTmp{end+1} = 'echo ''   ''already done, skipping';
        end
        fRun(I,E).cmd =  strjoin(cmdTmp,newline);
        
        cmd = [cmd cmdTmp];
        % [status,cmdout] = system(strjoin(cmdTmp,newline),'-echo');
    end
end

%% Extract and plot design matrix
if verbose>1
    mri = MRIread(fIn,1); trMri = mri.tr; nFrame = mri.nframes; clear mri
    dt = trStim;
    cmdX = cmdTmp;
    cmdX{contains(cmdTmp,'-input ')} = ['-nodata ' num2str(nFrame) ' ' num2str(trMri,'%f') ' \'];
    cmdX = [cmdX(1:end-2) {'-x1D_stop \'} cmdX(end-1:end)];
    cmdX{end} = [cmdX{end}];
    [status,cmdout] = system(strjoin(cmdX,newline));
    cmdX = {srcAfni};
    cmdX{end+1} = ['1dcat ' fMat];
    [status,cmdout] = system(strjoin(cmdX,newline));
    param.funDsgn.mat = str2num(cmdout);
    figure('WindowStyle','docked');
    imagesc(param.funDsgn.mat(:,4:end)); colormap gray
    ax = gca;
    ax.XTick = 1:size(param.funDsgn.mat(:,4:end),2);
    ax.XTickLabel = num2str(((ax.XTick)*dt)','%f');
    clim([-0.001 0.001])
end

%% Run afni command
cmd = strjoin(cmd,newline); % disp(strjoin(cmd,newline))
[status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end

% Extract and plot design matrix
if verbose>1
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
    clim([-0.001 0.001])
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
        
        fStat = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_stats.nii.gz');
        fMat = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_stats.xmat.1D');
        fStim = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_startTime.1D');
        fResp = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_resp.nii.gz');
        fFit = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_fit.nii.gz');
        fResid = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_resid.nii.gz');

        fSes(1,E).fResp = fResp;
        fSes(1,E).fFit = fFit;
        fSes(1,E).fResid = fResid;
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

            cmdTmp{end+1} = ['echo ''   ''' fResp];
            cmdTmp{end+1} = ['echo ''   ''' fStat];
            cmdTmp{end+1} = ['echo ''   ''' fMat];
            cmdTmp{end+1} = 'echo ''   ''done';
        else
            cmdTmp{end+1} = ['echo ''   ''' fResp];
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
    [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end


    % Extract and plot design matrix
    if verbose>1
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
        clim([-0.001 0.001])
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
nFrame = mri.nframes;
% end
% trMri = 3;
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
deconWin = mode(diff(startSeq(condSeq==1)));
dt = trStim; %trStim;
deconWin = floor(deconWin/dt)*dt;
b = 0;
c = deconWin-dt;
n = (c-b)/dt + 1;
% (c-b)/(n-1)
k = 1;
cmd{end+1} = ['-stim_times ' num2str(k) ' ' fStim ' ''TENTzero(' num2str(b) ',' num2str(c) ',' num2str(n) ')'' \'];
cmd{end+1} = ['-TR_times ' num2str(dt,'%f') ' \'];
cmd{end+1} = ['-iresp ' num2str(k) ' ' fResp ' \'];
cmd{end+1} = ['-fitts ' fFit ' \'];
cmd{end+1} = ['-errts ' fResid ' \'];
cmd{end+1} = '-bout \';
% cmdTmp{end+1} = ['-TR_times ' num2str(trStim) ' \'];
cmd{end+1} = ['-x1D ' fMat ' \'];
if verbose>1
    cmd{end+1} = ['-bucket ' fStat];
else
    cmd{end+1} = ['-bucket ' fStat ' 2>/dev/null'];
end