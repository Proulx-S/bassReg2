function [fRun,fSes,fSes_echoCat,param] = funcAna(fFunc,param,fAnat,fFmap,force)
global srcAfni srcFs
verbose = 0;
if ~exist('force','var'); force = []; end
if ~exist('fAnat','var'); fAnat = []; end
if ~exist('fFmap','var'); fFmap = []; end
if ~isfield(param,'verbose'); param.verbose = []; end
if isempty(force); force = 0; end
if isempty(param.verbose); param.verbose = 0; end

meFlag = isfield(fFunc,'fEchoRmsList') && ~isempty(fFunc.fEchoRmsList);

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
[fRun,fSes,param] = runAfni(fFunc.fList,param,fMask,force); % analysis performed on each echoe within that function
if meFlag % relaunch the function to perform analysis on cross-echo rms
    [fRun_echoRms,fSes_echoRms] = runAfni(fFunc.fEchoRmsList,param,fMask,force);
    % fRun = cat(2,fRun,fRun_echoRms);
    % fSes = cat(2,fSes,fSes_echoRms);
end

% tmp = MRIread(fSes(1).fFit);
% plot(squeeze(mean(mean(tmp.vol,1),2)))

%% Generate visualization commands
if ~isempty(fAnat)
    fT1w = fAnat{contains(fAnat,'proc-RMS_T1w.nii.gz')};
    fSatinIndex_echoRms = fAnat{contains(fAnat,'echo-rms_satIndex.nii.gz')};
    fSatinIndexPlus_echoRms = fAnat{contains(fAnat,'echo-rms_satIndexNumOnBack.nii.gz')};
    fSatinIndexAbs_echoRms = fAnat{contains(fAnat,'echo-rms_satIndexNumAbs.nii.gz')};
    fSatinIndex_echoCat = fAnat{contains(fAnat,'echo-cat_satIndex.nii.gz')};
    fSatinIndexPlus_echoCat = fAnat{contains(fAnat,'echo-cat_satIndexNumOnBack.nii.gz')};
    fSatinIndexAbs_echoCat = fAnat{contains(fAnat,'echo-cat_satIndexNumAbs.nii.gz')};
else
    fT1w = [];
    fSatinIndex_echoRms = [];
    fSatinIndexPlus_echoRms = [];
    fSatinIndexAbs_echoRms = [];
    fSatinIndex_echoCat = [];
    fSatinIndexPlus_echoCat = [];
    fSatinIndexAbs_echoCat = [];
end
if ~isempty(fFmap)
    fB1 = fFmap{contains(fFmap,'rec-FA_TB1SRGE.nii.gz')};
else
    fB1 = [];
end

%%% individual runs
for I = 1:size(fRun,1)
    E = 1; % multi-echo visualization not implemented
    fRun(I,E).fUnder = fFunc.fAvList{I,E};

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
    cmd = strjoin(cmd,newline); disp(cmd)
    fRun(I,E).cmdVisAfni = cmd;

    if verbose
        [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    end

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
        cmd{end+1} = [fB1 ':resample=cubic \'];
    end
    cmd{end+1} = [fRun(I).fResp ':resample=cubic \'];
    if ~isempty(fSatinIndexPlus_echoRms)
        cmd{end+1} = [fSatinIndexPlus_echoRms ':resample=cubic \'];
    end
    thresh = abs(norminv(0.95));
    if ~isempty(fSatinIndex_echoRms)
        cmd{end+1} = [fSatinIndex_echoRms ':colormap=heat:resample=cubic \'];
    end
    cmd{end+1} = [fRun(I).fStat ':colormap=heat:resample=cubic &'];
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
    fSes_echoCat.fFit = {fSes.fFit};
    fSes_echoCat.fResid = {fSes.fResid};

    fOut = strsplit(fSes(1).fStat,'_'); fOut{contains(fOut,'echo-')} = 'echo-cat'; fOut = replace(strjoin(fOut,'_'),'_stats.nii.gz','_fullF.nii.gz'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
    cmd = {srcAfni};
    cmd{end+1} = '3dTcat -overwrite \';
    cmd{end+1} = ['-prefix ' fOut ' \'];
    cmd{end+1} = [strjoin({fSes.fStat},'[0] ') '[0]'];
    cmd = strjoin(cmd,newline);
    [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    fSes_echoCat.cmd = cmd;

    fSes_echoCat.fStat = fOut;
    fSes_echoCat.fUnder = fFunc.fAvCatAvEchoCat;

    %%%%% catenate cross-echo rms
    fSes_echoCat.rms.fResp = fSes_echoRms.fResp;
    fSes_echoCat.rms.fFit = fSes_echoRms.fFit;
    fSes_echoCat.rms.fResid = fSes_echoRms.fResid;
    fSes_echoCat.rms.fStat = fSes_echoRms.fStat;

    fSes_echoCat.rms.fUnder = fFunc.fAvCatAvEchoRms;


    %%%%% visulaize with afni
    cmd = {srcAfni};
    if isfield(param,'layout') && ~isempty(param.layout) && exist(param.layout,'file')
        cmd{end+1} = ['afni -layout ' param.layout ' \'];
    else
        cmd{end+1} = 'afni \';
    end
    cmd{end+1} = ['-tbar ' num2str(size(fRun,1)) 'runs \'];
    cmd{end+1} = [fSes_echoCat.fUnder ' \'];
    cmd{end+1} = [fSes_echoCat.rms.fUnder ' \'];
    cmd{end+1} = [fSes_echoCat.fStat ' \'];
    cmd{end+1} = [fSes_echoCat.rms.fStat ' \'];
    cmd{end+1} = [strjoin([fSes_echoCat.fResp fSes_echoCat.rms.fResp],' ') ' &'];
    cmd = strjoin(cmd,newline); disp(cmd)
    fSes_echoCat.cmdVisAfni = cmd;

    if verbose
        [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    end

    %%%%% visulaize with freeview
    cmd = {srcFs};
    cmd{end+1} = 'freeview \';
    cmd{end+1} = '-timecourse \';
    % cmd{end+1} = '-subtitle \';
    cmd{end+1} = [fSes_echoCat.fUnder ':name=echo-cat_funcAv \'];
    cmd{end+1} = [fSes_echoCat.rms.fUnder ':name=echo-rms_funcAv \'];
    if ~isempty(fT1w)
        cmd{end+1} = [fT1w ':resample=cubic:name=T1w \'];
    end
    if ~isempty(fB1)
        cmd{end+1} = [fB1 ':colormap=jet:resample=cubic:name=B1map \'];
    end
    for E = 1:length(fSes_echoCat.fResp)
        tmp = strsplit(fSes_echoCat.fResp{E},'_');
        cmd{end+1} = [fSes_echoCat.fResp{E} ':name=' tmp{contains(tmp,'echo-')} '_resp:visible=0 \'];
    end
    cmd{end+1} = [fSes_echoCat.rms.fResp ':name=echo-rms_resp:visible=0 \'];
    if ~isempty(fSatinIndexPlus_echoRms)
        cmd{end+1} = [fSatinIndexPlus_echoRms ':resample=cubic:name=echo-rms_satIndexNumPlusBack \'];
    end
    if ~isempty(fSatinIndexPlus_echoCat)
        cmd{end+1} = [fSatinIndexPlus_echoCat ':resample=cubic:name=echo-cat_satIndexNumPlusBack \'];
    end
    if ~isempty(fSatinIndexAbs_echoRms)
        cmd{end+1} = [fSatinIndexAbs_echoRms ':resample=cubic:name=echo-rms_satIndexNumAbs \'];
    end
    if ~isempty(fSatinIndexAbs_echoCat)
        cmd{end+1} = [fSatinIndexAbs_echoCat ':resample=cubic:name=echo-cat_satIndexNumAbs \'];
    end
    thresh = abs(norminv(0.95));
    if ~isempty(fSatinIndex_echoRms)
        cmd{end+1} = [fSatinIndex_echoRms ':colormap=heat:resample=cubic:name=echo-rms_satIndex \'];
    end
    if ~isempty(fSatinIndex_echoCat)
        cmd{end+1} = [fSatinIndex_echoCat ':colormap=heat:resample=cubic:name=echo-cat_satIndex \'];
    end
    cmd{end+1} = [fSes_echoCat.fStat ':colormap=heat:name=echo-cat_fullF \'];
    cmd{end+1} = [fSes_echoCat.rms.fStat ':colormap=heat:name=echo-rms_fullF &'];
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
    cmd = strjoin(cmd,newline); disp(cmd)
    fSes(1,E).cmdVisAfni = cmd;

    if verbose
        [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    end

    %%%%% visualize with freeview
    cmd = {srcFs};
    cmd{end+1} = 'freeview \';
    cmd{end+1} = [fSes(1,E).fUnder ' \'];
    if ~isempty(fT1w)
        cmd{end+1} = [fT1w ':resample=cubic \'];
    end
    if ~isempty(fB1)
        cmd{end+1} = [fB1 ':resample=cubic \'];
    end
    cmd{end+1} = [fSes(1,E).fResp ':resample=cubic \'];
    if ~isempty(fSatinIndexPlus_echoRms)
        cmd{end+1} = [fSatinIndexPlus_echoRms ':resample=cubic \'];
    end
    thresh = abs(norminv(0.95));
    if ~isempty(fSatinIndex_echoRms)
        cmd{end+1} = [fSatinIndex_echoRms ':colormap=heat:resample=cubic \'];
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

k = param.funDsgn.k;
trStim = param.funDsgn.trStim;
durSeq =  param.funDsgn.durSeq;
condSeq = param.funDsgn.condSeq;
% if isfield(param.funDsgn,'trMri')
%     trMri = param.funDsgn.trMri; else; trMri = []; end
startSeq = cumsum(durSeq)-durSeq;
% deconWin = mode(diff(startSeq(condSeq==1)));
% dt = trStim;
% b = 0;
% c = deconWin-dt*1;
% n = (c-b)/dt + 1;


%% fixed-effect model on individual runs
cmd = {srcAfni};
for E = 1:size(fList,2)
    for I = 1:size(fList,1)
        cmdTmp = {srcAfni};
        fIn = fList{I,E};
        fOut = fIn;
        fStat = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_stats.nii.gz');
        fStim = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_startTime.1D');
        fResp = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_resp.nii.gz');
        fFit = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_fit.nii.gz');
        fResid = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_resid.nii.gz');

        cmdTmp{end+1} = '3dDeconvolve -overwrite \';
        % cmdTmp{end+1} = ['-force_TR ' num2str(trStim) ' \'];
        cmdTmp{end+1} = ['-input ' fIn ' \'];
        if ~isempty(fMask)
            cmdTmp{end+1} = ['-mask ' fMask ' \'];
        end
        cmdTmp{end+1} = '-polort A \';
        if ~exist('trMri','var') || isempty(trMri)
            trMri = MRIread(fIn,1); trMri = trMri.tr/1000;
        end
        % trMri = 3;
        cmdTmp{end+1} = ['-stim_times_subtract ' num2str(trMri*param.nDummy,'%f') ' \'];
        cmdTmp{end+1} = '-num_stimts 1 \';
        fido = fopen(fStim, 'w');
        fprintf(fido,'%.3f ',startSeq(condSeq==1));
        fclose(fido);
        deconWin = mode(diff(startSeq(condSeq==1)));
        deconWin = floor(deconWin/trMri)*trMri;
        dt = trMri*1.5; %trStim;
        b = 0;
        c = deconWin-dt;
        n = (c-b)/dt + 1;
        % (c-b)/(n-1)
        cmdTmp{end+1} = ['-stim_times ' num2str(k) ' ' fStim ' ''TENT(' num2str(b) ',' num2str(c) ',' num2str(n) ')'' \'];
        cmdTmp{end+1} = ['-iresp ' num2str(k) ' ' fResp ' \'];
        cmdTmp{end+1} = ['-fitts ' fFit ' \'];
        cmdTmp{end+1} = ['-errts ' fResid ' \'];
        cmdTmp{end+1} = ['-TR_times ' num2str(trStim) ' \'];
        cmdTmp{end+1} = ['-bucket ' fStat];

        fRun(I,E).fResp = fResp;
        fRun(I,E).fFit = fFit;
        fRun(I,E).fResid = fResid;
        fRun(I,E).fStat = fStat;
        fRun(I,E).cmd =  strjoin(cmdTmp,newline);
        cmd = [cmd cmdTmp(2:end)];

        % [status,cmdout] = system(strjoin(cmdTmp,newline),'-echo');

        %%% extract design matrix
        cmdX = {srcAfni};
        cmdX{end+1} = ['1dcat ' replace(fStat,'.nii.gz','.xmat.1D')];
        [status,cmdout] = system(strjoin(cmdX,newline));
        param.funDsgn.mat = str2num(cmdout);
        if verbose
            figure('WindowStyle','docked');
            imagesc(param.funDsgn.mat); colormap gray
        end

    end
end
if force   ||   ( ~exist(fStat,'file') || ~exist(fStim,'file') || ~exist(fResp,'file') || ~exist(fFit,'file') || ~exist(fResid,'file') )...
        && length(cmd)>1
    cmd = strjoin(cmd,newline); % disp(strjoin(cmd,newline))
    [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
else
    cmd = [];
end

%% fixed-effect model on concatenated runs (separate baselines)
if size(fList,1)>1
    for E = 1:size(fList,2)
        cmd = {srcAfni};
        fIn = fList{1,E};
        tmp = strsplit(fIn,'_'); tmp = char(tmp(contains(tmp,'run-')));
        fOut = replace(fIn,tmp,'run-cat'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        fStim = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_startTime.1D');
        fResp = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_resp.nii.gz');
        fFit = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_fit.nii.gz');
        fResid = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_resid.nii.gz');
        fStat = fullfile(fileparts(replace(fOut,'.nii.gz','')),'cond-visOn_stats.nii.gz');

        cmd = {srcAfni};
        cmd{end+1} = '3dDeconvolve -overwrite \';
        cmd{end+1} = ['-input ' strjoin(fList(:,E),' ') ' \'];
        cmd{end+1} = ['-mask ' fMask ' \'];
        cmd{end+1} = '-polort A \';
        if isempty(trMri)
            trMri = MRIread(fIn,1); trMri = trMri.tr/1000;
        end
        cmd{end+1} = ['-stim_times_subtract ' num2str(trMri*param.nDummy,'%f') ' \'];
        cmd{end+1} = '-num_stimts 1 \';
        fido = fopen(fStim, 'w');
        for i = 1:size(fList,1)
            fprintf(fido,'%.3f ',startSeq(condSeq==1)); fprintf(fido,newline);
        end
        fclose(fido);
        cmd{end+1} = ['-stim_times ' num2str(k) ' ' fStim ' ''TENT(' num2str(b) ',' num2str(c) ',' num2str(n) ')'' \'];
        cmd{end+1} = ['-iresp ' num2str(k) ' ' fResp ' \'];
        cmd{end+1} = ['-fitts ' fFit ' \'];
        cmd{end+1} = ['-errts ' fResid ' \'];
        cmd{end+1} = ['-TR_times ' num2str(trStim) ' \'];
        cmd{end+1} = ['-bucket ' fStat];

        fSes(1,E).fResp = fResp;
        fSes(1,E).fFit = fFit;
        fSes(1,E).fResid = fResid;
        fSes(1,E).fStat = fStat;
        fSes(1,E).cmd = strjoin(cmd,newline);

        if force || ~exist(fStat,'file') || ~exist(fStim,'file') || ~exist(fResp,'file') || ~exist(fFit,'file') || ~exist(fResid,'file')
            cmd = strjoin(cmd,newline); % disp(cmd)
            [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
        end
        cmdX = {srcAfni};
        cmdX{end+1} = ['1dcat ' replace(fStat,'.nii.gz','.xmat.1D')];
        [status,cmdout] = system(strjoin(cmdX,newline));
        param.funDsgn.matCat = str2num(cmdout);
        if param.verbose
            figure('WindowStyle','docked');
            imagesc(param.funDsgn.matCat); colormap gray
        end
    end

else
    fSes = [];
end


