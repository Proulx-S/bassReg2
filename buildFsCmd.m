function buildFsCmd(xSet,fT1w,fB1)
global srcFs srcAfni
if ~exist('fT1w','var'); fT1w = []; end
if ~exist('fB1' ,'var'); fB1  = []; end

cmd = {srcFs};
cmd{end+1} = 'freeview \';

%% Get set labels
labelList = cell(size(xSet))';
for S = 1:length(xSet)
    labelList{S} = xSet{S}.label;
end

%% Base image
S = find(contains(labelList,'f1b1') | contains(labelList,'f1b0'),1);
candidate = {'fAvCatAvEchoRms' 'fAvEchoRms' 'fEchoRms' 'fAvCatAv' 'fAv' 'f'};
candidate = candidate(ismember(candidate,fields(xSet{S}.finalFiles)));
fBase = xSet{S}.finalFiles.(candidate{1}){1};

cmd{end+1} = [fBase ':name=ref:visible=0 \'];

%% Underlays

%%% functional run averages
S = find(contains(labelList,'f1b1'),1);
labelSuf = {'runAv'};
if ~isempty(S)
    cmd = [cmd buildFuncRunAverage(xSet{S},labelSuf)];
end
S = find(contains(labelList,'f0b0'),1);
labelSuf = {'runAv'};
if ~isempty(S)
    cmd = [cmd buildFuncRunAverage(xSet{S},labelSuf)];
end
S = find(contains(labelList,'f1b0'),1);
labelSuf = {'runAv'};
if ~isempty(S)
    cmd = [cmd buildFuncRunAverage(xSet{S},labelSuf)];
end
S = find(contains(labelList,'f0b1'),1);
labelSuf = {'runAv'};
if ~isempty(S)
    cmd = [cmd buildFuncRunAverage(xSet{S},labelSuf)];
end


%%% T1w
if ~isempty(fT1w)
    cmd{end+1} = [fT1w ':resample=cubic:name=T1w:visible=1 \'];
end

%%% B1map
if ~isempty(fB1)
    cmd{end+1} = [fB1 ':resample=cubic:colormap=turbo:colorscale=60,120:name=B1map:visible=0 \'];
end

%%% avMap
S = find(contains(labelList,'avMap'),1);

candidate = {'fAvCatAvEchoCat' 'fAvEchoCat' 'fEchoCat'};
candidate = candidate(ismember(candidate,fields(xSet{S}.finalFiles)));
f = xSet{S}.finalFiles.(candidate{1}){1};
setLabel = replace(xSet{S}.label,'set-','');
cmd{end+1} = [f ':resample=cubic:name=' setLabel '_echo-cat_avMap:visible=1 \'];

candidate = {'fAvCatAvEchoRms' 'fAvEchoRms' 'fEchoRms'};
candidate = candidate(ismember(candidate,fields(xSet{S}.finalFiles)));
f = xSet{S}.finalFiles.(candidate{1}){1};
setLabel = replace(xSet{S}.label,'set-','');
cmd{end+1} = [f ':resample=cubic:name=' setLabel '_echo-rms_avMap:visible=0 \'];

%%% functional averages
S = find(contains(labelList,'f1b0'),1);
labelSuf = {'tsAv'};
if ~isempty(S)
    cmd = [cmd buildFuncAverage(xSet{S},labelSuf)];
end

S = find(contains(labelList,'f0b1'),1);
labelSuf = {'tsAv'};
if ~isempty(S)
    cmd = [cmd buildFuncAverage(xSet{S},labelSuf)];
end

S = find(contains(labelList,'f1b1'),1);
labelSuf = {'tsAv'};
if ~isempty(S)
    cmd = [cmd buildFuncAverage(xSet{S},labelSuf)];
end

S = find(contains(labelList,'f0b0'),1);
labelSuf = {'tsAv'};
if ~isempty(S)
    cmd = [cmd buildFuncAverage(xSet{S},labelSuf)];
end

S = find(contains(labelList,'bold'),1);
labelSuf = {'tsAv'};
if ~isempty(S)
    cmd = [cmd buildFuncAverage(xSet{S},labelSuf)];
end

%%% response timecourse
onBaseFlag = 0;
for S = 1:length(xSet)
    if any(ismember(xSet{S}.dataType,'task'))
        labelSuf = {'tsResp'};
        if ~isempty(S) && xSet{S}.nEcho>1
            cmd = [cmd buildRespTs(xSet{S},labelSuf,onBaseFlag)];
        end
    end
end
onBaseFlag = 1;
for S = 1:length(xSet)
    if any(ismember(xSet{S}.dataType,'task'))
        labelSuf = {'tsRespOnBase'};
        if ~isempty(S) && xSet{S}.nEcho>1
            cmd = [cmd buildRespTs(xSet{S},labelSuf,onBaseFlag)];
        end
    end
end



%%% satin
cmd = [cmd buildSatin(xSet)];

%%% pc
cmd = [cmd buildPC(xSet)];


%%% activation maps
cmd = [cmd buildActivationMaps(xSet)];


%%% mask
cmd{end+1} = [xSet{1}.finalFiles.manBrainMaskInv ' \'];



cmd{end} = replace(cmd{end},' \',' &');
clipboard('copy',strjoin(cmd,newline))
disp('copied to clipboard')


function cmd = buildActivationMaps(xSet)
cmd = {};

for S = 1:length(xSet)
    if any(ismember(xSet{S}.dataType,'task'))

        % Model-free
        modelField = 'modelFree';
        statField = 'respF';
        outLabel = 'respF';
        statMax = 10;
        cmd = [cmd buildStatWithThresh(xSet{S},modelField,statField,outLabel,[],[],statMax)];

        % Custom model
        modelField = 'customModel'; modelField1mFDR = [];%'modelFree';
        statField  = 'respT'      ;  statField1mFDR = [];%'respF_fdr';
        outLabel   = 'ampT';
        statMax = 10;
        cmd = [cmd buildStatWithThresh(xSet{S},modelField,statField,outLabel,modelField1mFDR,statField1mFDR,statMax)];

    end
end


function cmd = buildStatWithThresh(xSet,modelField,statField,outLabel,modelField1mFDR,statField1mFDR,threshMax)
global srcAfni
if ~exist('modelField1mFDR','var'); modelField1mFDR = []; end
if ~exist('statField1mFDR' ,'var');  statField1mFDR = []; end
if isempty(modelField1mFDR) && isempty(statField1mFDR); maskFlag = 0; else maskFlag = 1; end
if ~exist('threshMax','var'); threshMax = []; end

% if ~exist('modelFieldThresh','var'); modelFieldThresh = modelField; end
% if ~exist('statFieldThresh' ,'var');  statFieldThresh = statField; end
cmd = {};

setLabel = replace(xSet.label,'set-','');
if isfield(xSet.statFiles,modelField)
    if xSet.nEcho>1
        if maskFlag
            %% Load masks for thresholding
            fFdr = xSet.statFiles.(modelField1mFDR).(statField1mFDR);
            
            %%% each echo
            candidate = {'fSes' 'f'}; candidate = candidate(ismember(candidate,fields(fFdr)));
            if ~isempty(candidate)
                candidate = candidate{1};
                maskName = cell(size(fFdr.(candidate)));
                for E = xSet.nEcho:-1:1
                    if ~exist(fFdr.(candidate){E},'file'); dbstack; error('something wrong'); end
                    maskName{E} = [setLabel '_echo-' num2str(E) '_' outLabel 'fdrMaskFromModelFree'];
                    cmd{end+1} = [fFdr.(candidate){E} ':resample=cubic:name=' maskName{E} ':visible=0 \'];
                end
            end
            
            %%% echo rms
            candidate = {'fSesEchoRms' 'fEchoRms'}; candidate = candidate(ismember(candidate,fields(fFdr)));
            if ~isempty(candidate)
                candidate = candidate{1};
                maskName_rms = cell(size(fFdr.(candidate)));
                E = 1;
                if ~exist(fFdr.(candidate){E},'file'); dbstack; error('something wrong'); end
                maskName_rms{E} = [setLabel '_echo-rms_' outLabel 'fdrMaskFromModelFree'];
                cmd{end+1} = [fFdr.(candidate){E} ':resample=cubic:name=' maskName_rms{E} ':visible=0 \'];
            end
        end
        
        %% Load actual statistical maps
        f = xSet.statFiles.(modelField).(statField);

        %%% each echo
        candidate = {'fSes' 'f'}; candidate = candidate(ismember(candidate,fields(f)));
        if ~isempty(candidate)
            candidate = candidate{1};
        
            %%%% get colorscale limits from stats (fdr thresh and max)
            cmdTmp = {srcAfni};
            for E = 1:xSet.nEcho    
                cmdTmp{end+1} = ['fdrval -qinput ' f.(candidate){E} ' 0 0.05 2> /dev/null'];
                cmdTmp{end+1} = ['3dBrickStat -max ' f.(candidate){E} ' 2> /dev/null'];
            end
            [~,thresh] = system(strjoin(cmdTmp,newline));
            thresh = str2num(thresh); thresh = cat(2,thresh(1:2:end),thresh(2:2:end));
            thresh(diff(thresh,[],2)<0,2) = thresh(diff(thresh,[],2)<0,1); % make sure max>min
            if ~isempty(threshMax) % arbitrary max
                thresh(:,2) = threshMax;
            end
            if any(isnan(thresh(:)))
                keyboard
            end
            %%%% make the command
            for E = xSet.nEcho:-1:1
                if maskFlag
                    %%%%% threshold using a custom mask
                    cmd{end+1} = [f.(candidate){E} ':colormap=heat:heatscale=' num2str(thresh(E,1)) ',' num2str(thresh(E,2)) ':resample=cubic:name=' setLabel '_echo-' num2str(E) '_' outLabel ':visible=0:mask=' maskName{E} ' \'];
                else
                    %%%%% threshold using colorscale
                    cmd{end+1} = [f.(candidate){E} ':colormap=heat:heatscale=' num2str(thresh(E,1)) ',' num2str(thresh(E,2)) ':resample=cubic:name=' setLabel '_echo-' num2str(E) '_' outLabel ':visible=0 \'];
                end
            end
        end

        %%% echo cat
        candidate = {'fSesEchoCat' 'fEchoCat'}; candidate = candidate(ismember(candidate,fields(f)));
        if ~isempty(candidate)
            candidate = candidate{1};
            %%%% get colorscal limits from stats
            thresh = mean(thresh,1);
            %%%% make the command
            E = 1;
            cmd{end+1} = [f.(candidate){E} ':colormap=heat:heatscale=' num2str(thresh(E,1)) ',' num2str(thresh(E,2)) ':resample=cubic:name=' setLabel '_echo-cat_' outLabel ':visible=0 \'];
        end

        %%% echo rms
        candidate = {'fSesEchoRms' 'fEchoRms'}; candidate = candidate(ismember(candidate,fields(f)));
        if ~isempty(candidate)
            candidate = candidate{1};
            %%%% get colorscal limits from stats
            cmdTmp = {srcAfni};
            E = 1;
            cmdTmp{end+1} = ['fdrval -qinput ' f.(candidate){E} ' 0 0.05 2> /dev/null'];
            cmdTmp{end+1} = ['3dBrickStat -max ' f.(candidate){E} ' 2> /dev/null'];
            [~,thresh] = system(strjoin(cmdTmp,newline));
            thresh = str2num(thresh); thresh = cat(2,thresh(1:2:end),thresh(2:2:end));
            thresh(diff(thresh,[],2)<0,2) = thresh(diff(thresh,[],2)<0,1); % make sure max>min
            if ~isempty(threshMax) % arbitrary max
                thresh(:,2) = threshMax;
            end
            if any(isnan(thresh(:)))
                keyboard
            end
            %%%% make the command
            E = 1;
            if maskFlag
                %%%%% threshold using a custom mask
                cmd{end+1} = [f.(candidate){E} ':colormap=heat:heatscale=' num2str(thresh(E,1)) ',' num2str(thresh(E,2)) ':resample=cubic:name=' setLabel '_echo-rms_' outLabel ':visible=0:mask=' maskName_rms{E} ' \'];
            else
                %%%%% threshold using colorscale
                cmd{end+1} = [f.(candidate){E} ':colormap=heat:heatscale=' num2str(thresh(E,1)) ',' num2str(thresh(E,2)) ':resample=cubic:name=' setLabel '_echo-rms_' outLabel ':visible=0 \'];
            end
        end
    end
end




function cmd = buildPC(xSet)
cmd = {};
labelList = cell(size(xSet))';
for S = 1:length(xSet)
    labelList{S} = xSet{S}.label;
end

S = find(contains(labelList,'-pc'),1);
if isempty(S); return; end
for i = 1:size(xSet{S}.pcVel.f,3)
    f = char(xSet{S}.pcVel.f(:,:,i));
    venc = strsplit(f,'_'); venc = strsplit(venc{contains(venc,'proc-venc')},'-'); venc = replace(replace(venc{end},'venc',''),'m0','');
    cmd{end+1} = [f ':colormap=heat:heatscale=0.00000000001,' venc ':name=venc-' num2str(venc) '_CmPerSec:visible=0 \'];
end



function cmd = buildSatin(xSet)
cmd = {};

labelList = cell(size(xSet))';
for S = 1:length(xSet)
    labelList{S} = xSet{S}.label;
end

SList = [];
S = find(contains(labelList,'f1b1'),1); if ~isempty(S); SList(end+1) = S; end;
S = find(contains(labelList,'f0b0'),1); if ~isempty(S); SList(end+1) = S; end;
S = find(contains(labelList,'f1b0'),1); if ~isempty(S); SList(end+1) = S; end;
S = find(contains(labelList,'f0b1'),1); if ~isempty(S); SList(end+1) = S; end;

for S = SList
    if isfield(xSet{S},'satIndexAbs')
        if xSet{S}.nEcho>1
            f = xSet{S}.satIndexAbs.fEchoCat;
            setLabel = replace(xSet{S}.label,'set-','');
            cmd{end+1} = [f{1} ':colormap=heat:heatscale=0.00000000001,1:name=' setLabel '_echo-cat_satIndexAbs:visible=0 \'];

            f = xSet{S}.satIndexAbs.fEchoRms;
            setLabel = replace(xSet{S}.label,'set-','');
            cmd{end+1} = [f{1} ':colormap=heat:heatscale=0.00000000001,1:name=' setLabel '_echo-rms_satIndexAbs:visible=0 \'];

            f = xSet{S}.satIndexAbs.fS0;
            setLabel = replace(xSet{S}.label,'set-','');
            cmd{end+1} = [f{1} ':colormap=heat:heatscale=0.00000000001,1:name=' setLabel '_echo-S0_satIndexAbs:visible=0 \'];
        else
            dbstack; error('code that');
        end
    end
    if isfield(xSet{S},'satIndexDir')
        if xSet{S}.nEcho>1
            f = xSet{S}.satIndexDir.fEchoCat;
            setLabel = replace(xSet{S}.label,'set-','');
            cmd{end+1} = [f{1} ':colormap=heat:heatscale=0.00000000001,1:name=' setLabel '_echo-cat_satIndexDir:visible=0 \'];

            f = xSet{S}.satIndexDir.fEchoRms;
            setLabel = replace(xSet{S}.label,'set-','');
            cmd{end+1} = [f{1} ':colormap=heat:heatscale=0.00000000001,1:name=' setLabel '_echo-rms_satIndexDir:visible=0 \'];

            f = xSet{S}.satIndexDir.fS0;
            setLabel = replace(xSet{S}.label,'set-','');
            cmd{end+1} = [f{1} ':colormap=heat:heatscale=0.00000000001,1:name=' setLabel '_echo-S0_satIndexDir:visible=0 \'];
        else
            dbstack; error('code that');
        end
    end
end



function cmd = buildFuncAverage(xSet,labelSuf)
cmd = {};
setLabel = replace(xSet.label,'set-','');
if xSet.nEcho>1
    label = {'echo-cat'};
    candidate = {'fAvCatAvEchoCat' 'fAvEchoCat' 'fEchoCat'};
    candidate = candidate(ismember(candidate,fields(xSet.finalFiles)));
    if ~isempty(candidate)
        f = xSet.finalFiles.(candidate{1}){1};
        cmd{end+1} = [f ':resample=cubic:name=' setLabel '_' strjoin([label labelSuf],'_') ':visible=0 \'];
    end

    label = {'echo-rms'};
    candidate = {'fAvCatAvEchoRms' 'fAvEchoRms' 'fEchoRms'};
    candidate = candidate(ismember(candidate,fields(xSet.finalFiles)));
    if ~isempty(candidate)
        f = xSet.finalFiles.(candidate{1}){1};
        cmd{end+1} = [f ':resample=cubic:name=' setLabel '_' strjoin([label labelSuf],'_') ':visible=0 \'];
    end

    label = {'echo-S0'};
    candidate = {'fAvCatAvS0' 'fAvS0' 'fS0'};
    candidate = candidate(ismember(candidate,fields(xSet.finalFiles)));
    if ~isempty(candidate)
        f = xSet.finalFiles.(candidate{1}){1};
        cmd{end+1} = [f ':resample=cubic:name=' setLabel '_' strjoin([label labelSuf],'_') ':visible=0 \'];
    end

    label = {'echo-R2s'};
    candidate = {'fAvCatAvR2s' 'fAvR2s' 'fR2s'};
    candidate = candidate(ismember(candidate,fields(xSet.finalFiles)));
    if ~isempty(candidate)
        f = xSet.finalFiles.(candidate{1}){1};
        cmd{end+1} = [f ':resample=cubic:name=' setLabel '_' strjoin([label labelSuf],'_') ':visible=0 \'];
    end
else
    candidate = {'fAv' 'f'};
    candidate = candidate(ismember(candidate,fields(xSet.finalFiles)));
    if ~isempty(candidate)
        f = xSet.finalFiles.(candidate{1}){1};
        cmd{end+1} = [f ':resample=cubic:name=' strjoin([{setLabel} labelSuf],'_') ':visible=0 \'];
    end
end



function cmd = buildFuncRunAverage(xSet,labelSuf)
cmd = {};
setLabel = replace(xSet.label,'set-','');
if xSet.nEcho>1
    label = {'run-cat' 'echo-rms'};
    candidate = {'fCatAvEchoRms' 'fCatEchoRms' 'fAvEchoRms' 'fEchoRms'};
    candidate = candidate(ismember(candidate,fields(xSet.finalFiles)));
    if ~isempty(candidate)
        candidate = candidate{1};
        f = xSet.finalFiles.(candidate){1};
        cmd{end+1} = [f ':resample=cubic:name=' setLabel '_' strjoin([label labelSuf],'_') ':visible=0 \'];
    end

    label = {'run-cat' 'echo-S0'};
    candidate = {'fCatAvS0' 'fCatS0' 'fAvS0' 'fS0'};
    candidate = candidate(ismember(candidate,fields(xSet.finalFiles)));
    if ~isempty(candidate)
        f = xSet.finalFiles.(candidate{1}){1};
        cmd{end+1} = [f ':resample=cubic:name=' setLabel '_' strjoin([label labelSuf],'_') ':visible=0 \'];
    end

    % label = {'run-cat' 'echo-S0Perc'};
    % candidate = {'fCatAvS0Perc' 'fCatS0Perc' 'fAvS0Perc' 'fS0Perc'};
    % candidate = candidate(ismember(candidate,fields(xSet.finalFiles)));
    % if ~isempty(candidate)
    %     f = xSet.finalFiles.(candidate{1}){1};
    %     cmd{end+1} = [f ':resample=cubic:name=' setLabel '_' strjoin([label labelSuf],'_') ':visible=0 \'];
    % end

    label = {'run-cat' 'echo-R2s'};
    candidate = {'fCatAvR2s' 'fCatR2s' 'fAvR2s' 'fR2s'};
    candidate = candidate(ismember(candidate,fields(xSet.finalFiles)));
    if ~isempty(candidate)
        f = xSet.finalFiles.(candidate{1}){1};
        cmd{end+1} = [f ':resample=cubic:name=' setLabel '_' strjoin([label labelSuf],'_') ':visible=0 \'];
    end

    % label = {'run-cat' 'echo-R2sPerc'};
    % candidate = {'fCatAvR2sPerc' 'fCatR2sPerc' 'fAvR2sPerc' 'fR2sPerc'};
    % candidate = candidate(ismember(candidate,fields(xSet.finalFiles)));
    % if ~isempty(candidate)
    %     f = xSet.finalFiles.(candidate{1}){1};
    %     cmd{end+1} = [f ':resample=cubic:name=' setLabel '_' strjoin([label labelSuf],'_') ':visible=0 \'];
    % end
else
    dbstack; error('code that')
end



function cmd = buildRespTs(xSet,labelSuf,onBaseFlag)
if ~exist('onBaseFlag','var'); onBaseFlag = 0; end
cmd = {};
if xSet.nEcho>1
    
    % each echo
    if onBaseFlag
        f = xSet.statFiles.modelFree.respOnBase;
    else
        f = xSet.statFiles.modelFree.resp;
    end
    setLabel = replace(xSet.label,'set-','');
    candidate = {'fSes' 'f'}; candidate = candidate(ismember(candidate,fields(f))); candidate = candidate{1};
    if ~isempty(candidate)
        f = f.(candidate);
        for E = length(f):-1:1
            cmd{end+1} = [f{E} ':name=' setLabel '_echo-' num2str(E) '_' char(labelSuf) ':visible=0 \'];
        end
    end

    % echo-rms
    if onBaseFlag
        f = xSet.statFiles.modelFree.respOnBase;
    else
        f = xSet.statFiles.modelFree.resp;
    end
    candidate = {'fSesEchoRms' 'fEchoRms'}; candidate = candidate(ismember(candidate,fields(f))); candidate = candidate{1};
    if ~isempty(candidate)
        f = f.(candidate);
        cmd{end+1} = [char(f) ':name=' setLabel '_echo-rms_' char(labelSuf) ':visible=0 \'];
    end
    
    % S0
    if onBaseFlag
        f = xSet.statFiles.modelFree.respOnBase;
    else
        f = xSet.statFiles.modelFree.resp;
    end
    candidate = {'fSesS0' 'fS0'}; candidate = candidate(ismember(candidate,fields(f)));
    if ~isempty(candidate)
        candidate = candidate{1};
        f = f.(candidate);
        cmd{end+1} = [char(f) ':name=' setLabel '_echo-S0_' char(labelSuf) ':visible=0 \'];
    end

    % S0 %change
    if onBaseFlag
        f = xSet.statFiles.modelFree.respOnBase;
    else
        f = xSet.statFiles.modelFree.resp;
    end
    candidate = {'fSesS0Perc' 'fS0Perc'}; candidate = candidate(ismember(candidate,fields(f)));
    if ~isempty(candidate)
        candidate = candidate{1};
        f = f.(candidate);
        cmd{end+1} = [char(f) ':name=' setLabel '_echo-S0Perc_' char(labelSuf) ':visible=0 \'];
    end

    % T2*
    if onBaseFlag
        f = xSet.statFiles.modelFree.respOnBase;
    else
        f = xSet.statFiles.modelFree.resp;
    end
    candidate = {'fSesT2s' 'fT2s'}; candidate = candidate(ismember(candidate,fields(f)));
    if ~isempty(candidate)
        candidate = candidate{1};
        f = f.(candidate);
        cmd{end+1} = [char(f) ':name=' setLabel '_echo-T2s_' char(labelSuf) ':visible=0 \'];
    end

    % R2*
    if onBaseFlag
        f = xSet.statFiles.modelFree.respOnBase;
    else
        f = xSet.statFiles.modelFree.resp;
    end
    candidate = {'fSesR2s' 'fR2s'}; candidate = candidate(ismember(candidate,fields(f)));
    if ~isempty(candidate)
        candidate = candidate{1};
        f = f.(candidate);
        cmd{end+1} = [char(f) ':name=' setLabel '_echo-R2s_' char(labelSuf) ':visible=0 \'];
    end

    % R2* %change
    if onBaseFlag
        f = xSet.statFiles.modelFree.respOnBase;
    else
        f = xSet.statFiles.modelFree.resp;
    end
    candidate = {'fSesR2sPerc' 'fR2sPerc'}; candidate = candidate(ismember(candidate,fields(f)));
    if ~isempty(candidate)
        candidate = candidate{1};
        f = f.(candidate);
        cmd{end+1} = [char(f) ':name=' setLabel '_echo-R2sPerc_' char(labelSuf) ':visible=0 \'];
    end

else

    dbstack; error('code that')
    f = xSet.statFiles.fullSes_echoCat.fRespOnBase;
    cmd{end+1} = [f ':name=' char(labelSuf) ':visible=0 \'];

end







