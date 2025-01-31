function runSet = estimMotionWR2(runSet,param,baseFileList,fMask,force,verbose)
global srcAfni srcFs

ppLabel = 'withinRunMoco';

if exist('param','var') && ~isempty(param) && isfield(param,'explicitlyFixThroughPlane'); explicitlyFixThroughPlane = param.explicitlyFixThroughPlane; else; explicitlyFixThroughPlane = []; end
if exist('param','var') && ~isempty(param) && isfield(param,'afni3dAlineateArg'); afni3dAlineateArg = param.afni3dAlineateArg; else; afni3dAlineateArg = {}; end
if exist('param','var') && ~isempty(param) && isfield(param,'spSmFac'); spSmFac = param.spSmFac; else; spSmFac = []; end

if ~exist('baseFileList','var'); baseFileList = []; end


if ~exist('fMask','var'); fMask = []; end
if ~exist('force','var'); force = []; end
if ~exist('verbose','var'); verbose = []; end
if ~exist('spSmFac','var'); spSmFac = []; end
if isempty(explicitlyFixThroughPlane); explicitlyFixThroughPlane = 0; end
if isempty(afni3dAlineateArg); afni3dAlineateArg = {'-cost ls' '-interp quintic' '-final wsinc5'}; end
if isempty(force); force = 0; end
if isempty(verbose); verbose = 1; end
if isempty(spSmFac); spSmFac = 0; end

nRun = size(runSet.fEstimList,1);


%% Get base files
if isempty(baseFileList)
    switch param.baseType
        case 'first'
            runSet.fEstimListBase = cell(size(runSet.fEstimList));
            for r = 1:nRun
                runSet.fEstimListBase{r} = runSet.fEstimList{r};
                runSet.fEstimListBase{r}(strfind(runSet.fEstimListBase{r},'['):strfind(runSet.fEstimListBase{r},']')) = [];
            end
            runSet.fEstimListBase = strcat(runSet.fEstimListBase,'[0]');
        case 'av'
            runSet.fEstimListBase = cell(size(runSet.fEstimList));
            for r = 1:nRun
                runSet.fEstimListBase{r} = runSet.fEstimList{r};
                runSet.fEstimListBase{r}(strfind(runSet.fEstimListBase{r},'['):strfind(runSet.fEstimListBase{r},']')) = [];
            end
            [a,b,~] = fileparts(replace(runSet.fEstimListBase,'.nii.gz',''));
            runSet.fEstimListBase = fullfile(a,strcat('av_',b,'.nii.gz'));
        case 'mcAv'
            dbstack; error('code that')
            paramBase = param;
            paramBase.baseType = 'av';
            runSet = estimMotionWR(runSet,paramBase,[],[],[]);
    end
else
    dbstack; error('code that');
end


%% Run moco on each run
runSet.cmd            = cell(size(runSet.fEstimList));
runSet.fMocoList      = cell(size(runSet.fEstimList));
runSet.fMocoParamList = cell(size(runSet.fEstimList));
runSet.fMocoMatList   = cell(size(runSet.fEstimList));
runSet.fMaskList      = cell(size(runSet.fEstimList));
for r = 1:nRun
    disp([' run' num2str(r) '/' num2str(nRun)])

    %%% set filename
    fIn = runSet.fEstimList{r};
    fBase = runSet.fEstimListBase{r};
    fOut = strsplit(fIn,filesep); fOut{end} = ['mcWR_' fOut{end}]; fOut = strjoin(fOut,filesep);
    fOut(strfind(fOut,'['):strfind(fOut,']')) = [];
    % fOutWeights = strsplit(fIn,filesep); fOutWeights{end} = ['mcWR_' fOutWeights{end}]; fOutWeights{end} = strsplit(fOutWeights{end},'_'); fOutWeights{end}{end} = 'weights.nii.gz'; fOutWeights{end} = strjoin(fOutWeights{end},'_'); fOutWeights = strjoin(fOutWeights,filesep);
    % fOutPear = strsplit(fIn,filesep); fOutPear{end} = ['mcWR_' fOutPear{end}]; fOutPear{end} = strsplit(fOutPear{end},'_'); fOutPear{end}{end} = 'pearCor.nii.gz'; fOutPear{end} = strjoin(fOutPear{end},'_'); fOutPear = strjoin(fOutPear,filesep);
    fOutParam = replace(fOut,'.nii.gz','');
    % fOutAv = strsplit(fOut,filesep); fOutAv{end} = ['av_' fOutAv{end}]; fOutAv = strjoin(fOutAv,filesep);
    cmd = {srcAfni};
    %%% moco
    cmd{end+1} = '3dAllineate -overwrite \';
    cmd{end+1} = ['-base ' fBase ' \'];
    cmd{end+1} = ['-source ' fIn ' \'];
    cmd{end+1} = ['-prefix ' fOut ' \'];
    cmd{end+1} = ['-wtprefix ' replace(fOut,'_volTs.nii.gz','_volWeigths.nii.gz') ' \'];
    cmd{end+1} = ['-1Dparam_save ' fOutParam ' \'];
    cmd{end+1} = ['-1Dmatrix_save ' fOutParam ' \'];
    cmd{end+1} = ['-SavePear ' fOutParam '_pear.nii.gz \'];
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
    if spSmFac>0
        fineBlur = mean(runSet.vSize(1:2))*param.spSmFac;
        cmd{end+1} = ['-fineblur ' num2str(fineBlur) ' \'];
    end
    % cmd{end+1} = ['-wtprefix ' fOutWeights ' \'];
    % cmd{end+1} = ['-PearSave ' tempname ' \']; % not working (because single slice ??)
    cmd{end+1} = '-warp shift_rotate'; % cmd{end+1} = ['-warp shift_rotate -parfix 2 0 -parfix 4 0 -parfix 5 0'];
    % disp(strjoin(cmd,newline))

    % %%% detect smoothing
    % sm = strsplit(fIn,filesep); sm = strsplit(sm{end},'_'); ind = ~cellfun('isempty',regexp(sm,'^sm\d+$')); if any(ind); sm = sm{ind}; else sm = 'sm1'; end; sm = str2num(sm(3:end));
    % n = MRIread(fIn,1); n = n.nframes - 1;
    % nLim = [0 n] + [1 -1].*((sm+1)/2-1);

    % %%% average
    % cmd{end+1} = '3dTstat -overwrite \';
    % cmd{end+1} = ['-prefix ' fOutAv ' \'];
    % cmd{end+1} = '-mean \';
    % cmd{end+1} = [fOut '[' num2str(nLim(1)) '..' num2str(nLim(2)) ']'];

    runSet.fMocoList{r} = fOut;
    if force || ~exist(fOut,'file')
        %%% execute command
        if verbose>1
            [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
        else
            [status,cmdout] = system(strjoin(cmd,newline)); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
        end

        % %%% adjust motion estimates for smoothing effects
        % editMocoParam([fOutParam '.param.1D'],sm)
        % editMocoParam([fOutParam '.aff12.1D'],sm)


        disp(['  outputs: ' fOutParam])
        disp('  done')
    else
        disp(['  outputs: ' fOutParam])
        disp('  already done, skipping')
    end

    %%% output files
    runSet.cmd{r} = cmd';
    runSet.fMocoList{r} = fOut;
    % runSet.fMocoAvList{r} = fOutAv;
    runSet.fMocoParamList{r} = [fOutParam '.param.1D'];
    runSet.fMocoMatList{r} = [fOutParam '.aff12.1D'];

    runSet.fMaskList{r} = fMask;
end


%% Summarize moco data
forceThis   = force;
verboseThis = verbose;
% runSet.fMocoSmr = summarizeVolTs2(runSet.fMocoList,[],[],runSet.nDummy,[],forceThis,verboseThis);
runSet.fMocoSmr = summarizeVolTs4(runSet.fMocoList,0,runSet.dataType,forceThis,verboseThis);
runSet.param   = param;
runSet.ppLabel = ppLabel;


% forceThis   = force;
% verboseThis = verbose;
% for r = 1:nRun
%     ind = strfind(runSet.fEstimList{r},'['):strfind(runSet.fEstimList{r},']');
%     if ~isempty(ind)
%         ind([1 end]) = [];
%         tmp = strsplit(runSet.fEstimList{r}(ind),'..');
%         if strcmp(tmp{2},'$'); tmp{2} = runSet.nFrame(r); else tmp{2} = str2double(tmp{2}); end
%         tmp{1} = str2double(tmp{1});
%         runSet.nFrame(r) = diff([tmp{:}]);
%     end
% end
% summarizeVolTs(runSet.fMocoList,[],runSet.nFrame,[],runSet.dataType,forceThis,verboseThis)


% runSet.param   = param;
% runSet.ppLabel = ppLabel;

return

disp(['estimating within-run motion (second-pass moco to ' param.baseType '; accounting for smoothing)'])
fMask = runSet.manBrainMaskInv;
runSet.fMocoList = cell(size(runSet.fEstimList));
runSet.fMocoAvList = cell(size(runSet.fEstimList));
runSet.fMocoParamList = cell(size(runSet.fEstimList));
runSet.fMocoMatList = cell(size(runSet.fEstimList));
% files.fBase = cell(size(files.fEstim));
for I = 1:length(runSet.fEstimList)
    disp([' run' num2str(I) '/' num2str(length(runSet.fEstimList))])
    %%% set filename
    fIn = runSet.fEstimList{I};
    vsize = MRIread(fIn,1); vsize = mean([vsize.xsize vsize.ysize]);
    fBase = runSet.fEstimListBase{I};
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
    runSet.fMocoList{I} = fOut;
    runSet.fMocoAvList{I} = fOutAv;
    runSet.fMocoParamList{I} = [fOutParam '.param.1D'];
    runSet.fMocoMatList{I} = [fOutParam '.aff12.1D'];
end

%%% write means
if nRun>1
    cmd = {srcFs};
    fIn = runSet.fMocoAvList;
    fOut = replace(fIn{1},char(regexp(fIn{1},'run-\d+','match')),'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
    fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
    if force || ~exist(fOut,'file')
        cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
    end
    runSet.fMocoCatAv = fOut;

    fIn = fOut;
    fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
    if force || ~exist(fOut,'file')
        cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
    end
    runSet.fMocoAvListCatAv = fOut;
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
    cmd{end+1} = ['fslview -m single ' runSet.fMocoCatAv ' \'];
    if isfield(runSet,'manBrainMaskInv') && ~isempty(runSet.manBrainMaskInv)
        cmd{end+1} = [runSet.manBrainMaskInv ' &'];
    end
    cmd = strjoin(cmd,newline); % disp(cmd)
    if verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    runSet.qaFiles.fFslviewBR = cmd;
end
%%%% within-run motion
%%%%% rewrite files with [frame] because not working with fslview grrrrr
cmd = {srcAfni};
for i = 1:length(runSet.fEstimListBase)
    if strcmp(runSet.fEstimListBase{i}(end),']')
        fIn = runSet.fEstimListBase{i};
        fOut = fIn; [~,fOut] = regexp(fOut,'\[.*\]','match','split'); fOut = strsplit(fOut{1},filesep); fOut{end} = ['mcRef-first_' fOut{end}]; fOut = strjoin(fOut,filesep);

        cmd{end+1} = '3dcalc -overwrite \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = ['-a ' fIn ' \'];
        cmd{end+1} = '-expr ''a''';

        runSet.fEstimListBase{i} = fOut;
    end
end
if length(cmd)>1
    cmd = strjoin(cmd,newline); % disp(cmd)
    [status,cmdout] = system(cmd); if status; dbstack; error(cmdout); error('x'); end
end
%%%%% generate the command including the bases
cmd = {srcFs};
cmd{end+1} = ['fslview -m single ' strjoin([runSet.fEstimListBase runSet.fMocoList],' ') ' \'];
if isfield(runSet,'manBrainMaskInv') && ~isempty(runSet.manBrainMaskInv)
    cmd{end+1} = [runSet.manBrainMaskInv ' &'];
end
cmd = strjoin(cmd,newline); % disp(cmd)
if verbose
    [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
end
runSet.qaFiles.fFslviewWR = cmd;

%%%% motion between first, mid and last frames (acounting for smoothing) of each run
runSet.qaFiles.fFslviewWRfstMdLst = qaFstMdLst(runSet.fMocoList,force,verbose);
runSet = addMaskToCmd(runSet);

runSet.param = param;



