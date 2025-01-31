function runSet = estimMotionBS2(runSet,fBase,fMask,param,force,verbose)
global srcAfni srcFs

ppLabel = 'betweenSesMoco';


if exist('param','var') && ~isempty(param) && isfield(param,'explicitlyFixThroughPlane'); explicitlyFixThroughPlane = param.explicitlyFixThroughPlane; else; explicitlyFixThroughPlane = []; end
if exist('param','var') && ~isempty(param) && isfield(param,'afni3dAlineateArg'); afni3dAlineateArg = param.afni3dAlineateArg; else; afni3dAlineateArg = {}; end
if exist('param','var') && ~isempty(param) && isfield(param,'spSmFac'); spSmFac = param.spSmFac; else; spSmFac = []; end % fraction of voxel size

if ~exist('fMask','var'); fMask = []; end
if ~exist('force','var'); force = []; end
if ~exist('verbose','var'); verbose = []; end
if ~exist('spSmFac','var'); spSmFac = 0; end % fraction of voxel size
if isempty(explicitlyFixThroughPlane); explicitlyFixThroughPlane = 0; end
if isempty(afni3dAlineateArg); afni3dAlineateArg = {'-cost lpa+ZZ' '-interp quintic' '-final wsinc5'}; end
if isempty(force); force = 0; end
if isempty(verbose); verbose = 0; end





%% Get source and base files
if isempty(fBase)
    dbstack; error('code that');
    switch param.baseType
        case 'firstRun_avFrame'
            [a,b,c] = fileparts(runSet.fMocoList);
            fSourceList = fullfile(a,strcat('av_',cellstr(b),cellstr(c)));
            runSet = wr2br(runSet);
            runSet.fEstimList = cellstr(fSourceList);
            runSet.fBase = char(fSourceList(1)); clear fSourceList
        otherwise
            dbstack; error('code that')
        % case 'av'
        %     [a,b,~] = fileparts(replace(runSet.fEstimList,'.nii.gz',''));
        %     runSet.fEstimListBase = fullfile(a,strcat('av_',b,'.nii.gz'));
        % case 'mcAv'
        %     dbstack; error('code that')
        %     paramBase = param;
        %     paramBase.baseType = 'av';
        %     runSet = estimMotionWR(runSet,paramBase,[],[],[]);
    end
else
    switch param.sourceType
        % case 'frstRun_avFrame'
        %     % Define source
        %     runSet.fMocoSmr.runAv
        %     fSourceList = cellstr(fullfile(runSet.wd,'av_cat_mcBR_av_mcWR_setPlumb_volTs.nii.gz'));
        %     runSet = wr2br(runSet);
        %     runSet.fEstimList = fSourceList;
        %     runSet.fMaskList = {};
        %     runSet.fMask     = {};
        % 
        %     % Define base in source session set
        %     mriBase = MRIread(fBase);
        %     mriBaseSetPlumb = MRIread(runSet.fPlumbList{1},1);
        %     mriBaseSetPlumb.vol = mriBase.vol;
        % 
        %     fBaseSetPlumb = strsplit(char(fSourceList),filesep);
        %     fBaseSetPlumb{end} = ['mcBSbase_' fBaseSetPlumb{end}];
        %     fBaseSetPlumb = strjoin(fBaseSetPlumb,filesep);
        % 
        %     MRIwrite(mriBaseSetPlumb,fBaseSetPlumb);
        % 
        %     runSet.fBase         = fBase;
        %     runSet.fBaseSetPlumb = fBaseSetPlumb;
        case 'avRun_avFrame'
            % Define source
            fSourceList = cellstr(runSet.fMocoSmr.sesAv.runAv.fList);
            % fSourceList = cellstr(fullfile(runSet.wd,'av_cat_mcBR_av_mcWR_setPlumb_volTs.nii.gz'));
            runSet = wr2br(runSet);
            runSet.fEstimList = fSourceList;
            runSet.fMaskList = {};
            runSet.fMask     = {};

            % Define base in source session set
            mriBase = MRIread(fBase);
            mriBaseSetPlumb = MRIread(runSet.fPlumbList{1},1);
            mriBaseSetPlumb.vol = mriBase.vol;

            fBaseSetPlumb = strsplit(char(fSourceList),filesep);
            fBaseSetPlumb{end} = ['mcBSbase_' fBaseSetPlumb{end}];
            fBaseSetPlumb = strjoin(fBaseSetPlumb,filesep);

            MRIwrite(mriBaseSetPlumb,fBaseSetPlumb);
            
            runSet.fBase         = fBase;
            runSet.fBaseSetPlumb = fBaseSetPlumb;

        otherwise
            dbstack; error('X');
    end
end



%% Run moco from source to base session
disp('estimating between-session motion')
runSet.fMocoList = cell(size(runSet.fEstimList));
% runSet.manBrainMaskInv = fMask;

for I = 1:numel(runSet.fEstimList)
    switch param.sourceType
        case 'avRun_avFrame'
            disp([' source            : ' runSet.fEstimList{I} ])
            disp([' base              : ' char(runSet.fBase) ])
            disp([' base in source set: ' char(runSet.fBaseSetPlumb) ])
            disp([' mask              : ' char(fMask) ])
        otherwise
            dbstack; error('X');
            disp([' run' num2str(I) '/' num2str(length(runSet.fEstimList))])
    end
    %%% set filename
    fIn = runSet.fEstimList{I};
    % vsize = MRIread(fIn,1); vsize = mean([vsize.xsize vsize.ysize]);
    fOut = strsplit(fIn,filesep); fOut{end} = ['mcBS_' fOut{end}]; fOut = strjoin(fOut,filesep);
    disp([' out               : ' char(fOut) ])
    % fOutWeights = strsplit(fIn,filesep); fOutWeights{end} = ['mcBR_' fOutWeights{end}]; fOutWeights{end} = strsplit(fOutWeights{end},'_'); fOutWeights{end}{end} = 'weights.nii.gz'; fOutWeights{end} = strjoin(fOutWeights{end},'_'); fOutWeights = strjoin(fOutWeights,filesep);
    fOutParam = replace(fOut,'.nii.gz','');
    % fOutAv = strsplit(fOut,filesep); fOutAv{end} = ['av_' fOutAv{end}]; fOutAv = strjoin(fOutAv,filesep);
    if force || ~exist(fOut,'file')
        cmd = {srcAfni};
        %%% moco
        cmd{end+1} = '3dAllineate -overwrite \';
        cmd{end+1} = ['-base ' runSet.fBaseSetPlumb ' \'];
        cmd{end+1} = ['-source ' fIn ' \'];
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = ['-wtprefix ' replace(fOut,'_volTs.nii.gz','_volWeights.nii.gz') ' \'];
        cmd{end+1} = ['-1Dparam_save ' fOutParam ' \'];
        cmd{end+1} = ['-1Dmatrix_save ' fOutParam ' \'];
        cmd{end+1} = [strjoin(afni3dAlineateArg,' ') ' \'];
        if ~isempty(fMask)
            % disp(['  using mask: ' fMask])
            cmd{end+1} = ['-emask ' fMask ' \'];
        else
            disp('  not using mask')
        end
        if explicitlyFixThroughPlane
            cmd{end+1} = '-parfix 2 0 -parfix 4 0 -parfix 5 0 \';
            disp('  explicitly enforcing no through-plane motion')
        end
        % cmd{end+1} = '-nopad -conv 0 -nmatch 100% -onepass -nocmass \';
        cmd{end+1} = '-nopad -conv 0 -nmatch 100% -onepass -cmass+xz \';
        % cmd{end+1} = '-maxrot 1 -maxshf 0.5 \';
        if spSmFac>0 % fraction of voxel size
            fineBlur = mean(runSet.vSize(I,1:2))*spSmFac; % fraction of voxel size to mm
            cmd{end+1} = ['-fineblur ' num2str(fineBlur) ' \']; % mm
        end
        % cmd{end+1} = ['-wtprefix ' fOutWeights ' \'];
        cmd{end+1} = '-warp shift_rotate'; % cmd{end+1} = ['-warp shift_rotate -parfix 2 0 -parfix 4 0 -parfix 5 0'];
        if verbose>1; disp(strjoin(cmd,newline)); end

        % %%% detect smoothing
        % sm = strsplit(fIn,filesep); sm = strsplit(sm{end},'_'); ind = ~cellfun('isempty',regexp(sm,'^sm\d+$')); if any(ind); sm = sm{ind}; else sm = 'sm1'; end; sm = str2num(sm(3:end));
        % n = MRIread(fIn,1); n = n.nframes - 1;
        % nLim = [0 n] + [1 -1].*((sm+1)/2-1);

        %%% execute command
        if verbose>1
            [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
        else
            [status,cmdout] = system(strjoin(cmd,newline)); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
        end

        disp(['  outputs: ' fOutParam])
        disp('  done')
    else
        disp(['  outputs: ' fOutParam])
        disp('  already done, skipping')
    end

    %%% output files
    runSet.fMocoList{I} = fOut;
    runSet.fMaskList{I} = fMask;
    runSet = rmfield(runSet,'fMask');
    runSet = rmfield(runSet,'fMasks');
end



% %% Summarize moco data
% forceThis   = force;
% verboseThis = verbose;
% summarizeVolTs(runSet.fMocoList,[],runSet.nFrame,[],runSet.dataType,forceThis,verboseThis)


runSet.param   = param;
runSet.ppLabel = ppLabel;

return

%%% write means
cmd = {srcFs};
fIn = files.fMoco;
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
    if verbose
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
    else
        [status,cmdout] = system(strjoin(cmd,newline)); if status; dbstack; error(cmdout); error('x'); end
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

files.param = param;

