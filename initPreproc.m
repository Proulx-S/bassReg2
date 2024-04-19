function xSet = initPreproc(xSet,dOut,geomRef,param,force)
if ~exist('force','var');        force = []       ; end
if ~isfield(xSet,'dataType'); dataType = {}       ; else; dataType = xSet.dataType; end
if ~isfield(param,'verbose');  verbose = []       ; else; verbose = param.verbose; end
if isempty(force);               force = 0        ; end
if isempty(dataType);         dataType = {'volTs'}; elseif ischar(dataType); dataType = {dataType}; end
if isempty(verbose);           verbose = 0        ; end



forceRewriteAtPlumb = force;
forceRecomputeRMS = force;
forceResmoothing = force;
global srcFs srcAfni

% geomRef refers to the file to be used as the reference for the slice
% prescription in scanner space
% Hint: this should be a scan with a scanner-space slice presciption that
% matches that of the most important auxiliary scan, so they show nicely in
% freeview. This matters moslty when areg was used.




%% Massage inputs
%%% funcSet
%%%% sort by acquisition time
[~,b] = sort(getAcqTime(fullfile({xSet.files.folder},{xSet.files.name})));
xSet.files = xSet.files(b);
%%%% refactor
dIn = xSet.files(1).folder;
fOrig = fullfile(dIn,{xSet.files.name}');
label = xSet.label;
bidsList = replace({xSet.files.name}','.nii.gz',''); for R = 1:length(bidsList); bidsList{R} = strsplit(bidsList{R},'_'); end; bidsList = cat(1,bidsList{:})';
%%% dOut
dOut = fullfile(dOut,label); if ~exist(dOut,'dir'); mkdir(dOut); end
%%% geomRef
if ~exist('geomRef','var'); geomRef = []; end
if isempty(geomRef)
    if isfield(xSet,'fGeom') && ~isempty(xSet.fGeom)
        fRef = xSet.fGeom;
    else
        fRef = fOrig{1};
    end
else
    if ~ischar(geomRef)
        fRef = fOrig{geomRef};
    else
        fRef = geomRef;
    end
end
xSet.fGeom = fRef;
% param
if ~exist('param','var'); param = []; end
if ~isfield(param,'nDummy') || isempty(param.nDummy); param.nDummy = 0; end
if ~isfield(param,'verbose') || isempty(param.verbose); param.verbose = 1; end


%% Rewrite data at plumb without dummies
disp('Rewriting data at plumb')

%%% Get oblique and plumb files
disp([' writing oblique and plumb references (for ' num2str(length(fOrig)) ' files)'])
cmd = {srcAfni};
%%%% references
fIn = fRef;
nframes = MRIread(fIn,1); nframes = nframes.nframes; if nframes>8; nframes = 8; end
fPlumbRef = fullfile(dOut,'setPlumb_volRef.nii.gz');
if forceRewriteAtPlumb || ~exist(fPlumbRef,'file')
    cmd{end+1} = '3dcalc -overwrite \';
    cmd{end+1} = ['-a ' fIn '[0..' num2str(nframes-1) '] \'];
    cmd{end+1} = '-expr a \';
    if verbose
        cmd{end+1} = ['-prefix ' fPlumbRef];
        cmd{end+1} = ['3drefit -deoblique ' fPlumbRef];
    else
        cmd{end+1} = ['-prefix ' fPlumbRef ' > /dev/null 2>&1'];
        cmd{end+1} = ['3drefit -deoblique ' fPlumbRef ' > /dev/null 2>&1'];
    end
    % cmd{end+1} = '3dNwarpApply -overwrite \';
    % cmd{end+1} = ['-nwarp ''IDENT(' fIn ')'' \'];
    % cmd{end+1} = ['-prefix ' fPlumbRef ' \'];
    % cmd{end+1} = ['-source ' fIn '[0..' num2str(nframes-1) ']'];
end
fObliqueRef = fullfile(dOut,'setOblique_volRef.nii.gz');
if forceRewriteAtPlumb || ~exist(fObliqueRef,'file')
    cmd{end+1} = '3dcalc -overwrite \';
    cmd{end+1} = ['-a ' fIn '[0..' num2str(nframes-1) '] \'];
    cmd{end+1} = ['-expr a \'];
    if verbose
        cmd{end+1} = ['-prefix ' fObliqueRef ' > /dev/null 2>&1'];
    else
        cmd{end+1} = ['-prefix ' fObliqueRef];
    end
end
%%%% individual files
for R = 1:length(fOrig)
    fIn = fOrig{R};
    nFrame{R} = MRIread(fIn,1); nFrame{R} = nFrame{R}.nframes;
    if nFrame{R}>8; nframes = 8; else; nframes = nFrame{R}; end
    fOut = replace(fIn,dIn,dOut); fOut = replace(fOut,'.nii.gz',''); if ~exist(fOut,'dir'); mkdir(fOut); end;
    fOut = fullfile(fOut,'runPlumb_volRef.nii.gz');
    if forceRewriteAtPlumb || ~exist(fOut,'file')
        cmd{end+1} = '3dcalc -overwrite \';
        cmd{end+1} = ['-a ' fIn '[0..' num2str(nframes-1) '] \'];
        cmd{end+1} = ['-expr a \'];
        if verbose
            cmd{end+1} = ['-prefix ' fOut];
            cmd{end+1} = ['3drefit -deoblique ' fOut];
        else
            cmd{end+1} = ['-prefix ' fOut ' > /dev/null 2>&1'];
            cmd{end+1} = ['3drefit -deoblique ' fOut ' > /dev/null 2>&1'];
        end
        % cmd{end+1} = '3dNwarpApply -overwrite \';
        % cmd{end+1} = ['-nwarp ''IDENT(' fIn ')'' \'];
        % cmd{end+1} = ['-prefix ' fOut ' \'];
        % cmd{end+1} = ['-source ' fIn '[0..' num2str(nframes-1) ']'];
    end
    fOut = replace(fOut,'runPlumb_volRef.nii.gz','runOblique_volRef.nii.gz');
    if forceRewriteAtPlumb || ~exist(fOut,'file')
        cmd{end+1} = '3dcalc -overwrite \';
        cmd{end+1} = ['-a ' fIn '[0..' num2str(nframes-1) '] \'];
        cmd{end+1} = ['-expr a \'];
        if verbose
            cmd{end+1} = ['-prefix ' fOut ' > /dev/null 2>&1'];
        else
            cmd{end+1} = ['-prefix ' fOut];
        end
    end
end

nFrame = [nFrame{:}];
if any(diff(nFrame)); dbstack; error('something wierd'); end
nFrame = nFrame(1);

%%%% launch command
if length(cmd)>1
    cmd = strjoin(cmd,newline); % disp(cmd)
    [status,cmdout] = system(cmd); if status; dbstack; error(cmdout); error('x'); end
end



%%% Rewrite at plumb (and means for later visualization) --- 3drefit -deoblique
disp(' writing data (and means for later visualization)')
fPlumb = cell(size(fOrig));
if nFrame>1
    fPlumbAv = cell(size(fOrig));
end
cmd = {srcAfni};
for R = 1:numel(fOrig)
    % disp(['  file' num2str(i) '/' num2str(numel(fOrig))])

    fIn = fOrig{R};
    fOut = replace(replace(fOrig{R},'.nii.gz',''),dIn,dOut); if ~exist(fOut,'dir'); mkdir(fOut); end
    if nFrame>1
        fOut_plumb = fullfile(fOut,'setPlumb_volTs.nii.gz');
        fOut_av_plumb = fullfile(fOut,'av_setPlumb_volTs.nii.gz');
        fOut_av_oblique = fullfile(fOut,'av_runOblique_volTs.nii.gz'); clear fOut
    elseif any(ismember(dataType,'vol'))
        fOut_plumb = fullfile(fOut,'setPlumb_vol.nii.gz'); clear fOut
        % fOut_av_plumb = fullfile(fOut,'av_setPlumb_vol.nii.gz');
        % fOut_av_oblique = fullfile(fOut,'av_runOblique_vol.nii.gz'); clear fOut
    else
        dbstack; error(['don''t know which suffix to give' newline 'please set ''vol'' or ''volTs'' as dataType'])
    end

    cmd{end+1} = ['echo ''  ''file ' num2str(R) '/' num2str(numel(fOrig))];
    if forceRewriteAtPlumb || ~exist(fOut_plumb,'file')
        cmd{end+1} = ['cp ' fIn ' ' fOut_plumb];
        if verbose
            cmd{end+1} = ['3drefit -deoblique ' fOut_plumb];
        else
            cmd{end+1} = ['3drefit -deoblique ' fOut_plumb ' > /dev/null 2>&1'];
        end
    end
    fPlumb{R} = fOut_plumb;


    if nFrame>1 && (forceRewriteAtPlumb || ~exist(fOut_av_plumb,'file'))
        cmd{end+1} = '3dTstat -overwrite \';
        cmd{end+1} = '-mean \';
        cmd{end+1} = ['-prefix ' fOut_av_plumb ' \'];
        if verbose
            cmd{end+1} = [fOut_plumb ' > /dev/null 2>&1'];
        else
            cmd{end+1} = fOut_plumb;
        end
    end
    if nFrame>1
        fPlumbAv{R} = fOut_av_plumb;
    end
    
    if nFrame>1 && (~exist(fOut_av_oblique,'file') || forceRewriteAtPlumb)
        cmd{end+1} = '3dTstat -overwrite \';
        cmd{end+1} = '-mean \';
        cmd{end+1} = ['-prefix ' fOut_av_oblique ' \'];
        if verbose
            cmd{end+1} = [fIn ' > /dev/null 2>&1'];
        else
            cmd{end+1} = fIn;
        end
    end
    cmd{end+1} = ['echo ''   ''done'];
end

if length(cmd)>1
    cmd = strjoin(cmd,newline); % disp(cmd)
    [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
end



% %%% Rewrite at plumb (and means for later visualization) --- 3dWarpApply -nwarp IDENT
% disp(' writing data (and means for later visualization)')
% fPlumb = cell(size(fOrig));
% fPlumbAv = cell(size(fOrig));
% cmd = {srcAfni};
% for i = 1:numel(fOrig)
%     % disp(['  file' num2str(i) '/' num2str(numel(fOrig))])
% 
%     fIn = fOrig{i};
%     fOut = replace(replace(fOrig{i},'.nii.gz',''),dIn,dOut); if ~exist(fOut,'dir'); mkdir(fOut); end
%     if any(ismember(dataType,'volTs'))
%         fOut_plumb = fullfile(fOut,'setPlumb_volTs.nii.gz');
%         fOut_av_plumb = fullfile(fOut,'av_setPlumb_volTs.nii.gz');
%         fOut_av_oblique = fullfile(fOut,'av_runOblique_volTs.nii.gz'); clear fOut
%     elseif any(ismember(dataType,'vol'))
%         fOut_plumb = fullfile(fOut,'setPlumb_vol.nii.gz');
%         fOut_av_plumb = fullfile(fOut,'av_setPlumb_vol.nii.gz');
%         fOut_av_oblique = fullfile(fOut,'av_runOblique_vol.nii.gz'); clear fOut
%     else
%         dbstack; error(['don''t know which suffix to give' newline 'please set ''vol'' or ''volTs'' as dataType'])
%     end
% 
%     cmd{end+1} = ['echo ''  ''file ' num2str(i) '/' num2str(numel(fOrig))];
%     if forceRewriteAtPlumb || ...
%             any([~exist(fOut_plumb,'file') ~exist(fOut_av_plumb,'file') ~exist(fOut_av_oblique,'file')])
%         cmd{end+1} = '3dNwarpApply -overwrite -quiet \';
%         cmd{end+1} = ['-nwarp ''IDENT(' fIn ')'' \'];
%         cmd{end+1} = ['-prefix ' fOut_plumb ' \'];
%         cmd{end+1} = ['-source ' fIn];
% 
%         cmd{end+1} = '3dTstat -overwrite \';
%         cmd{end+1} = '-mean \';
%         cmd{end+1} = ['-prefix ' fOut_av_plumb ' \'];
%         cmd{end+1} = ['-source ' fOut_plumb];
% 
%         cmd{end+1} = '3dTstat -overwrite \';
%         cmd{end+1} = '-mean \';
%         cmd{end+1} = ['-prefix ' fOut_av_oblique ' \'];
%         cmd{end+1} = ['-source ' fIn];
% 
%         cmd{end+1} = ['echo ''   ''done'];
%     else
%         cmd{end+1} = ['echo ''   ''already done,skipping'];
%     end
%     fPlumb{i} = fOut_plumb;
%     fPlumbAv{i} = fOut_av_plumb;
% end
% 
% if length(cmd)>1
%     cmd = strjoin(cmd,newline); % disp(cmd)
%     [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
% end

% %%% Rewrite at plumb (and means for later visualization) --- MRIread MRIwrite
% mri_setPlumbRef = MRIread(fPlumbRef,1);
% fPlumb = cell(size(fOrig));
% fPlumbAv = cell(size(fOrig));
% for i = 1:numel(fOrig)
%     disp(['  file' num2str(i) '/' num2str(numel(fOrig))])
%     fIn = fOrig{i};
%     fOut = replace(replace(fOrig{i},'.nii.gz',''),dIn,dOut); if ~exist(fOut,'dir'); mkdir(fOut); end
%     if any(ismember(dataType,'volTs'))
%         fOut_plumb = fullfile(fOut,'setPlumb_volTs.nii.gz');
%         fOut_av_plumb = fullfile(fOut,'av_setPlumb_volTs.nii.gz');
%         fOut_av_oblique = fullfile(fOut,'av_runOblique_volTs.nii.gz'); clear fOut
%     elseif any(ismember(dataType,'vol'))
%         fOut_plumb = fullfile(fOut,'setPlumb_vol.nii.gz');
%         fOut_av_plumb = fullfile(fOut,'av_setPlumb_vol.nii.gz');
%         fOut_av_oblique = fullfile(fOut,'av_runOblique_vol.nii.gz'); clear fOut
%     else
%         dbstack; error(['don''t know which suffix to give' newline 'please set ''vol'' or ''volTs'' as dataType'])
%     end
%     if forceRewriteAtPlumb || ...
%             any([~exist(fOut_plumb,'file') ~exist(fOut_av_plumb,'file') ~exist(fOut_av_oblique,'file')])
%         mri_runOblique = MRIread(fIn);
%         mri_runOblique.vol = mri_runOblique.vol(:,:,:,param.nDummy+1:end);
%         mri_setPlumbRef.vol = mri_runOblique.vol;
%         MRIwrite(mri_setPlumbRef,fOut_plumb);
%         mri_setPlumbRef.vol = mean(mri_setPlumbRef.vol,4);
%         MRIwrite(mri_setPlumbRef,fOut_av_plumb);
%         mri_runOblique.vol = mean(mri_runOblique.vol,4);
%         MRIwrite(mri_runOblique,fOut_av_oblique);
%         disp('   done')
%     else
%         disp('   already done, skipping')
%     end
%     fPlumb{i} = fOut_plumb;
%     fPlumbAv{i} = fOut_av_plumb;
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare data for motion estimation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Multi-echo (cross-echo rms images for later motion/distortion estimation)

%%% Detect multi-echo data if not specified in dataType
nEcho = contains(bidsList(:,1),'echo-'); nEcho = unique(bidsList(nEcho,:));
if isempty(nEcho)
    nEcho = 1;
    dataType{end+1} = 'singleEcho';
else
    nEcho = length(nEcho);
    dataType{end+1} = 'multiEcho';
end
disp([dataType{end} ' data detected'])
% if ~any(ismember({'multiEcho' 'singleEcho'},dataType))
%     tmp = bidsList(contains(bidsList(:,1,1),'echo-'),:,:);
%     if ~isempty(tmp) && length(unique(tmp))>1
%     % if nnz(contains(fPlumb,'echo-'))==length(fPlumb)
%         dataType{end+1} = 'multiEcho';
%     else
%         dataType{end+1} = 'singleEcho';
%     end
%     disp([dataType{end} ' data detected'])
% else
%     disp([dataType{ismember(dataType,{'multiEcho' 'singleEcho'})} ' data specified'])
% end

nRun = contains(squeeze(bidsList(:,1,1)),'run-');
if ~any(nRun); nRun = 1; else; nRun = length(unique(bidsList(nRun,:))); end

if nRun==1 && all(ismember({'singleEcho' 'pc' 'vol'},dataType))
    nVenc = length(fPlumb)/nEcho/nRun;

    fPlumb = reshape(fPlumb,nRun,nEcho,nVenc);
    if nFrame>1
        fPlumbAv = reshape(fPlumbAv,nRun,nEcho,nVenc);   
    end
    fOrig = reshape(fOrig,nRun,nEcho,nVenc);
    fEstim = fPlumb(contains(fPlumb,'proc-venc0_') & contains(fPlumb,'part-mag'));
    if size(fEstim,1)~=size(fOrig,1); dbstack; error('double-check that'); end
    % bidsList;
    
    

elseif all(ismember({'singleEcho'},dataType))
    %%% Single echo data
    %%%% Write means (no rms) for later visualization
    if nRun>1
        cmd = {srcFs};
        fIn = fPlumbAv;
        if any(contains(squeeze(bidsList(:,1,1)),'run-'))
            fOut = replace(fIn{1},bidsList{contains(squeeze(bidsList(:,1,1)),'run-'),1,1},'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        else
            fOut = fIn{1};
            if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        end
        fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
        if forceRewriteAtPlumb || ~exist(fOut,'file')
            cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
        end
        fPlumbCatAv = fOut;
        fIn = fOut;
        fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
        if forceRewriteAtPlumb || ~exist(fOut,'file')
            cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
        end
        fPlumbAvCatAv = fOut;

        disp(' averaging')
        if length(cmd)>1
            cmd = strjoin(cmd,newline); % disp(cmd)
            [status,cmdout] = system(cmd); if status; dbstack; error(cmdout); error('x'); end
            disp('  done')
        else
            disp('  already done, skipping')
        end
    end
    
    %%% Set fEstim
    fEstim = fPlumb;
    fEstimAv = fPlumbAv;
    if nRun>1
        fEstimCatAv = fPlumbCatAv;
        fEstimAvCatAv = fPlumbAvCatAv;
    end

elseif all(ismember({'multiEcho'},dataType))
    
    %%% Multi-echo data
    
    %%%% Reshape lists with echos in different columns
    echoList = bidsList(contains(bidsList(:,1,1),'echo-'),:,:);
    echoUniqueList = unique(echoList);
    
    bidsList = permute(reshape(bidsList,size(bidsList,1),length(echoUniqueList),length(echoList)/length(echoUniqueList)),[1 3 2]);
    fPlumb = reshape(fPlumb,length(echoUniqueList),length(echoList)/length(echoUniqueList))';
    if nFrame>1
        fPlumbAv = reshape(fPlumbAv,length(echoUniqueList),length(echoList)/length(echoUniqueList))';
    end
    fOrig = reshape(fOrig,length(echoUniqueList),length(echoList)/length(echoUniqueList))';

    disp('---------------------------------------------------')
    disp('---------------------------------------------------')
    warning(['Make sure echoes are properly sorted below' ...
        newline '(runs across rows and echos across columns)'])
    if size(bidsList,2)>1
        disp(fullfile(squeeze(bidsList(contains(bidsList(:,1,1),'run-'),:,:)),squeeze(bidsList(contains(bidsList(:,1,1),'echo-'),:,:))))
    else
        disp(fullfile(permute((bidsList(contains(bidsList(:,1,1),'echo-'),:,:)),[1 3 2])));
    end
    disp('---------------------------------------------------')
    disp('---------------------------------------------------')

    %%%% Write means (before rms) for later visualization
    if nRun>1 || nFrame>1
        cmd = {srcFs};
        fPlumbCatAv = cell(1,size(fPlumbAv,2));
        fPlumbAvCatAv = cell(1,size(fPlumbAv,2));
        if nRun>1
            for E = 1:nEcho
                fIn = fPlumbAv(:,E);
                fOut = replace(fIn{1},bidsList{contains(squeeze(bidsList(:,1,1)),'run-'),1,1},'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
                fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
                if forceRecomputeRMS || ~exist(fOut,'file')
                    cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
                end
                fPlumbCatAv{E} = fOut;
                fIn = fOut;
                fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
                if forceRecomputeRMS || ~exist(fOut,'file')
                    cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
                end
                fPlumbAvCatAv{E} = fOut;
            end
            disp(' averaging')
            if length(cmd)>1
                cmd = strjoin(cmd,newline); % disp(cmd)
                [status,cmdout] = system(cmd); if status; dbstack; error(cmdout); error('x'); end
                disp('  done')
            else
                disp('  already done, skipping')
            end
            disp(' non need for averaging (single-run)')
        end
    end

    %%%% Compute and write rms
    disp('computing cross-echo rms')
    fEstim = cell(size(fPlumb,1),1);
    for R = 1:size(fPlumb,1)
        disp([' file' num2str(R) '/' num2str(size(fPlumb,1))])
        fOut = replace(fPlumb{R,1},'echo-1','echo-rms');
        if forceRecomputeRMS || ~exist(fOut,'file')
            [a,~] = fileparts(fOut); if ~exist(a,'dir'); mkdir(a); end

            mri = MRIread(fPlumb{R,1},1);
            mri.vol = zeros([mri.volsize mri.nframes]);
            for ii = 1:size(fPlumb,2)
                mriTmp = MRIread(fPlumb{R,ii});
                mri.vol = mri.vol + mriTmp.vol.^2;
            end
            mri.vol = sqrt(mri.vol./size(fPlumb,2));

            MRIwrite(mri,fOut);
            disp('  done')
        else
            disp('  already done, skipping')
        end
        fEstim{R} = fOut;
    end


    %%%% Write means (after rms) for later visualization
    cmd = {srcFs};
    if nFrame>1
        fEstimAv = cell(size(fEstim));
    end
    for R = 1:length(fEstim)
        fIn = fEstim{R};
        if nFrame>1
            fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
            if forceRecomputeRMS || ~exist(fOut,'file')
                cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
            end
            fEstimAv{R} = fOut;
        % else
        %     fEstimAv{R} = fIn;
        end
    end
    if nRun>1
        fIn = fEstimAv;
        fOut = replace(fIn{1},bidsList{contains(squeeze(bidsList(:,1,1)),'run-'),1,1},'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
        if forceRecomputeRMS || ~exist(fOut,'file')
            cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
        end
        fEstimCatAv = fOut;
        fIn = fOut;
        fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
        if forceRecomputeRMS || ~exist(fOut,'file')
            cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
        end
        fEstimAvCatAv = fOut;
    else
        fEstimCatAv = '';
        fEstimAvCatAv = '';
    end
    disp(' averaging')
    if length(cmd)>1
        cmd = strjoin(cmd,newline); % disp(cmd)
        [status,cmdout] = system(cmd); if status; dbstack; error(cmdout); error('x'); end
    else
        disp('  already done, skipping')
    end
end



%% Data on which to apply motion correction
fApply = fPlumb;



%% QA: Visualize motion
%%% between-run
if nRun>1 && ~isempty(fEstimCatAv)
    cmd = {srcFs};
    cmd{end+1} = ['fslview -m single ' fEstimCatAv ' &'];
    cmd = strjoin(cmd,newline); % disp(cmd)
    if param.verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    fFslviewBR = cmd;
% else
%     fFslviewBR = '';
end
%%% within-run
if nFrame>1
    cmd = {srcFs};
    cmd{end+1} = ['fslview -m single ' strjoin(fEstim,' ') ' &'];
    cmd = strjoin(cmd,newline); % disp(cmd)
    if param.verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    fFslviewWR = cmd;
    %%% within-run, first, middle and last frames
    fFslviewWRfstMdLst = qaFstMdLst(fEstim,forceRewriteAtPlumb || forceRecomputeRMS,param.verbose);
% else
%     fFslviewWR = '';
%     fFslviewWRfstMdLst = '';
end

if nRun>1
    xSet.initFiles.qaFiles.fFslviewBR = fFslviewBR; end
if nFrame>1
    xSet.initFiles.qaFiles.fFslviewWR = fFslviewWR;
    xSet.initFiles.qaFiles.fFslviewWRfstMdLst = fFslviewWRfstMdLst; end



%% Low SNR data (estimate motion on temporally smoothed data)


%%% Detect low SNR data it not specified in dataType
if any(ismember({'volTs'},dataType))
    if ~any(ismember({'lowSNR' 'highSNR'},dataType))
        dbstack; error(['automatic detection of lowSNR data not implemented' newline 'specify ''lowSNR'' or ''highSNR'' in dataType '])
    else
        disp([dataType{ismember(dataType,{'lowSNR' 'highSNR'})} ' data specified'])
    end
end


%%% Temporally smooth lowSNR data
if all(ismember({'volTs' 'lowSNR'},dataType))

    %%%% Visualize SNR before smoothing
    fEstimBeforeSm = fEstim{1};
    volTs = MRIread(fEstimBeforeSm);
    volAv = mean(volTs.vol,4);
    volEr = std(volTs.vol,[],4);
    if param.verbose
        hFig = figure('WindowStyle','docked');
    else
        hFig = figure('Visible','off','units','normalized','outerposition',[0 0 1 1]);
    end
    hTile = tiledlayout(2,3); hTile.TileSpacing = 'tight'; hTile.Padding = 'tight';
    ax = {};
    ax{end+1} = nexttile;
    imagesc(volAv(:,:,ceil(end/2)));
    ax{end}.DataAspectRatio = [1 1 1];
    ax{end}.XTick = []; ax{end}.YTick = [];
    ylabel(colorbar,'mean'); ax{end}.Colormap = gray;
    ylabel('without smoothing')
    ax{end+1} = nexttile;
    imagesc(volEr(:,:,ceil(end/2)));
    ax{end}.DataAspectRatio = [1 1 1];
    ax{end}.XTick = []; ax{end}.YTick = [];
    ylabel(colorbar,'std'); ax{end}.Colormap = gray;
    ax{end+1} = nexttile;
    imagesc(volAv(:,:,ceil(end/2))./volEr(:,:,ceil(end/2)));
    ax{end}.DataAspectRatio = [1 1 1];
    ax{end}.XTick = []; ax{end}.YTick = [];
    ylabel(colorbar,'tSNR'); ax{end}.Colormap = gray;
    if isempty(param) || ~isfield(param,'tSmWin_vol') || isempty(param.tSmWin_vol)
        dbstack; error('when dataType is lowSNR, must specify param.tSmWin_vol')
    end
    
    %%%% Smoothing
    disp('temporally smooth timeseries')
    cmd = {srcAfni};
    for R = 1:length(fEstim)
        fIn = fEstim{R};
        [d,fOut,~] = fileparts(replace(fIn,'.nii.gz',''));
        fOut = fullfile(d,[strjoin({['sm' num2str(param.tSmWin_vol)] fOut},'_') '.nii.gz']);
        cmd{end+1} = ['echo '' ''file' num2str(R) '/' num2str(length(fEstim))];
        if forceResmoothing || ~exist(fOut,'file')
            cmd{end+1} = '3dTsmooth -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = ['-hamming ' num2str(param.tSmWin_vol) ' \'];
            if verbose
                cmd{end+1} = fIn;
            else
                cmd{end+1} = [fIn ' > /dev/null 2>&1'];
            end
        else
            cmd{end+1} = ['echo ''  ''already smoothed, skipping'];
        end
        fEstim{R} = fOut;
    end
    cmd = strjoin(cmd,newline); % disp(cmd)
    [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end

    %%%% Visualize SNR after smoothing
    fEstimAfterSm = fEstim{1};
    volTs = MRIread(fEstimAfterSm);
    volAv = mean(volTs.vol,4);
    volEr = std(volTs.vol,[],4);
    ax{end+1} = nexttile;
    imagesc(volAv(:,:,ceil(end/2)));
    ax{end}.DataAspectRatio = [1 1 1];
    ax{end}.XTick = []; ax{end}.YTick = [];
    ylabel(colorbar,'mean'); ax{end}.Colormap = gray;
    ylabel(['with smoothing (' num2str(param.tSmWin_vol) 'volumes)'])
    ax{end+1} = nexttile;
    imagesc(volEr(:,:,ceil(end/2)));
    ax{end}.DataAspectRatio = [1 1 1];
    ax{end}.XTick = []; ax{end}.YTick = [];
    ylabel(colorbar,'std'); ax{end}.Colormap = gray;
    ax{end+1} = nexttile;
    imagesc(volAv(:,:,ceil(end/2))./volEr(:,:,ceil(end/2)));
    ax{end}.DataAspectRatio = [1 1 1];
    ax{end}.XTick = []; ax{end}.YTick = [];
    ylabel(colorbar,'tSNR'); ax{end}.Colormap = gray;
    for R = 1:3
        cLim = [min([ax{R}.CLim(1) ax{R+3}.CLim(1)]) max([ax{R}.CLim(2) ax{R+3}.CLim(2)])];
        ax{R}.CLim = cLim;
        ax{R+3}.CLim = cLim;
    end
    [~,b,~] = fileparts(fileparts(fEstimAfterSm));
    title(hTile,b,'interpreter','none')
    fFigSm = replace(fEstimAfterSm,'.nii.gz','.fig');
    set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    saveas(hFig,fFigSm)
    
    %%%% Visualize SNR before and after smoothing
    cmd = {srcFs};
    cmd{end+1} = ['fslview -m single ' fEstimBeforeSm ' ' fEstimAfterSm ' &'];
    cmd = strjoin(cmd,newline); % disp(cmd)
    fFslviewSm = cmd;
    if param.verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end

    %%%% Write means for later visualization
    cmd = {srcFs};
    fEstimAv = cell(size(fEstim));
    for R = 1:length(fEstim)
        fIn = fEstim{R};
        fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
        if forceResmoothing || ~exist(fOut,'file')
            cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
        end
        fEstimAv{R} = fOut;
    end
    fIn = fEstimAv;
    fOut = replace(fIn{1},bidsList{contains(squeeze(bidsList(:,1,1)),'run-'),1,1},'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
    fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
    if forceResmoothing || ~exist(fOut,'file')
        cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
    end
    fIn = fOut;
    fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
    if forceResmoothing || ~exist(fOut,'file')
        cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
    end
    disp(' averaging')
    if length(cmd)>1
        cmd = strjoin(cmd,newline); % disp(cmd)
        [status,cmdout] = system(cmd); if status; dbstack; error(cmdout); error('x'); end
    else
        disp('  already done, skipping')
    end

    %% QA: Visualize motion (with smoothing)
    %%% between-run
    cmd = {srcFs};
    cmd{end+1} = ['fslview -m single ' fEstimCatAv ' &'];
    cmd = strjoin(cmd,newline); % disp(cmd)
    if param.verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    fFslviewBRsm = cmd;
    %%% within-run
    cmd = {srcFs};
    cmd{end+1} = ['fslview -m single ' strjoin(fEstim,' ') ' &'];
    cmd = strjoin(cmd,newline); % disp(cmd)
    if param.verbose
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
    fFslviewWRsm = cmd;
    %%% within-run, first, middle and last frames
    fFslviewWRsmFstMdLst = qaFstMdLst(fEstim,forceRewriteAtPlumb||forceRecomputeRMS,param.verbose);

    xSet.initFiles.qaFilesSm.fFslviewBR = fFslviewBRsm;
    xSet.initFiles.qaFilesSm.fFslviewWR = fFslviewWRsm;
    xSet.initFiles.qaFilesSm.fFslviewWRfstMdLst = fFslviewWRsmFstMdLst;


end


%% Sumarize outputs
if ~exist('nRun','var') || isempty(nRun)
    nRun = size(fOrig,1); end
if ~exist('nEcho','var') || isempty(nEcho)
    nEcho = size(fOrig,2); end
if any(ismember(dataType,'pc')) && (~exist('nVenc','var') || isempty(nVenc))
    nVenc = size(fOrig,3); end

xSet.dataType = dataType;
xSet.nRun = nRun;
xSet.nEcho = nEcho;
if exist('nVenc','var')
    xSet.nVenc = nVenc; end
xSet.nFrame = nFrame;
xSet.initFiles.nRun = nRun;
xSet.initFiles.nEcho = nEcho;
if exist('nVenc','var')
    xSet.initFiles.nVenc = nVenc; end
xSet.initFiles.nFrame = nFrame;
xSet.initFiles.fOrig = fOrig;

xSet.initFiles.fPlumbRef = fPlumbRef;
xSet.initFiles.fPlumb = fPlumb;
% if exist('nVenc','var')
if nFrame>1
    xSet.initFiles.fPlumbAv = fPlumbAv; end
if nRun>1
    xSet.initFiles.fPlumbCatAv = fPlumbCatAv;
    xSet.initFiles.fPlumbAvCatAv = fPlumbAvCatAv; end

xSet.initFiles.fEstim = fEstim;
if nFrame>1
    xSet.initFiles.fEstimAv = fEstimAv; end
if nRun>1
    if iscell(fEstimCatAv)
        xSet.initFiles.fEstimCatAv = fEstimCatAv;
    else
        xSet.initFiles.fEstimCatAv = {fEstimCatAv};
    end
    if iscell(fEstimAvCatAv)
        xSet.initFiles.fEstimAvCatAv = fEstimAvCatAv;
    else
        xSet.initFiles.fEstimAvCatAv = {fEstimAvCatAv};
    end
end

xSet.initFiles.fApply = fApply;

xSet.initFiles.fObliqueRef = fObliqueRef;
xSet.initFiles.fPlumbRef = fPlumbRef;

%%% QA
if exist('fFigSm','var')
    xSet.initFiles.qaFilesSm.fFigSm = fFigSm;
end
if exist('fFslviewSm','var')
    xSet.initFiles.qaFilesSm.fFslviewSm = fFslviewSm;
end
