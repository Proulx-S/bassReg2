function runSet = writeObliqueAndPlumb(runSet,force,verbose)
global srcAfni
if ~exist('force','var'); force = []; end
if ~exist('verbose','var'); verbose = []; end
if isempty(force); force = 0; end
if isempty(verbose); verbose = 0; end

duporigin = 1;
cmd = {srcAfni};

if size(runSet.nDummy,1)==1
    runSet.nDummy = repmat(runSet.nDummy,size(runSet.fOrigList));
end
if size(runSet.nDummy,1)~=size(runSet.fOrigList,1)
    dbstack; error('somthing wrong')
end

%% Individual files
runSet.fPlumbList = cell(size(runSet.fOrigList));
for r = 1:length(runSet.fOrigList)
    cmd{end+1} = ['echo '' ''' num2str(r) '/' num2str(length(runSet.fOrigList))];
    fIn = runSet.fOrigList{r};
    [~,fOut,~] = fileparts(replace(fIn,'.nii.gz','')); fOut = fullfile(runSet.wd,fOut); if ~exist(fOut,'dir'); mkdir(fOut); end
    fOut = fullfile(fOut,'setPlumb_volTs.nii.gz');
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dcalc -overwrite \';
        cmd{end+1} = ['-a ' fIn '[' num2str(runSet.nDummy(r)) '..$] \'];
        cmd{end+1} = ['-expr a \'];
        if verbose
            cmd{end+1} = ['-prefix ' fOut];
        else
            cmd{end+1} = ['-prefix ' fOut ' > /dev/null 2>&1'];
        end
    end
    fIn = fOut;
    fPlumbRef = fullfile(runSet.wd,'setPlumb_volRef.nii.gz');
    if force || ~exist(fOut,'file')
        if duporigin
            cmd{end+1} = ['3drefit -duporigin ' fPlumbRef ' -deoblique ' fOut];
        else
            cmd{end+1} = ['3drefit -deoblique ' fOut];
        end
        if ~verbose
            cmd{end} = [cmd{end} ' > /dev/null 2>&1'];
        end
    end
    
    runSet.fPlumbList{r} = fOut;
end
runSet.nFrame = runSet.nFrame-runSet.nDummy;

%% Launch command
if length(cmd)>1
    if verbose
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
    else
        [status,cmdout] = system(strjoin(cmd,newline)); if status; dbstack; error(cmdout); error('x'); end
    end
end


return





%%% Rewrite at plumb (and means for later visualization) --- 3drefit -deoblique
disp(' writing data (and means for later visualization)')
fPlumb = cell(size(fOrigList));
if nFrame>1
    fPlumbAv = cell(size(fOrigList));
end
cmd = {srcAfni};
for r = 1:numel(fOrigList)
    % disp(['  file' num2str(i) '/' num2str(numel(fOrig))])

    fIn = fOrigList{r};
    fOut = replace(replace(fOrigList{r},'.nii.gz',''),dIn,dOut); if ~exist(fOut,'dir'); mkdir(fOut); end
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

    cmd{end+1} = ['echo ''  ''file ' num2str(r) '/' num2str(numel(fOrigList))];
    if force || ~exist(fOut_plumb,'file')
        cmd{end+1} = ['cp ' fIn ' ' fOut_plumb];
        if verbose
            if duporigin
                cmd{end+1} = ['3drefit -duporigin ' fPlumbRef ' -deoblique ' fOut_plumb];
            else
                cmd{end+1} = ['3drefit -deoblique ' fOut_plumb];
            end
        else
            if duporigin
                cmd{end+1} = ['3drefit -duporigin ' fPlumbRef ' -deoblique ' fOut_plumb ' > /dev/null 2>&1'];
            else
                cmd{end+1} = ['3drefit -deoblique ' fOut_plumb ' > /dev/null 2>&1'];
            end
        end
    end
    fPlumb{r} = fOut_plumb;


    if nFrame>1 && (force || ~exist(fOut_av_plumb,'file'))
        cmd{end+1} = '3dTstat -overwrite \';
        cmd{end+1} = '-mean \';
        cmd{end+1} = ['-prefix ' fOut_av_plumb ' \'];
        if verbose
            cmd{end+1} = fOut_plumb;
        else
            cmd{end+1} = [fOut_plumb ' > /dev/null 2>&1'];
        end
    end
    if nFrame>1
        fPlumbAv{r} = fOut_av_plumb;
    end
    
    if nFrame>1 && (~exist(fOut_av_oblique,'file') || force)
        cmd{end+1} = '3dTstat -overwrite \';
        cmd{end+1} = '-mean \';
        cmd{end+1} = ['-prefix ' fOut_av_oblique ' \'];
        if verbose
            cmd{end+1} = fIn;
        else
            cmd{end+1} = [fIn ' > /dev/null 2>&1'];
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

