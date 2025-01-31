function runSet = writeObliqueAndPlumbRef(runSet,force,verbose)
global srcAfni
if ~exist('force','var'); force = []; end
if ~exist('verbose','var'); verbose = []; end
if isempty(force); force = 0; end
if isempty(verbose); verbose = 0; end

duporigin = 1;
cmd = {srcAfni};

%% Geom
fIn = runSet.fGeom;
nFrame = MRIget(fIn,'nframes');
if nFrame>8; nFrameRef = 8; else nFrameRef = nFrame; end
fObliqueRef = fullfile(runSet.wd,'setOblique_volRef.nii.gz');
if force || ~exist(fObliqueRef,'file')
    cmd{end+1} = '3dcalc -overwrite \';
    cmd{end+1} = ['-a ' fIn '[0..' num2str(nFrameRef-1) '] \'];
    cmd{end+1} = ['-expr a \'];
    if verbose
        cmd{end+1} = ['-prefix ' fObliqueRef];
    else
        cmd{end+1} = ['-prefix ' fObliqueRef ' > /dev/null 2>&1'];
    end
end
fPlumbRef = fullfile(runSet.wd,'setPlumb_volRef.nii.gz');
if force || ~exist(fPlumbRef,'file')
    cmd{end+1} = ['cp ' fObliqueRef ' ' fPlumbRef];
    cmd{end+1} = ['3drefit -deoblique ' fPlumbRef];
end

%% Individual files
nFrameList = zeros([1 size(runSet.fOrigList)]);
vSizeList  = zeros([3 size(runSet.fOrigList)]);
for r = 1:numel(runSet.fOrigList)
    cmd{end+1} = ['echo '' ''' num2str(r) '/' num2str(numel(runSet.fOrigList))];
    fIn = runSet.fOrigList{r};
    [~,fOut,~] = fileparts(replace(fIn,'.nii.gz','')); fOut = fullfile(runSet.wd,fOut); if ~exist(fOut,'dir'); mkdir(fOut); end
    fOut = fullfile(fOut,'runOblique_volRef.nii.gz');
    nFrame = MRIget(fIn,'nframes'); if nFrame>8; nFrameRef = 8; else; nFrameRef = nFrame; end
    vSize = MRIget(fIn,'voxSize');
    if force || ~exist(fOut,'file')
        % mri = MRIread(fIn,1);
        % nFrame = mri.nframes; if nFrame>8; nFrameRef = 8; else; nFrameRef = nFrame; end
        % vSize  = [mri.xsize mri.ysize mri.zsize];
        cmd{end+1} = '3dcalc -overwrite \';
        cmd{end+1} = ['-a ' fIn '[0..' num2str(nFrameRef-1) '] \'];
        cmd{end+1} = ['-expr a \'];
        if verbose
            cmd{end+1} = ['-prefix ' fOut];
        else
            cmd{end+1} = ['-prefix ' fOut ' > /dev/null 2>&1'];
        end
    %     save(fullfile(fileparts(fOut),'imParam.mat'),'nFrame','vSize');
    % else
    %     load(fullfile(fileparts(fOut),'imParam.mat'),'nFrame','vSize');
    end
    nFrameList(1,r) = nFrame;
    vSizeList(:,r)  = vSize;

    fIn = fOut;
    fOut = replace(fOut,'runOblique_volRef.nii.gz','runPlumb_volRef.nii.gz');
    if force || ~exist(fOut,'file')
        cmd{end+1} = ['cp ' fIn ' ' fOut];
        if duporigin
            cmd{end+1} = ['3drefit -duporigin ' fPlumbRef ' -deoblique ' fOut];
        else
            cmd{end+1} = ['3drefit -deoblique ' fOut];
        end
        if ~verbose
            cmd{end} = [cmd{end} ' > /dev/null 2>&1'];
        end
    end
end

runSet.nFrame = permute(nFrameList,[2 3 1]); runSet.nFrame = runSet.nFrame(:,1,:);
runSet.vSize  = permute(vSizeList ,[2 3 1]); runSet.vSize  = runSet.vSize(:,1,:);

% runSet.nFrame = reshape(nFrameList,size(runSet.fList));                      runSet.nFrame = runSet.nFrame(:,1);
% runSet.vSize  = permute(reshape(vSizeList',[3 size(runSet.fList)]),[2 3 1]); runSet.vSize  = permute(runSet.vSize(:,1,:),[2 3 1]);
if any(diff(runSet.nFrame)); warning(['runs in ' runSet.label ' set don''t all have the same number of frame']); end



%% launch command
if length(cmd)>1
    if verbose
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
    else
        [status,cmdout] = system(strjoin(cmd,newline)); if status; dbstack; error(cmdout); error('x'); end
    end
end

