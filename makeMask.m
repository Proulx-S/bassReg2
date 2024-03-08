function files = makeMask(files,bidsDir,force)
global srcFs srcAfni
if ~exist('bidsDir','var'); bidsDir = []; end
if ~exist('force','var'); force = []; end
if isempty(force); force = 0; end

fIn = files.fEstimCatAv;
fOut = replace(fIn,'.nii.gz','-mask.nii.gz');
% save in bids deriv dir for persistence since this is manually drawn
if ~isempty(bidsDir)
    if ~exist(bidsDir,'dir'); mkdir(bidsDir); end
    files.manBrainMask = fullfile(bidsDir,'manBrainMask.nii.gz');
else
    files.manBrainMask = fullfile(fileparts(fOut),'manBrainMask.nii.gz');
end

%% Create mask manually
if force || ~exist(files.manBrainMask,'file')
    warning('draw mask, save with default name and close fslview')
    cmd = {srcFs};
    cmd{end+1} = ['fslview -m single ' fIn];
    cmd = strjoin(cmd,newline); % disp(cmd)
    [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
    movefile(fOut,files.manBrainMask)
end

%% Invert mask
files.manBrainMaskInv = replace(files.manBrainMask,'.nii.gz','Inv.nii.gz');
if force || ~exist(files.manBrainMaskInv,'file')
    cmd = {srcAfni};
    cmd{end+1} = '3dcalc -overwrite \';
    cmd{end+1} = ['-prefix ' files.manBrainMaskInv ' \'];
    cmd{end+1} = ['-a ' files.manBrainMask ' \'];
    cmd{end+1} = '-expr ''-(a-1)''';
    cmd = strjoin(cmd,newline); % disp(cmd)
    [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
end

%% Add mask to qa
files = addMaskToCmd(files);
