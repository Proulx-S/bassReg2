function surfReg(SUBJECTS_DIR,SUBJECT,MOV,outDir,contrast,force)
global srcFs
if ~exist('force','var');       force = [];     end
if ~exist('contrast','var'); contrast = [];     end
if isempty(force);              force = 0;      end
if isempty(contrast);        contrast = '--t1'; end % --bold, --dti, --t2, or --t1



% Grrrrrrrr, not enough white matter contrast!!!!



if ~isFS(fullfile(SUBJECTS_DIR,SUBJECT)); dbstack; error('bad fsDir'); end
if exist('outDir','var')
    if ~exist(outDir,'dir'); mkdir(outDir); end 
else
    dbstack; error('code that')
end






o   = fullfile(outDir,'toSurf.nii.gz');
reg = fullfile(outDir,'toSurf.dat');
lta = fullfile(outDir,'toSurf.lta');

cmd = {srcFs};
cmd{end+1,1} = ['SUBJECTS_DIR=' SUBJECTS_DIR];
cmd{end+1} = 'bbregister \';
cmd{end+1} = ['--s '   SUBJECT ' \'];
cmd{end+1} = ['--mov ' MOV ' \'];
% cmd{end+1} = ['--o '   o ' \'];
cmd{end+1} = ['--reg ' reg ' \'];
cmd{end+1} = ['--lta ' lta ' \'];
cmd{end+1} = '--init-header \';
cmd{end+1} = contrast;

if force || ~exist(o,'file') || ~exist(reg,'file') || ~exist(lta,'file')
    [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
end

tkregisterfv --mov /autofs/space/takoyaki_001/users/proulxs/vsmDriven/doIt_vsmDriven3/bids/sub-vsmDrivenP1/ses-1/anat/sub-vsmDrivenP1_ses-1_acq-avMap_run-5_echo-1_angio.nii.gz --reg /autofs/space/takoyaki_001/users/proulxs/vsmDriven/doIt_vsmDriven3/bids/sub-vsmDrivenP1/ses-1/derivatives/set-avMap/sub-vsmDrivenP1_ses-1_acq-avMap_run-5_echo-1_angio/toSurf.lta --surfs  --sd /autofs/space/takoyaki_001/users/proulxs/vsmDriven/doIt_vsmDriven3/bids/sub-vsmDrivenP1/ses-1/derivatives/set-fs

cmd = {srcFs};
cmd{end+1,1} = 'mri_vol2vol \';
cmd{end+1}   = ['--mov '  MOV ' \'];
cmd{end+1}   = ['--targ ' MOV ' \'];
cmd{end+1}   = ['--o   '  o   ' \'];
cmd{end+1}   = ['--lta '  lta ' \'];
cmd{end+1}   = '--no-resample';



cmd = {srcFs};
cmd{end+1,1} = 'freeview \';
cmd{end+1} = ['-v ' o ''];


cmd = {srcFs};
cmd{end+1,1} = 'freeview \';
cmd{end+1} = ['-v ' o                                             '                \'];
cmd{end+1} = ['-v ' fullfile(SUBJECTS_DIR,SUBJECT,'mri','nu.mgz') ':resample=cubic \'];
cmd{end+1} = ['-f ' fullfile(SUBJECTS_DIR,SUBJECT,'surf','lh.pial')];
clipboard('copy',strjoin(cmd,newline))

cmd = {srcFs};
cmd{end+1,1} = 'freeview \';
cmd{end+1} = ['-v ' MOV                                           ':resample=cubic \'];
cmd{end+1} = ['-v ' fullfile(SUBJECTS_DIR,SUBJECT,'mri','nu.mgz') ':resample=cubic \'];
cmd{end+1} = ['-f ' fullfile(SUBJECTS_DIR,SUBJECT,'surf','lh.pial')];
clipboard('copy',strjoin(cmd,newline))



cmd = {srcFs};
cmd{end+1,1} = 'freeview \';
% cmd{end+1} = ['-v ' o                                                           ':reg=' reg '  \'];
cmd{end+1} = ['-v ' MOV                                           ':reg=' reg ' \'];
% cmd{end+1} = ['-v ' fullfile(SUBJECTS_DIR,SUBJECT,'mri','nu.mgz') ':resample=cubic \'];
cmd{end+1} = ['-f ' fullfile(SUBJECTS_DIR,SUBJECT,'surf','lh.pial')];
clipboard('copy',strjoin(cmd,newline))

freeview
MOV
o
fullfile(SUBJECTS_DIR,SUBJECT,'mri','nu.mgz')
fullfile(SUBJECTS_DIR,SUBJECT,'surf','lh.pial')
fullfile(SUBJECTS_DIR,SUBJECT,'surf','lh.white')
fullfile(SUBJECTS_DIR,SUBJECT,'surf','rh.pial')
fullfile(SUBJECTS_DIR,SUBJECT,'surf','rh.white')




