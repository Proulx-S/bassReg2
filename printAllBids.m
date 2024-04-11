function fGeom = printAllBids(bidsDir,iGeom)

% get all data from bids directory
tmp = dir(fullfile(bidsDir,'anat','*.nii.gz'));
tmp = cat(1,tmp,dir(fullfile(bidsDir,'func','*.nii.gz')));
tmp = cat(1,tmp,dir(fullfile(bidsDir,'fmap','*.nii.gz')));
tmp = fullfile({tmp.folder},{tmp.name})';

% for multiecho data, keep only the first echo
tmpE = tmp(contains(tmp,'echo-'));
tmp = tmp(~contains(tmp,'echo-'));
tmp = cat(1,tmp,tmpE(contains(tmpE,'echo-1'))); clear tmpE

% sort according to acquisition time
[~,b] = sort(getAcqTime(tmp));
fTmp = tmp(b);

% prepare nice string matrix for printing
[~,tmp] = fileparts(replace(fTmp,'.nii.gz',''));
tmp = char(tmp);
tmp = cat(2,num2str((1:size(tmp,1))','%i-----'),tmp);

% identify seleted file (here geom file to be used as the reference slice prescription in scanner space)
if exist('iGeom','var') && ~isempty(iGeom)
    fGeom = fTmp{iGeom};
    tmp(iGeom,:) = replace(tmp(iGeom,:),'-----','--G--');
else
    fGeom = [];
end

% print
disp(tmp)
