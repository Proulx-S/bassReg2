function filesPcVel = PCphase2vel(files)

fList = files.f(:,:,contains(files.f(1,1,:),'part-phase'));
vencList = cell(size(fList));
for i = 1:numel(fList)
    vencList{i} = strsplit(fList{i},filesep); vencList{i} = strsplit(vencList{i}{end-2},'_'); vencList{i} = vencList{i}{contains(vencList{i},'proc-venc')};
    vencList{i} = strsplit(vencList{i},'-'); vencList{i} = str2double(replace(replace(vencList{i}{end},'venc',''),'m0',''));
end
vencList = cell2mat(vencList);

filesPcVel = files;
filesPcVel.f = cell(size(fList));

for i = 1:numel(fList)
    f = fList{i};
    fVel = replace(replace(f,'_vol.nii.gz','_vel.nii.gz'),'_volTs.nii.gz','_velTs.nii.gz');

    mri = MRIread(f);
    mri.vol = mri.vol./4096*vencList(i);
    
    % figure('WindowStyle','docked');
    % imagesc(mri.vol,[-10 10]);
    % ax = gca; ax.Colormap = turbo; ax.DataAspectRatio = [1 1 1]; ylabel(colorbar,'satin index')
    
    MRIwrite(mri,fVel);
    filesPcVel.f{i} = fVel;
end


