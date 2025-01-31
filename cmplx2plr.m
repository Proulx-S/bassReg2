function [fMag,fPhs] = cmplx2plr(fR,fI,force)
if ~exist('force','var'); force = []; end
if isempty(force);        force = 0 ; end

fPhs = cell(size(fR));
fMag = cell(size(fR));
for r = 1:length(fR)
    disp([num2str(r) '/' num2str(length(fR))])

    fPhs{r} = replace(fR{r},'part-real','part-phase');
    % fMag{r} = replace(fR{r},'part-real','part-mag');
    
    if force || ~exist(fPhs{r},'file')% || ~exist(fMag{r},'file')
        mriR = MRIread(fR{r});
        mriI = MRIread(fI{r});
        mriPhs = rmfield(mriR,'vol'); mriPhs.fspec = fPhs{r};
        % mriMag = rmfield(mriR,'vol'); mriMag.fspec = fMag{r};
        % [mriPhs.vol,mriMag.vol] = cart2pol(mriR.vol,mriI.vol);
        [mriPhs.vol,~] = cart2pol(mriR.vol,mriI.vol);
        mriPhs.vol = mriPhs.vol/pi*4096;

        if ~exist(fileparts(mriPhs.fspec),'dir'); mkdir(fileparts(mriPhs.fspec)); end
        MRIwrite(mriPhs,mriPhs.fspec);
        % MRIwrite(mriMag.vol,mriMag.fspec);
        disp(' done')
    else
        disp(' already done, skipping')
    end
end