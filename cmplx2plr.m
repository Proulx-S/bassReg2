function [fMag,fPhs] = cmplx2plr(fR,fI,force)
% function [fMag,fPhs] = cmplx2plr(fR,fI,force,siemensFlag)
if ~exist('force','var');         force = []; end
if isempty(force);                force = 0 ; end
% if ~exist('siemensFlag','var'); siemensFlag = []; end
% if isempty(siemensFlag);          siemensFlag = 0; end
fR = cellstr(fR);
fI = cellstr(fI);

fPhs = cell(size(fR));
fMag = cell(size(fR));
for r = 1:length(fR)
    disp([num2str(r) '/' num2str(length(fR))])

    fPhs{r} = replace(fR{r},'part-real','part-phase');
    fMag{r} = replace(fR{r},'part-real','part-mag');
    
    % avoid overwriting
    fMag{r} = strsplit(fMag{r},'_');
    fMag{r}{contains(fMag{r},'rec-vencDiff')} = 'rec-venc0';
    fMag{r} = strjoin(fMag{r},'_');
    fMag{r} = strsplit(fMag{r},filesep);
    fMag{r}{end} = ['cmplxIntrp_' fMag{r}{end}];
    fMag{r} = strjoin(fMag{r},filesep);

    % if siemensFlag
    %     fMag{r} = strsplit(fMag{r},'_');
    %     fMag{r}{contains(fMag{r},'rec-vencDiff')} = 'rec-venc0';
    %     fMag{r} = strjoin(fMag{r},'_');
    %     fPhs{r} = replace(fPhs{r},'av_','cmplxAv_');
    %     fMag{r} = replace(fMag{r},'av_','cmplxAv_');
    % end
    
    if force || ~exist(fPhs{r},'file') || ~exist(fMag{r},'file')
        mriR = MRIread(fR{r});
        mriI = MRIread(fI{r});
        mriPhs = rmfield(mriR,'vol'); mriPhs.fspec = fPhs{r};
        mriMag = rmfield(mriR,'vol'); mriMag.fspec = fMag{r};

        [mriPhs.vol,mriMag.vol] = cart2pol(mriR.vol,mriI.vol);
        mriPhs.vol = mriPhs.vol./pi.*4096;

        if ~exist(fileparts(mriPhs.fspec),'dir'); mkdir(fileparts(mriPhs.fspec)); end
        MRIwrite(mriPhs,mriPhs.fspec);
        if ~exist(fileparts(mriMag.fspec),'dir'); mkdir(fileparts(mriMag.fspec)); end
        MRIwrite(mriMag,mriMag.fspec);
        disp(' done')

        % f = 20;
        % figure('WindowStyle','docked');
        % imagesc(mriMag.vol(:,:,f));
        % ax = gca; ax.DataAspectRatio = [1 1 1]; ax.Colormap = gray; colorbar
        % mriTmp = MRIread(mriMag.fspec);
        % figure('WindowStyle','docked');
        % imagesc(mriTmp.vol(:,:,f));
        % ax = gca; ax.DataAspectRatio = [1 1 1]; ax.Colormap = gray; colorbar
        
    else
        disp(' already done, skipping')
    end
end