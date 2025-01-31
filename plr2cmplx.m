function [fPhsDiffR,fPhsDiffI] = plr2cmplx(fMag0,fPhsDiff,force)
if ~exist('force','var'); force = []; end
if isempty(force);        force = 0 ; end

fPhsDiffR = cell(size(fPhsDiff));
fPhsDiffI = cell(size(fPhsDiff));
for r = 1:length(fMag0)
    disp([num2str(r) '/' num2str(length(fMag0))])

    %% Just slap vencOff signal magnitude to the venc-induced phase shift
    fPhsDiffR{r} = replace(fPhsDiff{r},'part-phase','part-real');
    fPhsDiffI{r} = replace(fPhsDiff{r},'part-phase','part-imag');

    if force || ~exist(fPhsDiffR{r},'file') || ~exist(fPhsDiffI{r},'file')
        % Read vencOff signal magnitude (mag0) and vencOff-vencOn phase (phsDiff)
        mriMag0     = MRIread(fMag0{r}   );
        mriPhsDiff  = MRIread(fPhsDiff{r});
        mriPhsDiffR = rmfield(mriPhsDiff,'vol');
        mriPhsDiffI = rmfield(mriPhsDiff,'vol');

        % Slap mag0 to phsDiff and get real and imag
        [mriPhsDiffR.vol,mriPhsDiffI.vol] = pol2cart(mriPhsDiff.vol./4096.*pi,mriMag0.vol);

        % Write real and imag
        mriPhsDiffR.fspec = fPhsDiffR{r}; if ~exist(fileparts(mriPhsDiffR.fspec),'dir'); mkdir(fileparts(mriPhsDiffR.fspec)); end
        mriPhsDiffI.fspec = fPhsDiffI{r}; if ~exist(fileparts(mriPhsDiffI.fspec),'dir'); mkdir(fileparts(mriPhsDiffI.fspec)); end
        
        disp(['writing ' mriPhsDiffR.fspec])
        MRIwrite(mriPhsDiffR,mriPhsDiffR.fspec);
        disp(['writing ' mriPhsDiffI.fspec])
        MRIwrite(mriPhsDiffI,mriPhsDiffI.fspec);
    else
        disp([fPhsDiffR{r} 'already done, skipping'])
        disp([fPhsDiffI{r} 'already done, skipping'])
    end

    %% Other strategies
    % (1) Keep magnitude of venc-induced phase shift to 1
    % (2) Recover real and imaginary values of both vencOff and vencOn data
    %     by assuming 0-phase for the [vencOff vencOn] average
end