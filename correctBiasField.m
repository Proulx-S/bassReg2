function [fVolCorr,fVolTsCorr,fVol,fVolField] = correctBiasField(f, fMask, force)
    global src;
    if ~exist('fMask','var'); fMask = []; end
    if ~exist('force','var'); force = []; end
    if isempty(force);        force = 0 ; end

    if iscell(f)
        fVolCorr   = cell(size(f));
        fVolTsCorr = cell(size(f));
        fVolField  = cell(size(f));
        for r = 1:length(f)
            [fVolCorr{r},fVolTsCorr{r},fVolField{r}] = correctBiasField(f{r}, fMask, force);
        end
        return;
    end



    %% Detect if volTs
    if MRIget(f,'nFrame')==1
        fVolTs = [];
        fVol   = f;
    else
        fVolTs = f;
        fVol = strsplit(fVolTs,filesep);
        fVol{end} = ['av_' fVol{end}];
        fVol = strjoin(fVol,filesep);
        if ~exist(fVol,'file')
            dbstack; error('volTs not found')
        end
    end


    %% N4 correction
    cmd = {src.afni};
    cmd{end+1} = src.ants;
    fVolCorr  = strsplit(fVol,filesep); fVolCorr{end}  = ['N4_' fVolCorr{end}]; fVolCorr = strjoin(fVolCorr,filesep);
    if ~isempty(fVolTs)
        fVolTsCorr = strsplit(fVolTs,filesep); fVolTsCorr{end}  = ['N4_' fVolTsCorr{end}]; fVolTsCorr = strjoin(fVolTsCorr,filesep);
    else
        fVolTsCorr = [];
    end
    fVolField = replace(fVolCorr,{'_volTs.nii.gz' '_vol.nii.gz'},'_volN4field.nii.gz');
    if force || ~exist(fVolField,'file') || ~exist(fVolCorr,'file') || (~isempty(fVolTsCorr) && ~exist(fVolTsCorr,'file'))
        % mask out non-brain
        if isempty(fMask)
            cmd{end+1} = ['cp ' fVol ' ' fVolCorr];
        else
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-a ' fVol ' \'];
            cmd{end+1} = ['-b ' fMask ' \'];
            cmd{end+1} = ['-expr ''a*b'' \'];
            cmd{end+1} = ['-prefix ' fVolCorr];
        end

        % compute N4 correction
        cmd{end+1} = ['N4BiasFieldCorrection -d 2 \'];
        cmd{end+1} = ['-i ' fVolCorr ' \'];
        cmd{end+1} = ['-o [' fVolCorr ',' fVolField ']'];

        % manually apply field correction
        cmd{end+1} = '3dcalc -overwrite \';
        cmd{end+1} = ['-a ' fVol ' \'];
        cmd{end+1} = ['-b ' fVolField ' \'];
        cmd{end+1} = ['-expr ''a/b'' \'];
        cmd{end+1} = ['-prefix ' fVolCorr];

        % manually apply field correction to time series
        if ~isempty(fVolTs)
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-a ' fVolTs ' \'];
            cmd{end+1} = ['-b ' fVolField ' \'];
            cmd{end+1} = ['-expr ''a/b'' \'];
            cmd{end+1} = ['-prefix ' fVolTsCorr];
        end

    end
    [status,cmdout] = system(strjoin(cmd,newline),'-echo');

    % conform fVolField
    MRIconform(fVolField, f)
