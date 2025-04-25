function [fVolCorr,fApplyCorr,fVol,fVolField] = correctBiasField(f, fMask, fApply, fOblique, force, verbose)
    global src;
    if ~exist('fMask','var');       fMask = []; end
    if ~exist('fApply','var');     fApply = []; end
    if ~exist('fOblique','var'); fOblique = []; end
    if ~exist('force','var');       force = []; end
    if ~exist('verbose','var');   verbose = []; end
    if isempty(force);              force = 0 ; end
    if isempty(verbose);          verbose = 0 ; end
    fApply = cellstr(fApply);

    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % Need to detect when about to write output file in bids directory and write somewhere else.
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!thy
    for iii = 1:length(fApply)
        if exist(replace(fApply{iii},'.nii.gz','.json'),'file')
            dbstack; error('seems like you want to correctBiasField a bids file, but this will write the result in the same folder, and it is a very naughty thing to do');
        end
    end


        
    if iscell(f)
        dbstack; error('cell input not supported, code that');
        fVolCorr   = cell(size(f));
        fVolTsCorr = cell(size(f));
        fVolField  = cell(size(f));
        for r = 1:length(f)
            [fVolCorr{r},fVolTsCorr{r},fVolField{r}] = correctBiasField(f{r}, fMask, force, verbose);
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

    %% Detect and fix header issue
    cmd = {src.ants};
    cmd{end+1} = ['PrintHeader' fVol];
    [status,cmdout] = system(strjoin(cmd,newline));
    if status % detect header issue
        % check if the volume is oblique using AFNI's 3dinfo
        cmd2 = {src.afni};
        cmd2{end+1} = ['3dinfo -is_oblique ' fVol];
        [status2, isOblique] = system(strjoin(cmd2, newline)); isOblique = str2num(isOblique);
        if isOblique
            % there is no issue
            fVolOblique = [];
            nSlice = MRIget(fVol,'depth');
        else
            % there is an issue
            if isempty(fOblique)
                % cannot fix it without fOblique
                dbstack; error('Header issue with ANTs. Please provide an original file (not deoblique) as f or fOblique'); % Added error handling
            else
                % fix it using fOblique
                disp('header issue, using fOblique header')
                mri = MRIread(fOblique,1);
                mriTmp = MRIread(fVol);
                mri.vol = mriTmp.vol; clear mriTmp;
                fVolOblique = [tempname '.nii.gz'];
                MRIwrite(mri,fVolOblique);
                nSlice = mri.depth;
            end
        end
    else
        fVolOblique = [];
    end


    %% N4 correction
    cmd = {src.afni};
    cmd{end+1} = src.ants;
    fVolCorr  = strsplit(fVol,filesep); fVolCorr{end}  = ['N4_' fVolCorr{end}]; fVolCorr = strjoin(fVolCorr,filesep);
    % if ~isempty(fVolTs)
    %     fVolTsCorr = strsplit(fVolTs,filesep); fVolTsCorr{end}  = ['N4_' fVolTsCorr{end}]; fVolTsCorr = strjoin(fVolTsCorr,filesep);
    % else
    %     fVolTsCorr = [];
    % end
    fVolField = replace(fVolCorr,{'_volTs.nii.gz' '_vol.nii.gz'},'_volN4field.nii.gz');
    if isempty(fApply)
        fApplyCorr = [];
    else
        fApply     = cellstr(fApply);
        fApplyCorr = cell(size(fApply));
        for i = 1:length(fApply)
            fApplyCorr{i} = strsplit((fApply{i}),filesep);
            fApplyCorr{i}{end} = ['N4_' fApplyCorr{i}{end}];
            fApplyCorr{i} = strjoin(fApplyCorr{i},filesep);
        end
    end
    if force || ~exist(fVolField,'file') || ~exist(fVolCorr,'file')% || (~isempty(fVolTsCorr) && ~exist(fVolTsCorr,'file'))
        % mask out non-brain
        if isempty(fMask)
            cmd{end+1} = ['cp ' fVol ' ' fVolCorr];
        else
            cmd{end+1} = '3dcalc -overwrite \';
            if ~isempty(fVolOblique)
                cmd{end+1} = ['-a ' fVolOblique ' \'];
            else
                cmd{end+1} = ['-a ' fVol ' \'];
            end
            cmd{end+1} = ['-b ' fMask ' \'];
            cmd{end+1} = ['-expr ''a*b'' \'];
            cmd{end+1} = ['-prefix ' fVolCorr];
        end

        % compute N4 correction
        if nSlice>1
            cmd{end+1} = ['N4BiasFieldCorrection \'];
        else
            cmd{end+1} = ['N4BiasFieldCorrection -d 2 \'];
        end
        cmd{end+1} = ['-i ' fVolCorr ' \'];
        cmd{end+1} = ['-o [' fVolCorr ',' fVolField ']'];

        % if ~isempty(fVolOrig)
        %     fVol = fVolOrig; clear fVolOrig;
        % end

        % manually apply field correction
        cmd{end+1} = '3dcalc -overwrite \';
        cmd{end+1} = ['-a ' fVol ' \'];
        cmd{end+1} = ['-b ' fVolField ' \'];
        cmd{end+1} = ['-expr ''a/b'' \'];
        cmd{end+1} = ['-prefix ' fVolCorr];

        % % manually apply field correction to time series
        % if ~isempty(fVolTs)
        %     cmd{end+1} = '3dcalc -overwrite \';
        %     cmd{end+1} = ['-a ' fVolTs ' \'];
        %     cmd{end+1} = ['-b ' fVolField ' \'];
        %     cmd{end+1} = ['-expr ''a/b'' \'];
        %     cmd{end+1} = ['-prefix ' fVolTsCorr];
        % end
        %and to other files
        if ~isempty(fApply)
            for i = 1:length(fApply)
                cmd{end+1} = '3dcalc -overwrite \';
                cmd{end+1} = ['-a ' fApply{i} ' \'];
                cmd{end+1} = ['-b ' fVolField ' \'];
                cmd{end+1} = ['-expr ''a/b'' \'];
                cmd{end+1} = ['-prefix ' fApplyCorr{i}];
            end
        end




    end
    if length(cmd)>2
        if verbose
            [status,cmdout] = system(strjoin(cmd,newline),'-echo');
        else
            [status,cmdout] = system(strjoin(cmd,newline));
        end
    end

    % conform fVolField
    MRIconform(fVolField, f)
