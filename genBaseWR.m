function outFiles = genBaseWR(inFiles,param,force,verbose)
% baseType,baseTypeForGenBase,afni3dAlineateArg
global srcAfni
%% Init
if exist('param','var') && ~isempty(param) && isfield(param,'enforceFixThroughPlane'); enforceFixThroughPlane = param.enforceFixThroughPlane; else; enforceFixThroughPlane = []; end
if exist('param','var') && ~isempty(param) && isfield(param,'baseType'); baseType = param.baseType; else; baseType = ''; end
if exist('param','var') && ~isempty(param) && isfield(param,'baseTypeForGenBase'); baseTypeForGenBase = param.baseTypeForGenBase; else; baseTypeForGenBase = ''; end
if exist('param','var') && ~isempty(param) && isfield(param,'afni3dAlineateArg'); afni3dAlineateArg = param.afni3dAlineateArg; else; afni3dAlineateArg = ''; end
if ~exist('force','var'); force = []; end
if ~exist('verbose','var'); verbose = []; end
%% Defaults
if isempty(enforceFixThroughPlane); enforceFixThroughPlane = 0; end
if isempty(baseType); baseType = 'mcAv'; end
if isempty(baseTypeForGenBase)
    switch baseType
        case 'mcAv'; baseTypeForGenBase = 'first';
        case 'first'; baseTypeForGenBase = '';
        otherwise; dbstack; error('code that');
    end
end
if isempty(afni3dAlineateArg); afni3dAlineateArg = {'-cost ls' '-interp quintic' '-final wsinc5'}; end
if isempty(force); force = 0; end
if isempty(verbose); verbose = 0; end




if isfield(inFiles,'manBrainMaskInv') && ~isempty(inFiles.manBrainMaskInv)
    fMask = inFiles.manBrainMaskInv;
else
    fMask = '';
end

outFiles.fEstimList = inFiles.fEstimList;
outFiles.fEstimBaseList = cell(size(outFiles.fEstimList));
switch baseType
    case 'mcAv'
        %% Generate base as the average of a first-pass-motion-corrected time seires, accounting for smoothing
        disp(['generating base for motion estimation (first-pass moco to ' baseTypeForGenBase ' frame, accounting for smoothing)'])
        for I = 1:length(outFiles.fEstimList)
            disp([' run' num2str(I) '/' num2str(length(outFiles.fEstimList))])
            cmd = {srcAfni};
            %%% set filename
            fIn = outFiles.fEstimList{I};
            fOut = strsplit(fIn,filesep); fOut{end} = ['mcRef-' baseType '_' fOut{end}]; fOut = strjoin(fOut,filesep);
            if force || ~exist(fOut,'file')
                %%% detect smoothing
                sm = strsplit(fIn,filesep); sm = strsplit(sm{end},'_'); ind = ~cellfun('isempty',regexp(sm,'^sm\d+$')); if any(ind); sm = sm{ind}; else sm = 'sm1'; end; sm = str2num(sm(3:end));
                n = MRIread(fIn,1); n = n.nframes - 1;
                nLim = [0 n] + [1 -1].*((sm+1)/2-1);
                %%% moco
                cmd{end+1} = '3dAllineate -overwrite \';
                switch baseTypeForGenBase
                    case 'first'
                        cmd{end+1} = ['-base ' fIn '[' num2str(nLim(1)) '] \'];
                    otherwise
                        dbstack; error('code that')
                end
                cmd{end+1} = ['-source ' fIn '[' num2str(nLim(1)) '..' num2str(nLim(2)) '] \'];
                cmd{end+1} = ['-prefix ' fOut ' \'];
                cmd{end+1} = [strjoin(afni3dAlineateArg,' ') ' \'];
                if ~isempty(fMask)
                    disp(['  using mask: ' fMask])
                    cmd{end+1} = ['-emask ' fMask ' \'];
                else
                    disp('  not using mask')
                end
                if enforceFixThroughPlane
                    cmd{end+1} = '-parfix 2 0 -parfix 4 0 -parfix 5 0 \';
                    disp('  enforcing no through-plane motion')
                end
                cmd{end+1} = '-nopad \';
                cmd{end+1} = '-warp shift_rotate';
                %%% average
                cmd{end+1} = '3dTstat -overwrite \';
                cmd{end+1} = ['-prefix ' fOut ' \'];
                cmd{end+1} = '-mean \';
                cmd{end+1} = fOut;

                %%% run shell command
                cmd = strjoin(cmd,newline); % disp(cmd)
                if verbose
                    [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
                else
                    [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
                end
                disp('  done')
            else
                disp('  already done, skipping')
            end
            %%% set base filename
            outFiles.fEstimBaseList{I} = fOut;
        end


    case 'first'
        %% Set base as the first properly smoothed frame
        disp('setting base for motion estimation (first properly smoothed frame)')
        for I = 1:length(outFiles.fEstimList)
            %%% set filename
            fIn = outFiles.fEstimList{I};
            %%% detect smoothing
            sm = strsplit(fIn,filesep); sm = strsplit(sm{end},'_'); ind = ~cellfun('isempty',regexp(sm,'^sm\d+$')); if any(ind); sm = sm{ind}; else sm = 'sm1'; end; sm = str2num(sm(3:end));
            n = MRIread(fIn,1); n = n.nframes - 1;
            nLim = [0 n] + [1 -1].*((sm+1)/2-1);
            %%% set base filename
            outFiles.fEstimBaseList{I} = [fIn '[' num2str(nLim(1)) ']'];
        end
        disp(' done')

end

%% Outputs
if isfield(inFiles,'manBrainMask')
    outFiles.manBrainMask = inFiles.manBrainMask;
end
if isfield(inFiles,'manBrainMaskInv')
    outFiles.manBrainMaskInv = inFiles.manBrainMaskInv;
end

outFiles.param.baseType = baseType;
outFiles.param.baseTypeForGenBase = baseTypeForGenBase;