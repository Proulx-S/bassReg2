function avMapPreproc(avMapSet,force,verbose)
global srcFs srcAfni bidsDir

if ~exist('force','var');              force         = []                    ; end
if ~exist('verbose','var');            verbose       = []                    ; end
if ~isfield(avMapSet,'fRegBasePlumb'); fRegBasePlumb = []                    ; end
if ~isfield(avMapSet,'fMaskPlumb');    fMaskPlumb    = []                    ; end

if isempty(force);                     force         = 0                     ; end
if isempty(verbose);                   verbose       = 1                     ; end
if isfield(avMapSet,'fRegBasePlumb');  fRegBasePlumb = avMapSet.fRegBasePlumb; end %must already be at plumb
if isfield(avMapSet,'fMaskPlumb');     fMaskPlumb    = avMapSet.fMaskPlumb   ; end %must already be at plumb

wDirAnat    = avMapSet.wDirAnat;
bidsPrpList = avMapSet.bidsPrpList;

%% Initiate preproc
for i = 1:length(bidsPrpList)
    av = dir(fullfile(bidsDir,'anat',bidsPrpList{i}));
    % av = fullfile({av(:).folder},{av(:).name})';
    if ~isempty(av); break; end
end
avMapSet.files = av;
avMapSet.label
avMapSet.dataType = {'vol' 'highSNR'};

force = 0;
avMapSet = initPreproc(avMapSet,wDirAnat,[],[],[],force);

avMapSet.initFiles.fEstimAvCatAv


%% Catenate echoes and do rms
fAvEchoCat = fullfile(wDirAnat,'avMap.nii.gz');

voxSize = MRIread(av{1},1);
voxSize = [voxSize.xsize voxSize.ysize voxSize.zsize];
if force || ~exist(fAvEchoCat,'file')
    for i = 1:length(av)
        [~,tmp{i}] = system(['jq ''.EchoTime'' ' replace(av{i},'.nii.gz','.json')]); tmp{i} = str2num(replace(replace(tmp{i},'[0;39m',''),['[0m' newline],''));
        av{i} = MRIread(av{i});
    end
    TE = [tmp{:}]';
    av = [av{:}]; av(1).vol = cat(4,av.vol); av(2:end) = [];
    av.vol = cat(4,av.vol,rms(av.vol,4)); % catenate rms after the last echo
    if ~exist(fileparts(fAvEchoCat),'dir'); mkdir(fileparts(fAvEchoCat)); end
    MRIwrite(av,fAvEchoCat);
else
    for i = 1:length(av)
        [~,tmp{i}] = system(['jq ''.EchoTime'' ' replace(av{i},'.nii.gz','.json')]); tmp{i} = str2num(replace(replace(tmp{i},'[0;39m',''),['[0m' newline],''));
    end
    TE = [tmp{:}]';
end

%% Smooth
sm = mean(voxSize(1:2)) *4;
fAvEchoCatSm = strsplit(fAvEchoCat,filesep); fAvEchoCatSm{end} = ['sm-' replace(num2str(sm),'.','p') '_' fAvEchoCatSm{end}]; fAvEchoCatSm = strjoin(fAvEchoCatSm,filesep);
if force || ~exist(fAvEchoCatSm,'file')
    cmd = {srcAfni};
    cmd{end+1} = '3dmerge -overwrite -doall \';
    cmd{end+1} = ['-prefix ' fAvEchoCatSm ' \'];
    cmd{end+1} = ['-1blur_fwhm ' num2str(sm) ' \'];
    cmd{end+1} = fAvEchoCat;
    cmd = strjoin(cmd,newline);
    if verbose
        [status,cmdout] = system(cmd,'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    else
        [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    end
    disp(' done')
else
    disp(' already done, skipping')
end

%% Fit T2* and S0
nneg = 0;
fInList = {fAvEchoCat fAvEchoCatSm};
for i = 1:length(fInList)
    fIn = fInList{i};
    avS0 = replace(fIn,'.nii.gz','_S0.nii.gz');
    avT2s = replace(fIn,'.nii.gz','_T2s.nii.gz');
    if force || ~exist(avS0,'file') || ~exist(avT2s,'file')
        av = MRIread(fIn);
        [t2star,S0] = calc_t2s_vol(av.vol(:,:,:,1:end-1), TE*1000, nneg, verbose); % last echo is actually the rms average
        av.vol = S0; if force || ~exist(avS0,'file'); MRIwrite(av,avS0); end
        av.vol = t2star; if force || ~exist(avT2s,'file'); MRIwrite(av,avT2s); end
    end
end


%% Register
force = 1;
if ~isempty(fRegBasePlumb)
    [~,b,c] = fileparts(fRegBasePlumb);
    disp(['Estimating motion' newline 'from: ''avMap'' ' newline '  to: ' b c ' (' fRegBasePlumb ')'])
    
    %%% Initiate avMap (oblique to plumb)
    cmd = {srcAfni};
    fSource = [fAvEchoCat '[$]'];
    fSourcePlumb = strsplit(fAvEchoCat,filesep); fSourcePlumb{end} = ['plumb_' fSourcePlumb{end}]; fSourcePlumb = strjoin(fSourcePlumb,filesep);
    if force || ~exist(fSourcePlumb,'file')
        cmd{end+1} = '3dNwarpApply -overwrite \';
        cmd{end+1} = ['-nwarp ''IDENT(' fSource ')'' \'];
        cmd{end+1} = ['-prefix ' fSourcePlumb ' \'];
        cmd{end+1} = ['-source ' fSource];
    end
    

    %%% Estimate registration
    fBasePlumb = fRegBasePlumb;
    fSourcePlumb;
    fMaskPlumb;
    fOutPlumb = strsplit(fSourcePlumb,filesep); fOutPlumb{end} = ['reg2func_' fOutPlumb{end}]; fOutPlumb = strjoin(fOutPlumb,filesep);
    fOutPlumbParam = replace(fOutPlumb,'.nii.gz','');

    explicitlyFixThroughPlane = 1;
    spSmFac = 0;
    
    % cmd = {srcAfni};
    cmd{end+1} = '3dAllineate -overwrite \';
    cmd{end+1} = ['-base ' fBasePlumb ' \'];
    cmd{end+1} = ['-source ' fSourcePlumb ' \'];
    % cmd{end+1} = ['-prefix ' fOutPlumb ' \'];
    cmd{end+1} = ['-1Dparam_save ' fOutPlumbParam ' \'];
    cmd{end+1} = ['-1Dmatrix_save ' fOutPlumbParam ' \'];
    % afni3dAlineateArg = {'-cost ls' '-interp quintic' '-final wsinc5'};
    % afni3dAlineateArg = {'-cost ls+' '-interp quintic' '-final wsinc5'};
    % afni3dAlineateArg = {'-cost ls+ZZ' '-interp quintic' '-final wsinc5'};
    % afni3dAlineateArg = {'-cost lpa' '-interp quintic' '-final wsinc5'};
    % afni3dAlineateArg = {'-cost lpa+' '-interp quintic' '-final wsinc5'};
    % afni3dAlineateArg = {'-cost lpa+ZZ' '-interp quintic' '-final wsinc5'};
    % afni3dAlineateArg = {'-cost lpc' '-interp quintic' '-final wsinc5'};
    % afni3dAlineateArg = {'-cost lpc+' '-interp quintic' '-final wsinc5'};
    % afni3dAlineateArg = {'-cost lpc+ZZ' '-interp quintic' '-final wsinc5'};
    afni3dAlineateArg = {'-cost mi' '-interp quintic' '-final wsinc5'};
    % afni3dAlineateArg = {'-cost mi+' '-interp quintic' '-final wsinc5'};
    % afni3dAlineateArg = {'-cost mi+ZZ' '-interp quintic' '-final wsinc5'};
    % afni3dAlineateArg = {'-cost nmi' '-interp quintic' '-final wsinc5'};
    % afni3dAlineateArg = {'-cost nmi+' '-interp quintic' '-final wsinc5'};
    % afni3dAlineateArg = {'-cost nmi+ZZ' '-interp quintic' '-final wsinc5'};
    cmd{end+1} = [strjoin(afni3dAlineateArg,' ') ' \'];
    if ~isempty(fMaskPlumb)
        disp([' using mask: ' fMaskPlumb])
        cmd{end+1} = ['-emask ' fMaskPlumb ' \'];
    else
        disp(' not using mask')
    end
    if explicitlyFixThroughPlane
        cmd{end+1} = '-parfix 2 0 -parfix 4 0 -parfix 5 0 \';
        disp(' explicitly enforcing no through-plane motion')
    end
    cmd{end+1} = '-nopad -conv 0 -nmatch 100% -onepass -nocmass \';
    % cmd{end+1} = '-nopad -conv 0 -nmatch 100% -nocmass \';
    % cmd{end+1} = '-maxrot 2 -maxshf 1 \';
    if spSmFac>0
        cmd{end+1} = ['-fineblur ' num2str(vsize*spSmFac) ' \'];
    end
    cmd{end+1} = '-warp shift_rotate'; % cmd{end+1} = ['-warp shift_rotate -parfix 2 0 -parfix 4 0 -parfix 5 0'];
    
system(strjoin(cmd,newline),'-echo')
    %%% Apply registration
    fSource = fAvEchoCat;
    fSourcePlumb = strsplit(fAvEchoCat,filesep); fSourcePlumb{end} = ['plumb_' fSourcePlumb{end}]; fSourcePlumb = strjoin(fSourcePlumb,filesep);
    [a,b,c] = fileparts(fSource);
    fOutPlumb = [a,'preproc']; if ~exist(fOutPlumb,'dir'); mkdir(fOutPlumb); end
    fOutPlumb = fullfile(fOutPlumb,['reg2func_' b c]);
    %%%% Put to plumb
    if force || ~exist(fSourcePlumb,'file')
        cmd{end+1} = '3dNwarpApply -overwrite \';
        cmd{end+1} = ['-nwarp ''IDENT(' fSource ')'' \'];
        cmd{end+1} = ['-prefix ' fSourcePlumb ' \'];
        cmd{end+1} = ['-source ' fSource];
    end
    %%%% Then apply
    cmd{end+1} = '3dAllineate -overwrite \';
    cmd{end+1} = ['-source ' fSourcePlumb ' \'];
    cmd{end+1} = ['-prefix ' fOutPlumb ' \'];
    cmd{end+1} = ['-1Dparam_save ' fOutPlumbParam ' \'];
    cmd{end+1} = ['-1Dmatrix_save ' fOutPlumbParam ' \'];





    if verbose; disp(strjoin(cmd,newline)); end
    

    cmd{end+1} = ['-base ' fBasePlumb ' \'];
    cmd{end+1} = ['-source ' fSourcePlumb ' \'];
    cmd{end+1} = ['-prefix ' fOutPlumb ' \'];
    
    cmd{end+1} = 'freeview \';
    cmd{end+1} = [fBasePlumb ' \'];
    cmd{end+1} = [fSourcePlumb ' \'];
    cmd{end+1} = [fOutPlumb ' \'];
    cmd{end+1} = fMaskPlumb;
    
    clipboard('copy',strjoin(cmd,newline))






    
end

