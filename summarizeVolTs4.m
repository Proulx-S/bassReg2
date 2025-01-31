function smr = summarizeVolTs4(fList,nDummyNotRemoved,dataType,force,verbose)
global srcAfni
if ~exist('outTag','var');                   outTag = []; end
if ~exist('nDummyNotRemoved','var'); nDummyNotRemoved = []; end
if ~exist('force','var');                     force = []; end
if ~exist('verbose','var');                 verbose = []; end
if isempty(force);                     force = 0; end
if isempty(verbose);                 verbose = 0; end
if isempty(nDummyNotRemoved); nDummyNotRemoved = 0; end


if size(nDummyNotRemoved,1)==1
    nDummyNotRemoved = repmat(nDummyNotRemoved,size(fList));
end
% if size(nDummyNotRemove,2)~=size(fList,2)
%     nDummyNotRemove = repmat(nDummyNotRemove,size(fList,2));
% end

[fMeanList,fStdList,fSnrList,fFstMdLst,cat_fSummaryList,av_cat_fSummaryList,std_cat_fSummaryList]...
    = doIt(fList(:,1),[],nDummyNotRemoved,[],force,verbose);

smr.runAv.fList  = fMeanList;
smr.runSd.fList  = fStdList;
smr.runSnr.fList = fSnrList;
smr.runFML.fList = fFstMdLst;

smr.sesCat.runAv.fList  = cat_fSummaryList(:,:,1);
if size(cat_fSummaryList,3)==4
    smr.sesCat.runSd.fList  = cat_fSummaryList(:,:,2);
    smr.sesCat.runSnr.fList = cat_fSummaryList(:,:,3);
    smr.sesCat.runFML.fList = cat_fSummaryList(:,:,4);
else
    smr.sesCat.runSd.fList  = {};
    smr.sesCat.runSnr.fList = {};
    smr.sesCat.runFML.fList = {};
end

smr.sesAv.runAv.fList  = av_cat_fSummaryList(:,:,1);
if size(av_cat_fSummaryList,3)==4
    smr.sesAv.runSd.fList  = av_cat_fSummaryList(:,:,2);
    smr.sesAv.runSnr.fList = av_cat_fSummaryList(:,:,3);
    smr.sesAv.runFML.fList = av_cat_fSummaryList(:,:,4);
else
    smr.sesAv.runSd.fList  = {};
    smr.sesAv.runSnr.fList = {};
    smr.sesAv.runFML.fList = {};
end

smr.sesSd.runAv.fList  = std_cat_fSummaryList(:,:,1);
if size(std_cat_fSummaryList,3)==4
    smr.sesSd.runSd.fList  = std_cat_fSummaryList(:,:,2);
    smr.sesSd.runSnr.fList = std_cat_fSummaryList(:,:,3);
    smr.sesSd.runFML.fList = std_cat_fSummaryList(:,:,4);
else
    smr.sesSd.runSd.fList  = {};
    smr.sesSd.runSnr.fList = {};
    smr.sesSd.runFML.fList = {};
end

if ismember('PC',dataType)
    switch size(fList,2)
        case 1
            disp('only one variable, probably just the mag')
            disp(fList)
        case 2
            dbstack; error('amplitude and phase diff data? code that')
        case 3
            %%% summarize real
            disp('real')
            % [fMeanList,fStdList,fSnrList,fFstMdLst,cat_fSummaryList,av_cat_fSummaryList,std_cat_fSummaryList]...
            %     = doIt(fList(:,2),'vencDiffReal',nDummyNotRemoved,dataType,force,verbose);
            [fMeanList,fStdList,fSnrList,fFstMdLst,cat_fSummaryList,av_cat_fSummaryList,std_cat_fSummaryList]...
                = doIt(fList(:,2),[],nDummyNotRemoved,dataType,force,verbose);

            smr.runAv.fList  = cat(2,smr.runAv.fList ,fMeanList);
            smr.runSd.fList  = cat(2,smr.runSd.fList ,fStdList);
            smr.runSnr.fList = cat(2,smr.runSnr.fList,fSnrList);
            smr.runFML.fList = cat(2,smr.runFML.fList,fFstMdLst);

            smr.sesCat.runAv.fList  = cat(2,smr.sesCat.runAv.fList ,cat_fSummaryList(:,:,1));
            smr.sesCat.runSd.fList  = cat(2,smr.sesCat.runSd.fList ,cat_fSummaryList(:,:,2));
            smr.sesCat.runSnr.fList = cat(2,smr.sesCat.runSnr.fList,cat_fSummaryList(:,:,3));
            smr.sesCat.runFML.fList = cat(2,smr.sesCat.runFML.fList,cat_fSummaryList(:,:,4));

            smr.sesAv.runAv.fList  = cat(2,smr.sesAv.runAv.fList ,av_cat_fSummaryList(:,:,1));
            smr.sesAv.runSd.fList  = cat(2,smr.sesAv.runSd.fList ,av_cat_fSummaryList(:,:,2));
            smr.sesAv.runSnr.fList = cat(2,smr.sesAv.runSnr.fList,av_cat_fSummaryList(:,:,3));
            smr.sesAv.runFML.fList = cat(2,smr.sesAv.runFML.fList,av_cat_fSummaryList(:,:,4));

            smr.sesSd.runAv.fList  = cat(2,smr.sesSd.runAv.fList ,std_cat_fSummaryList(:,:,1));
            smr.sesSd.runSd.fList  = cat(2,smr.sesSd.runSd.fList ,std_cat_fSummaryList(:,:,2));
            smr.sesSd.runSnr.fList = cat(2,smr.sesSd.runSnr.fList,std_cat_fSummaryList(:,:,3));
            smr.sesSd.runFML.fList = cat(2,smr.sesSd.runFML.fList,std_cat_fSummaryList(:,:,4));

            %%% summarize imag
            disp('imag')
            % [fMeanList,fStdList,fSnrList,fFstMdLst,cat_fSummaryList,av_cat_fSummaryList,std_cat_fSummaryList]...
            %     = doIt(fList(:,3),'vencDiffImag',nDummyNotRemoved,dataType,force,verbose);
            [fMeanList,fStdList,fSnrList,fFstMdLst,cat_fSummaryList,av_cat_fSummaryList,std_cat_fSummaryList]...
                = doIt(fList(:,3),[],nDummyNotRemoved,dataType,force,verbose);

            smr.runAv.fList  = cat(2,smr.runAv.fList ,fMeanList);
            smr.runSd.fList  = cat(2,smr.runSd.fList ,fStdList);
            smr.runSnr.fList = cat(2,smr.runSnr.fList,fSnrList);
            smr.runFML.fList = cat(2,smr.runFML.fList,fFstMdLst);

            smr.sesCat.runAv.fList  = cat(2,smr.sesCat.runAv.fList ,cat_fSummaryList(:,:,1));
            smr.sesCat.runSd.fList  = cat(2,smr.sesCat.runSd.fList ,cat_fSummaryList(:,:,2));
            smr.sesCat.runSnr.fList = cat(2,smr.sesCat.runSnr.fList,cat_fSummaryList(:,:,3));
            smr.sesCat.runFML.fList = cat(2,smr.sesCat.runFML.fList,cat_fSummaryList(:,:,4));

            smr.sesAv.runAv.fList  = cat(2,smr.sesAv.runAv.fList ,av_cat_fSummaryList(:,:,1));
            smr.sesAv.runSd.fList  = cat(2,smr.sesAv.runSd.fList ,av_cat_fSummaryList(:,:,2));
            smr.sesAv.runSnr.fList = cat(2,smr.sesAv.runSnr.fList,av_cat_fSummaryList(:,:,3));
            smr.sesAv.runFML.fList = cat(2,smr.sesAv.runFML.fList,av_cat_fSummaryList(:,:,4));

            smr.sesSd.runAv.fList  = cat(2,smr.sesSd.runAv.fList ,std_cat_fSummaryList(:,:,1));
            smr.sesSd.runSd.fList  = cat(2,smr.sesSd.runSd.fList ,std_cat_fSummaryList(:,:,2));
            smr.sesSd.runSnr.fList = cat(2,smr.sesSd.runSnr.fList,std_cat_fSummaryList(:,:,3));
            smr.sesSd.runFML.fList = cat(2,smr.sesSd.runFML.fList,std_cat_fSummaryList(:,:,4));
        otherwise
            dbstack; error('more than one venc? code that')
    end
end


function [fMeanList,fStdList,fSnrList,fFstMdLst,cat_fSummaryList,av_cat_fSummaryList,std_cat_fSummaryList] = doIt(fList,outTag,nDummyNotRemove,dataType,force,verbose)
global srcAfni
if ~isempty(outTag) && ~strcmp(outTag(end),'_'); outTag(end+1) = '_'; end

%%%%%%%%%%%%%%%%
%% Within run %%
%%%%%%%%%%%%%%%%

% %% Deal with multivariate data
% if size(fList,2)>1
%     dbstack; error('this lookds like multi-echo data, code that')
% end

%% Summarize time
fMeanList = cell(size(fList));
fStdList  = cell(size(fList));
fSnrList  = cell(size(fList));
fFstMdLst = cell(size(fList));
for f = 1:numel(fList)
    % if ~all(ismember({'singleEcho' 'volTs'},dataType))
    %     dbstack; error('not singleEcho volTs')
    % end

    cmd = {srcAfni};
    fIn = fList{f};
    nFrame = MRIget(fIn,'nv');

    %mean
    if nFrame>1
        fOut = strsplit(fIn,filesep); fOut{end} = ['av_' outTag fOut{end}]; fOut = strjoin(fOut,filesep);
        if force || ~exist(fOut,'file')
            cmd{end+1} = '3dTstat -overwrite -mean \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fIn '[' num2str(nDummyNotRemove(f)) '..$]'];
        end
    else
        fOut = fIn;
    end
    fMeanList{f} = fOut;

    %std
    if nFrame>1
        fOut = strsplit(fIn,filesep); fOut{end} = ['std_' outTag fOut{end}]; fOut = strjoin(fOut,filesep);
        if force || ~exist(fOut,'file')
            cmd{end+1} = '3dTstat -overwrite -stdev \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fIn '[' num2str(nDummyNotRemove(f)) '..$]'];
        end
    else
        fOut = [];
    end
    fStdList{f} = fOut;

    %snr
    if nFrame>1
        fOut = strsplit(fIn,filesep); fOut{end} = ['snr_' outTag fOut{end}]; fOut = strjoin(fOut,filesep);
        if force || ~exist(fOut,'file')
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-a ' fMeanList{f} ' \'];
            cmd{end+1} = ['-b ' fStdList{f} ' \'];
            cmd{end+1} = ['-expr ''a/b'' \'];
            cmd{end+1} = ['-prefix ' fOut];
        end
    else
        fOut = [];
    end
    fSnrList{f} = fOut;

    % first, middle and last frames
    if nFrame>1
        fOut = strsplit(fIn,filesep); fOut{end} = ['fstMdLst_' outTag fOut{end}]; fOut = strjoin(fOut,filesep);
        if force || ~exist(fOut,'file')
            % nFrame = MRIget(fIn,'nt');
            cmd{end+1} = '3dcalc -overwrite \';
            % cmd{end+1} = ['-a ' fIn '[0,' num2str(round(nFrame/2-1)) ',' num2str(nFrame-1) '] \'];
            cmd{end+1} = ['-a ' fIn '[' num2str(nDummyNotRemove(f)) ',' num2str(round((nFrame-nDummyNotRemove(f))/2-1)) ',' num2str(nFrame-1) '] \'];
            % cmd{end+1} = ['-a ' fIn '[' num2str(nDummy(f)) ',' num2str(round((nFrame-nDummy(f))/2)) ',$] \'];
            cmd{end+1} = ['-expr a \'];
            cmd{end+1} = ['-prefix ' fOut];
        end
    else
        fOut = [];
    end
    fFstMdLst{f} = fOut;


    if length(cmd)>1
        if verbose
            disp([' summarizing volTs ' num2str(f) '/' num2str(numel(fList))])
        end
        if verbose>1
            [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status || contains(cmdout,'error','IgnoreCase',true); dbstack; error(cmdout); error('x'); end
        else
            [status,cmdout] = system(strjoin(cmd,newline)); if status || contains(cmdout,'error','IgnoreCase',true); dbstack; error(cmdout); error('x'); end
        end
    end
end

%%%%%%%%%%%%%%%%%
%% Between-run %%
%%%%%%%%%%%%%%%%%
cmd = {srcAfni};
if nFrame>1
    fSummaryList = cat(3,fMeanList,fStdList,fSnrList,fFstMdLst);
else
    fSummaryList = fMeanList;
end
cat_fSummaryList     = cell([1 size(fSummaryList,[2 3])]);
av_cat_fSummaryList  = cell([1 size(fSummaryList,[2 3])]);
std_cat_fSummaryList = cell([1 size(fSummaryList,[2 3])]);
for i = 1:prod(size(fSummaryList,[2 3]))
    %%% runCat
    fIn = fSummaryList(:,i);
    [a,b,~] = fileparts(replace(fIn{1},'.nii.gz',''));
    a = strsplit(a,'_'); a{contains(a,'run-')} = 'run-cat'; a = strjoin(a,'_');
    fOut = fullfile(a,[replace(b,'av_','') '.nii.gz']); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
    % fOut = fullfile(fileparts(a),['cat_' b '.nii.gz']);
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dTcat -overwrite \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = strjoin(fIn,' ');
    end
    cat_fSummaryList{1,i} = fOut;
    fCat = fOut;
    

    %%% runAv
    fIn = fCat;
    if size(fSummaryList,1)>1
        fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
        if force || ~exist(fOut,'file')
            cmd{end+1} = '3dTstat -overwrite -mean \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = fIn;
        end
    else
        fOut = fIn;
    end
    av_cat_fSummaryList{1,i} = fOut;

    %%% runStd
    fIn = fCat;
    if size(fSummaryList,1)>1
        fOut = strsplit(fIn,filesep); fOut{end} = ['std_' fOut{end}]; fOut = strjoin(fOut,filesep);
        if force || ~exist(fOut,'file')
            cmd{end+1} = '3dTstat -overwrite -stdev \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = fIn;
        end
    else
        fOut = [];
    end
    std_cat_fSummaryList{1,i} = fOut;
end

%%% launch command
if length(cmd)>1
    if verbose
        if any(ismember(dataType,'volTs'))
            disp(' summarizing volTs across runs')
        else
            disp(' summarizing vol across runs')
        end
    end
    if verbose>1
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status || contains(cmdout,'error','IgnoreCase',true); dbstack; error(cmdout); error('x'); end
    else
        [status,cmdout] = system(strjoin(cmd,newline)); if status || contains(cmdout,'error','IgnoreCase',true); dbstack; error(cmdout); error('x'); end
    end
end


% %% Output
% smr.runAv.fList  = fMeanList;
% smr.runSd.fList  = fStdList;
% smr.runSnr.fList = fSnrList;
% smr.runFML.fList = fFstMdLst;
% 
% smr.sesCat.runAv.fList  = cat_fSummaryList(:,:,1);
% smr.sesCat.runSd.fList  = cat_fSummaryList(:,:,2);
% smr.sesCat.runSnr.fList = cat_fSummaryList(:,:,3);
% smr.sesCat.runFML.fList = cat_fSummaryList(:,:,4);
% 
% smr.sesAv.runAv.fList  = av_cat_fSummaryList(:,:,1);
% smr.sesAv.runSd.fList  = av_cat_fSummaryList(:,:,2);
% smr.sesAv.runSnr.fList = av_cat_fSummaryList(:,:,3);
% smr.sesAv.runFML.fList = av_cat_fSummaryList(:,:,4);
% 
% smr.sesSd.runAv.fList  = std_cat_fSummaryList(:,:,1);
% smr.sesSd.runSd.fList  = std_cat_fSummaryList(:,:,2);
% smr.sesSd.runSnr.fList = std_cat_fSummaryList(:,:,3);
% smr.sesSd.runFML.fList = std_cat_fSummaryList(:,:,4);