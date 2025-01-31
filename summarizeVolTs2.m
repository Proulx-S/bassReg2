function smr = summarizeVolTs2(fList,fOutList,nFrame,nDummy,dataType,force,verbose)
global srcAfni
if ~exist('fOutList','var'); fOutList = []; end
if ~exist('nDummy','var');     nDummy = []; end
if ~exist('force','var'); force = []; end
if ~exist('verbose','var'); verbose = []; end
if isempty(force);     force = 0; end
if isempty(verbose); verbose = 0; end
if isempty(nDummy);   nDummy = 0; end

if size(nDummy,1)==1
    nDummy = repmat(nDummy,size(fList));
end



% % % %%%% Reshape lists with echos in different columns
% % % echoList = bidsList(contains(bidsList(1,:),'echo-'),:,:);
% % % echoUniqueList = unique(echoList);
% % % 
% % % bidsList = permute(reshape(bidsList,size(bidsList,1),length(echoUniqueList),length(echoList)/length(echoUniqueList)),[1 3 2]);
% % % fPlumb = reshape(fPlumb,length(echoUniqueList),length(echoList)/length(echoUniqueList))';
% % % if nFrame>1
% % %     fPlumbAv = reshape(fPlumbAv,length(echoUniqueList),length(echoList)/length(echoUniqueList))';
% % % end
% % % fOrigList = reshape(fOrigList,length(echoUniqueList),length(echoList)/length(echoUniqueList))';
% % % 
% % % disp('---------------------------------------------------')
% % % disp('---------------------------------------------------')
% % % warning(['Make sure echoes are properly sorted below' ...
% % %     newline '(runs across rows and echos across columns)'])
% % % if size(bidsList,2)>1
% % %     disp(fullfile(squeeze(bidsList(contains(bidsList(:,1,1),'run-'),:,:)),squeeze(bidsList(contains(bidsList(:,1,1),'echo-'),:,:))))
% % % else
% % %     disp(fullfile(permute((bidsList(contains(bidsList(:,1,1),'echo-'),:,:)),[1 3 2])));
% % % end
% % % disp('---------------------------------------------------')
% % % disp('---------------------------------------------------')

%%%%%%%%%%%%%%%%
%% Within run %%
%%%%%%%%%%%%%%%%

%% Deal with multivariate data
if size(fList,2)>1
    dbstack; error('this lookds like multi-echo data, code that')
end

%% Summarize time
fMeanList = cell(size(fList));
fStdList  = cell(size(fList));
fSnrList  = cell(size(fList));
fFstMdLst = cell(size(fList));
for f = 1:size(fList,1)
    % if ~all(ismember({'singleEcho' 'volTs'},dataType))
    %     dbstack; error('not singleEcho volTs')
    % end

    cmd = {srcAfni};
    fIn = fList{f};

    %mean
    if ~isempty(fOutList)
        fOut = strsplit(fOutList{f},filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
    else
        fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
    end
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dTstat -overwrite -mean \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = [fIn '[' num2str(nDummy(f)) '..$]'];
    end
    fMeanList{f} = fOut;

    %std
    if ~isempty(fOutList)
        fOut = strsplit(fOutList{f},filesep); fOut{end} = ['std_' fOut{end}]; fOut = strjoin(fOut,filesep);
    else
        fOut = strsplit(fIn,filesep); fOut{end} = ['std_' fOut{end}]; fOut = strjoin(fOut,filesep);
    end
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dTstat -overwrite -stdev \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = [fIn '[' num2str(nDummy(f)) '..$]'];
    end
    fStdList{f} = fOut;

    %snr
    if ~isempty(fOutList)
        fOut = strsplit(fOutList{f},filesep); fOut{end} = ['snr_' fOut{end}]; fOut = strjoin(fOut,filesep);
    else
        fOut = strsplit(fIn,filesep); fOut{end} = ['snr_' fOut{end}]; fOut = strjoin(fOut,filesep);
    end
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dcalc -overwrite \';
        cmd{end+1} = ['-a ' fMeanList{f} ' \'];
        cmd{end+1} = ['-b ' fStdList{f} ' \'];
        cmd{end+1} = ['-expr ''a/b'' \'];
        cmd{end+1} = ['-prefix ' fOut];
    end
    fSnrList{f} = fOut;

    % first, middle and last frames
    if ~isempty(fOutList)
        fOut = strsplit(fOutList{f},filesep); fOut{end} = ['fstMdLst_' fOut{end}]; fOut = strjoin(fOut,filesep);
    else
        fOut = strsplit(fIn,filesep); fOut{end} = ['fstMdLst_' fOut{end}]; fOut = strjoin(fOut,filesep);
    end

    if force || ~exist(fOut,'file')
        nFrame = MRIget(fIn,'nt');
        cmd{end+1} = '3dcalc -overwrite \';
        cmd{end+1} = ['-a ' fIn '[' num2str(nDummy(f)) ',' num2str(round((nFrame-nDummy(f))/2)) ',$] \'];
        cmd{end+1} = ['-expr a \'];
        cmd{end+1} = ['-prefix ' fOut];
    end
    fFstMdLst{f} = fOut;


    if length(cmd)>1
        if verbose
            disp([' summarizing volTs ' num2str(f) '/' num2str(size(fList,1))])
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
fSummaryList = cat(2,fMeanList,fStdList,fSnrList,fFstMdLst);
cat_fSummaryList     = cell(1,size(fSummaryList,2));
av_cat_fSummaryList  = cell(1,size(fSummaryList,2));
std_cat_fSummaryList = cell(1,size(fSummaryList,2));
for i = 1:size(fSummaryList,2)
    %%% runCat
    fIn = fSummaryList(:,i);
    [a,b,~] = fileparts(replace(fIn{1},'.nii.gz',''));
    fOut = fullfile(fileparts(a),['cat_' b '.nii.gz']);
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dTcat -overwrite \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = strjoin(fIn,' ');
    end
    cat_fSummaryList{i} = fOut;
    fCat = fOut;
    

    %%% runAv
    fIn = fCat;
    fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dTstat -overwrite -mean \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = fIn;
    end
    av_cat_fSummaryList{i} = fOut;

    %%% runStd
    fIn = fCat;
    fOut = strsplit(fIn,filesep); fOut{end} = ['std_' fOut{end}]; fOut = strjoin(fOut,filesep);
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dTstat -overwrite -stdev \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = fIn;
    end
    std_cat_fSummaryList{i} = fOut;
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


%% Output
smr.runAv.fList  = fMeanList;
smr.runSd.fList  = fStdList;
smr.runSnr.fList = fSnrList;
smr.runFML.fList = fFstMdLst;

smr.sesCat.runAv.fList  = cat_fSummaryList(1);
smr.sesCat.runSd.fList  = cat_fSummaryList(2);
smr.sesCat.runSnr.fList = cat_fSummaryList(3);
smr.sesCat.runFML.fList = cat_fSummaryList(4);

smr.sesAv.runAv.fList  = av_cat_fSummaryList(1);
smr.sesAv.runSd.fList  = av_cat_fSummaryList(2);
smr.sesAv.runSnr.fList = av_cat_fSummaryList(3);
smr.sesAv.runFML.fList = av_cat_fSummaryList(4);

smr.sesSd.runAv.fList  = std_cat_fSummaryList(1);
smr.sesSd.runSd.fList  = std_cat_fSummaryList(2);
smr.sesSd.runSnr.fList = std_cat_fSummaryList(3);
smr.sesSd.runFML.fList = std_cat_fSummaryList(4);