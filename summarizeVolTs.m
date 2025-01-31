function summarizeVolTs(fList,fOutList,nFrame,nDummy,dataType,force,verbose)
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


%% Individual run
if any(ismember(dataType,'volTs'))
    fMeanList = cell(size(fList));
    fStdList  = cell(size(fList));
    fSnrList  = cell(size(fList));
    for f = 1:size(fList,1)
        if ~all(ismember({'singleEcho' 'volTs'},dataType))
            dbstack; error('not singleEcho volTs')
        end

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
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-a ' fIn '[' num2str(nDummy(f)) ',' num2str(round((nFrame(f)-nDummy(f))/2)) ',$] \'];
            cmd{end+1} = ['-expr a \'];
            cmd{end+1} = ['-prefix ' fOut];
        end
        fFstMdLst = fOut;


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


    %     continue
    % 
    % 
    %     %%%% Write means (before rms) for later visualization
    %     if nRun>1
    %         cmd = {srcFs};
    %         if nFrame>1
    %             fPlumbCat = cell(1,size(fPlumbAv,2));
    %             fPlumbAvCat = cell(1,size(fPlumbAv,2));
    %         else
    %             fPlumbCat = cell(1,size(fPlumb,2));
    %             fPlumbAvCat = cell(1,size(fPlumb,2));
    %         end
    %         for E = 1:nEcho
    %             if nFrame>1
    %                 fIn = fPlumbAv(:,E);
    %             else
    %                 fIn = fPlumb(:,E);
    %             end
    %             fOut = replace(fIn{1},bidsList{contains(squeeze(bidsList(:,1,1)),'run-'),1,1},'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
    %             fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
    %             if forceRecomputeRMS || ~exist(fOut,'file')
    %                 cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
    %             end
    %             fPlumbCat{E} = fOut;
    %             fIn = fOut;
    %             fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
    %             if forceRecomputeRMS || ~exist(fOut,'file')
    %                 cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
    %             end
    %             fPlumbAvCat{E} = fOut;
    %         end
    %         disp(' averaging')
    %         if length(cmd)>1
    %             cmd = strjoin(cmd,newline); % disp(cmd)
    %             [status,cmdout] = system(cmd); if status; dbstack; error(cmdout); error('x'); end
    %             disp('  done')
    %         else
    %             disp('  already done, skipping')
    %         end
    %         disp(' non need for averaging (single-run)')
    %     end
    % 
    %     %%%% Compute and write rms
    %     disp('computing cross-echo rms')
    %     fEstim = cell(size(fPlumb,1),1);
    %     for R = 1:size(fPlumb,1)
    %         disp([' file' num2str(R) '/' num2str(size(fPlumb,1))])
    %         fOut = replace(fPlumb{R,1},'echo-1','echo-rms');
    %         if forceRecomputeRMS || ~exist(fOut,'file')
    %             [a,~] = fileparts(fOut); if ~exist(a,'dir'); mkdir(a); end
    % 
    %             mri = MRIread(fPlumb{R,1},1);
    %             mri.vol = zeros([mri.volsize mri.nframes]);
    %             for ii = 1:size(fPlumb,2)
    %                 mriTmp = MRIread(fPlumb{R,ii});
    %                 mri.vol = mri.vol + mriTmp.vol.^2;
    %             end
    %             mri.vol = sqrt(mri.vol./size(fPlumb,2));
    % 
    %             MRIwrite(mri,fOut);
    %             disp('  done')
    %         else
    %             disp('  already done, skipping')
    %         end
    %         fEstim{R} = fOut;
    %     end
    % 
    % 
    %     %%%% Write means (after rms) for later visualization
    %     cmd = {srcFs};
    %     if nFrame>1
    %         fEstimAv = cell(size(fEstim));
    %     end
    %     for R = 1:length(fEstim)
    %         fIn = fEstim{R};
    %         if nFrame>1
    %             fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
    %             if forceRecomputeRMS || ~exist(fOut,'file')
    %                 cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
    %             end
    %             fEstimAv{R} = fOut;
    %             % else
    %             %     fEstimAv{R} = fIn;
    %         end
    %     end
    %     if nRun>1
    %         if nFrame>1
    %             fIn = fEstimAv;
    %         else
    %             fIn = fEstim;
    %         end
    %         fOut = replace(fIn{1},bidsList{contains(squeeze(bidsList(:,1,1)),'run-'),1,1},'run-catAv'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
    %         fOut = strsplit(fOut,filesep); fOut{end} = replace(fOut{end},'av_',''); fOut = strjoin(fOut,filesep);
    %         if forceRecomputeRMS || ~exist(fOut,'file')
    %             cmd{end+1} = ['mri_concat --o ' fOut ' ' strjoin(fIn,' ')];
    %         end
    %         fEstimCat = fOut;
    %         fIn = fOut;
    %         fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
    %         if forceRecomputeRMS || ~exist(fOut,'file')
    %             cmd{end+1} = ['mri_concat --mean --o ' fOut ' ' fIn];
    %         end
    %         fEstimAvCat = fOut;
    %     else
    %         fEstimCat = '';
    %         fEstimAvCat = '';
    %     end
    %     disp(' averaging')
    %     if length(cmd)>1
    %         cmd = strjoin(cmd,newline); % disp(cmd)
    %         [status,cmdout] = system(strjoin(cmd,newline)); if status; dbstack; error(cmdout); error('x'); end
    %     else
    %         disp('  already done, skipping')
    %     end
    end
else
    fMeanList = fList;
end


%% Full session
cmd = {srcAfni};
if any(ismember(dataType,'volTs'))
    fSummaryList = cat(2,fMeanList,fStdList,fSnrList);
else
    fSummaryList = cat(2,fMeanList);
end
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
    fCat = fOut;

    %%% runAv
    fIn = fCat;
    fOut = strsplit(fIn,filesep); fOut{end} = ['av_' fOut{end}]; fOut = strjoin(fOut,filesep);
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dTstat -overwrite -mean \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = fIn;
    end

    %%% runStd
    if any(ismember(dataType,'vol'))
        fIn = fCat;
        fOut = strsplit(fIn,filesep); fOut{end} = ['std_' fOut{end}]; fOut = strjoin(fOut,filesep);
        if force || ~exist(fOut,'file')
            cmd{end+1} = '3dTstat -overwrite -stdev \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = fIn;
        end
    end
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