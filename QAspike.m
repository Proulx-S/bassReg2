function fNonVesselMask = QAspike(fVolTs,fMask,force)
    global src
    if ~exist('force','var'); force = []; end
    if isempty(force);        force = 0 ; end

    %% Bias Field Correction
    [fVolCorr,fVolTsCorr,fVol,fVolField] = correctBiasField(fVolTs, fMask, force);

    %% Compute spikiness
    fSpkns = computeSpkns(fVolTsCorr,force);

    %% Identify and mask out vessels and non-brain (is vfMRI because vessels show large physiological variations)
    [fVesselMask,fNonVesselMask] = computeVesselness(fVolCorr,fMask,force);

    %% Summarize each volume spikiness into a spikiness timeseries
    [spknsTimeseries,fSpknsTimeseries] = summarizeSpkns(fSpkns,fNonVesselMask,force)

    figure('WindowStyle','docked');
    plot(spknsTimeseries)
    

    
    %% Label spiky volumes










    % return
    %     % QAdendo - Interactive dendrogram clustering visualization
    %     %
    %     % Inputs:
    %     %   ax1 - Axes handle for correlation matrix 
    %     %   ax2 - Axes handle for dendrogram
    %     %   Z - Linkage matrix from hierarchical clustering
    %     %   rho - Correlation matrix
    %     %
    %     % Controls:
    %     %   Left/Right Arrow - Decrease/Increase number of clusters
    %     %   m - Manually enter number of clusters
    %     %   d - Done, return final k value
    %     %
    %     % Returns:
    %     %   k - Final number of clusters selected
        
    %     disp('--------------------------------')
    %     disp('Adjust frame clustering:')
    %     disp('up/dowm arrows to increase/decrease k number of clusters')
    %     disp('m to enter k value')
    %     disp('f to identify a file')
    %     disp('d when done')

    %     hFigTmp = open(char(fig));

    %     hFig = figure('WindowStyle','docked');
    %     hFig.UserData = hFigTmp.UserData;
    %     ht   = tiledlayout(1,2,'TileSpacing','tight','Padding','tight');
    %     ax1  = copyobj(findobj(hFigTmp.Children,'Type','Axes'),ht); close(hFigTmp);
    %     ax1.DataAspectRatioMode = 'auto'; ax1.PlotBoxAspectRatio = [1 1 1]; tmp = cell(size(ax1.YTickLabel)); for r = 1:length(ax1.YTickLabel); tmp{r} = strsplit(ax1.YTickLabel{r},'_'); tmp{r}(~contains(tmp{r},'ses-')) = []; tmp{r} = strjoin([tmp{r} {['run-' num2str(r)]}],'_'); end; ax1.YTickLabel = tmp;
    %     drawnow

    %     % Add dendrogram
    %     rho = get(findobj(ax1,'Type','Image'),'CData');
    %     Z = linkage(squareform(1-rho), 'average'); % Convert correlation to distance 
    %     ax2 = nexttile;
    %     [H, T, perm] = dendrogram(ax2,Z, 0, 'Reorder',1:length(rho),'Orientation','right');
    %     ax2.YTick = []; ax2.YDir = 'reverse'; linkaxes([ax1 ax2],'y');
        

    %     hFigF = figure('WindowStyle','docked');
    %     ax1 = copyobj(ax1,hFigF); drawnow
    %     ax2 = copyobj(ax2,hFigF); drawnow
    %     ax2.Position([2 4]) = ax1.Position([2 4]);
    %     ax2.Position(3) = 0.25;
    %     hFigF.UserData = hFig.UserData; close (hFig);
    %     hFig = hFigF; clear hFigF

    %     % Set up figure and key press callback
    %     set(hFig, 'KeyPressFcn', @(src,event)keypress_callback(src,event,ax1,ax2,Z,rho));

    %     % Initialize k
    %     k = 2;
    %     kI = cluster(Z,'MaxClust',k);
    %     [H, T, perm] = dendrogram(ax2,Z, 0, 'Reorder',1:length(rho),'ClusterIndices',kI,'Orientation','right','ShowMarkers',true);
    %     ax2.YTick = []; ax2.YDir = 'reverse'; linkaxes([ax1 ax2],'y');

    %     % Store k in figure userdata for access in callback
    %     hFig.UserData.k  = k;
    %     hFig.UserData.kI = kI;
    %     hFig.UserData.done = false;

    %     % Display file names
    %     disp('---Timeseries---')
    %     disp(char(hFig.UserData.fileNames))
    %     disp('---Run averages---')
    %     disp(char(replace(hFig.UserData.fileNames,'preproc_volTs.nii.gz','av_preproc_volTs.nii.gz')))
    %     tmp = strsplit(hFig.UserData.fileNames{1,1},'_'); tmp{contains(tmp,'run-')} = 'run-cat'; tmp = strjoin(tmp,'_');
    %     disp('---Session catenation---')
    %     disp(tmp)


    %     % Wait for user to be done
    %     while ~hFig.UserData.done
    %         pause(0.1);
    %     end



    %     % Get final k value
    %     k  = hFig.UserData.k;
    %     kI = hFig.UserData.kI;

    %     % Order clusters by size
    %     [kIu, ~, ~] = unique(kI);
    %     [~, b] = sort(histcounts(kI, 1:max(kI)+1), 'descend');
    %     newLabels = 1:length(kIu);
    %     newLabels = newLabels(b);
    %     kI = newLabels(kI);
        
    %     % Output to file
    %     fClust = cell(size(hFig.UserData.fileNames));
    %     for r = 1:length(hFig.UserData.fileNames)
    %         fClust{r} = replace(hFig.UserData.fileNames{r},'_volTs.nii.gz','_volTsClstIdx.1D');
    %         writematrix(kI(hFig.UserData.fileInd==r), fClust{r},'FileType','text');
    %     end

    %     disp('--------------------------------')


    %     fClust = spikiness(hFig,kI)
    % end






    function fSpkns = computeSpkns(f,force)
        global src
        if ~ischar(f); dbstack; error('code that'); end
        if ~exist('force','var'); force = []; end
        if isempty(force);        force = 0 ; end
        cmd = {src.afni};
        
        %%% Compute spikiness
        fSpkns = replace(f,'_volTs.nii.gz','_volTsSpkns.nii.gz');
        if force || ~exist(fSpkns,'file')
            cmd{end+1} = ['3dDespike -overwrite -NEW \'];
            cmd{end+1} = ['-ssave ' fSpkns ' \'];
            cmd{end+1} = f;
        end

        %%% Execute
        if length(cmd)>1
            [status,cmdout] = system(strjoin(cmd,newline),'-echo');
        end
    



    function [spknsMeanStd,fSpknsMeanStd] = summarizeSpkns(fSpkns,fMask,force)
        global src
        if ~ischar(fSpkns); dbstack; error('code that'); end
        if ~exist('force','var'); force = []; end
        if isempty(force);        force = 0 ; end
        cmd = {src.afni};
        

        %%% Summarize across space
        fSpknsMeanStd = replace(fSpkns,'_volTsSpkns.nii.gz','_volTsSpknsSpatialMean+Std.1D');
        if force || ~exist(fSpknsMeanStd,'file')
            cmd{end+1} = ['3dmaskave -q -sigma -mask ' fMask ' ' fSpkns ' > ' fSpknsMeanStd];
        end
        
        %%% Execute
        if length(cmd)>1
            [status,cmdout] = system(strjoin(cmd,newline),'-echo');
        end

        %%% Read back
        spknsMeanStd = readmatrix(fSpknsMeanStd,'FileType','text');
    



    


        







        % function keypress_callback(src,event,ax1,ax2,Z,rho)
        %     k = src.UserData.k;
            
        %     switch event.Key
        %         case 'downarrow'  % Decrease k
        %             k = max(k-1, 1);
        %         case 'uparrow' % Increase k 
        %             k = k+1;
        %         case 'm' % Manual entry
        %             k_str = inputdlg('Enter number of clusters:','Manual Entry',1,{num2str(k)});
        %             if ~isempty(k_str)
        %                 k = max(1,round(str2double(k_str{1})));
        %             end
        %         case 'f' % identify file
        %             disp('Select the run you want to inspect')
        %             [a,b] = ginput(1);
        %             disp(src.UserData.fileNames{src.UserData.fileInd(round(a))})
        %             disp(src.UserData.fileNames{src.UserData.fileInd(round(b))})
        %         case 'd' % Done
        %             src.UserData.done = true;
        %             return;
        %         otherwise
        %             return;
        %     end
            
        %     % Update dendrogram
        %     cla(ax2)
        %     kI = cluster(Z,'MaxClust',k); 
            
        %     % Randomly permute the cluster labels while maintaining groups
        %     [uniqueK, ~, ic] = unique(kI);
        %     newLabels = randperm(length(uniqueK));
        %     kI = newLabels(ic);
            
        %     % Update dendrogram
        %     [H, T, perm] = dendrogram(ax2,Z, 0, 'Reorder',1:length(rho),'ClusterIndices',kI,'Orientation','right','ShowMarkers',true);
        %     ax2.YTick = [];
        %     ax2.YDir = 'reverse';
        %     linkaxes([ax1 ax2],'y');
            
        %     % Store new k value
        %     src.UserData.k  = k;
        %     src.UserData.kI = kI;
            
        %     % Update title to show current k
        %     title(ax2,['k = ' num2str(k)]);
        % end
