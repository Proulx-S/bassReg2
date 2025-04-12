function [kI,k,hFig,fClust] = QAdendogram(fig)
    global src
        % QAdendo - Interactive dendrogram clustering visualization
        %
        % Inputs:
        %   ax1 - Axes handle for correlation matrix 
        %   ax2 - Axes handle for dendrogram
        %   Z - Linkage matrix from hierarchical clustering
        %   rho - Correlation matrix
        %
        % Controls:
        %   Left/Right Arrow - Decrease/Increase number of clusters
        %   m - Manually enter number of clusters
        %   d - Done, return final k value
        %
        % Returns:
        %   k - Final number of clusters selected
        
        disp('--------------------------------')
        disp('Adjust frame clustering:')
        disp('up/dowm arrows to increase/decrease k number of clusters')
        disp('m to enter k value')
        disp('f to identify a file')
        disp('d when done')

        hFigTmp = open(char(fig));

        hFig = figure('WindowStyle','docked');
        hFig.UserData = hFigTmp.UserData;
        ht   = tiledlayout(1,2,'TileSpacing','tight','Padding','tight');
        ax1  = copyobj(findobj(hFigTmp.Children,'Type','Axes'),ht); close(hFigTmp);
        ax1.DataAspectRatioMode = 'auto'; ax1.PlotBoxAspectRatio = [1 1 1]; tmp = cell(size(ax1.YTickLabel)); for r = 1:length(ax1.YTickLabel); tmp{r} = strsplit(ax1.YTickLabel{r},'_'); tmp{r}(~contains(tmp{r},'ses-')) = []; tmp{r} = strjoin([tmp{r} {['run-' num2str(r)]}],'_'); end; ax1.YTickLabel = tmp;
        drawnow

        % Add dendrogram
        rho = get(findobj(ax1,'Type','Image'),'CData');
        Z = linkage(squareform(1-rho), 'average'); % Convert correlation to distance 
        ax2 = nexttile;
        [H, T, perm] = dendrogram(ax2,Z, 0, 'Reorder',1:length(rho),'Orientation','right');
        ax2.YTick = []; ax2.YDir = 'reverse'; linkaxes([ax1 ax2],'y');
        

        hFigF = figure('WindowStyle','docked');
        ax1 = copyobj(ax1,hFigF); drawnow
        ax2 = copyobj(ax2,hFigF); drawnow
        ax2.Position([2 4]) = ax1.Position([2 4]);
        ax2.Position(3) = 0.25;
        hFigF.UserData = hFig.UserData; close (hFig);
        hFig = hFigF; clear hFigF

        % Set up figure and key press callback
        set(hFig, 'KeyPressFcn', @(src,event)keypress_callback(src,event,ax1,ax2,Z,rho));

        % Initialize k
        k = 2;
        kI = cluster(Z,'MaxClust',k);
        [H, T, perm] = dendrogram(ax2,Z, 0, 'Reorder',1:length(rho),'ClusterIndices',kI,'Orientation','right','ShowMarkers',true);
        ax2.YTick = []; ax2.YDir = 'reverse'; linkaxes([ax1 ax2],'y');

        % Store k in figure userdata for access in callback
        hFig.UserData.k  = k;
        hFig.UserData.kI = kI;
        hFig.UserData.done = false;

        % Display file names
        disp('---Timeseries---')
        disp(char(hFig.UserData.fileNames))
        disp('---Run averages---')
        disp(char(replace(hFig.UserData.fileNames,'preproc_volTs.nii.gz','av_preproc_volTs.nii.gz')))
        tmp = strsplit(hFig.UserData.fileNames{1,1},'_'); tmp{contains(tmp,'run-')} = 'run-cat'; tmp = strjoin(tmp,'_');
        disp('---Session catenation---')
        disp(tmp)


        % Wait for user to be done
        while ~hFig.UserData.done
            pause(0.1);
        end



        % Get final k value
        k  = hFig.UserData.k;
        kI = hFig.UserData.kI;

        % Order clusters by size
        [kIu, ~, ~] = unique(kI);
        [~, b] = sort(histcounts(kI, 1:max(kI)+1), 'descend');
        newLabels = 1:length(kIu);
        newLabels = newLabels(b);
        kI = newLabels(kI);
        
        % Output to file
        fClust = cell(size(hFig.UserData.fileNames));
        for r = 1:length(hFig.UserData.fileNames)
            fClust{r} = replace(hFig.UserData.fileNames{r},'_volTs.nii.gz','_volTsClstIdx.1D');
            writematrix(kI(hFig.UserData.fileInd==r), fClust{r},'FileType','text');
        end

        disp('--------------------------------')


        fClust = spikiness(hFig,kI)
    end






        function fClust = spikiness(hFig,kI)
            force = 1;
            if ~exist('kI','var'); kI = []; end
            if isempty( kI);       kI = ones(size(squeeze(hFig.UserData.frameNumber)))'; end
            %%% Get spikiness
            global src
            %%% Estimate spikiness
            cmd = {src.afni};
            frames = {};
            framesStr = {};
            fSpkns = {};
            fOutMeanStd = {};
            fOutMedian  = {};
            fOutMax     = {};
            for r = 1:length(hFig.UserData.fileNames)
                frames{end+1} = find(kI(hFig.UserData.fileInd==r)==1);
                framesStr{end+1} = ['[' strjoin(arrayfun(@num2str, frames{end}-1, 'UniformOutput', false), ',') ']'];
                
                fIn  = [hFig.UserData.fileNames{r} framesStr{end}];
                fSpkns{end+1} = replace(hFig.UserData.fileNames{r},'_volTs.nii.gz','_volTsSpkns.nii.gz');
                if force || ~exist(fSpkns{end},'file')
                    cmd{end+1} = ['3dDespike -overwrite -NEW \'];
                    cmd{end+1} = ['-ssave ' fSpkns{end} ' \']
                    cmd{end+1} = fIn;
                end

                fOutMeanStd{end+1} = replace(fSpkns{end},'_volTsSpkns.nii.gz','_volTsSpknsSpatialMean+Std.1D');
                fOutMedian{end+1}  = replace(fSpkns{end},'_volTsSpkns.nii.gz','_volTsSpknsSpatialMedian.1D'  );
                fOutMax{end+1}     = replace(fSpkns{end},'_volTsSpkns.nii.gz','_volTsSpknsSpatialMax.1D'     );
                if force || ~exist(fOutMeanStd{end},'file')
                    cmd{end+1} = ['3dmaskave -q -sigma  ' fSpkns{end} ' > ' fOutMeanStd{end}];
                end
                if force || ~exist(fOutMedian{end},'file')
                    cmd{end+1} = ['3dmaskave -q -median ' fSpkns{end} ' > ' fOutMedian{end} ];
                end
                if force || ~exist(fOutMax{end},'file')  % Corrected line
                    cmd{end+1} = ['3dmaskave -q -max    ' fSpkns{end} ' > ' fOutMax{end}    ];
                end
            end
            if length(cmd)>1
                [status,cmdout] = system(strjoin(cmd,newline),'-echo');
            end
    
            meanStd = cell(size(hFig.UserData.fileNames))';
            median  = cell(size(hFig.UserData.fileNames))';
            max     = cell(size(hFig.UserData.fileNames))';
            nFrames = 0;
            for r = 1:length(hFig.UserData.fileNames)
                frames{r} = frames{r} + nFrames;
                meanStd{r} = readmatrix(fOutMeanStd{r},'FileType','text')';
                median{r}  = readmatrix(fOutMedian{r},'FileType','text')';
                max{r}     = readmatrix(fOutMax{r},'FileType','text')';
                nFrames = nFrames + nnz(hFig.UserData.fileInd==r);
            end
            frames  = cat(2,frames{:});
            meanStd = cat(2,meanStd{:});
            median  = cat(2,median{:});
            max     = cat(2,max{:});

            
            axSpkn = axes(hFig,'Position',[sum(hFig.Children(1).Position([1 3])) hFig.Children(1).Position(2) 1-sum(hFig.Children(1).Position([1 3])) hFig.Children(1).Position(4)]);
            % plot(cat(1,meanStd,median,max),frames)
            plot(cat(1,meanStd),frames)
            ylim([1 nFrames])
            axSpkn.YDir = 'reverse';
            % title(legend({'mean','std','median','max'},'Location','best'),'spikiness')
            title(legend({'mean','std'},'Location','best'),'spikiness')
            linkaxes(findobj(hFig.Children,'Type','Axes'),'y'); 
        end

        







        function keypress_callback(src,event,ax1,ax2,Z,rho)
            k = src.UserData.k;
            
            switch event.Key
                case 'downarrow'  % Decrease k
                    k = max(k-1, 1);
                case 'uparrow' % Increase k 
                    k = k+1;
                case 'm' % Manual entry
                    k_str = inputdlg('Enter number of clusters:','Manual Entry',1,{num2str(k)});
                    if ~isempty(k_str)
                        k = max(1,round(str2double(k_str{1})));
                    end
                case 'f' % identify file
                    disp('Select the run you want to inspect')
                    [a,b] = ginput(1);
                    disp(src.UserData.fileNames{src.UserData.fileInd(round(a))})
                    disp(src.UserData.fileNames{src.UserData.fileInd(round(b))})
                case 'd' % Done
                    src.UserData.done = true;
                    return;
                otherwise
                    return;
            end
            
            % Update dendrogram
            cla(ax2)
            kI = cluster(Z,'MaxClust',k); 
            
            % Randomly permute the cluster labels while maintaining groups
            [uniqueK, ~, ic] = unique(kI);
            newLabels = randperm(length(uniqueK));
            kI = newLabels(ic);
            
            % Update dendrogram
            [H, T, perm] = dendrogram(ax2,Z, 0, 'Reorder',1:length(rho),'ClusterIndices',kI,'Orientation','right','ShowMarkers',true);
            ax2.YTick = [];
            ax2.YDir = 'reverse';
            linkaxes([ax1 ax2],'y');
            
            % Store new k value
            src.UserData.k  = k;
            src.UserData.kI = kI;
            
            % Update title to show current k
            title(ax2,['k = ' num2str(k)]);
        end
