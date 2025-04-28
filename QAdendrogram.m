function [kI,k,hFig,fClust,mainClust] = QAdendrogram(fig,force,verbose)
    global src
    if ~exist('force','var');     force = []; end
    if ~exist('verbose','var'); verbose = []; end
    if isempty(force);            force = 0; end
    if isempty(verbose);        verbose = 0; end
    warning('off', 'stats:linkage:NonMonotonicTree');


    %% Make or load dendrogram figure
    fFigDendro = replace(fig,'.fig','_dendro.fig');
    if ~exist(fFigDendro,'file') || force>1
        hFig = makeDendroFig(fig);
    else
        if verbose
            disp(['loading ' fFigDendro]);
            hFig = openfig(fFigDendro);
        else
            load(replace(fFigDendro,'_dendro.fig','_dendro.mat'),'kI','k','fClust','mainClust')
            hFig = [];
            return;
        end        
    end
    ax2 = findobj(hFig.Children,'type','axes');
    ttl = [ax2.Title]; ttl = {ttl.String};
    ax1 = ax2(~cellfun('isempty',ttl));
    ax2 = ax2(cellfun('isempty',ttl));

    % Display file names
    disp('---timeseries---')
    disp(char(hFig.UserData.fileNames))
    disp('---timeseries averages---')
    disp(char(hFig.UserData.fAvList))
    disp('---catenated averages---')
    disp(char(hFig.UserData.fCatAv))
    disp('---grand average---')
    disp(char(hFig.UserData.fAvCatAv))
    

    % Load saved k and kI and update dendrogram
    if force && verbose && exist(replace(fFigDendro,'_dendro.fig','_dendro.mat'),'file')
        load(replace(fFigDendro,'_dendro.fig','_dendro.mat'),'k')
        hFig.UserData.k  = k;
        hFig.UserData.kI = cluster(hFig.UserData.Z,'MaxClust',hFig.UserData.k);
        updateDendrogram(ax2, hFig.UserData.Z, hFig.UserData.k);
    end


    %% Define clustering interactively
    if verbose && force
        done = false;
        while ~done
            disp('Options:')
            disp(['Current k = ' num2str(hFig.UserData.k)])
            disp('  - Enter a number for desired clusters')
            disp('  - a/z keys to increase/decrease k by 1')
            disp('  - d when done')
            k = input('', 's');
            if strcmp(k, 'd')
                done = true;
            elseif strcmp(k, 'a') % Up arrow
                % Increase k by 1
                updateDendrogram(ax2, hFig.UserData.Z, hFig.UserData.k + 1);
            elseif strcmp(k, 'z') % Down arrow
                % Decrease k by 1
                updateDendrogram(ax2, hFig.UserData.Z, max(1, hFig.UserData.k - 1));
            elseif ~isnan(str2double(k)) && ~imag(str2double(k))
                % Try to convert to number
                updateDendrogram(ax2, hFig.UserData.Z, str2double(k));
            else
                % Change cluster coloring
                updateDendrogram(ax2, hFig.UserData.Z, hFig.UserData.k);
            end
        end
    end



    %% Finalize outputs
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

    % Output each run main cluster index
    mainClust = nan(size(hFig.UserData.fileNames));
    for r = 1:length(hFig.UserData.fileNames)
        mainClust(r) = mode(kI(hFig.UserData.fileInd==r));
    end

    % Save figure
    if force || ~exist(fFigDendro,'file')
        disp(['saving ' fFigDendro])
        saveas(hFig,fFigDendro)
        save(replace(fFigDendro,'_dendro.fig','_dendro.mat'),'kI','k','fClust','mainClust')
    end
    disp('--------------------------------')




    function hFig = makeDendroFig(fig)
        %% Figure setup
        % Open xCorr figure
        hFig = openfig(char(fig),'invisible');
        ax0 = findobj(hFig.Children,'Type','Axes');
        delete(findobj(hFig.Children,'Type','Colorbar'));
        drawnow

        % Copy to temporary figure to use tiledlayout
        hFigHt = figure('WindowStyle','docked');
        ht  = tiledlayout(1,2,'TileSpacing','tight','Padding','tight');
        ax1 = copyobj(findobj(hFig.Children,'Type','Axes'),ht);
        hIm = findobj(ax1.Children,'Type','Image');
        ax1.YTick = [];
        ax1.XTick = [];
        i = 1:size(hIm.CData,1);
        yTick              = [];
        task               = cell(size(hFig.UserData.fileNames));
        ses                = cell(size(hFig.UserData.fileNames));
        hFig.UserData.task = cell(size(hFig.UserData.fileNames));
        for f = 1:length(hFig.UserData.fileNames)
            yTick = [yTick; mean(i(hFig.UserData.fileInd==f))];
            bids{f} = strsplit(hFig.UserData.fileNames{f},filesep);
            bids{f} = strsplit(strjoin(bids{f},'_'),'_');
            task{f} = bids{f}{contains(bids{f},'task')};
            hFig.UserData.task{f} = task{f};
            task{f} = strsplit(task{f},'-');
            task{f} = task{f}{end};
            ses{f}  = {bids{f}{contains(bids{f},'ses')}};
            ses{f}  = ses{f}{1};
        end
        ax1.YTick      = yTick;
        ax1.YTickLabel = task;
        ax1.XTick      = yTick;
        ax1.XTickLabel = ses;
        drawnow

        % Add dendrogram
        rho = get(findobj(ax1,'Type','Image'),'CData');
        hFig.UserData.Z = linkage(squareform(1-rho), 'average'); % Convert correlation to distance 
        ax2 = nexttile(ht);
        [H, T, perm] = dendrogram(ax2,hFig.UserData.Z, 0, 'Reorder',1:length(rho),'Orientation','right');
        ax2.YDir = 'reverse';
        ax2.YTick = []; ax2.XTick = [];
        drawnow

        % Copy back to main figure without tiledlayout
        figure(hFig)
        delete(hFig.Children)
        ax1 = copyobj(ax1,hFig); drawnow
        ax2 = copyobj(ax2,hFig); drawnow
        delete(hFigHt)
        figure(hFig)
        ax2.Position([2 4]) = ax1.Position([2 4]);
        linkaxes([ax1 ax2],'y');
        hFig.Visible = 'on';
        drawnow

        hFig.UserData.k = 1;
        hFig.UserData.kI = ones(size(length(rho),1));

        % % Display file names
        % disp('---timeseries---')
        % disp(char(hFig.UserData.fileNames))
        % disp('---timeseries averages---')
        % disp(char(hFig.UserData.fAvList))
        % disp('---catenated averages---')
        % disp(char(hFig.UserData.fCatAv))
        % disp('---grand average---')
        % disp(char(hFig.UserData.fAvCatAv))

    
    function updateDendrogram(ax,Z,k)
        warning('off', 'stats:linkage:NonMonotonicTree');
        % Update clustering if k has changed
        if k~=ax.Parent.UserData.k
            ax.Parent.UserData.kI = cluster(Z,'MaxClust',k);
            ax.Parent.UserData.k  = k;
        end
        % Randomize coloring using kIrand
        kU     = unique(ax.Parent.UserData.kI);
        kMap   = containers.Map(kU, kU(randperm(length(kU))));
        kIrand = ax.Parent.UserData.kI;
        for i = 1:length(kU); kIrand(ax.Parent.UserData.kI==kU(i)) = kMap(kU(i)); end
        % Update dendrogram
        dendrogram(ax,Z, 0, 'Reorder',1:size(Z,1)+1,'Orientation','right','ClusterIndices',kIrand,'ShowCut',true);
        set(ax,'YDir','reverse','YTick',[],'XTick',[]);
        drawnow

        % clusterAssignments = cluster(Z,Cutoff=18,Criterion="inconsistent");
        % dendrogram(ax,Z, 0, 'Reorder',1:size(Z,1)+1,'Orientation','right','ClusterIndices',clusterAssignments,'ShowCut',true);
        % % set(ax,'YDir','reverse','YTick',[],'XTick',[]);
        % set(ax,'YDir','reverse','YTick',[]);


        