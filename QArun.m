function QArun(rSet,fMask,derivDir,force,verbose)
    if ~exist('force','var');     force = []; end
    if ~exist('verbose','var'); verbose = []; end
    if isempty(force);            force = 0; end
    if isempty(verbose);        verbose = 0; end



    imMask = MRIread(fMask); imMask = logical(imMask.vol);


    hFig    = cell([size(rSet.fPreprocList,1) 1]);
    ht      = cell([size(rSet.fPreprocList,1) 1]);
    axCorr  = cell([size(rSet.fPreprocList,1) 1]);
    axCorr2 = cell([size(rSet.fPreprocList,1) 1]);
    axSpkns = cell([size(rSet.fPreprocList,1) 1]);
    for R = 1:size(rSet.fPreprocList,1)
        %% Define files
        f       = rSet.fPreprocList{R,1,1};
        fOrig   = rSet.fOrigList{R,1,1};
        fMcWR   = replace(rSet.fTransList(R,:,:),'.aff12.','.param.');
        fMcWR   = fMcWR{find(contains(fMcWR,'mcWR'),1,"first")};
        fCensor = strsplit(replace(f,'.nii.gz', '.csv'),filesep); fCensor{end} = ['censor_' fCensor{end}]; fCensor = strjoin(fCensor,filesep);
        fQA     = strsplit(f,filesep); fQA{end} = ['QA_' replace(fQA{end},'.nii.gz','.fig')]; fQA = strjoin(fQA,filesep);

        fCensorDeriv = strsplit(fCensor,filesep); fCensorDeriv = strjoin(fCensorDeriv(find(ismember(fCensorDeriv,rSet.label)):length(fCensorDeriv)),filesep); fCensorDeriv = fullfile(derivDir,fCensorDeriv);
        if ~exist(fileparts(fCensorDeriv),'dir'); mkdir(fileparts(fCensorDeriv)); end

        disp('-------')
        disp(['Current file:' newline f])
        disp('-------')
    

        %% Restore bidsDerive file if exists
        if exist(fCensorDeriv,'file')
            copyfile(fCensorDeriv,fCensor)
        end

  
        


        if force || ~exist(fCensor,'file') || ~exist(fQA,'file')




            if verbose
                hFig{R} = figure('WindowStyle','docked');
            else
                hFig{R} = figure('WindowStyle','docked','Visible','off');
            end
            ht{R} = tiledlayout(1,8); ht{R}.TileSpacing = "none"; ht{R}.Padding = 'none';
    
    
            %% Correlate frames
            im = MRIread(f);
            im = permute(im.vol,[4 1 2 3]);
            rho = corr(permute(im(:,imMask),[2 1]));

            axCorr{R} = nexttile([1 5]);
            imagesc(rho,[0 1]);
            
            axCorr{R}.DataAspectRatio = [1 1 1];
            % axCorr{R}.YAxisLocation = 'right';
            ylabel(colorbar('Location', 'westoutside'),'cross-frame Pearson''s correlation')
            axCorr{R}.XTick = [];
            ylabel('frame indices')
            [~,b,~] = fileparts(fileparts(f));
            ylabel(ht{R},b,'Interpreter','none')
            drawnow

            axCorr2{R} = nexttile([1 1]);
            plot(mean(rho),1:rSet.nFrame(R))
            axCorr2{R}.YDir = 'reverse';
            xlabel('mean correlation');
            axCorr2{R}.YTickLabel = [];
            ylim(axCorr{R}.YLim)
            grid on
            drawnow
            


            %% Compute spikiness
            if contains(rSet.label,'vfMRI')
                %%% Bias Field Correction        
                [fVolCorr,fVolTsCorr,fVol,fVolField] = correctBiasField(f, fMask, fOrig, force, 0);
                
                %%% Bias Field Correction
                [fVesselMask,fNonVesselMask] = computeVesselness(fVolCorr,fMask,force,0);
                disp('!!!!!!!')
                disp(['inspect mask: ' fNonVesselMask])
                fVolTs_spkns = fVolTsCorr;
                fMask_spkns  = fNonVesselMask;
            else
                fVolTs_spkns = f;
                fMask_spkns  = fMask;
            end
            

            %%% Spikiness ts        
            spkns         = MRIread(fVolTs_spkns);             spkns = permute(spkns.vol,[4 1 2 3]);
            spknsMask = MRIread(fMask_spkns); spknsMask = logical(spknsMask.vol);
            spknsMean = mean(spkns(:,spknsMask)   ,2);
            spknsStd  = std( spkns(:,spknsMask),[],2); clear spkns


            drawnow
            try
                axSpkns{R} = nexttile([1 1]);
            catch
                disp('??????')
                keyboard
                axSpkns{R} = nexttile([1 1]);
            end
            plot(spknsMean - mean(spknsMean),1:rSet.nFrame(R))
            axSpkns{R}.YDir = 'reverse';
            ylim(axCorr{R}.YLim)
            xlabel('magnitude spikiness');
            axSpkns{R}.YTickLabel = [];
            grid on
            hold on
            plot(spknsStd - mean(spknsStd),1:rSet.nFrame(R))
            legend(axSpkns{R},{'spatial mean' 'spatial std'},'box','off')
            
            
            %% Show movement paramters
            %%% read motion parameters
            fid = fopen(fMcWR); fgetl(fid);
            mcLabel = fgetl(fid); mcLabel = strsplit(replace(mcLabel,'#',''),' '); mcLabel(cellfun('isempty',mcLabel)) = [];
            fclose(fid);
            mcWR = readmatrix(fMcWR,'FileType','text');
            mcWR(:,contains(mcLabel,'$')) = [];
            mcLabel(:,contains(mcLabel,'$')) = [];
            
            %%% plot motion parameters
            axMC{R} = nexttile([1 1]);
            plot(mcWR - mean(mcWR,1),1:rSet.nFrame(R))
            axMC{R}.YDir = 'reverse';
            ylim(axCorr{R}.YLim)
            xlabel('motion (mm/deg)');
            axMC{R}.YTickLabel = [];
            grid on
            legend(axMC{R},mcLabel,'box','off')
            




            %%% Display file names
            disp('!!!!!!!!!!!')
            if contains(rSet.label,'vfMRI')
                disp(strjoin({
                    f
                    fVolCorr
                    fNonVesselMask
                    fCensor
                },newline))
            else
                disp(strjoin({
                    f
                    fCensor
                },newline))
            end
            disp('!!!!!!!!!!!')
            

            
            %%% create censor file for manual censor point identification
            if force>1 || ~exist(fCensor, 'file')
                writematrix([(1:rSet.nFrame(R))' ones(rSet.nFrame(R),1)], fCensor, 'Delimiter', ',');
            end
            disp('!!!!!!!!!!!')
            disp(['enter censor points (0) in: ' newline fCensor])
            disp('then type "done"')
            disp('!!!!!!!!!!!')
            while ~strcmpi(input(' ', 's'), 'done')
                disp('type "done" to continue')
            end


            
            %%% visualize censored points
            cnsr = readmatrix(fCensor, 'Delimiter', ',');
            ax = findobj(hFig{R}.Children.Children,'Type','axes');
            lgnd = findobj(hFig{R}.Children.Children,'Type','legend');
            set(lgnd,'AutoUpdate','off')
            if ~all(cnsr(:,2))
                for i = 1:length(ax)
                    yline(ax(i),find(~cnsr(:,2)),'r')
                end
            end




            %% Store metadata for later use
            hFig{R}.UserData.fileNames.f = f;
            hFig{R}.UserData.fileNames.fOrig = fOrig;
            hFig{R}.UserData.fileNames.fMcWR = fMcWR;
            hFig{R}.UserData.nDummy      = rSet.nFrameOrig(R) - rSet.nFrame(R);
            hFig{R}.UserData.fileInd     = R;
            hFig{R}.UserData.censoredPoints = cnsr;
            
            %% Save figure
            if ~verbose
                set(hFig{R}, 'CreateFcn', @(obj, ~) set(obj, 'Visible', 'on'));
            end
            if force || ~exist(fQA,'file')
                saveas(hFig{R},fQA);
            end


            
            %% Save censor to permanent bids derivatives directory
            copyfile(fCensor,fCensorDeriv)
            



        else



            if verbose


                %% Load and upadte censor
                hFig{R} = open(fQA);
                disp('!!!!!!!!!!!')
                disp(['enter censor points (0) in: ' newline fCensor])
                disp('then type "done"')
                disp('!!!!!!!!!!!')
                while ~strcmpi(input(' ', 's'), 'done')
                    disp('type "done" to continue')
                end

                cnsr = readmatrix(fCensor, 'Delimiter', ',');
                ax = findobj(hFig{R}.Children.Children,'Type','axes');
                lgnd = findobj(hFig{R}.Children.Children,'Type','legend');
                set(lgnd,'AutoUpdate','off')

                for i = 1:length(ax)
                    delete(findobj(ax(i).Children,'Type','ConstantLine'));
                end

                if ~all(cnsr(:,2))
                    for i = 1:length(ax)
                        yline(ax(i),find(~cnsr(:,2)),'r')
                    end
                end

                hFig{R}.UserData.censoredPoints = cnsr;
                saveas(hFig{R},fQA);
            

                %% Save censor to permanent bids derivatives directory
                if exist(fCensorDeriv,'file')    
                    disp('!!!!!!!!!!')
                    disp('!!!!!!!!!!')
                    disp('!!!!!!!!!!')
                    disp(['Overwrite file: ' newline fCensorDeriv newline 'with' newline fCensor newline '? (y/n)'])
                    if strcmpi(input(' ','s'),'y')
                        disp('overwriting')
                        copyfile(fCensor,fCensorDeriv)
                    else
                        disp('skipping')
                    end
                else
                    copyfile(fCensor,fCensorDeriv)
                end


            end






        end




        
    end






