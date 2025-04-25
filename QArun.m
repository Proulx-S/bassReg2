function QArun(rSet,fMask,derivDir,force,verbose)
    if ~exist('force','var');     force = []; end
    if ~exist('verbose','var'); verbose = []; end
    if isempty(force);            force = 0; end
    if isempty(verbose);        verbose = 0; end
    % if force;                   verbose = 1; end


    % Load mask and conform mask
    mriConform = MRIread(rSet.fPreprocList{1},1);

    imMask = MRIread(fMask);
    if ~contains(fMask,'setPlumb_'); dbstack; error('sdfg'); end; fMask = replace(fMask,'setPlumb_','setOblique_');
    mriConform.vol = logical(imMask.vol); MRIwrite(mriConform,fMask);
    imMask = mriConform.vol;

    fMaskInv = replace(fMask,'.nii.gz','Inv.nii.gz');
    mriConform.vol = -(mriConform.vol-1);
    MRIwrite(mriConform,fMaskInv);

    disp(fMaskInv)
    


    %% Correct bias field for each run
    forceThis   = force;
    verboseThis = verbose;
    f          = cell(size(rSet.fPreprocList(:,1,1)));
    fOrig      = cell(size(rSet.fPreprocList(:,1,1)));
    fVolCorr   = cell(size(rSet.fPreprocList(:,1,1)));
    fVolTsCorr = cell(size(rSet.fPreprocList(:,1,1)));
    fVol       = cell(size(rSet.fPreprocList(:,1,1)));
    fVolField  = cell(size(rSet.fPreprocList(:,1,1)));
    
    for R = 1:size(rSet.fPreprocList,1)
        %% Define files
        f{R}     = rSet.fPreprocList{R,1,1};
        fOrig{R} = rSet.fOrigList{R,1,1};
        if contains(rSet.label,'vfMRI')
            [fVolCorr{R},fVolTsCorr{R},fVol{R},fVolField{R}] = correctBiasField(f{R}, fMask, f{R}, fOrig{R}, forceThis, verboseThis);
        end
    end

    %% Segment on voxels pooled across runs
    forceThis = force;
    verboseThis = verbose;
    if contains(rSet.label,'vfMRI')

        % test if all run uses the same acquisition condition
        [~,tmp,~] = fileparts(fileparts(fVolCorr));
        for R = 1:size(rSet.fPreprocList,1);
            tmp{R} = strsplit(tmp{R},'_');
            tmp{R}(contains(tmp{R},'run-' )) = [];
            tmp{R}(contains(tmp{R},'task-')) = [];
            tmp{R} = strjoin(tmp{R},'_');
        end
        [tmpU,~,c] = unique(tmp);

        if length(tmpU)==length(tmp) && length(tmp)~=1; dbstack; error('asedfarfge'); end
        
        % if not all the same, split into subgroup and computeVesselness on each
        fSegMask    = cell(size(tmpU));
        fNonSegMask = cell(size(tmpU));
        fVolCorr2   = cell(size(tmpU));
        fSegFig     = cell(size(tmpU));
        for i = 1:length(tmpU)
            [fSegMask{i},fNonSegMask{i},fSegFig{i},fVolCorr2{i}] = computeVesselness(fVolCorr(c==i,:),fMask,forceThis,verboseThis);
        end

        % catenate back all runs and sort to original
        fSegMask    = cat(1,fSegMask{:});
        fNonSegMask = cat(1,fNonSegMask{:});
        fVolCorr2   = cat(1,fVolCorr2{:});
        [~,b] = ismember(fVolCorr(:,1),fVolCorr2)
        fVolCorr2(b)
        fSegMask    = fSegMask(b);
        fNonSegMask = fNonSegMask(b);
        
        % [fSegMask,fNonSegMask,fSegFig,fVolCorr2] = computeVesselness(fVolCorr,fMask,forceThis,verboseThis);
    end





    hFig    = cell([size(rSet.fPreprocList,1) 1]);
    ht      = cell([size(rSet.fPreprocList,1) 1]);
    axCorr  = cell([size(rSet.fPreprocList,1) 1]);
    axCorr2 = cell([size(rSet.fPreprocList,1) 1]);
    axSpkns = cell([size(rSet.fPreprocList,1) 1]);
    for R = 1:size(rSet.fPreprocList,1)
        %% Define files
        % f       = rSet.fPreprocList{R,1,1};
        % fOrig   = rSet.fOrigList{R,1,1};
        fMcWR   = replace(rSet.fTransList(R,:,:),'.aff12.','.param.');
        fMcWR   = fMcWR{find(contains(fMcWR,'mcWR'),1,"first")};
        fCensor = strsplit(replace(f{R},'.nii.gz', '.csv'),filesep); fCensor{end} = ['censor_' fCensor{end}]; fCensor = strjoin(fCensor,filesep);
        fQA     = strsplit(f{R},filesep); fQA{end} = ['QA_' replace(fQA{end},'.nii.gz','.fig')]; fQA = strjoin(fQA,filesep);

        fCensorDeriv = strsplit(fCensor,filesep); fCensorDeriv = strjoin(fCensorDeriv(find(ismember(fCensorDeriv,rSet.label)):length(fCensorDeriv)),filesep); fCensorDeriv = fullfile(derivDir,fCensorDeriv);
        if ~exist(fileparts(fCensorDeriv),'dir'); mkdir(fileparts(fCensorDeriv)); end

        disp('-------')
        disp(['Current file:' newline f{R} newline fCensor])
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
            ht{R} = tiledlayout(1,9); ht{R}.TileSpacing = "none"; ht{R}.Padding = 'none';
    
    
            %% Correlate frames
            im = MRIread(f{R});
            im = permute(im.vol,[4 1 2 3]);
            rho = corr(permute(im(:,imMask),[2 1]));

            axCorr{R} = nexttile([1 5]);
            imagesc(rho,[0 1]);
            
            axCorr{R}.DataAspectRatio = [1 1 1];
            axCorr{R}.YAxisLocation = 'right';
            ylabel(colorbar('Location', 'westoutside'),'cross-frame Pearson''s correlation')
            axCorr{R}.XTick = [];
            ylabel('frame indices')
            [~,b,~] = fileparts(fileparts(f{R}));
            xlabel(axCorr{R},b,'Interpreter','none')
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
                % forceThis = force;
                % verboseThis = verbose;
                % %%% Bias Field Correction        
                % [fVolCorr,fVolTsCorr,fVol,fVolField] = correctBiasField(f, fMask, f, fOrig, forceThis, verboseThis);
                
                % forceThis = force;
                % verboseThis = verbose;
                % %%% Bias Field Correction
                % [fSegMask,fNonSegMask,fSegFig] = computeVesselness(fVolCorr{R},fMask,forceThis,verboseThis);
                fNonVesselMask{R,1} = fNonSegMask{R,contains(fNonSegMask(R,:),'_nonvesselSegMask.nii.gz')};
                fVesselMask{R,1}    = fSegMask{   R,contains(fNonSegMask(R,:),'_nonvesselSegMask.nii.gz')};
                disp('!!!!!!!')
                disp(['inspect mask: ' fNonVesselMask{R,1}])
                fVolTs_spkns = char(fVolTsCorr{R});
                fMask_spkns  = fNonVesselMask{R,1};
            else
                fVolTs_spkns = char(f{R});
                fMask_spkns  = fMask;
            end
            

            %%% Spikiness ts        
            spkns     = MRIread(fVolTs_spkns); spkns     = permute(spkns.vol,[4 1 2 3]);
            spknsMask = MRIread(fMask_spkns);  spknsMask = logical(spknsMask.vol);
            spknsMean = mean(spkns(:,spknsMask)   ,2);
            spknsStd  = std( spkns(:,spknsMask),[],2); clear spkns


            drawnow
            try
                axSpkns{R} = nexttile(ht{R},[1 1]);
            catch
                disp('??????')
                keyboard
                axSpkns{R} = nexttile(ht{R},[1 1]);
            end
            % plot(spknsMean - mean(spknsMean),1:rSet.nFrame(R))
            plot(zscore(spknsMean),1:rSet.nFrame(R))
            axSpkns{R}.YDir = 'reverse';
            ylim(axCorr{R}.YLim)
            xlabel('magnitude spikiness');
            axSpkns{R}.YTickLabel = [];
            grid on
            hold on
            % plot(spknsStd - mean(spknsStd),1:rSet.nFrame(R))
            plot(zscore(spknsStd),1:rSet.nFrame(R))
            xlim([-5 5])
            axSpkns{R}.XTick = [-2.5 0 2.5];
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
            axMC{R} = nexttile(ht{R},[1 1]);
            plot(mcWR - mean(mcWR,1),1:rSet.nFrame(R))
            axMC{R}.YDir = 'reverse';
            ylim(axCorr{R}.YLim)
            xlabel('motion (mm/deg)');
            axMC{R}.YTickLabel = [];
            grid on
            xlim([-1 1])
            axMC{R}.XTick = [-0.5 0 0.5];
            legend(axMC{R},mcLabel,'box','off')
            




            %%% Display file names
            disp('!!!!!!!!!!!')
            if contains(rSet.label,'vfMRI')
                disp(strjoin({
                    f{R}
                    fVolCorr{R}
                    fNonVesselMask{R,1}
                    fCensor
                },newline))
            else
                disp(strjoin({
                    f{R}
                    fCensor
                },newline))
            end
            disp('!!!!!!!!!!!')
            

            
            %%% create censor file for manual censor point identification
                if force>1 || ~exist(fCensor, 'file')
                writematrix([(1:rSet.nFrame(R))' ones(rSet.nFrame(R),1)], fCensor, 'Delimiter', ',');
            end

            if verbose
                disp('!!!!!!!!!!!')
                disp(['enter censor points (0) in: ' newline fCensor])
                disp('then type "done"')
                disp('!!!!!!!!!!!')
                while ~strcmpi(input(' ', 's'), 'done')
                    disp('type "done" to continue')
                end
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
            hFig{R}.UserData.fileNames.f     = f{R};
            hFig{R}.UserData.fileNames.fOrig = fOrig{R};
            hFig{R}.UserData.fileNames.fMcWR = fMcWR;
            hFig{R}.UserData.nDummy          = rSet.nFrameOrig(R) - rSet.nFrame(R);
            hFig{R}.UserData.fileInd         = R;
            hFig{R}.UserData.censoredPoints  = cnsr;
            if exist('fSegFig','var')
                hFig{R}.UserData.fSegFig         = fSegFig;
            else
                hFig{R}.UserData.fSegFig         = [];
            end
            
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

                if verbose>1
                    if ~isempty(hFig{R}.UserData.fSegFig) && R==1
                        if iscell(hFig{R}.UserData.fSegFig)
                            for iii = 1:length(hFig{R}.UserData.fSegFig)
                                open(hFig{R}.UserData.fSegFig{iii})
                            end
                            pause(1)
                        end
                    end
                    figure(hFig{R})

                end


                if verbose>2
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
                end


    

                % hFig{R}.UserData.censoredPoints = cnsr;
                % saveas(hFig{R},fQA);
            

                %% Save censor to permanent bids derivatives directory
                if exist(fCensorDeriv,'file')
                    if~strcmp(fileread(fCensor),fileread(fCensorDeriv))
                        disp('!!!!!!!!!!')
                        disp('!!!!!!!!!!')
                        disp('!!!!!!!!!!')
                        disp('censor points seem to have been updated')
                        disp(['Overwrite file: ' newline fCensorDeriv newline 'with' newline fCensor newline '? (y/n)'])
                        if strcmpi(input(' ','s'),'y')
                            disp('overwriting')
                            copyfile(fCensor,fCensorDeriv)
                        else
                            disp('not saving changes (yolo)')
                        end
                    end
                else
                    copyfile(fCensor,fCensorDeriv)
                end
            end
        end
    end






