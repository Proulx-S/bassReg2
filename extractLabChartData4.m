function phys = extractLabChartData4(physFileDb,mriFiles,physOutDir,force)
if ~exist('force','var') || isempty(force); force = 0; end
if ~exist('powerThrough','var'); powerThrough = []; end
if isempty(powerThrough); powerThrough = 0; end
powerThrough = 1;
funcFileAcqLabel = {'vfMRI' 'vfMRIpc' 'bold'};

%% Initiate info
info = doIt; if nargin==0; return; end
info.dbFile =  physFileDb;
info.outFile = fullfile(physOutDir,'minCurated.mat');

% If output file exists, load and exit
if exist(info.outFile,'file') && ~force
    load(info.outFile);
    if strcmp(phys.info.f1.Visible,'on')
        phys.info.f1.Visible = 'off';
        phys.info.f2.Visible = 'off';
        save(info.outFile,'phys');
    end
    return
end
% if exist(info.outFile,'file') && ~force
%     load(info.outFile);
%     phys.info.f1.Visible = 'off';
%     phys.info.f2.Visible = 'off';
%     save(info.outFile,'phys');
%     return
% end




%% Manual identification of runs on physio chanels
fManId = fullfile(fileparts(physFileDb),'manual','runAndChanId.mat');
if force>1 || ~exist(fManId,'file')
    if ~exist(fileparts(fManId),'dir'); mkdir(fileparts(fManId)); end
    info = manId(info,physFileDb);
    save(fManId,'info')
else
    load(fManId,'info')
end


%% Synchronize physio to mri
% get mriFiles in proper format and get runCond to runRun
if iscell(mriFiles)
    if ~ischar([mriFiles{:}])
        if isa([mriFiles{:}],'runCond')
            rCond = [mriFiles{:}]';
            rCond = rCond(ismember({rCond.acq},funcFileAcqLabel));
            runRun = {};
            for rc = 1:length(rCond)
                for rr = 1:size(rCond(rc).fList,1)
                    runRun{end+1,1} = {};
                    fieldList = fields(rCond(rc));
                    for fi = 1:length(fieldList)
                        if ismember(fieldList{fi},{'fList' 'date' 'bhvr'})
                            if ~isempty(rCond(rc).(fieldList{fi}))
                                runRun{end}.(fieldList{fi}) = rCond(rc).(fieldList{fi})(rr,1);
                            else
                                runRun{end}.(fieldList{fi}) = [];
                            end
                        elseif ismember(fieldList{fi},{'sub' 'ses' 'task' 'acq' 'dsgn'})
                            runRun{end}.(fieldList{fi}) = rCond(rc).(fieldList{fi});
                        % else
                        %     dbstack; error('extra field we don''t know what to do with...')
                        end
                    end
                end
            end
            runRun = [runRun{:}]';
            mriFiles = [runRun.fList]';
        elseif isstruct([mriFiles{:}])
            rCond = [mriFiles{:}]';

            if isfield(rCond,'labelAcq')
                rCond = rCond(ismember({rCond.labelAcq},funcFileAcqLabel));
            else
                if any(~contains({rCond.acq},funcFileAcqLabel)); warning('runCond may contain non functional files'); end
                % runCond = runCond(ismember({runCond.acq},funcFileAcqLabel));
            end
            
            runRun = {};
            for rc = 1:length(rCond)
                for rr = 1:size(rCond(rc).fList,1)
                    runRun{end+1,1} = {};
                    fieldList = fields(rCond(rc));
                    for fi = 1:length(fieldList)
                        if ismember(fieldList{fi},{'fList' 'date' 'bhvr'})
                            if ~isempty(rCond(rc).(fieldList{fi}))
                                runRun{end}.(fieldList{fi}) = rCond(rc).(fieldList{fi})(rr,1);
                            else
                                runRun{end}.(fieldList{fi}) = [];
                            end
                        elseif ismember(fieldList{fi},{'sub' 'ses' 'label' 'labelAcq' 'stim' 'acq' 'dsgn'})
                            runRun{end}.(fieldList{fi}) = rCond(rc).(fieldList{fi});
                        else
                            dbstack; error('extra field we don''t know what to do with...')
                        end
                    end
                end
            end
            runRun = [runRun{:}]';
            mriFiles = [runRun.fList]';
        else
            dbstack; error('code that');
        end
    end
elseif ischar(mriFiles)
else
    dbstack; error('code that');
end
clear rCond
% sort mri according to acquisition time
mriTimes = getAcqTime(mriFiles);
[mriTimes,b] = sort(mriTimes);
mriFiles = mriFiles(b);
if exist('runRun','var')
    runRun = runRun(b);
end
info.mriDate = cat(1,runRun.date);


% get sync param
forceThis = force;
info.maxLagSec = 20;
[info,phys] = doIt(physFileDb,info,mriFiles,forceThis,powerThrough);
if exist('runRun','var')
    % for r = 1:length(runRun)
    %     runRun(r).runPhysioTimes = phys.mriRunTimes(r,:);
    % end
    phys.mriRuns = runRun;
end
phys.fOrig = physFileDb;

f2 = info.f1;
title(findobj(f2.Children,'type','axes'),['physio to mri alignment' newline ' alignment NOT DONE (mri runs should roughly align with trigger trace)' newline 'estimated shift: ' num2str(info.FsShift,'%e') 'sec' newline 'estimated scaling: 1+' num2str(info.FsScale-1,'%e')])
drawnow

% apply sync (scale)
% physOrig = phys;
phys.samplerateOrig = phys.samplerate;
phys.samplerate = phys.samplerate .* phys.FsScale;
% blocktimes = zeros(size(phys.blocktimes));
% % blocktimesDate = repmat(datetime(0,1,1,0,0,0),size(blocktimes));
% blocktimesDate = cell(size(blocktimes));

% if numel(phys.blocktimes)>1; dbstack; error('double-check that'); end;
% phys.blocktimes = datetime(phys.blocktimes,'ConvertFrom','datenum');
phys.blocktimes(:,info.segmentInd) = phys.blocktimes(:,info.segmentInd) + duration(0,0,phys.FsShift);

% for blk = 1:size(phys.blocktimes,2)
%     blocktimesDate{1,blk} = datetime( phys.blocktimes(1,blk) + phys.FsShift ,'ConvertFrom','datenum');
%     % % anonymize date
%     % blocktimes(1,blk) = sum([hour(phys.blocktimes(1,blk))*60*60 minute(phys.blocktimes(1,blk))*60 second(phys.blocktimes(1,blk))]); % put in seconds
%     % H  = floor(blocktimes(1,blk)/60/60); blocktimes(1,blk) = blocktimes(1,blk) - H *60*60;
%     % MI = floor(blocktimes(1,blk)/60);    blocktimes(1,blk) = blocktimes(1,blk) - MI*60;
%     % S  = blocktimes(1,blk);
%     % phys.blocktimesOrig(1,blk) = datetime(0,1,1,H,MI,S);
%     % % apply sync (shift)
%     % blocktimes(1,blk) = sum([H*60*60 MI*60 S]); % put in seconds
%     % blocktimes(1,blk) = blocktimes(1,blk) + phys.FsShift;
%     % H  = floor(blocktimes(1,blk)/60/60); blocktimes(1,blk) = blocktimes(1,blk) - H *60*60;
%     % MI = floor(blocktimes(1,blk)/60);    blocktimes(1,blk) = blocktimes(1,blk) - MI*60;
%     % S  = blocktimes(1,blk);
%     % blocktimesDate(1,blk) = datetime(0,1,1,H,MI,S);
% end
% phys.blocktimes = cat(2,blocktimesDate{:});
phys.mriSyncApplied = 1;

%plot corrected aligment
forceThis = force;
[infoX,~] = doIt(physFileDb,info,mriFiles,forceThis,powerThrough);
f3 = gcf;
title(['physio to mri alignment' newline ' alignment DONE (mri runs should perfectly align with trigger trace)' newline 'estimated residual shift: ' num2str(infoX.FsShift,'%e') 'sec' newline 'estimated residual scaling: 1+' num2str(infoX.FsScale-1,'%e')])
drawnow


%% Mark run start and end times
dt = sum([hour(phys.blocktimes(:,info.segmentInd))*60*60 minute(phys.blocktimes(:,info.segmentInd))*60 second(phys.blocktimes(:,info.segmentInd))]);
runTimes = phys.mriRunTimes - dt;


% %% Plot traces
% chanList = phys.titles1(ismember(phys.titles1,{'cardiac' 'resp'}));
% for c = 1:length(chanList)
%     curChan = chanList{c};
%     chanInd = ismember(phys.titles1,curChan);
%     if ~any(chanInd); dbstack; error('X'); end
%     curData = phys.data(phys.datastart(chanInd,info.segmentInd):phys.dataend(chanInd,info.segmentInd));
%     Fs = phys.samplerate(chanInd,info.segmentInd);
%     n = length(curData);
%     t = linspace(0,(n-1)./Fs,n);
%     % t = 0:1/Fs:(n-1)./Fs;
%     for r = 1:size(runTimes,1)
%         ind = ( runTimes(r,1)<=t & t<runTimes(r,2) );
% 
%         f4(r,c) = figure('WindowStyle','docked');
%         ht = tiledlayout(2,1); ht.Padding = 'tight'; ht.TileSpacing = 'tight';
%         ax1(r) = nexttile;
%         plot(t(ind),curData(ind));
%         axis tight
%         xlabel('time (sec)')
%         drawnow
% 
%         ax2(r) = nexttile;
%         % df = info.W*2/5;
%         % nfft = Fs/2 / df;
%         nfft = length(curData(ind));
%         nfft = 2^(nextpow2(nfft)-1);
%         W = info.W*2;
%         T = length(curData(ind))./Fs;
%         TW = T.*W;
%         K = round(TW*2-1);
%         TW = (K+1)/2;
%         W = TW/T;
%         disp(['computing ' curChan ' mt spectum with ' num2str(K) ' tapers; run ' num2str(r) '/' num2str(size(runTimes,1))])
%         tic
%         [pw,f] = pmtm(curData(ind)-mean(curData(ind)),TW,nfft,Fs);
%         toc
%         plot(f,pw);
%         axis tight
%         xlim([0 10])
%         ax2(end).YScale = 'log';
%         grid on
%         yLim = ylim(ax2(end)); xLim = xlim(ax2(end));
%         hLine = line(ax2(end),1+W.*[-1 1],exp(mean(log(yLim))).*[1 1]);
%         hLine.Color = 'g';
%         hLine.LineWidth = 5;
%         uistack(hLine,'bottom')
% 
%         xlabel('Hz')
%         ylabel('PSD')
% 
%         title(ht,[curChan '; run' num2str(r)])
%         drawnow
%     end
% end


%% Plot traces
% datetime(runTimes,'ConvertFrom','datenum')
% chanList = phys.titles1(ismember(phys.titles1,{'cardiac' 'resp'}));
chanList = {'cardiac' 'resp'};
for r = 1:size(runTimes,1)
    f4(r) = figure('WindowStyle','docked');
    ht = tiledlayout(2,2); ht.Padding = 'tight'; ht.TileSpacing = 'tight';
    for c = 1:length(chanList)
        chanInd = ismember(phys.titles1,chanList{c});
        if ~any(chanInd); continue; end

        curData = phys.data(phys.datastart(chanInd,info.segmentInd):phys.dataend(chanInd,info.segmentInd));
        Fs = phys.samplerate(chanInd,info.segmentInd);
        n = length(curData);
        t = linspace(0,(n-1)./Fs,n);
        % t = 0:1/Fs:(n-1)./Fs;
        ind = ( runTimes(r,1)<=t & t<runTimes(r,2) );

        % f4(r,c) = figure('WindowStyle','docked');
        % ht = tiledlayout(2,1); ht.Padding = 'tight'; ht.TileSpacing = 'tight';
        ax1(r,c) = nexttile;
        plot(t(ind),curData(ind));
        axis tight
        grid on
        grid minor
        xlabel('time (sec)')
        ylabel('signal amplitude (a.u.)')
        title([chanList{c} ' timecourse'])
        drawnow

        ax2(r,c) = nexttile;
        % df = info.W*2/5;
        % nfft = Fs/2 / df;
        nfft = length(curData(ind));
        nfft = 2^(nextpow2(nfft)-1);
        W = info.W*2;
        T = length(curData(ind))./Fs;
        TW = T.*W;
        K = round(TW*2-1);
        TW = (K+1)/2;
        W = TW/T;
        disp(['computing ' chanList{c} ' mt spectum with ' num2str(K) ' tapers; run ' num2str(r) '/' num2str(size(runTimes,1))])
        tic
        [pw,f] = pmtm(curData(ind)-mean(curData(ind)),TW,nfft,Fs);
        toc
        plot(f,pw);
        axis tight
        xlim([0 10])
        ax2(r,c).YScale = 'log';
        grid on
        yLim = ylim(ax2(r,c)); % xLim = xlim(ax2(r,c));
        hLine = line(ax2(r,c),1+W.*[-1 1],exp(mean(log(yLim))).*[1 1]);
        hLine.Color = 'g';
        hLine.LineWidth = 5;
        uistack(hLine,'bottom')

        xlabel('Hz')
        ylabel('PSD')
        title([chanList{c} ' spectrum'])

        title(ht,['run' num2str(r)])
        drawnow

    end
end
set(ax1,'XMinorTick','on','XMinorGrid','on')
set(ax2,'XMinorTick','on','XMinorGrid','on')
for r = 1:size(runTimes,1)
    linkaxes(ax1(r,:),'x')
    linkaxes(ax2(r,:),'x')
end


%% save
disp(['saving minimally curated physio data' newline info.outFile])
phys.info = info;
phys.info.f1.Visible = 'off';
phys.info.f2.Visible = 'off';
save(info.outFile,'phys')
fileName1 = replace(info.outFile,'.mat','_chanId.png');
fileName2 = replace(info.outFile,'.mat','_phys2mri.png');
disp(['saving physio data curation figure' newline fileName1 newline fileName2])
copyfile(info.chanManId,fileName1)
saveas(f2,fileName2)

for r = 1:length(f4)
    fileName4 = replace(info.outFile,'.mat',[ '_run' num2str(r) '.fig']);
    disp(['saving physio data curation figure (' num2str(r) '/' num2str(length(f4)) ')' newline fileName4])
    saveas(f4(r),fileName4)
end




function [info,physio] = doIt(physioFile,info,mriFileList,figFlag,powerThrough)
% Automatically extract cardiac, respiration and scanner triggers from
% LabChart data in matlab format, sorting out runs at the same time. The
% data type corresponding to each channel is automatically identified based
% on power spectra, and the appropriate data segment is selected as the
% longest, assuming all mri runs were acquired during a single long
% LabChart data segment. This latter behavior is overridden by specifying
% the segmentInd variable as input.

% info.chanLabel must be a cell array of strings, each string being one of
% {'trigger' 'cardiac' 'resp'}.


% By Sebastien Proulx (jsproulx@mgh.harvard.edu)
% 2023-08-21
if ~exist('powerThrough','var'); powerThrough = []; end
if isempty(powerThrough); powerThrough = 0; end


if ~exist('info','var');     info.chanLabel = [];
    info.dbFile     = [];
    info.outFile    = [];
    info.segmentInd = [];
    info.runTimes = [];
    info.chanKnown = [];
    info.W = []; end
if ~exist('figFlag','var');         figFlag = []; end
% if ~exist('segmentInd','var');   segmentInd = []; end
% if ~exist('chanLabel','var');     chanLabel = {}; end
if ~exist('mriFileList','var'); mriFileList = {}; end
if isempty(figFlag)
    if ~isempty(info.segmentInd); figFlag = 1;
    else;                         figFlag = 2; end
end
if nargin==0; info.chanLabel = {'trigger' 'cardiac' 'resp'}; return; end
if isempty(info.chanLabel) || figFlag; doFFT = true; else doFFT = false; end
if isempty(info.chanKnown)
    info.chanKnown = ~(  isempty(info.chanLabel)  ||  ~any(ismember(info.chanLabel,{'trigger' 'cardiac' 'resp'}))  );
end
if isempty(info.W); info.W = 0.02; end

%% Load physio
physio = load(physioFile);
physio.blocktimes = datetime(physio.blocktimes,'ConvertFrom','datenum');
if isempty(info.chanLabel); info.chanLabel = cellstr(num2str((1:size(physio.datastart,1))','chan%i'))'; end

%% Load MRI info
if ~isempty(mriFileList)
    mri = MRIread(mriFileList{1},1);
    volTr = mri.tr/1000;
    nframes = mri.nframes;
    runDur = volTr*nframes;
    rfTr = ['jq ''.RepetitionTimeExcitation'' ' replace(mriFileList{1},'.nii.gz','.json')]; [~,rfTr] = system(rfTr); rfTr = str2num(replace(replace(rfTr,'[0;39m',''),['[0m' newline],''));
    if isempty(rfTr); rfTr = volTr; end
end

%% Identify segment
if isempty(info.segmentInd)
    if figFlag>1; figure('WindowStyle','docked'); end
    startTime = nan(size(physio.datastart,2),size(physio.datastart,1));
    endTime = nan(size(physio.datastart,2),size(physio.datastart,1));
    for chanInd = 1:size(physio.datastart,1)
        for segmentIndTmp = 1:size(physio.datastart,2)
            data = physio.data(physio.datastart(chanInd,segmentIndTmp):physio.dataend(chanInd,segmentIndTmp));
            Fs = physio.samplerate(chanInd,segmentIndTmp);
            FsTarg = 50;
            ds = round(Fs/FsTarg);
            FsDs = Fs/ds;
            if segmentIndTmp==1
                tPhys = 0:1/Fs:length(data)/Fs-1/Fs;
                if figFlag>1; hPlotTmp(chanInd) = plot(tPhys,data); hold on; end
                % if figFlag>1; hPlotTmp(chanInd) = plot(downsample(tPhys,ds),interp1(tPhys,data,downsample(tPhys,ds))); hold on; end
            else
                tPhys = tLast + 1/Fs + (0:1/Fs:length(data)/Fs-1/Fs);
                if figFlag>1; plot(tPhys,data,'Color',hPlotTmp(chanInd).Color); end
            end
            tLast = tPhys(end);
            startTime(segmentIndTmp,chanInd) = tPhys(1);
            endTime(segmentIndTmp,chanInd) = tPhys(end);
        end
    end
    x = [startTime(:,1); endTime(end,1)];
    if figFlag>1
        plot(repmat(x',[2 1]),repmat(ylim',[1 length(x)]),'k','LineWidth',2);
        yLim = ylim;
        text(startTime(:,1),repmat(yLim(1),[length(startTime(:,1)) 1]),replace(fullfile('Segment',cellstr(num2str((1:size(physio.datastart,2))'))),filesep,' '),'Rotation',90,'VerticalAlignment','top','FontSize',20)
        xlabel('time (sec)')
        title('Select the appropriate segment and runs')
    end
    %%% If unspecified, choose the longest segment
    if isempty(info.segmentInd)
        [~,info.segmentInd] = max(endTime(:,1) - startTime(:,1));
    end
    if figFlag>1
        x = [startTime(info.segmentInd,1) endTime(info.segmentInd,1)];
        x = x([1 2 2 1 1]);
        y = ylim;
        y = y([1 1 2 2 1]);
        hPatch = patch(x,y,'k');
        uistack(hPatch,'bottom')
        set(hPatch,'EdgeColor','none')
        % set(hPatch,'FaceColor',[1 1 1].*0.9,'EdgeColor','none')
        hPatch.FaceAlpha = 0.1;

        hLeg = legend([hPlotTmp hPatch],[replace(fullfile('Channel',cellstr(num2str((1:size(physio.datastart,1))'))),filesep,' '); {'selected segment'}],'AutoUpdate','off');
    end
    dataTime = [];
    sampleRate = [];
    volTr = [];
    nframes = [];
    dataTrigTimes = [];
    dataChan = [];
    dataRun = [];
    info.chanLabel = {'trigger' 'cardiac' 'resp'};

    yLim = ylim;
    if length(unique(physio.com(:,5)))~=length(physio.com(:,5))
        tmp = unique(physio.com(:,5));
        ind = zeros(size(tmp));
        for i = 1:length(tmp)
            ind(i) = find(tmp(i)==physio.com(:,5),1,'last');
        end
        physio.com = physio.com(ind,:);
    end
    hText = text(tPhys(physio.com(:,3)),ones(size(physio.com(:,3))).*yLim(1),physio.comtext,'Rotation',90);


    % manually identify rough run limits
    xlim([startTime(info.segmentInd,1) endTime(info.segmentInd,1)])
    info.runTimes = {};
    drawnow
    commandwindow();
    disp('zoom as desired, then press enter')
    input('')
    % input('zoom as desired, then press enter')


    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    disp('for channel identification and initial timing estimate')
    disp('click the begining and end of a single functional run')
    disp('when done press enter')
    commandwindow();
    [info.runTimes{end+1},~,button] = ginput;
    
    info.runTimes{end} = info.runTimes{end}(button==1);
    info.runTimes{end} = info.runTimes{end}([end-1 end]);
    x = info.runTimes{end}([1 2 2 1 1]);
    hP = patch(x,y,'k','FaceAlpha',0.2,'EdgeColor','non');
    hLeg.String{end} = ['run' num2str(length(info.runTimes))];
    uistack(hP,'bottom')
    % uistack(hP,'up')
    % uistack(hP,'up')
    drawnow
    
    if ~info.chanKnown
        info.chanLabel = cellstr(num2str((1:size(physio.datastart,1))','chan%i'))';
    end


    return
end

%% Identify channels
if isempty(info.segmentInd)
    info.segmentInd = size(physio.datastart,2);
end
if figFlag>0
    figure('WindowStyle','docked');
    hTile = tiledlayout(2,1);
    hTile.Padding = 'tight';
    hTile.TileSpacing = 'tight';
    ax1 = nexttile;
    ax2 = nexttile;
end
n = physio.dataend(:,info.segmentInd) - physio.datastart(:,info.segmentInd) + 1;
if any(diff(n)); stack = dbstack; error([stack(end).file ' line ' num2str(stack(end).line)]); end
% physioSpec = nan(size(physio.datastart,1),n(1));
physioSpec = cell(size(physio.datastart));
f = cell(size(physio.datastart));
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')
if figFlag>0; clear hPlot1 hPlot2; end
if figFlag>0; axes(ax1); end
for chanInd = 1:size(physio.datastart,1)
    % if figFlag>0; axes(ax1); end
    Fs = physio.samplerate(chanInd,info.segmentInd);
    curData = physio.data(physio.datastart(chanInd,info.segmentInd):physio.dataend(chanInd,info.segmentInd));
    tPhys = (0:length(curData)-1)/Fs;
    % tPhys = 0:1/Fs:length(curData)/Fs-1/Fs;
    % if ~info.chanKnown
    ind = tPhys>=info.runTimes{1}(1) & tPhys<=info.runTimes{1}(2);
    tPhys = tPhys(ind);
    curData = curData(ind);
    % end


    % FsTarg = 50;
    % ds = round(n(1)/2^nextpow2(n(1)./(Fs/FsTarg)));
    % ds = 1;
    % FsDs = Fs/ds;
    %
    % t = 0:1/Fs:length(curData)/Fs-1/Fs;
    % curData = curData - polyval(polyfit(t',curData',4),t')';
    % curData = curData - mean(curData,2);
    % curData = curData./max(curData);
    if figFlag>0
        % axes(ax1);
        % hPlot1(chanInd) = plot(downsample(t,ds),downsample(curData,ds)); hold on
        hPlot1(chanInd) = plot(tPhys,curData); hold on
        uistack(hPlot1(chanInd),'bottom')
        axis tight
    end
    % % Fs = Fs/length(curData);
    % % f = linspace(0,Fs,length(curData));

    if doFFT

        % df = Fs/2 / size(curData,2);
        df = info.W/5;
        nfft = Fs/2 / df;
        nfft = 2^(nextpow2(nfft));
        % nfft = 2^(nextpow2(size(curData,2))-1);
        % if info.chanKnown
        %     K = 10;
        %     TW = (K+1)/2;
        %     T = length(curData)./Fs;
        %     W = TW/T;
        % else
        W = info.W;
        T = length(curData)./Fs;
        TW = T.*W;
        K = round(TW*2-1);
        TW = (K+1)/2;
        W = TW/T;
        % end



        disp(['computing mt spectum with ' num2str(K) ' tapers for channel identification ' num2str(chanInd) '/' num2str(size(physio.datastart,1))])
        tic
        % [physioSpec{chanInd},f{chanInd}] = pmtm(downsample(curData,ds),TW,nfft,FsDs);
        [physioSpec{chanInd},f{chanInd}] = pmtm(curData-mean(curData),TW,nfft,Fs);
        toc
    end
end

if doFFT
    physioSpec = cat(2,physioSpec{:})';
    f = cat(2,f{:})';
end

% physioSpec(chanInd,:) = abs(fft(curData));
if figFlag>0 && doFFT
    hPlot2 = plot(ax2,f',physioSpec'); hold on
    ax2.YScale = 'log';
    for chanInd = 1:size(physio.datastart,1)
        uistack(hPlot2(chanInd),'bottom');
    end
end


if figFlag>0
    drawnow
    % axes(ax1)
    yLim = ylim;
    % hText = text(ax1,t(physio.com(:,3)),ones(size(physio.com(:,3))).*yLim(1),physio.comtext,'Rotation',90);
    xlabel(ax1,'time(sec)')
    % axes(ax2)
    if doFFT
        xlabel(ax2,'Hz')
    end
    title('run 1')
    legend(hPlot1,info.chanLabel,'AutoUpdate','off');
    xlim(ax2,[0 10]);
    grid(ax2,'on')
    grid(ax1,'on')

    yLim = ylim(ax2); xLim = xlim(ax2);
    hLine = line(ax2,1+W.*[-1 1],exp(mean(log(yLim))).*[1 1]);
    hLine.Color = 'g';
    hLine.LineWidth = 5;
    uistack(hLine,'bottom')

    drawnow
end

if ~info.chanKnown
    info.chanLabel = repmat({''},size(physio.datastart))';
    info.chanLabel(1:3) = {'trigger' 'cardiac' 'resp'};
    return
end

% %%% specify frequency ranges
% if figFlag>0; yLim = ylim(ax2); end; hold(ax2,'on');
% trFreq = [-0.5 0.5]+1/rfTr; if figFlag>0; hPlotFreq = plot(ax2,trFreq,[1 1].*mean(yLim),'k'); end
% trPower = max(physioSpec(:,f>trFreq(1) & f<trFreq(2)),[],2);
% if figFlag>0; xlim(ax2,[0 trFreq(2)]); end
% runFreq = [0 1/runDur]; if figFlag>0; plot(ax2,runFreq,[1 1].*mean(ylim),'k'); end
% runPower = max(physioSpec(:,f>runFreq(1) & f<runFreq(2)),[],2);
% cardFreq = [0.75 1.25]; if figFlag>0; plot(ax2,cardFreq,[1 1].*mean(yLim),'k'); end
% % cardPower = max(physioSpec(:,f>cardFreq(1) & f<cardFreq(2)),[],2);
% cardPower = max(physioSpec(:,f>cardFreq(1) & f<cardFreq(2)),[],2)./sum(physioSpec,2);
% respFreq = [0.1 0.5]; if figFlag>0; plot(ax2,respFreq,[1 1].*mean(yLim),'k'); end
% % respPower = max(physioSpec(:,f>respFreq(1) & f<respFreq(2)),[],2);
% respPower = max(physioSpec(:,f>respFreq(1) & f<respFreq(2)),[],2)./sum(physioSpec,2);
% %%% identify based on maximum power at corresponding frequencies
% %%%% trigger
% [~,b1] = max(trPower); [~,b2] = max(runPower);
% trigChanInd = b1;
% if b1~=b2
%     warning(['The trigger channel does not show highest power at the lowest frequency.' newline 'Respiration trace may be bad.'])
% end
% %%%% cardiac
% [~,cardChanInd] = max(cardPower);
% % %%%% respiration
% % [~,respChanInd] = max(respPower);
% % %%% set labels
% % if ~all(diff([trigChanInd cardChanInd respChanInd]))
% %     stack = dbstack; error([stack(end).file ' line ' num2str(stack(end).line)]);
% % end
% % chanLabel = {'trigger' 'cardiac' 'resp'}';
% % physio.titles1 = cell(size(chanLabel));
% % physio.titles1([trigChanInd cardChanInd respChanInd]) = chanLabel;

if size(physio.datastart,1) ~= length(info.chanLabel); dbstack; error(['info.chanLabel must be length ' num2str(size(physio.datastart,1))]); end
% info.chanLabel = {'trigger' 'cardiac'}';
physio.titles1 = info.chanLabel';
% physio.titles1 = cell(size(info.chanLabel));
% physio.titles1([trigChanInd cardChanInd]) = info.chanLabel;


%% Adjust timing to MRI
chanInd = ismember(physio.titles1,'trigger');
Fs = physio.samplerate(chanInd,info.segmentInd);
tPhysStart = physio.blocktimes(info.segmentInd);
if isfield(info,'FsScale') && ~isempty(info.FsScale)
    Fs = Fs*info.FsScale;
end
if isfield(info,'FsShift') && ~isempty(info.FsShift)
    tPhysStart = tPhysStart + duration(0,0,info.FsShift);
end
tPhysStartSecs = seconds(duration(hour(tPhysStart),minute(tPhysStart),second(tPhysStart)));
physDat = physio.data(physio.datastart(chanInd,info.segmentInd):physio.dataend(chanInd,info.segmentInd));
nPhys = (physio.dataend(1,info.segmentInd) - physio.datastart(1,info.segmentInd)) + 1;
tPhys = (0:nPhys-1)/Fs;

% tPhys = linspace(0,(nPhys-1)/Fs,nPhys);
% iSegStart = physio.datastart(1,info.segmentInd);
% iSegEnd   = physio.dataend(1,info.segmentInd);
% tPhysStart = tPhysStart + (iSegStart-1)/Fs; % time [s] of physio recording segment start (0=midnigh)
% t = (iSegStart-1:iSegEnd-1)/Fs + tPhysStart; % time [s] of physio recording segment (0=midnigh)
% t = 0:1/Fs:(iSegEnd-iSegStart)/Fs; % time [s] (0=midnigh)



% t = tPhysSegment:1/Fs:(length(physDat)-1+)/Fs; % 0 = time of physio session recording start
% get mri start and end times
mriTimes = info.mriDate + getAcqTime(mriFileList); % time of mri recording start (0=midnigh)
% [mriTimes,b] = sort(mriTimes);
% mriFileList = mriFileList(b);
tMRI    = zeros(size(mriTimes,1),2);
tr      = zeros(size(mriTimes,1),1);
nframes = zeros(size(mriTimes,1),1);
for i = 1:length(mriTimes)
    disp([num2str(i) '/' num2str(length(mriTimes))])
    tMRI(i,1) = sum([hour(mriTimes(i))*60*60 minute(mriTimes(i))*60 second(mriTimes(i))]); % time of mri [s] run recording start (0=midnight) !!! not the same as t
    if ~isfield(info,'mri') || isempty(info.mri) || length(info.mri)<i
        info.mri(i,1) = MRIread(mriFileList{i},1);
    end
    tr(i)      = info.mri(i,1).tr/1000;
    nframes(i) = info.mri(i,1).nframes;
    tMRI(i,2) = tMRI(i,1) + tr(i) * nframes(i);
end



info.f1 = figure('WindowStyle','docked');
plot(tPhys,physDat); hold on
% xlabel(['seconds since physio segment; physio recording started at ' char(datetime(physio.blocktimes(info.segmentInd),'ConvertFrom','datenum','Format','HH:mm:ss.SSS'))]);
% tmpTime = tPhysStart;
% H  = floor(tmpTime/60/60); tmpTime = tmpTime - H *60*60;
% MI = floor(tmpTime/60);    tmpTime = tmpTime - MI*60;
% S  = tmpTime;
xlabel([...
    'seconds since physio segment' newline ...
    'physio recording started at ' char(datetime(physio.blocktimes(1),'Format','HH:mm:ss.SSS')) newline ...
    'physio recording segment started at ' char(datetime(tPhysStart,'Format','HH:mm:ss.SSS'))]);
% xlabel([...
%     'seconds since physio segment' newline ...
%     'physio recording started at ' char(datetime(physio.blocktimes(1),'ConvertFrom','datenum','Format','HH:mm:ss.SSS')) newline ...
%     'physio recording segment started at ' char(datetime(tPhysStart,'Format','HH:mm:ss.SSS'))]);
xx = tMRI(:,[1 2 2 1 1])' -  tPhysStartSecs;
% xx = tMRI(:,[1 2 2 1 1])' - tPhysStart;
yy = ylim; yy = repmat(yy([1 1 2 2 1])',[1 length(mriTimes)]);
hP = patch(xx,yy,'k');
uistack(hP,'bottom')
legend({'mri runs' 'trigger trace'},'AutoUpdate','off')



% cross correlate mri and physio to get a rough time shift due to non-synchronous
% MRI and physio
mriDat = zeros(size(physDat));
mriDatMax = max(physDat);
for i = 1:length(mriTimes)
    mriDat(tPhys>=tMRI(i,1)-tPhysStartSecs & tPhys<tMRI(i,2)-tPhysStartSecs) = mriDatMax;
end
% plot(tPhys,mriDat,'r')
[r,mriLag] = xcorr(physDat,mriDat,round(info.maxLagSec*Fs));
[~,b] = max(r);
mriLagSec = (mriLag(b)-1)/Fs;
info.f2 = figure('WindowStyle','docked');
plot(mriLag/Fs,r);
hold on
h = xline(mriLagSec,'r');
xlabel(['delay between mri and physio times (sec)' newline 'physio leads for positive values'])
ylabel('cross-correlation between trigger recording and MRI square wave function')
grid on
legend(h,'estimated time shift (cross-correlation max)','AutoUpdate','off')

% mriLagSec = (mriLag(b)-1) / Fs; % not sure about that minus one

% scale physio Fs to MRI tr
if volTr==rfTr
    dtTrig = volTr;
else
    dtTrig = rfTr;
end

figure(info.f1);
x = tPhys(  tPhys>tMRI(1  ,1)-tPhysStartSecs+mriLagSec-1 & tPhys<tMRI(1  ,1)-tPhysStartSecs+mriLagSec+1);
y = physDat(tPhys>tMRI(1  ,1)-tPhysStartSecs+mriLagSec-1 & tPhys<tMRI(1  ,1)-tPhysStartSecs+mriLagSec+1);
tFirst = x(find(y>max(abs(diff(y)))/2,1,'first')-1);
% plot(x,y,'r')
xline(tFirst,'r')
xlim([-5 5] + tFirst)
if ~powerThrough
    disp('is run onset found correctly? dbcont')
    keyboard
end

x = tPhys(  tPhys>tMRI(end,end)-tPhysStartSecs+mriLagSec-0.5 & tPhys<tMRI(end,end)-tPhysStartSecs+mriLagSec+0.5);
y = physDat(tPhys>tMRI(end,end)-tPhysStartSecs+mriLagSec-0.5 & tPhys<tMRI(end,end)-tPhysStartSecs+mriLagSec+0.5);
tLast = x(find(y>max(abs(diff(y)))/2,1,'last')) + dtTrig;
% plot(x,y,'r')
xline(tLast,'r')
xlim([-5 5] + tLast)
if ~powerThrough
    disp('is run offset found correctly? dbcont')
    keyboard
end
xlim auto


physDur = tLast - tFirst;
mriDur = tMRI(end,2)-tMRI(1,1);
% if ~isfield(info,'FsScale') || isempty(info.FsScale)
FsScale = physDur/mriDur;
% end


% Redefing lag based on the more precise trigger time estimations
mriLagSec = mean([tFirst - (tMRI(1  ,1)-tPhysStartSecs)
                  tLast  - (tMRI(end,2)-tPhysStartSecs)]);


% get time shift due to above scaling (shift at physio time=0)
% tFirst*FsScale
% FsShift = tFirst - tFirst/FsScale;
T = mean([tFirst tLast]);
FsShift = T/FsScale - T;
% FsShift = 0 + extraShift;
% % pBefore = mean([tFirst tLast]) * Fs;
% % pAfter  = mean([tFirst tLast]) * (Fs*FsScale);
% % FsShift = (pAfter - pBefore) / Fs;

% Combine both time shift into the physio time
disp(['FsShift=' num2str(-FsShift) 'sec'])
% FsScale = 1;
% FsShift = 0;
info.FsScale = FsScale;
info.FsShift = -FsShift - mriLagSec;

disp(['  scale=1+' num2str(physDur/mriDur-1)])
disp(['    lag=' num2str(-mriLagSec) 'sec'])
disp([' mriLag=' num2str(-mriLagSec) 'sec'])

physio.FsShift   = info.FsShift;
physio.FsScale   = info.FsScale;
info.physioTime = 'tPhys = sum([hour(physio.blocktimes(info.segmentInd))*60*60 minute(physio.blocktimes(info.segmentInd))*60 second(physio.blocktimes(info.segmentInd))]);';
info.physio2mri = ['shift physio time axis: tPhys + physio.FsShift' newline ...
    'scale physio time axis: physio.samplerate * physio.FsScale'];
physio.physioTime = info.physioTime;
physio.physio2mri = info.physio2mri;

% % Plot data
% 
% % datestr(physio.blocktimes(info.segmentInd))
% figure('WindowStyle','docked');
% plot(tPhysSegment,physDat); hold on
% xlabel(['seconds since midnight; physio recording started at ' char(datetime(physio.blocktimes(info.segmentInd),'ConvertFrom','datenum','Format','HH:mm:ss.SSS'))]);
% xx = tMRI(:,[1 2 2 1 1])';
% % xx = tMRI(:,[1 2 2 1 1])' - tPhysStart;
% yy = ylim; yy = repmat(yy([1 1 2 2 1])',[1 length(mriTimes)]);
% hP = patch(xx,yy,'k');
% uistack(hP,'bottom')
% legend({'mri runs' 'trigger trace'})
%
%
%
% Fs = physio.samplerate(chanInd,info.segmentInd);
% Fs = Fs*info.FsScale;
% physDat = physio.data(physio.datastart(chanInd,info.segmentInd):physio.dataend(chanInd,info.segmentInd));
% t = 0:1/Fs:length(physDat)/Fs-1/Fs;
% figure('WindowStyle','docked');
% plot(t,physDat); hold on
% xlabel(['seconds since physio recording starts (' erase(datestr(physio.blocktimes(info.segmentInd)),' ') ')']);
%
% tMRI    = zeros(size(mriTimes,1),2);
% for i = 1:length(mriTimes)
%     tMRI(i,1) = sum([hour(mriTimes(i))*60*60 minute(mriTimes(i))*60 second(mriTimes(i))]);
%     tMRI(i,2) = tMRI(i,1) + tr(i) * nframes(i);
% end
%
% % apply shift
% tPhys = tPhys + info.FsShift;
%
% xx = tMRI(:,[1 2 2 1 1])' - tPhys;
% yy = ylim; yy = repmat(yy([1 1 2 2 1])',[1 length(mriTimes)]);
% hP = patch(xx,yy,'k');
% uistack(hP,'bottom')
% legend({'mri runs' 'trigger trace'})
%
%
%
% disp(['scale=1+' num2str(info.FsScale-1)])
% disp(['  lag=' num2str(info.FsShift)])


%% Output
physio.mriRunFiles = mriFileList;
physio.mriRunTimes = tMRI;
physio.mriRunDates = [mriTimes mriTimes + duration(0,0,diff(tMRI,[],2))];


return

%%% Find trigger times
trigInd = ismember(physio.titles1,'trigger');
trigData = physio.data(physio.datastart(trigInd,info.segmentInd):physio.dataend(trigInd,info.segmentInd));
thresh = min(trigData)+range(trigData)/2;
trigTimes = diff(trigData>thresh); trigTimes(trigTimes<0) = 0; trigTimes = find(trigTimes)+1; trigTimes = tPhys(trigTimes);
%%% Find candidate run start times
runStart = [trigTimes(1) trigTimes(find(diff(trigTimes)>volTr*2)+1)];
%%% Keep runs with expected number of frames
runEnd = nan(size(runStart));
runFrames = nan(size(runStart));
for i = 1:length(runStart)
    indStart = find(runStart(i)==trigTimes);
    if indStart+nframes<=length(trigTimes) && any(diff(trigTimes(indStart:end))>volTr*1.5)
        indEnd = indStart + find(diff(trigTimes(indStart:end))>volTr*1.5,1) - 1;
    else
        indEnd = length(trigTimes);
    end
    runFrames(i) = indEnd - indStart + 1;
    runEnd(i) = trigTimes(indEnd);
end
runStart = runStart(runFrames==nframes);
runEnd = runEnd(runFrames==nframes);
if figFlag>0
    axes(ax1);
    yLim = ylim;
    for i = 1:length(runStart)
        hPatch(i) = patch([runStart(i) runEnd(i) runEnd(i) runStart(i) runStart(i)],yLim([1 1 2 2 1]),'k');
    end
    uistack(hPatch,'bottom')
    set(hPatch,'FaceColor',[1 1 1].*0.9,'EdgeColor','none')

    axes(ax1); legend([hPlot1 hPatch(1)],[physio.titles1; {'runs'}]','location','northwest');
    axes(ax2); legend([hPlot2 hPlotFreq],[physio.titles1; {'freq of interest'}],'location','northwest');
end
%%% Output data
data = cell(size(runStart,2),length(chanInd));
dataTime = cell(size(runStart,2),length(chanInd));
dataTrigTimes = cell(size(runStart,2),length(chanInd));
dataRun = nan(size(runStart,2),1);
dataChan = cell(1,length(chanInd));
timePad = 20; % in sec
for chanInd = 1:length(physio.titles1)
    Fs = physio.samplerate(chanInd,info.segmentInd);
    curData = physio.data(physio.datastart(chanInd,info.segmentInd):physio.dataend(chanInd,info.segmentInd));
    tPhys = 0:1/Fs:length(curData)*1/Fs-1/Fs;
    for runInd = 1:length(runStart)
        [~,indStart] = min(abs(tPhys-(runStart(runInd)-timePad)));
        [~,indEnd] = min(abs(tPhys-(runEnd(runInd)+timePad)));
        data{runInd,chanInd} = curData(indStart:indEnd);
        dataTime{runInd,chanInd} = tPhys(indStart:indEnd) - runStart(runInd);
        curTrigTimes = trigTimes - runStart(runInd); curTrigTimes = curTrigTimes(curTrigTimes>dataTime{runInd,chanInd}(1) & curTrigTimes<dataTime{runInd,chanInd}(end));
        dataTrigTimes{runInd,chanInd} = curTrigTimes;
        dataChan(chanInd) = physio.titles1(chanInd);
        dataRun(runInd) = runInd;
    end
end
sampleRate = physio.samplerate(:,info.segmentInd)';



function info = manId(info,physFileDb)

%% Identify runs on physio trace
info2 = doIt(physFileDb);
info.chanLabel  = info2.chanLabel;
info.segmentInd = info2.segmentInd;
info.runTimes   = info2.runTimes;
info.chanKnown  = 0;
info.W          = info2.W; clear info2

chanLabel = info.chanLabel;

f1 = gcf;


%% Identify channels
info = doIt(physFileDb,info);
done = 0;
while ~done
    i = zeros(1,3);
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    disp('identify each channel by number')
    for chanInd = 1:3
        disp(['input channel number for ' info.chanLabel{chanInd} ':'])
        tmp = input('')
        % tmp = input(['input channel number for ' info.chanLabel{chanInd} ':' newline]);
        if ~isempty(tmp)
            i(chanInd) = tmp;
        end
    end
    for chanInd = 1:length(i)
        if i(chanInd)
            chanLabel(i(chanInd)) = info.chanLabel(chanInd);
        end
    end
    hLeg = gcf; hLeg = findobj(hLeg.Children.Children,'type','Legend');
    hLeg.String = chanLabel;
    disp('satisfied?');
    done = input(['satisfied?' newline]);
    disp(newline)
end
% chanLabel = {'cardiac' 'trigger' 'chan3' 'chan4'};
info.chanLabel = chanLabel;
info.chanKnown = 1;

f2 = gcf;
title(f2.Children,'channel identification')
drawnow

info.runManId = fullfile(fileparts(physFileDb),'manual','runId.png');
saveas(f1,info.runManId);
info.chanManId = fullfile(fileparts(physFileDb),'manual','chanId.png');
saveas(f2,info.chanManId)




