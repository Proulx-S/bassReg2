function xSet = bidsCount(xSet)

%% Number of echo
xSet.nEcho = length(unique(xSet.files_bidsList(contains(xSet.files_bidsList(:,1),'echo-'),:)));
if ~xSet.nEcho; xSet.nEcho = 1; end

%% Number of runs
xSet.nRun  = length(unique(xSet.files_bidsList(contains(xSet.files_bidsList(:,1),'run-'),:)));
if ~xSet.nRun;  xSet.nRun  = 1; end

%% Number of venc
% Special case:
% nVenc    is the number of venc acquired, excluding the null venc
% nVencDat is the number of data file in the structure
% The reason for that is that data is currently stored
% mag-venc0, mag-venc30m0, phase-venc30m0, ..., mag-vencNm0, phase-vencNm0
vencInd = contains(xSet.files_bidsList(:,1),'proc-venc');
partInd = contains(xSet.files_bidsList(:,1),'part-');
if any(vencInd)
    xSet.nVencDat = length(unique(fullfile(xSet.files_bidsList(vencInd,:),xSet.files_bidsList(partInd,:))));
    xSet.nVenc = length(unique(xSet.files_bidsList(vencInd,~ismember(xSet.files_bidsList(vencInd,:),'proc-venc0') & ismember(xSet.files_bidsList(partInd,:),'part-phase'))));
else
    xSet.nVenc    = 0;
    xSet.nVencDat = 1;
end


disp([xSet.label ': ' 'found ' num2str(xSet.nRun) 'run/' num2str(xSet.nEcho) 'echo/' num2str(xSet.nVenc) 'venc'])
