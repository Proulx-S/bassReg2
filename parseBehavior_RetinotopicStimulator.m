function d = parseBehavior_RetinotopicStimulator(fd,mriTime)
% need to find a way to estimate when was the first trigger received...
if ischar(fd)
    f = fd; clear fd
    %% Parse file
    M = strsplit(fileread(f),'----------------------------------------------------------------------------')';
    for i = 1:length(M); if ~contains(M{i},'/bin/retinotopicStimulator'); M{i} = []; end; end
    M(cellfun('isempty',M)) = [];
    aborted      = true(size(M));
    for i = 1:length(M)
        tmp = strsplit(M{i},newline)';
        if any(contains(tmp,'stim exit status:'))
            aborted(i) = contains(tmp{contains(tmp,'stim exit status:')},'aborted');
        end
    end
    % M(aborted) = [];

    d = struct();
    for iM = 1:length(M)
        % if aborted(iM)
        %     %time
        %     d(iM,1).time = [];
        %     %COMMAND
        %     d(iM,1).command = '';
        %     %par
        %     d(iM,1).par = '';
        %     %note
        %     d(iM,1).note = '';
        %     %index
        %     d(iM,1).iM = iM;
        %     d(iM,1).Mi = M{iM};
        %     %file
        %     d(iM,1).f = f;
        % else
            Mi = strsplit(M{iM},newline)';
            %time
            tmp = strsplit(Mi{3},'_');
            d(iM,1).time = datetime([strjoin(tmp(1:3),'-') 'T' strjoin(tmp(4:6),':')]);
            %performance
            if any(contains(Mi,'performance:'))
                tmp = Mi{contains(Mi,'performance:')}; tmp = strsplit(tmp,':');
                d(iM,1).performance = tmp{2};
            else
                d(iM,1).performance = [];
            end
            %false positive
            if any(contains(Mi,'percent false positives:'))
                tmp = Mi{contains(Mi,'percent false positives:')}; tmp = strsplit(tmp,':');
                d(iM,1).percFalsePositive = tmp{2};
            else
                d(iM,1).percFalsePositive = [];
            end
            %correct
            if any(contains(Mi,'percent correct:'))
                tmp = Mi{contains(Mi,'percent correct:')}; tmp = strsplit(tmp,':');
                d(iM,1).percTruePositive = tmp{2};
            else
                d(iM,1).percTruePositive = [];
            end
            %COMMAND
            tmp = Mi{contains(Mi,'COMMAND')}; tmp = strsplit(tmp,':');
            d(iM,1).command = tmp{2};
            %par
            tmp = Mi{contains(Mi,'paradigm file:')}; tmp = strsplit(tmp,':');
            d(iM,1).par = tmp{2};
            %note
            if any(contains(Mi,'note:'))
                tmp = Mi{contains(Mi,'note:')}; tmp = strsplit(tmp,':');
                d(iM,1).note = tmp{2};
            else
                d(iM,1).note = [];
            end
            %index
            d(iM,1).iM = iM;
            d(iM,1).Mi = M{iM};
            %file
            d(iM,1).f = f;
            %aborted
            d(iM,1).aborted = aborted(iM);
        % end
    end
    % d(contains({d.note}','test')) = [];
    d(cellfun('isempty',{d.par}')) = [];
elseif isstruct(fd)
    d = fd;
end

%% Get frame rate file to adjust for actual trigger time
[a,b,~] = fileparts(f);
tridDir = dir(a);
tridDir = tridDir([tridDir.isdir]);
tridDir = tridDir(ismember({tridDir.name},{'log' b}));
tridDir = fullfile({tridDir.folder},{tridDir.name});
if numel(tridDir)>1; dbstack; error('X'); end
tridDir = char(tridDir);
for i = 1:length(d)
    frameFile = dir(fullfile(tridDir,['frameRate' char(datetime(d(i).time,"Format",'uuuu_MM_dd_HH_mm_ss'))]));
    if isempty(frameFile)
        d(i).frameFile    = [];
        d(i).trigTime    = [];
        d(i).timeAdj     = [];
        d(i).hasTrigTime = false;
    else
        frameFile = fullfile(frameFile.folder,frameFile.name);
        trigTime = strsplit(fileread(frameFile),newline)';
        trigTime = strsplit(trigTime{1},' ');
        d(i).frameFile    = frameFile;
        d(i).trigTime    = seconds(str2double(trigTime{end})/1000);
        d(i).timeAdj     = d(i).time + d(i).trigTime;
        d(i).hasTrigTime = true;
    end
end



%% Get the ones that most closely match mriTime
% d = matchWithMri(d,mriTime,0);
d = matchWithMri(d([d.hasTrigTime]),mriTime,1);

function d = matchWithMri(d,mriTime,trigAdjFlag)
di       = zeros(size(mriTime));
dTimeAdj = repmat(duration,size(di));
dTime    = repmat(duration,size(di));
for i = 1:length(mriTime)
    td    = [d.time]'    - mriTime(i);
    tdAdj = [d.timeAdj]' - mriTime(i);
    if trigAdjFlag
        % [~,b] = sort(tdAdj); b = max(b(tdAdj<=0));
        [~,b] = min(abs(tdAdj));
    else
        [~,b] = sort(td); b = max(b(td<=0));
    end
    di(i) = b;
    dTimeAdj(i) = tdAdj(b);
    dTime(i) = td(b);
end
d = d(di);
for di = 1:length(d)
    d(di).timeAdjDelay = dTimeAdj(di);
    d(di).timeDelay = dTime(di);
    d(di).mriTime = mriTime(di);
end

