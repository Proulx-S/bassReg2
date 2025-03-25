function assertBehavior_RetinotopicStimulator2(rCond)


% make sure the stimulus condition matches between the
% MRI filename and the behavior file info
mriFile   = {};
logFile  = {};
frameFile = {};
par       = {};
task      = {};
acq       = {};
rcInd     = {};
t         = {};
tAdj      = {};
tDelay    = {};
tAdjDelay = {};
note      = {};
mriTime   = {};
for rc = 1:length(rCond)
    if isempty(rCond{rc}.fList); continue; end
    mriFile{end+1,1} = rCond{rc}.fList(:,1);
    task{end+1,1}    = repmat({rCond{rc}.task},size(rCond{rc}.fList,1),1);
    acq{end+1,1}     = repmat({rCond{rc}.acq} ,size(rCond{rc}.fList,1),1);
    rcInd{end+1,1}   = repmat({rc},size(rCond{rc}.fList,1),1);
    if isempty(rCond{rc}.bhvr)
        logFile{end+1,1}  = 'x';
        frameFile{end+1,1} = 'x';
        par{end+1,1}       = 'x';
        t{end+1,1}         = 'x';
        tDelay{end+1,1}    = 'x';
        tAdj{end+1,1}      = 'x';
        tAdjDelay{end+1,1} = 'x';
        note{end+1,1}      = 'x';
        mriTime{end+1,1}   = 'x';
    else
        logFile{end+1,1} = {rCond{rc}.bhvr.f}';
        frameFile{end+1,1}  = {rCond{rc}.bhvr.frameFile}';
        par{end+1,1}       = {rCond{rc}.bhvr.par}';
        t{end+1,1}         = {rCond{rc}.bhvr.time}';
        tDelay{end+1,1}    = {rCond{rc}.bhvr.timeDelay}';
        tAdj{end+1,1}      = {rCond{rc}.bhvr.timeAdj}';
        tAdjDelay{end+1,1} = {rCond{rc}.bhvr.timeAdjDelay}';
        note{end+1,1}      = {rCond{rc}.bhvr.note}';
        mriTime{end+1,1}   = {rCond{rc}.bhvr.mriTime}';
    end
end
mriFile   = cat(1,mriFile{:});
frameFile = cat(1,frameFile{:});
logFile   = cat(1,logFile{:});
par       = cat(1,par{:});
task      = cat(1,task{:});
acq       = cat(1,acq{:});
rcInd     = cat(1,rcInd{:});
t         = cat(1,t{:});
tDelay    = cat(1,tDelay{:});
tAdj      = cat(1,tAdj{:});
tAdjDelay = cat(1,tAdjDelay{:});
note      = cat(1,note{:});
mriTime   = cat(1,mriTime{:});

if isempty(mriFile); disp('no MRI to assert behavior for'); return; end

% par = replace(replace(par,' 			conf/vsmDriven_',''),'.par','');
par = replace(par,' 			conf/vsmDriven_','');

% ifOk = cell(size(par));
% for rc = 1:length(par)
%     if strcmp(par{rc},'x')
%         ifOk{rc} = true;
%     else
%         ifOk{rc} = ...
%             ~isempty(strfind(par{rc},task{rc}(1:4))) & ...
%             ~isempty(strfind(par{rc},task{rc}(7:9)));
%     end
% end

b = zeros(size(mriTime)); [~,b(cellfun('isclass',mriTime,'datetime'))] = sort([mriTime{cellfun('isclass',mriTime,'datetime')}]); b(b==0) = length(b)-nnz(b==0)+1:length(b);
% disp([{'stimFile' 'stimNote' 'mriStimLabel' 'stimTime' '->trigAdjusted' 'mriTimeDelay' '->TrigAdjusted' 'looksOk' 'condInd'}
%         par(b)     note(b)    label(b)       t(b)       tAdj(b)          tDelay(b)      tAdjDelay(b)     ifOk(b)   rcInd(b)]);
% disp([{'stimFile' 'stimNote' 'mriStimLabel' 'stimTime' 'mriTimeDelay' 'looksOk' 'condInd'}
%         par        note       label          tAdj       tAdjDelay      ifOk      rcInd]);
% table( par(b),    note(b),   acq(b),       task(b),       t(b),      tAdjDelay(b), ifOk(b),  rcInd(b),'VariableNames',...
%     { 'stimFile' 'stimNote' 'mriAcqLabel' 'mriTaskLabel' 'stimTime' 'mriDelayAdj' 'looksOk' 'condInd'})
% disp(table( par(b),    note(b),   acq(b),       task(b),       t(b),      tAdjDelay(b),  rcInd(b),'VariableNames',...
%     { 'stimFile' 'stimNote' 'mriAcqLabel' 'mriTaskLabel' 'stimTime' 'mriDelayAdj'  'condInd'}))
disp(table( note(b),   par(b),    task(b),       acq(b),       t(b),      tAdjDelay(b), rcInd(b)  ,'VariableNames',...
           {'stimNote' 'stimFile' 'mriTaskLabel' 'mriAcqLabel' 'stimTime' 'mriDelayAdj' 'condInd'}))

% if any(~[ifOk{:}]) || length(unique([t{cellfun('isclass',t,'datetime')}]))~=nnz(cellfun('isclass',t,'datetime'))
%     tmp = [rCond{:}];
%     sub = char(unique({tmp.sub}));
%     ses = char(unique({tmp.ses}));
%     warning(['sub-' num2str(sub) '_ses-' num2str(ses) ', some behavioral data seem to not match mri data']);
%     tmpInd = find(~[ifOk{:}]);
%     for i = 1:length(tmpInd)
%         disp(cat(1,mriFile(tmpInd(i)),logFile(tmpInd(i)),frameFile(tmpInd(i))))
%     end
% end
