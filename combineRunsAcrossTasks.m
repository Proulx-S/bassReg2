function [fList,fMaskList,nDummy,taskList2,acqTime] = combineRunsAcrossTasks(rCond)
    
    taskList = fields(rCond); taskList = taskList(contains(taskList,'task_'));

    fList     = cell(size(taskList));
    fMaskList = cell(size(taskList));
    nDummy    = cell(size(taskList));
    acqTime   = cell(size(taskList));
    taskList2 = cell(size(taskList));
    for T = 1:length(taskList)
        fList{T}     = rCond.(taskList{T}).fPreprocList(:,1);
        ind = find(squeeze(all(~cellfun('isempty',rCond.(taskList{T}).fPreprocMaskList),1)),1,'last');
        fMaskList{T} = rCond.(taskList{T}).fPreprocMaskList(:,ind);
        nDummy{T}    = rCond.(taskList{T}).nFrameOrig - rCond.(taskList{T}).nFrame;
        taskList2{T} = repmat(taskList(T),size(fList{T}));
        acqTime{T}   = rCond.(taskList{T}).acqTime;
    end
    fList     = cat(1,fList{:}    );
    fMaskList = cat(1,fMaskList{:});
    nDummy    = cat(1,nDummy{:});
    taskList2 = cat(1,taskList2{:});
    acqTime   = cat(1,acqTime{:});
    [acqTime,b] = sort(acqTime);
    fList     = fList(b);
    fMaskList = fMaskList(b);
    nDummy    = nDummy(b);
    taskList2 = taskList2(b);