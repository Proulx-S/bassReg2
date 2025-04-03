function [aSet,subListU,QA] = runSet_combSes(rSet,subList,sesList)

%% Combine cells from multiple sessions into a single cell per subject
subListU = unique(subList);
aSet     = cell(size(subListU));
for S = 1:length(subListU)
    cSet = rSet(ismember(subList,subListU{S}));
    cSes = sesList(ismember(subList,subListU{S})); for s = 1:length(cSes); cSes{s} = repmat(cSes(s),size(cSet{s})); end
    cSet = cat(2,cSet{:}); cSet = [cSet{:}];
    cSes = cat(2,cSes{:});

    setList  = {cSet.label};
    setListU = unique(setList);
    
    aSet{S} = cell(size(setListU));
    for s = 1:length(setListU)
        aSet{S}{s} = cSet(ismember(setList,setListU{s}))';
    end    
end

%% Get just filenames for QA
QA.label        = cell(size(aSet));
QA.fOrigList    = cell(size(aSet));
QA.fPreprocList = cell(size(aSet));
QA.fMaskList    = cell(size(aSet));
QA.nDummy       = cell(size(aSet));
for S = 1:length(aSet)
    QA.label{S}        = cell(size(aSet{S}));
    QA.fOrigList{S}    = cell(size(aSet{S}));
    QA.fPreprocList{S} = cell(size(aSet{S}));
    QA.fMaskList{S}    = cell(size(aSet{S}));
    QA.nDummy{S}       = cell(size(aSet{S}));
    for A = 1:length(aSet{S})
        QA.label{S}{A}        = {aSet{S}{A}.label}';
        finalFiles = [aSet{S}{A}.finalFiles];
        QA.fOrigList{S}{A}    = cat(1,finalFiles.fOrigList);
        QA.fPreprocList{S}{A} = cat(1,finalFiles.fPreprocList);
        QA.fMaskList{S}{A} = {};
        for s = 1:length(finalFiles)
            QA.fMaskList{S}{A} = cat(1,QA.fMaskList{S}{A},repmat(cellstr(aSet{S}{A}(s).fMasks.fMaskInv),size(finalFiles(s).fPreprocList)));
        end
        QA.nDummy{S}{A}       = cat(1,finalFiles.nDummy);
    end    
end
