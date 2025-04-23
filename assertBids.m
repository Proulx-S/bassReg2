function [S,str] = assertBids(rCond)

fList1 = dir(fullfile(rCond{1}.dirs.bids,'func','*.nii.gz')); fList1 = fullfile({fList1.folder},{fList1.name})';
fList2 = [rCond{:}];
for i = 1:numel(fList2)
    fList2(i).bhvr = repmat(fList2(i).bhvr,[1 size(fList2(i).fList,[2 3 4])]);
end
fBhvr  = {fList2.bhvr}';
fList2 = {fList2.fList}';
fBhvr  = fBhvr(~cellfun('isempty',fList2));
fList2 = fList2(~cellfun('isempty',fList2));
for r = 1:length(fList2)
    fList2{r} = fList2{r}(:);
    if isempty(fBhvr{r})
        fBhvr{r} = repmat({'noTask'},size(fList2{r}));
    else
        fBhvr{r} = {fBhvr{r}.par}';
    end
end
rcGrp  = {}; for i = 1:length(fList2); rcGrp{end+1} = num2str(i.*ones(size(fList2{i}))); end
fBhvr  = cat(1,fBhvr{:}); if ~iscell(fBhvr); fBhvr = {}; end
fList2 = cat(1,fList2{:}); if ~iscell(fList2); fList2 = {}; end
rcGrp  = cat(1,rcGrp{:} );% if ~iscell(rcGrp);  rcGrp  = {}; end
fBhvr  = replace(strtrim(fBhvr),'conf/','');


disp('db func files accouted for')
[~,b] = fileparts(fList2);
disp(char(strcat(rcGrp,'---',b,'---',fBhvr)))
disp('db func files NOT accouted for')
[~,b] = fileparts(fList1(~ismember(fList1,fList2)));
if isempty(b)
    disp('none')
else
    disp(char(b))
end



%% Identify cases where run set seems badly defined
fList2 = [rCond{:}];
fList2 = {fList2.fList};
fList2 = fList2(~cellfun('isempty',fList2));
S   = [];
str = {};
for s = 1:numel(fList2)
    if size(fList2{s},1) > 1
        tmp = fList2{s}(:,1);
        for r = 1:length(tmp)
            tmp{r} = strsplit(tmp{r},'_');
            tmp{r}(contains(tmp{r},'run-')) = [];
            tmp{r} = strjoin(tmp{r},'_');
        end
        [tmp,ia,ic] = unique(tmp);
        if length(tmp) > 1
            disp(' ')
            disp(' ')
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            disp('WARNING: run set seems badly defined')
            disp(['SET-' num2str(s) '--' strjoin({['sub-' rCond{s}.sub]
            ['ses-' rCond{s}.ses]
            ['acq-' rCond{s}.acq]
            ['task-' rCond{s}.task]},'_') newline '--------------------------------'])
            [~,b] = fileparts(replace(fList2{s}(:,:),'.nii.gz',''));
            disp(char(b))
            disp(' ')
            disp(['should be split in ' num2str(length(tmp)) 'sets (?):'])
            for i = 1:length(tmp)
                disp('--------------------------------')
                [~,b] = fileparts(replace(fList2{s}(ic==i,:),'.nii.gz',''));
                disp(char(b))
            end
            disp('--------------------------------')
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            disp(' ')
            disp(' ')

            S(end+1)   = s;
            str{end+1} = strjoin({['sub-' rCond{s}.sub]
            ['ses-' rCond{s}.ses]
            ['acq-' rCond{s}.acq]
            ['task-' rCond{s}.task]},'_');
        end
    end
end


