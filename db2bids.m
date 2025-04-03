function [dirs,dirsOrig] = db2bids(sesDb,sub,ses,info,force)
if ~exist('force','var'); force = []; end
if isempty(force);        force = 0 ; end


%% mri --> bids
bidsSubDir = {'anat' 'func' 'fmap'}';
for bs = 1:length(bidsSubDir)
    from = fullfile(sesDb,'bids',bidsSubDir{bs});
    to   = fullfile(info.prcDir,'bids',['sub-' sub],['ses-' ses],bidsSubDir{bs});
    disp(['-' bidsSubDir{bs}])
    disp(from);
    disp('to');
    disp(to)
    if exist(from,'dir') && ~isempty(dir(from))
        if force || isempty(dir(to))
            if exist(to,'dir'); system(['rm -r ' to]); end
            copyfile(from,to);
            disp('done')
        else
            disp('already done, skipping it')
        end
    else
        disp('bids data not found')
    end
end
dirsOrig.bids = fullfile(sesDb,'bids');
dirsOrig.bidsDeriv = fullfile(sesDb,'bids','derivatives');
dirs.bids     = fullfile(info.prcDir,'bids',['sub-' sub],['ses-' ses]);
dirs.bidsDeriv     = fullfile(info.prcDir,'bids','derivatives',['sub-' sub],['ses-' ses]);

renameBids(dirs.bids,sub,ses);



%% behav --> source
from = fullfile(sesDb,'stim');
to   = fullfile(info.prcDir,'source','retinoPoli',['sub-' sub],['ses-' ses]);
disp('-bhvr')
disp(from); disp('to'); disp(to)
if exist(from,'dir') && ~isempty(dir(from))
    if force || isempty(dir(to))
        if exist(to,'dir'); system(['rm -r ' to]); end
        copyfile(from,to);
        disp('done')
    else
        disp('already done, skipping it')
    end
    dirsOrig.bhvr = from;
    dirs.bhvr     = to;
else
    dirsOrig.bhvr = [];
    dirs.bhvr     = [];
    disp('no data available')
end

%% physio --> source
from = fullfile(sesDb,'physio');
to   = fullfile(info.prcDir,'source','labChart',['sub-' sub],['ses-' ses]);
disp('-phs')
disp(fullfile(from,'*.mat')); disp('to'); disp(to)
if exist(from,'dir') && ~isempty(dir(fullfile(from,'*.mat')))
    if force || isempty(dir(to))
        if exist(to,'dir'); system(['rm -r ' to]); end
        copyfile(fullfile(from,'*.mat'),to);
        disp('done')
    else
        disp('already done, skipping it')
    end
    dirsOrig.phs = from;
    dirs.phs     = to;
else
    dirsOrig.phs = [];
    dirs.phs     = [];
    disp('no data available')
end


function renameBids(bidsDir,sub,ses)
    folder = {};
    fOrig    = {};
    fNew   = {};
    fList  = dir(bidsDir);
    for f = 1:length(fList)
        if ~fList(f).isdir
            folder{end+1,1} = fList(f).folder;
            fOrig{end+1,1}    = fList(f).name;
            fNew{end+1,1}   = strsplit(fOrig{end},'_');
            fNew{end}{contains(fNew{end},'sub-')} = ['sub-' sub];
            fNew{end}{contains(fNew{end},'ses-')} = ['ses-' ses];
            fNew{end} = strjoin(fNew{end},'_');
            movefile(fullfile(folder{end},fOrig{end}),fullfile(folder{end},fNew{end}));
        end
    end

    % Here we might want to also update run- for cases where they are not an uninterupted series.
    % This however would run the risk of a mismatch with the acquisition notes, so let's not do that.
