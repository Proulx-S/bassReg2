function cmd = qaFstMdLst(fList,force,verbose)
global srcAfni srcFs
if ~exist('force','var'); force = []; end
if ~exist('verbose','var'); verbose = []; end
if isempty(force); force = 0; end
if isempty(verbose); verbose = 0; end

%% Catenate first, middle and last frames
fTmpList = cell(size(fList));
cmd = {srcAfni};
for i = 1:length(fList)
    fIn = fList{i};
    fOut = fIn; fOut = strsplit(fIn,filesep); fOut{end} = ['fstMdLst_' fOut{end}]; fOut = strjoin(fOut,filesep);
    
    if ~exist(fOut,'file') || force
        sm = strsplit(fIn,filesep); sm = strsplit(sm{end},'_'); ind = ~cellfun('isempty',regexp(sm,'^sm\d+$')); if any(ind); sm = sm{ind}; else sm = 'sm1'; end; sm = str2num(sm(3:end));
        n = MRIread(fIn,1); n = n.nframes - 1;
        nLim = [0 n] + [1 -1].*((sm+1)/2-1);
        nLim = [nLim(1) round(mean(nLim)) nLim(2)];
        nLim = num2str(nLim,'%i,'); nLim(end) = [];

        cmd{end+1} = '3dcalc -overwrite \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = ['-a ' fIn '[' nLim '] \'];
        cmd{end+1} = '-expr ''a''';
    end

    fTmpList{i} = fOut;
end
if length(cmd)>1
    cmd = strjoin(cmd,newline); % disp(cmd)
    [status,cmdout] = system(cmd); if status; dbstack; error(cmdout); error('x'); end
end

%%%%% generate fslview command
cmd = {srcFs};
if isfield(fList,'manBrainMaskInv') && ~isempty(fList.manBrainMaskInv)
    cmd{end+1} = ['fslview -m single ' strjoin(fTmpList,' ') ' \'];
    cmd{end+1} = [fList.manBrainMaskInv ' &'];
else
    cmd{end+1} = ['fslview -m single ' strjoin(fTmpList,' ') ' &'];
end
cmd = strjoin(cmd,newline); % disp(cmd)
if verbose
    [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
end

