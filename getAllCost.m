function [cost,costLabel] = getAllCost(fList)
global srcAfni

nFrame = nan(size(fList));
for r = 1:length(fList)
    mri = MRIload3(fList{r});
    nFrame(r) = mri.nframes;
end

cost = cell(length(fList),length(fList));
for r1 = 1:length(fList)
    for r2 = 1:length(fList)
        tic
        if r2>r1; break; end
        disp(['cost run' num2str(r2) '->run' num2str(r1)])
        % costList{r1,r2} = cell(nFrame(r1),nFrame(r2));
        cost{r1,r2} = nan(nFrame(r1),nFrame(r2),14);
        for f1 = 1:nFrame(r1)
            disp(['cost run' num2str(r2) '->run' num2str(r1) 'frame' num2str(f1)])
            cmd = {srcAfni};
            cmd{end+1} = '3dAllineate -overwrite \';
            cmd{end+1} = ['-base   ' fList{r1} '[' num2str(f1-1) '] \'];
            cmd{end+1} = ['-source ' fList{r2} ' \'];
            cmd{end+1} = '-allcostX';
            [~,cmdout] = system(strjoin(cmd,newline));
            [cost{r1,r2}(f1,:,:),costLabel] = parseCost(cmdout);
        end
        toc
    end
end


function [cost,costLabel] = parseCost(cmdout)
cmdout = strsplit(cmdout,'++ allcost output:'); cmdout(1) = [];
cost      = nan(length(cmdout),14);
costLabel = cell(1            ,14);
for i = 1:length(cmdout)
    cmdout{i} = strsplit(cmdout{i},newline)';
    cmdout{i} = cmdout{i}(2:15);
    for ii = 1:length(cmdout{i})
        cmdout{i}{ii} = strsplit(cmdout{i}{ii},'=');
        cost(i,ii) = str2double(cmdout{i}{ii}{2});
        if i==1
            costLabel{ii} = char(replace(cmdout{i}{ii}(1),' ',''));
        end
    end
end