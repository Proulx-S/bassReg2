function cmd = mergeVisCmd(cmd1,tag1,cmd2,tag2)

% parse command into files and their associated tags
[cmd1f,cmd1t,cmdh] = parseCmd(cmd1);
[cmd2f,cmd2t] = parseCmd(cmd2);

% remove duplicates from the second command
dup1Ind = ismember(cmd1f,cmd2f);
dup2Ind = ismember(cmd2f,cmd1f);
cmdf = cmd1f(dup1Ind); cmdt = cmd1t(dup1Ind);
cmd1f(dup1Ind) = []; cmd1t(dup1Ind) = [];
cmd2f(dup2Ind) = []; cmd2t(dup2Ind) = [];

% add set to name tag
cmd1t = addNameTag(cmd1t,tag1);
cmd2t = addNameTag(cmd2t,tag2);

% reconstruct command
cmdf = [cmd1f; cmd2f; cmdf];
cmdt = [cmd1t; cmd2t; cmdt];
for i = 1:length(cmdt)
    cmd{i,1} = [cmdf{i} ':' strjoin(cmdt{i},':') ' \'];
end
cmd = [cmdh; cmd];
cmd{end} = replace(cmd{end},' \',' &');
cmd = strjoin(cmd,newline);
disp(cmd)


function cmdt = addNameTag(cmdt,tag)
for i = 1:length(cmdt)
    if isempty(cmdt{i})
        nameInd = 0;
    else
        nameInd = contains(cmdt{i},'name=');
    end
    if any(nameInd)
        % if name tag already present, append to it
        tmp = strsplit(cmdt{i}{nameInd},'=');
        tmp{2} = [tag '_' tmp{2}];
        cmdt{i}{nameInd} = strjoin(tmp,'=');
    else
        % if no name tag, create one
        cmdt{i} = [cmdt{i} {['name=' tag]}];
    end
end





function [cmdf,cmdt,cmdh] = parseCmd(cmd)
cmd = strsplit(cmd,newline)';
tmp = char(cmd);
cmdh = cmd(~ismember(tmp(:,1),'/'));
cmd = cmd(ismember(tmp(:,1),'/'));
for i = 1:length(cmd)
    cmd{i} = strsplit(cmd{i},':');
    cmdf{i,1} = cmd{i}{1};
    if length(cmd{i})>1
        cmdt{i,1} = replace(replace(cmd{i}(2:end),' \',''),' &','');
    else
        cmdt{i,1} = [];
    end
end

