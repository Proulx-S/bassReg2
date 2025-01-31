function [dummyParam, dummyMatrix] = getDummyMcParamFile(cmd,frame1flag)
if ~exist('frame1flag','var'); frame1flag = []; end
if isempty(frame1flag);        frame1flag = 0; end

%% Alter command to produce dummy param.1D and aff12.1D mc param files that include the header and all parameters
% consider changing base to source

cmd(contains(cmd,'-prefix')) = [];

if frame1flag
    tmp = cmd{contains(cmd,'-source')}; tmp = strsplit(tmp,' \'); tmp{1} = strsplit(tmp{1},'['); tmp{1} = tmp{1}{1}; tmp{1} = [tmp{1} '[0]']; tmp = strjoin(tmp,' \');
    cmd{contains(cmd,'-source')} = tmp;
end

tmp = cmd{contains(cmd,'-1Dparam_save')};
tmp = strsplit(tmp,' ');
tmp{2} = [tmp{2} '.dummy'];
dummyParam = [tmp{2} '.param.1D'];
tmp = strjoin(tmp,' ');
cmd{contains(cmd,'-1Dparam_save')} = tmp;

tmp = cmd{contains(cmd,'-1Dmatrix_save')};
tmp = strsplit(tmp,' ');
tmp{2} = [tmp{2} '.dummy'];
dummyMatrix = [tmp{2} '.aff12.1D'];
tmp = strjoin(tmp,' ');
cmd{contains(cmd,'-1Dmatrix_save')} = tmp;

cmd(contains(cmd,'-warp')) = [];

cmd{end} = replace(cmd{end},' \','');

[status,cmdout] = system(strjoin(cmd,newline)); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end


%% Rewrite param.1D
fid = fopen(dummyParam);
tline = {fgetl(fid)};
while ischar(tline{end})
    tline{end+1,1} = fgetl(fid);
end
fclose(fid);

tline{3} = strsplit(tline{3},' ');
for i = 1:length(tline{3})
    if isempty(tline{3}{i}); continue; end
    tline{3}{i} = num2str(round(str2double(tline{3}{i})),'%0.6f');
end
tline{3} = strjoin(tline{3},' ');

tmp = tempname;
fid = fopen(tmp, 'w');
for i = 1:numel(tline)
    if tline{i+1} == -1
        fprintf(fid,'%s', tline{i});
        break
    else
        fprintf(fid,'%s\n', tline{i});
    end
end
fclose(fid);
movefile(tmp,dummyParam)


%% Rewrite aff12.1D
fid = fopen(dummyMatrix);
tline = {fgetl(fid)};
while ischar(tline{end})
    tline{end+1,1} = fgetl(fid);
end
fclose(fid);

tline{2} = strsplit(tline{2},' ');
for i = 1:length(tline{2})
    if isempty(tline{2}{i}); continue; end
    tline{2}{i} = num2str(round(str2double(tline{2}{i})),'%0.6f');
end
tline{2} = strjoin(tline{2},' ');

tmp = tempname;
fid = fopen(tmp, 'w');
for i = 1:numel(tline)
    if tline{i+1} == -1
        fprintf(fid,'%s', tline{i});
        break
    else
        fprintf(fid,'%s\n', tline{i});
    end
end
fclose(fid);
movefile(tmp,dummyMatrix)
