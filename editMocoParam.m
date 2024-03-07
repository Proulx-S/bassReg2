function editMocoParam(curFile,tSmWin_vol)

if ~exist("tSmWin_vol",'var'); tSmWin_vol = []; end
if isempty(tSmWin_vol)
    tSmWin_vol = strsplit(fIn,filesep); tSmWin_vol = strsplit(tSmWin_vol{end},'_'); ind = ~cellfun('isempty',regexp(tSmWin_vol,'^sm\d+$')); if any(ind); tSmWin_vol = tSmWin_vol{ind}; else; tSmWin_vol = 'sm1'; end; tSmWin_vol = str2double(tSmWin_vol(3:end));
end

fid = fopen(curFile,'r');
%read all lines
tline = {fgetl(fid)}; while tline{end}(1)~=-1; tline{end+1} = fgetl(fid); end; tline = tline(1:end-1)';
%find index of first data line
i = 1; while strcmp(tline{i}(1),'#'); i = i+1; end; i = i-1;
%replace the first tSmWin_vol/2 lines
tline((1:floor(tSmWin_vol/2))+i) = tline(floor(tSmWin_vol/2)+i+1);
%replace the last tSmWin_vol/2 lines
tline(end-floor(tSmWin_vol/2)+1:end) = tline(end-floor(tSmWin_vol/2));
fclose(fid);

%rewrite the file
tmpFile = tempname;
fid = fopen(tmpFile, 'w');
for i = 1:numel(tline)
    if i == numel(tline)
        fprintf(fid,'%s', tline{i});
        break
    else
        fprintf(fid,'%s\n', tline{i});
    end
end
fclose(fid);
movefile(tmpFile,curFile);