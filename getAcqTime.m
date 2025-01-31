function [fileTime,runDur] = getAcqTime(fileList)

% fileTime = repmat(datetime,size(fileList));
fileTime = repmat(duration,size(fileList));
for I = 1:length(fileList)
    jsnFile = replace(fileList{I},'.nii.gz','.json');
    cmd = ['jq ''.AcquisitionTime'' ' jsnFile];
    [~,cmdout] = system(cmd);
    % tmpTime = strsplit(cmdout,':'); tmpTime{1} = tmpTime{1}(end-1:end); tmpTime{3} = strsplit(tmpTime{3},'.'); tmpTime{3}{2} = strsplit(tmpTime{3}{2},'"'); tmpTime{3}{2} = tmpTime{3}{2}{1}; tmpTime{3} = strjoin(tmpTime{3},'.'); tmpTime = strjoin(tmpTime,':');
    % fileTime(I) = datetime(['0001-01-01 ' tmpTime],'InputFormat','yyyy-MM-dd HH:mm:ss.SSSSSS');
    tmpTime = strsplit(cmdout,':'); tmpTime{1} = tmpTime{1}(end-1:end); tmpTime{3} = strsplit(tmpTime{3},'.'); tmpTime{3}{2} = strsplit(tmpTime{3}{2},'"'); tmpTime{3}{2} = tmpTime{3}{2}{1}; tmpTime{3} = strjoin(tmpTime{3},'.'); tmpTime = str2double(tmpTime);
    fileTime(I) = duration(tmpTime(1),tmpTime(2),tmpTime(3));
end

