function fileTime = getAcqTime(fileList)

fileTime = repmat(datetime,size(fileList));
for I = 1:length(fileList)
    jsnFile = replace(fileList{I},'.nii.gz','.json');
    cmd = ['jq ''.AcquisitionTime'' ' jsnFile];
    [status,cmdout] = system(cmd);
    tmpTime = strsplit(cmdout,':'); tmpTime{1} = tmpTime{1}(end-1:end); tmpTime{3} = strsplit(tmpTime{3},'.'); tmpTime{3}{2} = strsplit(tmpTime{3}{2},'"'); tmpTime{3}{2} = tmpTime{3}{2}{1}; tmpTime{3} = strjoin(tmpTime{3},'.'); tmpTime = strjoin(tmpTime,':');
    fileTime(I) = datetime(tmpTime,'InputFormat','HH:mm:ss.SSSSSS');
end

