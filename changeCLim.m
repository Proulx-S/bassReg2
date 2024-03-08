function cmd = changeCLim(cmd,cLim)

if ~isempty(regexp(cmd,'.nii.gz -b \d+,\d+','match'))
    % When cLim is already set, change it
    dbstack; error('code that, or make sure to never have cLim already set')
else
    % When cLim is not already set, set it
    cmd = strjoin(strsplit(cmd,'.nii.gz'),['.nii.gz -b ' num2str(cLim(1)) ',' num2str(cLim(2))]);
    cmd = strjoin(strsplit(cmd,['manBrainMaskInv.nii.gz -b ' num2str(cLim(1)) ',' num2str(cLim(2))]),'manBrainMaskInv.nii.gz');
end