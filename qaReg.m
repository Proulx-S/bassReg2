function cmd = qaReg(xSet,fileLabel)
global srcFs
if ~exist('fileLabel','var'); fileLabel = ''     ; end
if isempty(fileLabel);        fileLabel = 'final'; end


cmd = {srcFs};
cmd{end+1} = 'freeview \';

for S = 1:length(xSet)

    switch fileLabel

        case 'final'
            fieldList = fields(xSet{S}.finalFiles);
            
            fieldList(contains(fieldList,'rainMask')) = [];
            fieldList(contains(fieldList,'EchoCat')) = [];

            if any(contains(fieldList,'EchoRms'))
                fieldList(~contains(fieldList,'EchoRms')) = [];
            end
            if any(contains(fieldList,'fCatAv'))
                fieldList = fieldList(contains(fieldList,'fCatAv'));
            elseif any(contains(fieldList,'fAv'))
                fieldList = fieldList(contains(fieldList,'fAv'));
            elseif any(ismember(fieldList,'f')) || any(ismember(fieldList,'fEchoRms'))
                ismember(fieldList,'f') | ismember(fieldList,'fEchoRms')
                fieldList = fieldList(ismember(fieldList,'f') | ismember(fieldList,'fEchoRms'));
                % fieldList = fieldList(contains(fieldList,'fAv'));
            end

            if length(fieldList)~=1; dbstack; error('can''t find what to show'); end

            if xSet{S}.nVenc
                cmd{end+1} = [char(xSet{S}.finalFiles.(char(fieldList)){1}) ':resample=cubic:name=' replace(xSet{S}.label,'set-','') ' \'];
            else
                cmd{end+1} = [char(xSet{S}.finalFiles.(char(fieldList))) ':resample=cubic:name=' replace(xSet{S}.label,'set-','') ' \'];
            end



        case 'preproc'

            fieldList = fields(xSet{S}.preprocFiles);
            if any(contains(fieldList,'EchoRms'))
                fieldList(~contains(fieldList,'EchoRms')) = [];
            end
            if any(contains(fieldList,'fCorrectedCatAv'))
                fieldList = fieldList(contains(fieldList,'fCorrectedCatAv'));
            elseif any(contains(fieldList,'fCorrectedAv'))
                fieldList = fieldList(contains(fieldList,'fCorrectedAv'));
            elseif any(contains(fieldList,'fCorrected'))
                fieldList = fieldList(contains(fieldList,'fCorrected'));
            end

            if length(fieldList)>1; dbstack; error('can''t find what to show'); end

            if xSet{S}.nVenc
                cmd{end+1} = [xSet{S}.preprocFiles.(char(fieldList)){1} ':resample=cubic:name=' replace(xSet{S}.label,'set-','') ' \'];
            else
                cmd{end+1} = [char(xSet{S}.preprocFiles.(char(fieldList))) ':resample=cubic:name=' replace(xSet{S}.label,'set-','') ' \'];
            end


        otherwise
            dbstack; error('please specify valid fileLabel')
    end
end



switch fileLabel

    case 'final'
        cmd{end+1} = [xSet{1}.finalFiles.manBrainMaskInv ':name=mask &'];

    case 'preproc'
        cmd{end+1} = [xSet{1}.preprocFiles.manBrainMaskInv ':name=mask &'];

end

cmd = strjoin(cmd,newline);
disp('QA registration')
disp(cmd);
clipboard('copy',cmd);
disp('command copied to clipboard');
