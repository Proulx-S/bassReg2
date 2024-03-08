function files = addMaskToCmd(files)
fieldList1 = fields(files); fieldList1 = fieldList1(contains(fieldList1,'qaFiles'));
if ~isempty(fieldList1) && isfield(files,'manBrainMaskInv') && ~isempty(files.manBrainMaskInv)
    for i1 = 1:length(fieldList1)
        fieldList2 = fields(files.(fieldList1{i1})); fieldList2 = fieldList2(contains(fieldList2,'Fslview'));
        for i2 = 1:length(fieldList2)
            if ~contains(files.(fieldList1{i1}).(fieldList2{i2}),files.manBrainMaskInv)
                files.(fieldList1{i1}).(fieldList2{i2}) = replace(files.(fieldList1{i1}).(fieldList2{i2}),' &',[' \' newline files.manBrainMaskInv ' &']);
            end
        end
    end
end

% if isfield(files,'qaFiles')
%     fieldList2 = fields(files.qaFiles); fieldList2 = fieldList2(contains(fieldList2,'Fslview'));
%     for i1 = 1:length(fieldList2)
%         if ~contains(files.qaFiles.(fieldList2{i1}),files.manBrainMaskInv)
%             files.qaFiles.(fieldList2{i1}) = replace(files.qaFiles.(fieldList2{i1}),' &',[' \' newline files.manBrainMaskInv ' &']);
%         end
%     end
% end