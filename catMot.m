function motFiles = catMot(xSet,candidate)

motFiles = {};
candidate = candidate(ismember(candidate,fields(xSet)));
if isempty(candidate); return; end

for C = 1:length(candidate)
    motFiles{C} = xSet.(candidate{C});
end



% C = 1;
% motFiles = xSet.(candidate{C});
% if length(candidate)==1; return; end
% 
% for C = 2:length(candidate)
%     motFiles.fMocoMat = cat(10,motFiles.fMocoMat,xSet.(candidate{C}).fMocoMat);
% end
