function brMocoFiles = wr2br(wrMocoFiles)

fieldList = fields(wrMocoFiles);
fieldList = fieldList(contains(fieldList,{'fEstim' 'fMoco'}));
brMocoFiles = rmfield(wrMocoFiles,fieldList);
brMocoFiles.nFrame = ones(size(brMocoFiles.nFrame));
brMocoFiles.dataType(ismember(brMocoFiles.dataType,'volTs')) = {'vol'};