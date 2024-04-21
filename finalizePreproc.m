function fFinal = finalizePreproc(fPreproc,fInit,force)
if ~exist('force','var'); force = []; end
if isempty(force);        force = 0; end
if fPreproc.updatedFlag;  force = 1; end % override force=0 when preproc was updated

%% Define all filenames
fieldList1 = {'fCorrected' 'fCorrectedAv'};
fieldList2 = {'f'          'fAv'};

fieldList2 = fieldList2(ismember(fieldList1,fields(fPreproc)));
fieldList1 = fieldList1(ismember(fieldList1,fields(fPreproc)));
for i = 1:length(fieldList1)
    fFinal.(fieldList2{i}) = cell(size(fPreproc.fOrig));
    for I = 1:numel(fPreproc.fOrig)
        fIn = fPreproc.(fieldList1{i}){I};
        fOut = replace(fIn,'setPlumb_',''); fOut = strsplit(fOut,filesep); fOut = strjoin([fOut(1:end-1) {'preproc'} fOut(end)],filesep); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        fFinal.(fieldList2{i}){I} = fOut;

        fInListA{I}{i} = fIn;
        fOutListA{I}{i} = fOut;
    end
end




fieldList1 = {'fCorrectedEchoRms' 'fCorrectedAvEchoRms' 'fCorrectedAvEchoCat' 'fCorrectedEchoRms' 'fCorrectedAvEchoRms'};
fieldList2 = {'fEchoRms'          'fAvEchoRms'          'fAvEchoCat'          'fEchoRms'          'fAvEchoRms'};
clear fInListA2 fOutListA2
for i = 1:length(fieldList1)
    if ~isfield(fPreproc,fieldList1{i}); continue; end
    fFinal.(fieldList2{i}) = cell(size(fPreproc.(fieldList1{i})));
    for I = 1:size(fPreproc.(fieldList1{i}),1)
        fIn = fPreproc.(fieldList1{i}){I};
        fOut = replace(fIn,'setPlumb_',''); fOut = strsplit(fOut,filesep); fOut = strjoin([fOut(1:end-1) {'preproc'} fOut(end)],filesep); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        fFinal.(fieldList2{i}){I} = fOut;

        fInListA2{I}{i} = fIn;
        fOutListA2{I}{i} = fOut;
    end
end
if exist('fInListA2','var') && exist('fOutListA2','var')
    for I = 1:length(fInListA2)
        ind = cellfun('isempty',fInListA2{I});
        fInListA2{I}(ind) = [];
        fOutListA2{I}(ind) = [];
    end

    fInListA = cat(2,fInListA,fInListA2);
    fOutListA = cat(2,fOutListA,fOutListA2);
end


fieldList1 = {'fCorrectedCatAv' 'fCorrectedAvCatAv'};
fieldList2 = {'fCatAv'          'fAvCatAv'};
fInListB = {};
fOutListB = {};
for i = 1:length(fieldList1)
    if ~isfield(fPreproc,fieldList1{i}); continue; end
    for E = 1:size(fPreproc.(fieldList1{i}),2)
        fIn = fPreproc.(fieldList1{i}){E};
        fOut = replace(fIn,'setPlumb_',''); fOut = strsplit(fOut,filesep); fOut = strjoin([fOut(1:end-1) {'preproc'} fOut(end)],filesep); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        fFinal.(fieldList2{i}){1,E} = fOut;

        fInListB{end+1} = fIn;
        fOutListB{end+1} = fOut;
    end
end


fieldList1 = {'fCorrectedAvCatAvEchoCat' 'fCorrectedCatAvEchoRms' 'fCorrectedAvCatAvEchoRms'};
fieldList2 = {'fAvCatAvEchoCat'          'fCatAvEchoRms'          'fAvCatAvEchoRms'         };
% fieldList1 = {'fCorrectedCatAvEchoRms' 'fCorrectedAvCatAvEchoRms' 'fCorrectedCatAvEchoCat' 'fCorrectedAvCatAvEchoCat'};
% fieldList2 = {'fCatAvEchoRms' 'fAvCatAvEchoRms' 'fCatAvEchoCat' 'fAvCatAvEchoCat'};
fInListB2 = {};
fOutListB2 = {};
for i = 1:length(fieldList1)
    if ~isfield(fPreproc,fieldList1{i}); continue; end
        fIn = char(fPreproc.(fieldList1{i}));
        fOut = replace(fIn,'setPlumb_',''); fOut = strsplit(fOut,filesep); fOut = strjoin([fOut(1:end-1) {'preproc'} fOut(end)],filesep); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        fFinal.(fieldList2{i}) = {fOut};

        fInListB2{end+1} = fIn;
        fOutListB2{end+1} = fOut;
end

fInListB = cat(2,fInListB,fInListB2);
fOutListB = cat(2,fOutListB,fOutListB2);



fInListB{end+1} = fPreproc.manBrainMask;
fOutListB{end+1} = replace(fInListB{end},'.nii.gz','Oblique.nii.gz');
fFinal.manBrainMask = fOutListB{end};

fInListB{end+1} = fPreproc.manBrainMaskInv;
fOutListB{end+1} = replace(fInListB{end},'.nii.gz','Oblique.nii.gz');
fFinal.manBrainMaskInv = fOutListB{end};

%% Rewrite all to oblique
obliqueRef = MRIread(fInit.fObliqueRef,1);
disp('writing back to set oblique (for visualization with other images)')
for I = 1:length(fInListA)
    disp([' file ' num2str(I) '/' num2str(length(fInListA))])
    for i = 1:length(fInListA{I})
        if force || ~exist(fOutListA{I}{i},'file')
            mriTmp     = MRIread(fInListA{I}{i});
            mri        = obliqueRef;
            mri.vol    = mriTmp.vol;
            mri.xsize  = mriTmp.xsize;
            mri.ysize  = mriTmp.ysize;
            mri.zsize  = mriTmp.zsize;
            mri.volres = mriTmp.volres;
            mri.tr     = mriTmp.tr;
            mri.te     = mriTmp.te;
            mri.ti     = mriTmp.ti;
            MRIwrite(mri,fOutListA{I}{i});
        end
        disp(['  in : ' fInListA{I}{i}])
        disp(['  out: ' fOutListA{I}{i}])
    end
    disp(' done')
end
disp(' runCat')
for i = 1:length(fInListB)
    if force || ~exist(fOutListB{i},'file')
        mriTmp = MRIread(fInListB{i});
        obliqueRef.vol = mriTmp.vol;
        MRIwrite(obliqueRef,fOutListB{i});
    end
    disp(['  in : ' fInListB{i}])
    disp(['  out: ' fOutListB{i}])
end
disp(' done')

