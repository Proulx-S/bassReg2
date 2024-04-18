function xSet = fitS0andT2s(xSet,nneg,force,verbose)
if ~exist('nneg','var');    nneg    = []; end
if ~exist('force','var');   force   = []; end
if ~exist('verbose','var'); verbose = []; end
if isempty(nneg);           nneg    = 0 ; end
if isempty(force);          force   = 0 ; end
if isempty(verbose);        verbose = 0 ; end


for S = 1:length(xSet)
    if xSet{S}.nEcho>2
        disp(['Fitting T2* to ' xSet{S}.label])
        % Pick the most averaged data
        fieldList = {'fAvCatAv' 'fAv' 'f'};
        fieldList = fieldList(ismember(fieldList,fields(xSet{S}.finalFiles)));

        % Get file names
        nii = xSet{S}.finalFiles.(fieldList{1});
        jsn = xSet{S}.initFiles.fOrig(1,:); jsn = replace(jsn,'.nii.gz','.json');

        fOut = xSet{S}.finalFiles.(fieldList{1}){1};
        fOut = strsplit(fOut,filesep);
        dOut = strjoin(fOut(1:end-1),filesep);
        dOut = strsplit(dOut,'_'); tmp = strsplit(dOut{contains(dOut,'echo-')},'-');
        tmp{2} = 'R'; dOut{contains(dOut,'echo-')} = strjoin(tmp,'-'); dOut = strjoin(dOut,'_');
        if ~exist(dOut,'dir'); mkdir(dOut); end

        fOut = fOut{end};
        fOutS0  = replace(replace(fOut,'_volTs.nii.gz','_S0.nii.gz') ,'_vol.nii.gz','_S0.nii.gz' ); fOutS0  = fullfile(dOut,fOutS0);
        fOutT2s = replace(replace(fOut,'_volTs.nii.gz','_T2s.nii.gz'),'_vol.nii.gz','_T2s.nii.gz'); fOutT2s = fullfile(dOut,fOutT2s);


        if force || ~exist(fOutS0,'file') || ~exist(fOutT2s,'file')
            
            % Get echo time
            TE = cell(size(jsn));
            for i = 1:length(TE)
                [~,TE{i}] = system(['jq ''.EchoTime'' ' replace(jsn{i},'.nii.gz','.json')]); TE{i} = str2num(replace(replace(TE{i},'[0;39m',''),['[0m' newline],''));
            end
            TE = [TE{:}]';

            % Fit
            mriS0  = MRIread(nii{1},1);
            mriT2s = mriS0;
            imME = cell(size(nii));
            for E = 1:size(nii,2)
                imME{E} = MRIread(nii{E});
                imME{E} = imME{E}.vol;
            end
            imME = cat(4,imME{:});
            [mriT2s.vol,mriS0.vol] = calc_t2s_vol(imME, TE*1000, nneg, verbose);
            disp(' writting')
            disp(['  T2*: ' fOutT2s])
            MRIwrite(mriT2s,fOutT2s);
            disp(['  S0 : ' fOutS0])
            MRIwrite(mriS0,fOutS0);
            disp(' done')
        else
            disp(['  T2*: ' fOutT2s])
            disp(['  S0 : ' fOutS0])
            disp(' already done, skipping')
        end
        xSet{S}.finalFiles.([fieldList{1} 'S0']) = fOutS0;
        xSet{S}.finalFiles.([fieldList{1} 'T2s']) = fOutT2s;
    end
end
end



%
%
% %% S0 and T2*
% TE = fullfile({avMapSet.files.folder},{avMapSet.files.name})';
% for i = 1:length(TE)
%     [~,TE{i}] = system(['jq ''.EchoTime'' ' replace(TE{i},'.nii.gz','.json')]); TE{i} = str2num(replace(replace(TE{i},'[0;39m',''),['[0m' newline],''));
% end
% TE = cell2mat(TE);
%
% nneg = 0;
% nii = avMapSet.finalFiles.fAvEchoCat;
% avMapSet.finalFiles.fAvEchoS0 = cell(size(nii,1),1,1);
% avMapSet.finalFiles.fAvEchoT2s = cell(size(nii,1),1,1);
% for i = 1:size(nii,1)
%     fIn = nii{i};
%     avS0 = replace(fIn,'_vol.nii.gz','_S0.nii.gz');
%     avT2s = replace(fIn,'_vol.nii.gz','_T2s.nii.gz');
%     if force || ~exist(avS0,'file') || ~exist(avT2s,'file')
%         av = MRIread(fIn);
%         [t2star,S0] = calc_t2s_vol(av.vol(:,:,:,1:end), TE*1000, nneg, verbose); % last echo is actually the rms average
%         av.vol = S0; if force || ~exist(avS0,'file'); MRIwrite(av,avS0); end
%         av.vol = t2star; if force || ~exist(avT2s,'file'); MRIwrite(av,avT2s); end
%     end
%     avMapSet.finalFiles.fAvEchoS0{i} = avS0;
%     avMapSet.finalFiles.fAvEchoT2s{i} = avT2s;
% end
%


function [t2star,S0] = calc_t2s_vol(mef, te, non_neg,verbose)
% mef: multiecho flash volume MxNxPxE, where E = #echoes
% te: Echo times [Ex1] (ms)
% non_neg (optional): set non_neg=1 if you want to use non negative least
% squares. (default=0).
%
% Shared by Divya Varadarajan, 2022-11-04
% Minor modification (verbose option) by Sebastien Proulx, 2022-11-07

if nargin <4
    verbose=0;
end
if nargin <3
    non_neg=0;
end
sz = size(mef);
if verbose
    disp([num2str(sz(3)) ' slices detected.'])
    disp(['Calculating T2 star...'])
    if non_neg == 0
        disp('Using least squares.')
    else
        disp('Using non-negative least squares. Slow calculation - it can take a few minutes)')
    end
end
for nslice = 1:sz(3)
    % log
    lmef = reshape(log(squeeze(double(mef(:,:,nslice,:)))),[],length(te))';

    % Design matrix
    A = [ones(length(te),1) -te];

    if non_neg == 0
        % Least squares solution
        xout = A\lmef;
    else
        %non negative least squares
        for nv = 1:size(lmef,2)
            xout(:,nv) = lsqnonneg(A,lmef(:,nv));
        end
    end

    S0(:,:,nslice) = reshape(exp(xout(1,:)),sz(1:2));
    r2star(:,:,nslice) = reshape(xout(2,:),sz(1:2));
end

t2star = 1./r2star; t2star(r2star==0)=0;
% t2star(t2star<0)=0;
% t2star(t2star>1000)=0;

if verbose
    disp(['Calculation complete.'])
end
end
