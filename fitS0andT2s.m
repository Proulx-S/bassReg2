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
        candidate = {'fSes' 'fAvCatAv' 'fAv' 'f'};
        fieldList = {};

        candidateTmp = candidate(ismember(candidate,fields(xSet{S}.finalFiles)));
        candidateTmp = candidateTmp(1);
        for iii = 1:length(candidateTmp)
            fieldList{end+1} = ['finalFiles.' candidateTmp{iii}];
        end
        if isfield(xSet{S},'statFiles')
            if isfield(xSet{S}.statFiles,'modelFree')
                if isfield(xSet{S}.statFiles.modelFree,'base')
                    candidateTmp = candidate(ismember(candidate,fields(xSet{S}.statFiles.modelFree.base)));
                    candidateTmp = candidateTmp(1);
                    for iii = 1:length(candidateTmp)
                        fieldList{end+1} = ['statFiles.modelFree.base.' candidateTmp{iii}];
                    end
                end
                if isfield(xSet{S}.statFiles.modelFree,'respOnBase')
                    candidateTmp = candidate(ismember(candidate,fields(xSet{S}.statFiles.modelFree.respOnBase)));
                    candidateTmp = candidateTmp(1);
                    for iii = 1:length(candidateTmp)
                        fieldList{end+1} = ['statFiles.modelFree.respOnBase.' candidateTmp{iii}];
                    end
                end
            end

            if isfield(xSet{S}.statFiles,'customModel')
                if isfield(xSet{S}.statFiles.customModel,'base')
                    candidateTmp = candidate(ismember(candidate,fields(xSet{S}.statFiles.customModel.base)));
                    candidateTmp = candidateTmp(1);
                    for iii = 1:length(candidateTmp)
                        fieldList{end+1} = ['statFiles.customModel.base.' candidateTmp{iii}];
                    end
                end
                if isfield(xSet{S}.statFiles.customModel,'respAmpOnBase')
                    candidateTmp = candidate(ismember(candidate,fields(xSet{S}.statFiles.customModel.respAmpOnBase)));
                    candidateTmp = candidateTmp(1);
                    for iii = 1:length(candidateTmp)
                        fieldList{end+1} = ['statFiles.customModel.respAmpOnBase.' candidateTmp{iii}];
                    end
                end
            end
        end





        for ii = 1:length(fieldList)

            % Get file names
            eval(['nii = xSet{' num2str(S) '}.' fieldList{ii} ';']);
            jsn = xSet{S}.initFiles.fOrig(1,:); jsn = replace(jsn,'.nii.gz','.json');
            nRun = size(nii,1);
            if nRun>1; dbstack; error('code that'); end
            nEcho = size(nii,2);
            niiFit = nii;
            for E = 1:nEcho
                niiFit{E} = replace(nii{E},'.nii.gz','RelaxFit.nii.gz');
            end


            % fOut = xSet{S}.finalFiles.(fieldList{ii}){1};
            fOut = strsplit(nii{1},filesep);
            dOut = strjoin(fOut(1:end-1),filesep);
            dOut = strsplit(dOut,'_'); tmp = strsplit(dOut{contains(dOut,'echo-')},'-');
            tmp{2} = 'R'; dOut{contains(dOut,'echo-')} = strjoin(tmp,'-'); dOut = strjoin(dOut,'_');
            if ~exist(dOut,'dir'); mkdir(dOut); end


            fOut = fOut{end};
            fOutS0  = replace(replace(replace(fOut,'_volTs.nii.gz','_S0.nii.gz') ,'_vol.nii.gz','_S0.nii.gz' ),'_respOnBase.nii.gz','_respOnBaseS0.nii.gz' ); fOutS0  = fullfile(dOut,fOutS0);
            fOutS0Perc = replace(replace(replace(fOut,'_volTs.nii.gz','_S0Perc.nii.gz'),'_vol.nii.gz','_S0Perc.nii.gz'),'_respOnBase.nii.gz','_respOnBaseS0Perc.nii.gz'); fOutS0Perc = fullfile(dOut,fOutS0Perc);
            fOutT2s = replace(replace(replace(fOut,'_volTs.nii.gz','_T2s.nii.gz'),'_vol.nii.gz','_T2s.nii.gz'),'_respOnBase.nii.gz','_respOnBaseT2s.nii.gz'); fOutT2s = fullfile(dOut,fOutT2s);
            fOutR2s = replace(replace(replace(fOut,'_volTs.nii.gz','_R2s.nii.gz'),'_vol.nii.gz','_R2s.nii.gz'),'_respOnBase.nii.gz','_respOnBaseR2s.nii.gz'); fOutR2s = fullfile(dOut,fOutR2s);
            fOutR2sPerc = replace(replace(replace(fOut,'_volTs.nii.gz','_R2sPerc.nii.gz'),'_vol.nii.gz','_R2sPerc.nii.gz'),'_respOnBase.nii.gz','_respOnBaseR2sPerc.nii.gz'); fOutR2sPerc = fullfile(dOut,fOutR2sPerc);


            niiFitFlag = false; for E = 1:nEcho; niiFitFlag = niiFitFlag | ~exist(niiFit{E},'file'); end
            if force || ~exist(fOutS0,'file') || ~exist(fOutS0Perc,'file') || ~exist(fOutT2s,'file') || ~exist(fOutR2s,'file') || ~exist(fOutR2sPerc,'file')  || niiFitFlag

                % Get echo time
                TE = cell(size(jsn));
                for i = 1:length(TE)
                    [~,TE{i}] = system(['jq ''.EchoTime'' ' replace(jsn{i},'.nii.gz','.json')]); TE{i} = str2num(replace(replace(TE{i},'[0;39m',''),['[0m' newline],''));
                end
                TE = [TE{:}]';

                % Fit
                mriS0  = MRIread(nii{1},1);
                mriT2s = mriS0;
                mriR2s = mriS0;
                imMEmri = cell(size(nii));
                imME = cell(size(nii));
                for E = 1:size(nii,2)
                    imMEmri{E} = MRIread(nii{E});
                    imME{E} = imMEmri{E}.vol;
                end
                imME = permute(cat(5,imME{:}),[1 2 3 5 4]);


                mriT2s.vol = zeros(size(imME,[1:3 5]));
                mriR2s.vol = zeros(size(imME,[1:3 5]));
                mriS0.vol = zeros(size(imME,[1:3 5]));
                imMEfit = zeros(size(imME));
                for t = 1:size(imME,5)
                    [mriT2s.vol(:,:,:,t),mriR2s.vol(:,:,:,t),mriS0.vol(:,:,:,t),imMEfit(:,:,:,:,t)] = calc_t2s_vol(imME(:,:,:,:,t), TE*1000, nneg, verbose);
                end
                mriT2s.vol(isnan(mriT2s.vol)) = 0;
                mriS0.vol(isnan(mriS0.vol)) = 0;
                mriR2s.vol(isnan(mriR2s.vol)) = 0;

                % Normalize
                if size(mriR2s.vol,4)>1
                    mriR2sPerc = mriR2s;
                    mriR2sPerc.vol = (mriR2s.vol-mriR2s.vol(:,:,:,1))./mriR2s.vol(:,:,:,1)*100;
                    mriS0Perc = mriS0;
                    mriS0Perc.vol = (mriS0.vol-mriS0.vol(:,:,:,1))./mriS0.vol(:,:,:,1)*100;
                end

                % Write
                disp(' writting')
                disp(['  S0         : ' fOutS0])
                MRIwrite(mriS0,fOutS0);
                if size(mriS0.vol,4)>1
                    disp(['  S0 %change: ' fOutS0Perc])
                    MRIwrite(mriS0Perc,fOutS0Perc);
                end
                disp(['  T2*        : ' fOutT2s])
                MRIwrite(mriT2s,fOutT2s);
                disp(['  R2*        : ' fOutR2s])
                MRIwrite(mriR2s,fOutR2s);
                if size(mriR2s.vol,4)>1
                    disp(['  R2* %change: ' fOutR2sPerc])
                    MRIwrite(mriR2sPerc,fOutR2sPerc);
                end
                disp(['  fit        : ' newline '   ' strjoin(niiFit,[newline '   '])])
                for E = 1:nEcho
                    imMEmri{E}.vol = imMEfit(:,:,:,E);
                    MRIwrite(imMEmri{E},niiFit{E});
                end
                disp(' done')
            else
                disp(['  S0 : ' fOutS0])
                disp(['  T2*: ' fOutT2s])
                disp(['  R2*: ' fOutR2s])
                disp(['  fit: ' niiFit])
                disp(' already done, skipping')
            end

            eval(['xSet{' num2str(S) '}.' fieldList{ii} 'S0 = {fOutS0};']);
            eval(['xSet{' num2str(S) '}.' fieldList{ii} 'S0Perc = {fOutS0Perc};']);
            eval(['xSet{' num2str(S) '}.' fieldList{ii} 'T2s = {fOutT2s};']);
            eval(['xSet{' num2str(S) '}.' fieldList{ii} 'R2s = {fOutR2s};']);
            eval(['xSet{' num2str(S) '}.' fieldList{ii} 'R2sPerc = {fOutR2sPerc};']);
            for E = 1:nEcho
                eval(['xSet{' num2str(S) '}.' fieldList{ii} 'ExpFit(' num2str(E) ') = niiFit(' num2str(E) ');']);
            end
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


function [t2star,r2star,S0,mefFit] = calc_t2s_vol(mef, te, non_neg,verbose)
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
sz = size(mef,1:5);
if verbose
    disp([num2str(sz(3)) ' slices detected.'])
    disp(['Calculating T2 star...'])
    if non_neg == 0
        disp('Using least squares.')
    else
        disp('Using non-negative least squares. Slow calculation - it can take a few minutes)')
    end
end
S0 = zeros(sz(1:3));
r2star = zeros(sz(1:3));
mefFit = zeros(sz(1:4));
for nslice = 1:sz(3)
    % log
    % lmef = reshape(log(squeeze(double(mef(:,:,nslice,:)))),[],length(te))';
    lmef = reshape(log(double(mef(:,:,nslice,:))),[],length(te))';

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
    mefFit(:,:,nslice,:) = reshape(exp(A*xout),[sz(1:2) 1 sz(4)]);
end

t2star = 1./r2star;
% t2star(t2star<0)=0;
% t2star(t2star>1000)=0;

if verbose
    disp(['Calculation complete.'])
end

