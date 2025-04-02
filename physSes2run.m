function physRun = physSes2run(physSes)

if nargin==0
    physRun.vec = [];
    physRun.t0 = [];
    physRun.chanLabel = {};
    physRun.Fs = [];
    physRun.mri = [];
    physRun.info = [];
    physRun.vecInfo = 'time/freq x chan x taper/mode x run';
    physRun.isempty = 1;
    return
end

s = physSes.info.segmentInd;
tPhys0 = sum([hour(physSes.blocktimes(s))*60*60 minute(physSes.blocktimes(s))*60 second(physSes.blocktimes(s))]);
chanList = {'cardiac' 'resp' 'trigger'}';
for c = 1:length(chanList)
    chanInd = ismember(physSes.titles1,chanList{c});
    if any(chanInd)
        data = physSes.data(physSes.datastart(chanInd,s):physSes.dataend(chanInd,s));
        n = size(data,2);
        Fs = physSes.samplerate(chanInd,s);
        t = linspace(0,(n-1)/Fs,n);
        for r = 1:size(physSes.mriRunTimes,1)
            % extract data
            tMri = physSes.mriRunTimes(r,:) - tPhys0;
            physRun(r,1).vec(:,c) = data(t>=tMri(1) & t<tMri(2));
            physRun(r,1).t0(1,c)  = t(find(t>=tMri(1),1,'first'));

            % extract other info
            physRun(r,1).chanLabel(1,c) = chanList(c);
            physRun(r,1).Fs(1,c) = physSes.samplerate(c,s);
            if c==1
                physRun(r,1).mri  = physSes.info.mri(r);
                physRun(r,1).info = rmfield(physSes.info,'mri');
            end
        end
    else
        %% fill non-recorded channels
        for r = 1:size(physSes.mriRunTimes,1)
            % extract data
            tMri = physSes.mriRunTimes(r,:) - tPhys0;
            physRun(r,1).vec(:,c) = nan(nnz(t>=tMri(1) & t<tMri(2)),1);
            physRun(r,1).t0(1,c)  = nan;

            % extract other info
            physRun(r,1).chanLabel(1,c) = chanList(c);
            physRun(r,1).Fs(1,c) = nan;
            if c==1
                physRun(r,1).mri = physSes.info.mri(r);
                physRun(r,1).info = rmfield(physSes.info,'mri');
            end
        end
    end
end
[physRun.vecInfo] = deal(strjoin({'time/freq' 'chan' 'taper/mode' 'run'},' x '));
[physRun.isempty] = deal(0);

