%% PREPROCESS THE TRIGGERS
cd ~/Desktop/
load('forYves.mat');
d = forYves;

%%
d.triggerFiles{1} = d.triggerFiles{1}';
gt = cell2mat(d.triggerFiles);
timings = find(gt(2:end)>gt(1:(end-1)))+1;
firsttrigAudio = timings(1);
lasttrigAudio = timings(end);
timings = timings/44100;
timings = timings-timings(1);

T = 2; % threshold croosing for detecting triggers
anchan = audioTriggerChannel;
antimings = find(anchan(2:end)>T &...
    anchan(1:(end-1))<T)+1;
antimingspl = antimings;
firsttrigAn = antimings(1);
lasttrigAn = antimings(end);
antimings = antimings/srEphys;
% fill in the cut trigs
tooshortint = find(diff(antimings)<0.008);
for intnum = 1:length(tooshortint)
    anchan(antimingspl(tooshortint(intnum)):antimingspl(tooshortint(intnum)+1))=max(anchan);
end
srAudio=44100

antimings = find(anchan(2:end)>T &...
    anchan(1:(end-1))<T)+1;
antimingspl = antimings;
firsttrigAn = antimings(1);
lasttrigAn = antimings(end);
antimings = antimings/srEphys;
%% DTW
% INITIALIZE
dur = 2.5; % duration of the processed segment (in s)
% segment duration
splnbephys = round(srEphys*dur);
splnbaudio = round(srAudio*dur);
Cephys = 0; Caudio = 0; Clst = [0;0];
jumpcount = 0; modC = 0;
while (Caudio+firsttrigAudio)<lasttrigAudio
    % check if we are at the very end of the open ephys trigs
    splnbaudio = min([splnbaudio,lasttrigAudio-(Caudio+firsttrigAudio)]);
    jumpcount = jumpcount +1;
    % cut and downsample segments in open ephys (x) and audio-matlab (y) data
    x = double(anchan(firsttrigAn+Cephys+(1:splnbephys)));
    x = x(1:10:end);
    x = x/std(x); x = x-mean(x);
    y = gt(firsttrigAudio+Caudio+(1:splnbaudio));
    y = y(1:15:end);
    y = y/std(y); y = y-mean(y);
    % align with DTW
    [di(jumpcount),ix,iy] = dtw(x,y);
    % lasttrig = find(x(1:end-1)<1.5&x(2:end)>1.5,20,'last'); lasttrig = lasttrig(1);
    if splnbaudio==round(srAudio*dur)
        % take the last adjusted trig from open ephys triggers (the next
        % segment will start from it) and adjust for the downsampling
        lasttrig = find(x(1:end-1)<1.5&x(2:end)>1.5,25,'first'); lasttrig = lasttrig(end);
        Cephys = Cephys+(ix(find(ix == lasttrig,1,'first'))-1)*10+1;
        Caudio = Caudio+(iy(find(ix == lasttrig,1,'first'))-1)*15+1; 
    else
        % when at the end of the openephys triggers
        lasttrig = find(x(1:end-1)<1.5&x(2:end)>1.5,1,'last');
        Caudio = Caudio+(length(x)-1)*15+1;
    end
    audiotrigs = find(y(1:end-5)<0&y(2:end-4)>0);
    at = 1;
    ephystrigs = [];
    while audiotrigs(at)<=(iy(find(ix == lasttrig)))
        % make the correspondance between openphys and audio triggers one
        % by one
        audiotrig = audiotrigs(at);
        ephystrigs = [ephystrigs [ (iy(find(iy == audiotrig,1,'first'))-1)*15+1;...
            (ix(find(iy == audiotrig,1,'first'))-1 )*10+1 ]  ];
        at = at+1;
        if at < length(audiotrigs); break; end
    end
    % zad
    Clst = [ Clst Clst(:,end)+ephystrigs(:,2:end) ];
    % plot the result of the DTW every 50 segments
    if mod(jumpcount,50)==0 || splnbaudio~=round(srAudio*dur)
        modC = modC+1;
        figure('name',num2str((Caudio+firsttrigAudio+splnbaudio+200000)));
        dtw(x,y);
    end
end

MatchedTrigs = [Clst+[firsttrigAudio;firsttrigAn]];


