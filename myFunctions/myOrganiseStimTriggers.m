function [order1Idx, order2Idx] = myOrganiseStimTriggers(stimuli_names, triggers, dur)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finds patterns in presentation of stimuli of given duration and returns the indexes  
if length(stimuli_names) ~= length(triggers)
    error('Unequeal number of stimuli names and triggers')
    
end

durStimuli = myGetDuration(stimuli_names);
durStimuli = cellfun(@(x) str2num(x)/1000, durStimuli);
[durIdx,~, dur]  = grp2idx(durStimuli);
[stimIdx,namesTCI] = grp2idx(stimuli_names);

pattern1 = stimIdx(1:length(namesTCI));
pattern2 = stimIdx(length(namesTCI)+1:(length(namesTCI)+length(namesTCI)));

pattern1onsets = strfind(stimIdx',pattern1');
pattern2onsets = strfind(stimIdx',pattern2');

sz = min(length(pattern1onsets),length(pattern2onsets));

pattern1onsets=pattern1onsets(1:sz);
pattern2onsets=pattern2onsets(1:sz);

A = pattern1onsets+pattern1-1;
A = reshape(A,[],1) ;

B = pattern2onsets+pattern1-1;
B = reshape(B,[],1) ;

order1 = [];
order2 = [];
order1Idx = [];
order2Idx = [];

indexesDur = find(durIdx == j); %get indexes of all stim of a given dur

x1 = length(intersect(pattern1,indexesDur));
x2 = length(intersect(pattern2,indexesDur));

id1=intersect(indexesDur, A);
id2=intersect(indexesDur, B);

order1Idx = triggers(id1);
order1Idx = reshape(order1Idx, [x1, length(order1Idx)/x1]);

order2Idx = triggers(id2);
order2Idx = reshape(order2Idx, [x2, length(order2Idx)/x2]);

%         order1 = NaN(size(order1Idx,1), size(order1Idx,2) , length(timebins));
%         order2 = NaN(size(order2Idx,1), size(order2Idx,2) , length(timebins));

%for tests
%         afterT = floor(dur(j)*44100);
order1 = NaN(size(order1Idx,1), size(order1Idx,2) , afterT/10);
order2 = NaN(size(order2Idx,1), size(order2Idx,2) , afterT/10);

v = stimuli_names(id1);
v(2,:) = stimuli_names(id2);
v = v(:,1:x1);
v=[v(1,:) v(2,:)];
cool = grp2idx(v);
cool2 = cool(x1+1:end);
back=order2Idx;

for k=1:size(order2Idx,2)
    for m=1:size(cool2,1)
        order2Idx(m,k) = back(find(cool2==m),k);
    end         
end
        

end