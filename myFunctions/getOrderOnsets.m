function [pattern1onsets, pattern2onsets] = getOrderOnsets(stimuli_names,triggers)

% [x,y]  = grp2idx(stimuli_names);
if strcmp(stimuli_names{1}, 'silence - buffer init') | contains(stimuli_names{1}, 'target presentation')

    stimuli_names(1) = [];
    triggers(1)=[];


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
end
        