function onset_target_in_stimuli = myGetOnsetInTarget(stimuli,target)
% it returns the onset of the sound in some target stimuli based on the
% indicated duration and part. 
[~,duration] =  myGetDuration(target);

% onset in tareget
% find where last _partX_ occurs to get the 'core name'  
expression = "_part";
[~,endIndex] = cellfun(@(x) regexp(x,expression), target,'UniformOutput', false);
part_in_target = str2num(target{1}(endIndex{1}(end)+1:end-4));

target_onset_in_source = (part_in_target-1) * duration;

% stimuli
% find where last _partX_ occurs to get the 'core name' 
[~,duration] =  myGetDuration(stimuli);
expression = "_part";
[~,endIndex] = cellfun(@(x) regexp(x,expression), stimuli,'UniformOutput', false);
part_in_stimuli= str2num(stimuli{1}(endIndex{1}(end)+1:end-4));

stimuli_onset_in_source = (part_in_stimuli-1) * duration;

onset_target_in_stimuli = target_onset_in_source - stimuli_onset_in_source;
end

