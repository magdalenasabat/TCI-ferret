function [onset_in_source] = myGetOnsetInsource(name)
% it returns the onset of the sound in some source stimuli based on the
% indicated duration and part. 
[~,duration] =  myGetDuration(name);

% onset in source
% find where last _partX_ occurs to get the 'core name'  
expression = "_part";
[~,endIndex] = cellfun(@(x) regexp(x,expression), name,'UniformOutput', false);
part_in_source = str2num(name{1}(endIndex{1}(end)+1:end-4));

onset_in_source = (part_in_source-1) * duration;
        
        
end

