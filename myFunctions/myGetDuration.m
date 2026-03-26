function [durs_str, durs] = myGetDuration(stimuliNames)

% find where _XXXms_ 
expression = "_[0-9]+\.*[0-9]+ms_";
[startIndex,endIndex] = cellfun(@(x) regexp(x,expression), stimuliNames,'UniformOutput', false);

%get the durations
durs_str = cellfun(@(x,startIdx,endIdx) x(startIdx(end)+1:endIdx(end)-3), stimuliNames,startIndex,endIndex,'UniformOutput', false);

%return as a number oooor nt
durs = cellfun(@str2num, durs_str);



end