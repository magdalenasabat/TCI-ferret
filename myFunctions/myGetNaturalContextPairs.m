function [natural_pairs, onset_instim] = myGetNaturalContextPairs(names)
% For each name in names returns possible natural context names
%
% assumes my naming convention : 
% <category>-<stim_name>-<duration>-<part_number>
%
% -<duration>-<part_number> can be repeated if the sounds were processed
% more than one time (legacy)
%
% IN:
% names             - Nx1 cell array of strings
%
% OUT:
% natural_pairs     - Nx1 cell array of vectors with indexes of natural
% comprisons
%
% magdalena, 15/01/2024

durs_sec = cellfun(@(x) str2num(x)/1000, myGetDuration(names));
[durIdx,~] = grp2idx(durs_sec);
durations = unique(durs_sec);
% n_seg_durs = numel(durations);

%get the indexs of where the last indication of duration starts in all
%names
expression = '_\d+(\.\d+)?ms_part\d+\.wav';
[namesStartIndexs, ~] = cellfun(@(x) regexp(x,expression), names);




 %%

natural_pairs  = cell(1,length(names));
onset_instim  = cell(1,length(names));

for i=1:length(names) % for each stimuli
%     fprintf(num2str(i))
    
    % get the name and duration
    target=names(i);
    duration = durs_sec(i);
    
    % find where first _XXXms_ occurs to get the 'core name' ~ the source file name  
    expression = '_\d+(\.\d+)?ms';
    [startIndex, endIndex] =  regexp(target,expression);
    
    % get the name and duration of the source file
    source_name = target{1}(1:startIndex{1}(end)-1);
    source_duration = str2num(target{1}(startIndex{1}(1)+1:endIndex{1}(1)-3));
    
    % when we include stretched and compressed sounds the normal sound will
    % match to them (and we don't want that), we check that that's not the
    % case by checking whether the end of the matched name ends where the
    % last indication of duration starts 
    

    namesEndIndexs = zeros(length(names),1);
    for ix = 1:length(names)
        namesEndIndexs(ix) = max([strfind(names{ix}, source_name)+length(source_name)-1,0]);
    end

%     namesEndIndexs = cellfun(@(x) max([strfind(x, source_name)+length(source_name)-1,0]), names);
    
    % get the names that match, get only sounds that have longer duration
    idx_longer_names = find((namesEndIndexs == namesStartIndexs-1) & (durs_sec > duration )) ;
    longer_names = names(idx_longer_names)  ;
    longer_durations = durs_sec(idx_longer_names);
    
    % onset in source
    % find where last _partX_ occurs to get part in source
    expression = "_part";
    [~,endIndex] = cellfun(@(x) regexp(x,expression), target,'UniformOutput', false);
    part_in_source = str2num(target{1}(endIndex{1}(end)+1:end-4));

    name_onset_insource = myGetOnsetInsource(target)/1000; 
    longer_names_onset_insource = cellfun(@(x) myGetOnsetInsource({x})/1000,  longer_names);

    % source files might have onsets >= longest duration , overwrite with 0
    longer_names_onset_insource(longer_names_onset_insource >= max(durations)) = name_onset_insource;

    valid_longer_names = (longer_names_onset_insource < name_onset_insource) & ... onset before the begining of the target
                         ((longer_names_onset_insource + longer_durations) > (name_onset_insource + duration) ); % offset after the end of the target
    
    % check if there is no mixup with indexes
    assert(all(strcmp( names(idx_longer_names(valid_longer_names)),longer_names(valid_longer_names)) ))
    
    valid_pairs = longer_names(valid_longer_names);
    % calculate onset of target in pair
    valid_paris_onsit_in_stim = cellfun(@(x) myGetOnsetInTarget({x},target)/1000,  valid_pairs);
    
    natural_pairs(i) = {[idx_longer_names(valid_longer_names)]};
    onset_instim(i) = {valid_paris_onsit_in_stim } ; % shift by the ramp                    
end
   


end

