function categories = myGetCategory(stimuli_names)
%returns the category of the target : speech, music or ferrets
%
% stimuli_names must be a cell array of char vecs
%
% author: Magdalena Sabat email:magdalena.sabat@ens.fr
% created: 27/10/2023


categories = cellfun(@checkCategory,  stimuli_names,'UniformOutput',false);

function category = checkCategory(stim_name)

    speech = 'speech';
    music =  'music';
    ferrets =  'ferrets';

    if strcmp(stim_name(1:length(speech)), speech)
        category = speech;
    elseif strcmp(stim_name(1:length(music)), music)
        category = music;
    elseif strcmp(stim_name(1:length(ferrets)), ferrets)
        category = ferrets;
    else
        fprintf('no cat')
        category = '';
    end

end

end