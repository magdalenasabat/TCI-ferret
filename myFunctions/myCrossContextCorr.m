function L = myCrossContextCorr(SAR, sourceFolder, indiv,sessN, taskType, trialType,signalType, chann_numbers,lag_t, binWidthMs, args )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Calculates same and cross context correlation for
% random_random and random_natural comparisons needed for TCI
%
% IN:
% SAR struct with length == length(chann_numbers)data
% SAR.random_random_context - size(number of durations x [])
% SAR.random_random_context.order1 - stimuli x rep x time for each duration
% SAR.random_random_context.order2 - stimuli x rep x time for each duration
% SAR.random_random_context.durMs - dur in ms 
% SAR.random_random_context.durN - N of stimuli
% SAR.random_random_context.n_total_segs - 
% SAR.random_random_context.surrounding_durations/onsets - infor about what
% segments surrounded
% SAR.random_natural_context - srim x rep x time  
%
% lag_t - length(time) in sec
% binWidthMs - scalar  width of the bin in miliseconds
%
% indiv - char subject name
% sessN  - scalar sessN 
% taskType  string 
% trialType - string
% signalType
% chann_number
% create_dir 
%
% OUT:
% CCC.random_random 
% CCC.random_natural
% 
% magdalena.sabat@ens.psl.eu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments 
    SAR struct 
    sourceFolder char
    indiv char
    sessN (1,1) 
    taskType char
    trialType char
    signalType char
    chann_numbers {mustBeDoubleorCell}
    lag_t double
    binWidthMs double
    args.smoothing logical = false  
    args.compute_natural_context logical = true
    args.compute_rep_splits logical = true
    args.create_dir logical = true
end

sr=1000/binWidthMs;
unique_segs = [SAR(1).random_random_context.durMs];
fprintf('calculating correlations  \n')

L.random_random_context = myPrepareL( sourceFolder, indiv,sessN, taskType, trialType, signalType, chann_numbers,unique_segs,sr,lag_t,[SAR(1).random_random_context.n_total_segs]);
L.output_directory =  L.random_random_context.output_directory;
L.random_random_context.figure_directory = strcat([L.random_random_context.figure_directory filesep 'random_random_context']);
L.random_random_context.output_directory = strcat([L.random_random_context.output_directory filesep 'random_random_context']);

if ~isfolder(L.random_random_context.output_directory)
    mkdir(L.random_random_context.figure_directory)  
    mkdir(L.random_random_context.output_directory)  
end  

if args.compute_rep_splits
    L.random_random_context.odd = myPrepareL( sourceFolder, indiv,sessN, taskType, trialType, signalType, chann_numbers,unique_segs,sr,lag_t,[SAR(1).random_random_context.n_total_segs]);
    L.random_random_context.odd.figure_directory = strcat([L.random_random_context.odd.figure_directory filesep 'random_random_context' filesep 'reps-odd']);
    L.random_random_context.odd.output_directory = strcat([L.random_random_context.odd.output_directory filesep 'random_random_context' filesep 'reps-odd']);

    L.random_random_context.even = myPrepareL( sourceFolder, indiv,sessN, taskType, trialType, signalType, chann_numbers,unique_segs,sr,lag_t,[SAR(1).random_random_context.n_total_segs]);
    L.random_random_context.even.figure_directory = strcat([L.random_random_context.even.figure_directory filesep 'random_random_context' filesep 'reps-even']);
    L.random_random_context.even.output_directory = strcat([L.random_random_context.even.output_directory filesep 'random_random_context' filesep 'reps-even']);
    if args.create_dir & ~isfolder(L.random_random_context.odd.output_directory)
        mkdir(L.random_random_context.odd.output_directory)  
        mkdir(L.random_random_context.odd.figure_directory)
        
        mkdir(L.random_random_context.even.figure_directory)
        mkdir(L.random_random_context.even.output_directory)
    end 
end

if isfield(SAR, 'random_natural_context')
    fields = {'random_random_context','random_natural_context'};
    unique_segs_natural = [SAR(1).random_natural_context.durMs];
    unique_segs_natural = unique_segs_natural(cellfun(@(x) ~isempty(x),{SAR(1).random_natural_context.order1}));
    
    n_total_segs_natural = [SAR(1).random_natural_context.n_total_segs];
    n_total_segs_natural = n_total_segs_natural(cellfun(@(x) ~isempty(x),{SAR(1).random_natural_context.order1}));
    
    L.random_natural_context = myPrepareL( sourceFolder, indiv,sessN, taskType, trialType, signalType, chann_numbers,unique_segs_natural,sr,lag_t,n_total_segs_natural);
    L.random_natural_context.figure_directory = strcat([L.random_natural_context.figure_directory filesep 'random_natural_context']);
    L.random_natural_context.output_directory = strcat([L.random_natural_context.output_directory filesep 'random_natural_context']);
    
    if args.create_dir &  ~isfolder(L.random_natural_context.output_directory)
        mkdir(L.random_natural_context.figure_directory)  
        mkdir(L.random_natural_context.output_directory)  
    end
    if args.compute_rep_splits
        L.random_natural_context.odd = myPrepareL( sourceFolder, indiv,sessN, taskType, trialType, signalType, chann_numbers,unique_segs_natural,sr,lag_t,n_total_segs_natural);    
        L.random_natural_context.odd.figure_directory = strcat([L.random_natural_context.odd.figure_directory filesep 'random_natural_context' filesep 'reps-odd']);
        L.random_natural_context.odd.output_directory = strcat([L.random_natural_context.odd.output_directory filesep 'random_natural_context' filesep 'reps-odd']);
        
        L.random_natural_context.even = myPrepareL( sourceFolder, indiv,sessN, taskType, trialType, signalType, chann_numbers,unique_segs_natural,sr,lag_t,n_total_segs_natural)  ;
        L.random_natural_context.even.figure_directory = strcat([L.random_natural_context.even.figure_directory filesep 'random_natural_context' filesep 'reps-even']);
        L.random_natural_context.even.output_directory = strcat([L.random_natural_context.even.output_directory filesep 'random_natural_context' filesep 'reps-even']);
        if args.create_dir &  ~isfolder(L.random_natural_context.odd.output_directory)
            mkdir(L.random_natural_context.odd.output_directory)  
            mkdir(L.random_natural_context.odd.figure_directory)

            mkdir(L.random_natural_context.even.figure_directory)
            mkdir(L.random_natural_context.even.output_directory)
        end 
    end
else
    fields = {'random_random_context'};
end


reps = {'odd' , 'even'};
n_rep = size(SAR(1).random_random_context(1).order1,2);
for chann = 1:length(chann_numbers)
    for dur=1:size(unique_segs,2) % for each duration

        % all reps random_random & random_natural
        for field=fields
            if isempty(SAR(chann).(field{:})(dur).order1);continue;end
            o1 = SAR(chann).(field{:})(dur).order1; 
            o2 = SAR(chann).(field{:})(dur).order2; 
            
            tm = logical(tril(ones(size(o1 ,2)),-1));

            x = []; y = []; z = [];
            for t=1:size(o1  ,3)
                x(t) = mean(corr(squeeze(o1 (:,:,t)),squeeze(o2 (:,:,t))),'all','omitnan');

                r = corr(squeeze(o1 (:,:,t)));
                y(t) = mean(r(tm),'all','omitnan');

                r = corr(squeeze(o2 (:,:,t)));
                z(t) = mean(r(tm),'all','omitnan');
            end

            L.(field{:}).diff_context(:,dur,chann)  = x;
                
            L.(field{:}).same_context(:,dur,chann) = (y/2 + z/2);    
            L.(field{:}).same_context_err(:,dur,chann) = (y/2 - z/2).^2;
            
            L.(field{:}).same_context_err(:,dur,chann) = (y/2 - z/2).^2;
            if strcmp(field{:},'random_random_context')
                if isfield(SAR(chann).(field{:})(dur),'surrounding_onsets')
                    L.(field{:}).surrounding_segments{dur} = [SAR(chann).(field{:})(dur).surrounding_onsets; SAR(chann).(field{:})(dur).surrounding_durations];
                end
            end
        end

        % odd & even reps random_random & random_natural
        for field=fields
            for rep=1:length(reps) 
                if isempty(SAR(chann).(field{:})(dur).order1);continue;end
                o1 = SAR(chann).(field{:})(dur).order1(:,rep:2:n_rep,:); 
                o2 = SAR(chann).(field{:})(dur).order2(:,rep:2:n_rep,:); 
                tm = logical(tril(ones(size(o1 ,2)),-1));

                x = []; y = [];z=[];
                for t=1:size(o1  ,3)
                    x(t) = mean(corr(squeeze(o1 (:,:,t)),squeeze(o2 (:,:,t))),'all','omitnan');

                    r = corr(squeeze(o1 (:,:,t)));
                    y(t) = mean(r(tm),'all','omitnan');

                    r = corr(squeeze(o2 (:,:,t)));
                    z(t) = mean(r(tm),'all','omitnan');
                end
                L.(field{:}).(reps{rep}).diff_context(:,dur,chann)   = x;
                
                L.(field{:}).(reps{rep}).same_context(:,dur,chann)  = (y/2 + z/2);   
                L.(field{:}).(reps{rep}).same_context_err(:,dur,chann)  = (y/2 - z/2).^2;
            end
        end


    end
end



end



function L = myPrepareL( sourceFolder,indiv,sessN, taskType, trialType, signalType, chann_numbers,unique_segs,sr,lag_t,n_total_segs)
    L = [];
    L.param_string = ['boundary-any_lag_win-0-1'];
    if isa(chann_numbers, 'double')
        L.channels = chann_numbers;
        L.chnames = cellstr(num2str(L.channels'));
    else
        L.channels = 1:length(chann_numbers);
        L.chnames = chann_numbers;
    end
    L.boundary = 'any';

    L.rampwin = 'hann';
    L.rampdur = 0.015625 ;

    L.indiv = indiv;
    L.sessN = sessN;
    L.taskType = taskType;
    L.trialType = trialType;
    L.signalType = signalType;
    
    L.unique_segs = unique_segs;
    L.sr = sr;
    L.lag_t = lag_t;
    L.n_total_segs = n_total_segs';
    L.figure_directory = strcat([sourceFolder filesep 'analysis' filesep 'figures' filesep ...
                'sub-' lower(indiv) filesep 'sess-' sprintf('%03d',sessN) filesep taskType filesep trialType filesep ...
                'TCI' filesep 'stim-aligned' filesep signalType filesep 'binWidth-' num2str(1000/sr) 'ms' ]);
             
    L.output_directory = strcat([sourceFolder filesep 'analysis' filesep 'data-work' filesep ...
                'sub-' lower(indiv) filesep 'sess-' sprintf('%03d',sessN) filesep taskType filesep trialType filesep ...
                'TCI' filesep 'stim-aligned' filesep signalType filesep 'binWidth-' num2str(1000/sr) 'ms'  ]);
             

end
% Custom validation function
function mustBeDoubleorCell(chann_numbers)

    if ~isa(chann_numbers, 'double') && ~isa(chann_numbers, 'cell')
        error('chann_numbers argument must be a double or a cell array.');
    end
end

