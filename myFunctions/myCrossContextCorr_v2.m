function L = myCrossContextCorr_v2(SAR, sourceFolder, indiv,sessN, taskType, trialType,signalType, chann_numbers,lag_t, binWidthMs, args )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Calculates same and cross context correlation for
% random_random and random_natural comparisons needed for TCI
% version that uses cross_context_helper from Sam
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
    args.simfunc char = 'corr'
end

% same as Sam's

make_samedur_comparisons = true; % this is to cheat on Sam's helper function
make_diffdur_comparisons = false;
diffdur_order_pairs = [];
diffdur_seg_pairs = [];
Y_embed_seg = [];

interleave_samedur = false;
interleave_diffdur = false;
n_orders=2;

switch args.simfunc
    case 'corr'
        simfunc = @nanfastcorr;
    case 'nse'
        simfunc = @nan_nse_match_columns;
    case 'mae'
        simfunc = @(a,b)nanmean(-abs(a-b),1);
    otherwise
        error('No matching similarity function')
end
    
tranweightfn = @(x)x;
trancorrfn = @(x)x;
samerep = false; 

reps = {'odd' , 'even'};
n_reps = size(SAR(1).random_random_context(1).order1,2);

%% Pairs of orders for same duration comparisons

% for same-duration comparisons, must always compare different orders
if make_samedur_comparisons
    samedur_order_pairs = [];
    for o1 = 1:n_orders
        for o2 = o1+1:n_orders
            new_pair = [o1, o2]';
            samedur_order_pairs = cat(2, samedur_order_pairs, new_pair);
        end
    end
    clear o1 o2 new_pair;
else
    samedur_order_pairs = [];
end
 %% Pairs of repetitions for same vs. different context comparisons
    
% pairs of repetitions to use for diff context comparisons
% when comparing different contexts/order
% we need to consider all possible pairs of repetitions
% since order1-rep1 vs. order2-rep2 and order1-rep2 vs. order2-rep1
% are both valid
% we can optionally also include pairs from the same
% repetition, though it is arguably more conservative to not do this
diffcontext_rep_pairs = [];
for r1 = 1:n_reps
    if samerep
        second_reps = 1:n_reps;
    else
        second_reps = [1:r1-1,r1+1:n_reps];
    end
    for r2 = second_reps
        new_pair = [r1, r2]';
        diffcontext_rep_pairs = cat(2, diffcontext_rep_pairs, new_pair);
    end
end
clear r1 r2 second_reps new_pair;

% for reliability w also need this for odd/even reps
diffcontext_rep_pairs_perrep = [];
for r1 = 1:floor(n_reps/2)
    if samerep
        second_reps = 1:floor(n_reps/2);
    else
        second_reps = [1:r1-1,r1+1:floor(n_reps/2)];
    end
    for r2 = second_reps
        new_pair = [r1, r2]';
        diffcontext_rep_pairs_perrep = cat(2, diffcontext_rep_pairs_perrep, new_pair);
    end
end

% pairs of repetitions to use for same context comparisons
% here we need to consider unique rep pairs obviously
% since we have exactly the same segment/context
samecontext_rep_pairs = [];
for r1 = 1:n_reps
    for r2 = r1+1:n_reps
        new_pair = [r1, r2]';
        samecontext_rep_pairs = cat(2, samecontext_rep_pairs, new_pair);
    end
end
clear r1 r2 new_pair;

% for reliability w also need this for odd/even reps
samecontext_rep_pairs_perrep = [];
for r1 = 1:floor(n_reps/2)
    for r2 = r1+1:floor(n_reps/2)
        new_pair = [r1, r2]';
        samecontext_rep_pairs_perrep = cat(2, samecontext_rep_pairs_perrep, new_pair);
    end
end
clear r1 r2 new_pair;

 %%

            
sr=1000/binWidthMs;
fprintf([num2str(binWidthMs) '\n'])
fprintf([num2str(sr) '\n'])
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



for chann = 1:length(chann_numbers)
    for dur=1:size(unique_segs,2) % for each duration
        
        % all reps random_random & random_natural
        for field=fields
            Y_seg = [];
            samedur_seg_pairs=[];
        
            if isempty(SAR(chann).(field{:})(dur).order1);continue;end
            o1 = SAR(chann).(field{:})(dur).order1; 
            o2 = SAR(chann).(field{:})(dur).order2;          
          
%             Y_seg = stim x time x n_orders x n_reps
            Y_seg(:,:,1,:) = permute(o1,[1,3,2]);
            Y_seg(:,:,2,:) = permute(o2,[1,3,2]);
                 

            samedur_seg_pairs = [[1:size(Y_seg,1)]' [1:size(Y_seg,1)]'];
            

            [diff_context, same_context, same_context_err, ~] = ...
                cross_context_corr_helper(...
                Y_seg, Y_embed_seg, ...
                samedur_order_pairs, diffdur_order_pairs, ...
                samecontext_rep_pairs, diffcontext_rep_pairs, ...
                samedur_seg_pairs, diffdur_seg_pairs, ...
                make_samedur_comparisons, make_diffdur_comparisons, ...
                interleave_samedur, interleave_diffdur, ...
                simfunc, tranweightfn, trancorrfn);

            L.(field{:}).diff_context(:,dur,chann)  = diff_context;
            L.(field{:}).same_context(:,dur,chann) = same_context;    
            L.(field{:}).same_context_err(:,dur,chann) = same_context_err;
                    
        end

        % odd & even reps random_random & random_natural
        for field=fields
            for rep=1:length(reps) 
                Y_seg=[];
                samedur_seg_pairs=[];
                
                if isempty(SAR(chann).(field{:})(dur).order1);continue;end
                o1 = SAR(chann).(field{:})(dur).order1(:,rep:2:n_reps,:); 
                o2 = SAR(chann).(field{:})(dur).order2(:,rep:2:n_reps,:); 
                    
                % Y_seg = stim x time x n_orders x n_reps
                Y_seg(:,:,1,:) = permute(o1,[1,3,2]);
                Y_seg(:,:,2,:) = permute(o2,[1,3,2]);
                
                samedur_seg_pairs = [[1:size(Y_seg,1)]' [1:size(Y_seg,1)]'];
                
                [diff_context, same_context, same_context_err, ~] = ...
                    cross_context_corr_helper(...
                    Y_seg, Y_embed_seg, ...
                    samedur_order_pairs, diffdur_order_pairs, ...
                    samecontext_rep_pairs_perrep, diffcontext_rep_pairs_perrep, ...
                    samedur_seg_pairs, diffdur_seg_pairs, ...
                    make_samedur_comparisons, make_diffdur_comparisons, ...
                    interleave_samedur, interleave_diffdur, ...
                    simfunc, tranweightfn, trancorrfn);

                L.(field{:}).(reps{rep}).diff_context(:,dur,chann)  = diff_context;
                L.(field{:}).(reps{rep}).same_context(:,dur,chann) = same_context;    
                L.(field{:}).(reps{rep}).same_context_err(:,dur,chann) = same_context_err;
            
            
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
                'TCI' filesep 'stim-aligned-v2' filesep signalType filesep 'binWidth-' num2str(1000/sr) 'ms' ]);
             
    L.output_directory = strcat([sourceFolder filesep 'analysis' filesep 'data-work' filesep ...
                'sub-' lower(indiv) filesep 'sess-' sprintf('%03d',sessN) filesep taskType filesep trialType filesep ...
                'TCI' filesep 'stim-aligned-v2' filesep signalType filesep 'binWidth-' num2str(1000/sr) 'ms'  ]);
             

end
% Custom validation function
function mustBeDoubleorCell(chann_numbers)

    if ~isa(chann_numbers, 'double') && ~isa(chann_numbers, 'cell')
        error('chann_numbers argument must be a double or a cell array.');
    end
end

