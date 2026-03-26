function L = myCrossContextCorr_v3(SAR, sourceFolder, indiv,sessN, taskType, trialType,signalType, chann_numbers, binWidthMs, args )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Calculates same and cross context correlation for
% random_random and random_natural comparisons needed for TCI
% version that uses cross_context_helper from Sam
%
% IN:
% SAR struct with length == length(chann_numbers)data
% SAR.samedur - size(number of durations x [])
% SAR.samedur.order1 - stimuli x rep x time for each duration
% SAR.samedur.order2 - stimuli x rep x time for each duration
% SAR.samedur.durMs - dur in ms 
% SAR.samedur.durN - N of stimuli
% SAR.samedur.n_total_segs - 
% SAR.samedur.surrounding_durations/onsets - infor about what
% segments surrounded
% SAR.diffdur - srim x rep x time  
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
    binWidthMs double
    args.smoothing logical = false  
    args.compute_natural_context logical = true
    args.create_dir logical = true
    args.simfunc char = 'corr'
    args.prefix char = ''
end
prefix= args.prefix;
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
n_reps = size(SAR(1).samedur(1).order1,2);


%output L will be intrpolated so that binsize = 1ms regardless of input
lag_t_original = SAR.samedur(1).lag_t;
lag_t=lag_t_original(1):0.001:lag_t_original(end);


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
unique_segs = [SAR(1).samedur.durMs];
n_stim=[SAR(1).samedur.durN];
fprintf('calculating correlations  \n')

L.samedur = myPrepareL( sourceFolder, indiv,sessN, taskType, trialType, signalType, chann_numbers,unique_segs,sr,lag_t,[SAR(1).samedur.n_total_segs]);
L.samedur.param_string = strcat(['boundary-samedur_lag_win-' num2str(lag_t(1)) '-' num2str(lag_t(end))]);
L.output_directory =  L.samedur.output_directory;
L.samedur.figure_directory = strcat([L.samedur.figure_directory prefix filesep L.samedur.param_string]);
L.samedur.output_directory = strcat([L.samedur.output_directory prefix filesep L.samedur.param_string]);

if ~isfolder(L.samedur.output_directory)
    mkdir(L.samedur.figure_directory)  
    mkdir(L.samedur.output_directory)  
end  



if isfield(SAR, 'diffdur')
    fields = {'samedur','diffdur'};
    unique_segs_natural = [SAR(1).diffdur.durMs];
    unique_segs_natural = unique_segs_natural(cellfun(@(x) ~isempty(x),{SAR(1).diffdur.order1}));
    
    n_total_segs_natural = [SAR(1).diffdur.n_total_segs];
    n_total_segs_natural = n_total_segs_natural(cellfun(@(x) ~isempty(x),{SAR(1).diffdur.order1}));
    
    
    L.diffdur = myPrepareL( sourceFolder, indiv,sessN, taskType, trialType, signalType, chann_numbers,unique_segs_natural,sr,lag_t,n_total_segs_natural);
    L.diffdur.param_string = strcat(['boundary-noleftright_lag_win-' num2str(lag_t(1)) '-' num2str(lag_t(end)) '_interleave_diffdur-1' ]);
    L.diffdur.figure_directory = strcat([L.diffdur.figure_directory prefix filesep L.diffdur.param_string]);
    L.diffdur.output_directory = strcat([L.diffdur.output_directory prefix filesep L.diffdur.param_string]);
   
    if args.create_dir &  ~isfolder(L.diffdur.output_directory)
        mkdir(L.diffdur.figure_directory)  
        mkdir(L.diffdur.output_directory)  
    end
   
else
    fields = {'samedur'};
end



for chann = 1:length(chann_numbers)
    for dur=1:size(unique_segs,2) % for each duration
        
        % all reps random_random & random_natural
        for field=fields
            Y_seg = [];
            samedur_seg_pairs=[];
        
            if isempty(SAR.(field{:})(dur).order1);continue;end
            o1 = squeeze(SAR.(field{:})(dur).order1(:,:,:,chann)); 
            o2 = squeeze(SAR.(field{:})(dur).order2(:,:,:,chann));          
          
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

            L.(field{:}).diff_context(:,dur,chann)  = interp1(lag_t_original,diff_context,lag_t); 
            L.(field{:}).same_context(:,dur,chann) = interp1(lag_t_original,same_context,lag_t);
            L.(field{:}).same_context_err(:,dur,chann) = interp1(lag_t_original,same_context_err,lag_t);
                    
        end

        % odd & even reps random_random & random_natural
        if n_reps>=4
            for field=fields
                if isempty(SAR.(field{:})(dur).order1);continue;end
                for rep=1:length(reps) 
                    Y_seg=[];
                    samedur_seg_pairs=[];


                    o1 = squeeze(SAR.(field{:})(dur).order1(:,rep:2:n_reps,:,chann)); 
                    o2 = squeeze(SAR.(field{:})(dur).order2(:,rep:2:n_reps,:,chann)); 

                    %make a note of how many stims go into this
                     L.(field{:}).splits_n_total_segs(dur,rep,1) = size(o1,1);

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




                    %changed august 2024 -> model fitting is faster that wat

                    L.(field{:}).splits_diff_context(:,dur,chann,rep,1)  =interp1(lag_t_original,diff_context,lag_t) ;
                    L.(field{:}).splits_same_context(:,dur,chann,rep,1) = interp1(lag_t_original,same_context,lag_t) ;    
                    L.(field{:}).splits_same_context_err(:,dur,chann,rep,1) = interp1(lag_t_original,same_context_err,lag_t) ;


                end

            end
        end
        
         % half of the stimuli random_random & random_natural
        for field=fields
            %if empty skip (random context doesn't have data for all
            %durations
            if isempty(SAR.(field{:})(dur).order1);continue;end

            % if not enough stimuli to split ski[
            n_stim_this = SAR.(field{:})(dur).durN;
            if n_stim_this<4;continue;end

            for stim=1:2
         
                Y_seg=[];
                samedur_seg_pairs=[];
                
                o1 = squeeze(SAR.(field{:})(dur).order1(stim:2:n_stim_this,:,:,chann)); 
                o2 = squeeze(SAR.(field{:})(dur).order2(stim:2:n_stim_this,:,:,chann)); 
                
                 %make a note of how many stims go into this
                L.(field{:}).splits_n_total_segs(dur,stim,2) = size(o1,1);
                
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
                
                 [diff_context, same_context, same_context_err, ~] = ...
                    cross_context_corr_helper(...
                    Y_seg, Y_embed_seg, ...
                    samedur_order_pairs, diffdur_order_pairs, ...
                    samecontext_rep_pairs, diffcontext_rep_pairs, ...
                    samedur_seg_pairs, diffdur_seg_pairs, ...
                    make_samedur_comparisons, make_diffdur_comparisons, ...
                    interleave_samedur, interleave_diffdur, ...
                    simfunc, tranweightfn, trancorrfn);

                
                %changed august 2024 -> model fitting is faster that way
                % we add two split - by number of reps and by stimuli
                L.(field{:}).splits_diff_context(:,dur,chann,stim,2)  = interp1(lag_t_original,diff_context,lag_t) ;
                L.(field{:}).splits_same_context(:,dur,chann,stim,2) = interp1(lag_t_original,same_context,lag_t) ;    
                L.(field{:}).splits_same_context_err(:,dur,chann,stim,2) = interp1(lag_t_original,same_context_err,lag_t) ;
            
            
            end
        
        end
        
        
      
        end
    end
L=myAddPooled(L);
end


function L=myAddPooled(L)
    if isfield(L, 'diffdur')

        L_this = L.samedur ;
        % compute weighted average of the two contexts
        valid_segs = L.diffdur.n_total_segs > 1 ;
        L_this.diff_context(:,valid_segs,:)     = [ L.diffdur.diff_context(:,valid_segs,:).* L.diffdur.n_total_segs(valid_segs)' + ...
                                                    L.samedur.diff_context(:,valid_segs,:).* L.samedur.n_total_segs(valid_segs)']./ ...
                                                    (L.diffdur.n_total_segs(valid_segs) + L.samedur.n_total_segs(valid_segs))';

        L_this.same_context(:,valid_segs,:)     = [ L.diffdur.same_context(:,valid_segs,:).* L.diffdur.n_total_segs(valid_segs)' + ...
                                                    L.samedur.same_context(:,valid_segs,:).* L.samedur.n_total_segs(valid_segs)']./ ...
                                                    (L.diffdur.n_total_segs(valid_segs) + L.samedur.n_total_segs(valid_segs))';

        L_this.same_context_err(:,valid_segs,:) = [ L.diffdur.same_context_err(:,valid_segs,:).* L.diffdur.n_total_segs(valid_segs)' + ...
                                                    L.samedur.same_context_err(:,valid_segs,:).* L.samedur.n_total_segs(valid_segs)']./ ...
                                                   (L.diffdur.n_total_segs(valid_segs) + L.samedur.n_total_segs(valid_segs))';
        
       
                                               
        L_this.n_total_segs(valid_segs) = L_this.n_total_segs(valid_segs) + L.diffdur.n_total_segs(valid_segs);
        
        
        
        % same for splits
        for s=1:size(L.diffdur.splits_n_total_segs ,3)
            for p=1:size(L.diffdur.splits_n_total_segs ,2)
                valid_segs = L.diffdur.splits_n_total_segs(:,p,s) > 1 ;
                 L_this.splits_diff_context(:,valid_segs,:,p,s)     = [ L.diffdur.splits_diff_context(:,valid_segs,:,p,s).* L.diffdur.splits_n_total_segs(valid_segs,p,s)' + ...
                                                            L.samedur.splits_diff_context(:,valid_segs,:,p,s).* L.samedur.splits_n_total_segs(valid_segs,p,s)']./ ...
                                                            (L.diffdur.splits_n_total_segs(valid_segs,p,s) + L.samedur.splits_n_total_segs(valid_segs,p,s))';

                L_this.splits_same_context(:,valid_segs,:,p,s)     = [ L.diffdur.splits_same_context(:,valid_segs,:,p,s).* L.diffdur.splits_n_total_segs(valid_segs,p,s)' + ...
                                                            L.samedur.splits_same_context(:,valid_segs,:,p,s).* L.samedur.splits_n_total_segs(valid_segs,p,s)']./ ...
                                                            (L.diffdur.splits_n_total_segs(valid_segs,p,s) + L.samedur.splits_n_total_segs(valid_segs,p,s))';

                L_this.splits_same_context_err(:,valid_segs,:,p,s) = [ L.diffdur.splits_same_context_err(:,valid_segs,:,p,s).* L.diffdur.splits_n_total_segs(valid_segs,p,s)' + ...
                                                            L.samedur.splits_same_context_err(:,valid_segs,:,p,s).* L.samedur.splits_n_total_segs(valid_segs,p,s)']./ ...
                                                           (L.diffdur.splits_n_total_segs(valid_segs,p,s) + L.samedur.splits_n_total_segs(valid_segs,p,s))';
                
                 L_this.splits_n_total_segs(valid_segs,p,s) = L.diffdur.splits_n_total_segs(valid_segs,p,s)  + L.samedur.splits_n_total_segs(valid_segs,p,s);
    
            end
        end
        
              
        %add new path name for pooled comparison
        L_this.param_string = strcat(['boundary-samedur-noleftright_lag_win-' num2str(L_this.lag_t(1)) '-' num2str(L_this.lag_t(end)) '_interleave_diffdur-1']);
        L_this.output_directory = [fileparts(L.samedur.output_directory) filesep L_this.param_string  ];
        L_this.figure_directory = [fileparts(L.samedur.output_directory) filesep L_this.param_string ];
        L.alldur = L_this;

    else
        warning('No natural comparison, nothing to pool over.')
    end
end






function L = myPrepareL( sourceFolder,indiv,sessN, taskType, trialType, signalType, chann_numbers,unique_segs,sr,lag_t,n_total_segs)
    L = [];
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
                'TCI' filesep 'stim-aligned-v3' filesep signalType filesep 'binWidth-' num2str(1000/sr) 'ms' ]);
             
    L.output_directory = strcat([sourceFolder filesep 'analysis' filesep 'data-work' filesep ...
                'sub-' lower(indiv) filesep 'sess-' sprintf('%03d',sessN) filesep taskType filesep trialType filesep ...
                'TCI' filesep 'stim-aligned-v3' filesep signalType filesep 'binWidth-' num2str(1000/sr) 'ms'  ]);
             

end
% Custom validation function
function mustBeDoubleorCell(chann_numbers)

    if ~isa(chann_numbers, 'double') && ~isa(chann_numbers, 'cell')
        error('chann_numbers argument must be a double or a cell array.');
    end
end

