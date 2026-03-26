function L = myCCC(data, S, t, chann_numbers , indiv, sessN, args)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Calculates correlation for
% random_random and random_natural ('diffdur') 
%
% 
% magdalena.sabat@ens.psl.eu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments 
    data double % stim x rep x time x channel x context
    S struct % stim info
    t double % time for data
    chann_numbers double % channel numbers
    indiv char % subject name
    sessN double % session number 
    args.diffdur double = 0   % if passed it should be stim x rep x time x channel x context
    args.simfunc char = 'corr'
    args.interpolate logical = true

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
n_reps = size(data,2);


%output L will be intrpolated so that binsize = 1ms regardless of input
lag_t_original = t;
if args.interpolate; 
    lag_t=lag_t_original(1):0.001:lag_t_original(end);
else
    lag_t=lag_t_original;
end

% make D
D{1}=data;
D{2}=args.diffdur;


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

binWidthMs = mean(diff(t))*1000 ;           
sr=1000/binWidthMs;

unique_segs = S.unique_segs{1};
n_stim=S.durN{1};


fprintf('calculating correlations  \n')

% L.samedur = myPrepareL( sourceFolder, indiv,sessN, taskType, trialType, signalType, chann_numbers,unique_segs,sr,lag_t,[SAR(1).samedur.n_total_segs]);
L.samedur = myPrepareL( chann_numbers,unique_segs,sr,lag_t,n_stim, indiv, sessN);


if ~isfolder(L.samedur.output_directory)
    mkdir(L.samedur.figure_directory)  
    mkdir(L.samedur.output_directory)  
end  



if length(size(args.diffdur))>1
    fields = {'samedur','diffdur'};
    unique_segs_natural = S.unique_segs{2};
    unique_segs_natural = unique_segs_natural; %(S.durN{2}~=0);
    
    n_total_segs_natural = S.durN{2}; %(S.durN{2}~=0);
    
    
%     L.diffdur = myPrepareL( sourceFolder, indiv,sessN, taskType, trialType, signalType, chann_numbers,unique_segs_natural,sr,lag_t,n_total_segs_natural);
    L.diffdur =myPrepareL( chann_numbers,unique_segs_natural,sr,lag_t,n_total_segs_natural, indiv, sessN);
   
%     if args.create_dir &  ~isfolder(L.diffdur.output_directory)
%         mkdir(L.diffdur.figure_directory)  
%         mkdir(L.diffdur.output_directory)  
%     end
   
else
    fields = {'samedur'};
end



for chann = 1:length(chann_numbers)
    for dur=1:size(unique_segs,2) % for each duration
        durMs=unique_segs(dur);
        f=1;
        % all reps random_random & random_natural
        for field=fields
            Y_seg = [];
            samedur_seg_pairs=[];
        
%             if isempty(SAR.(field{:})(dur).order1); continue; end
            o1 = D{f}(S.segs_dur{f}==durMs,:,:,chann,1);
            o2 = D{f}(S.segs_dur{f}==durMs,:,:,chann,2);        
          
%             Y_seg = stim x time x n_orders x n_reps ; data = stim x rep x time x channel x context
            Y_seg(:,:,1,:) = permute(o1,[1,3,2]);
            Y_seg(:,:,2,:) = permute(o2,[1,3,2]);
                 

            samedur_seg_pairs = [[1:size(Y_seg,1)]' [1:size(Y_seg,1)]'];
            
%             if isempty(samedur_seg_pairs); continue; end

            [diff_context, same_context, same_context_err, ~] = ...
                cross_context_corr_helper(...
                Y_seg, Y_embed_seg, ...
                samedur_order_pairs, diffdur_order_pairs, ...
                samecontext_rep_pairs, diffcontext_rep_pairs, ...
                samedur_seg_pairs, diffdur_seg_pairs, ...
                make_samedur_comparisons, make_diffdur_comparisons, ...
                interleave_samedur, interleave_diffdur, ...
                simfunc, tranweightfn, trancorrfn);
            
            if args.interpolate
                L.(field{:}).diff_context(:,dur,chann)  = interp1(lag_t_original,diff_context,lag_t); 
                L.(field{:}).same_context(:,dur,chann) = interp1(lag_t_original,same_context,lag_t);
                L.(field{:}).same_context_err(:,dur,chann) = interp1(lag_t_original,same_context_err,lag_t);
            else
                L.(field{:}).diff_context(:,dur,chann)  = diff_context; 
                L.(field{:}).same_context(:,dur,chann) = same_context;
                L.(field{:}).same_context_err(:,dur,chann) = same_context_err;
            end
            f=f+1;
                    
        end

       
          % odd & even reps random_random & random_natural
        if n_reps>=4
            f=1;
            for field=fields
                for rep=1:length(reps) 
                    Y_seg=[];
                    samedur_seg_pairs=[];

                    o1 = D{f}(S.segs_dur{f}==durMs,rep:2:n_reps,:,chann,1);
                    o2 = D{f}(S.segs_dur{f}==durMs,rep:2:n_reps,:,chann,2); 
            

                    %make a note of how many stims go into this
                     L.(field{:}).splits_n_total_segs(dur,rep,1) = size(o1,1);

                    % Y_seg = stim x time x n_orders x n_reps
                    Y_seg(:,:,1,:) = permute(o1,[1,3,2]);
                    Y_seg(:,:,2,:) = permute(o2,[1,3,2]);

                    samedur_seg_pairs = [[1:size(Y_seg,1)]' [1:size(Y_seg,1)]'];
%                     if isempty(samedur_seg_pairs); continue; end

                    [diff_context, same_context, same_context_err, ~] = ...
                        cross_context_corr_helper(...
                        Y_seg, Y_embed_seg, ...
                        samedur_order_pairs, diffdur_order_pairs, ...
                        samecontext_rep_pairs_perrep, diffcontext_rep_pairs_perrep, ...
                        samedur_seg_pairs, diffdur_seg_pairs, ...
                        make_samedur_comparisons, make_diffdur_comparisons, ...
                        interleave_samedur, interleave_diffdur, ...
                        simfunc, tranweightfn, trancorrfn);



                    L.(field{:}).splits_diff_context(:,dur,chann,rep,1)  =interp1(lag_t_original,diff_context,lag_t) ;
                    L.(field{:}).splits_same_context(:,dur,chann,rep,1) = interp1(lag_t_original,same_context,lag_t) ;    
                    L.(field{:}).splits_same_context_err(:,dur,chann,rep,1) = interp1(lag_t_original,same_context_err,lag_t) ;


                end
                f=f+1;

            end
        end
        
         % half of the stimuli random_random & random_natural
        f=1;
        for field=fields
            
            
            % if not enough stimuli to split skip
            n_stim_this = sum(S.segs_dur{f}==durMs);
            if n_stim_this<4;continue;end

            for stim=1:2
         
                Y_seg=[];
                samedur_seg_pairs=[];
                
                o1 = D{f}(S.segs_dur{f}==durMs,:,:,chann,1); 
                o1=o1(stim:2:n_stim_this,:,:);
                o2 = D{f}(S.segs_dur{f}==durMs,:,:,chann,2); 
                o2=o2(stim:2:n_stim_this,:,:);
                
                 %make a note of how many stims go into this
                L.(field{:}).splits_n_total_segs(dur,stim,2) = size(o1,1);
                
                % Y_seg = stim x time x n_orders x n_reps
                Y_seg(:,:,1,:) = permute(o1,[1,3,2]);
                Y_seg(:,:,2,:) = permute(o2,[1,3,2]);
                               
                samedur_seg_pairs = [[1:size(Y_seg,1)]' [1:size(Y_seg,1)]'];
%                 if isempty(samedur_seg_pairs); continue; end
                
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
            f=f+1;
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
   
        L.alldur = L_this;

    else
        warning('No natural comparison, nothing to pool over.')
    end
end






function L = myPrepareL( chann_numbers,unique_segs,sr,lag_t,n_total_segs, indiv,sessN)
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
   
    
    L.unique_segs = unique_segs;
    L.sr = sr;
    L.lag_t = lag_t;
    L.n_total_segs = n_total_segs';
    L.figure_directory = strcat(['~/' ]);
             
    L.output_directory = strcat(['~/'  ]);
             

end
% Custom validation function
function mustBeDoubleorCell(chann_numbers)

    if ~isa(chann_numbers, 'double') && ~isa(chann_numbers, 'cell')
        error('chann_numbers argument must be a double or a cell array.');
    end
end

