function [random_random_context,random_natural_context,lag_t,binWidthMs] = myStimAlignTCI(data, dtime, sr, M, args)
%MYSTIMALIGNTCI creates a SAR matrix for TCI analysis
%  
% data  - a vector of data to align, given it's always a whole number and
% it's always positive it should be unit16
% time  - a vector of time for each point in data 
% sr - sampling rate of data
% M - struct with fields : 
% M.stimId -  id of stimulus
% M.durId -  id of duration of each stimulus
% M.order - id of random context
% M.rep     - no of rep 
% M.catId - id of category of stimulus 
% M.triggers_smp - triggers in samples (same as SR as data) 
% M.stimuli_param with fields 
%   M.stimuli_param.natural_pairs - ids (refers to stimId) of natural pairs
%   for each stim
%   M.stimuli_param.names - char array name for each stim 
%   M.stimuli_param.durations - duration for each id in durId 
%   M.stimuli_param.categories - category name for each if in m.category
%
% OUT:
%
% 16/01/2024
% Magdalena, magdalena.sabat@ens.psl.eu
%
arguments 
    data 
    dtime {mustBeEqualSize(data,dtime)} 
    sr (1,1) 
    M struct
    args.method string = 'hist'  
    args.binWidthSec (1,1) double = 0.01
    args.plots (1,1) logical = false
    args.compute_natural_context (1,1) logical = true
    args.beforeT double = 0.1;
    args.afterT double = 0.6; 
end
args.binWidthMs = args.binWidthSec *1000;
plots = args.plots;
newSr = 1/args.binWidthSec;

beforeT = args.beforeT;
afterT = args.afterT;

dt_data = mean(diff(dtime));
dt= round(args.binWidthSec / mean(diff(dtime)),0);

% round triggers
% M.triggers_sec = round(M.triggers_sec,abs(floor(log10(dt_data))));

timebins = -beforeT:dt_data:afterT;
lags = (-beforeT:(dt*dt_data):afterT);

n_reps = max(M.rep);

lag_t = lags(1:end-1);
binWidthMs = dt*dt_data*1000;  
   
random_random_context.order1 = [];
random_random_context.order2 = [];
random_random_context.durMs = [];
random_random_context.durN = [];
    
random_random_context.n_total_segs = NaN(length(M.stimuli_param.durations),1);



random_natural_context.order1 = [];
random_natural_context.order2 = [];
random_natural_context.durMs = [];
random_natural_context.durN = [];

random_natural_context.n_total_segs = NaN(length(M.stimuli_param.durations),1);


if plots; figure; end
% for each duration
for dur=unique(M.durId)
   
    fprintf( strcat([ ' *** Duration: '  num2str(M.stimuli_param.durations(dur))  ' sec\n' ] ));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % random context
%     fprintf( strcat([ ' random context ... \n' ] ));
    
    % order1
    mask = M.durId(:) == dur & M.order(:) == 1  ;
    trigs1 = M.triggers_sec(mask);
    stimId1 = M.stimId(mask);
    stim_unique = unique(stimId1);
    rep1  =M.rep(mask);
    
    switch args.method
        case 'hist'
            O = cell2mat(arrayfun( @(x) sum(reshape(data( find(dtime >= (x-beforeT),length(timebins)-1) ), dt, []),1) , trigs1','UniformOutput',false )) ;
            P1 = []; %reshape(O,[],n_reps,size(O,2));
            for sp=1:length(stim_unique)
                for rp=unique(rep1)
                    P1(sp,rp,:) =  O(rep1==rp & stimId1==stim_unique(sp),:);
                end
            end
        case 'resample'
            O = cell2mat(arrayfun( @(x) myResample(data( find(dtime >= (x-beforeT),length(timebins)-1) ), sr, newSr) , trigs1','UniformOutput',false )) ;
            P1 = reshape(O,[],n_reps,size(O,2));
            
    end
   

    % order2
    mask = M.durId(:) == dur & M.order(:) == 2  ;
    trigs_org = M.triggers_sec(mask);
    stimId2_org = M.stimId(mask);
    stimId2=[];
    rep2  =M.rep(mask);
   
    %%% very important step - order 2 is shuffled differently thant order 1
    %%% but is made of the same stimuli - in order to correlate them
    %%% directly we need to reorganise the matrix so that at index 1 we
    %%% have the response to the same stimulus in both P1 and P2
    
    reorganize = reshape(grp2idx([stimId1 stimId2_org]) , [], 2);
    trigs2 = zeros(size(trigs_org));
    for m=unique(reorganize)'     
        trigs2(find(reorganize(:,1)==m)) = trigs_org(find(reorganize(:,2)==m));
        stimId2(find(reorganize(:,1)==m)) = stimId2_org(find(reorganize(:,2)==m));
    end 
    


    switch args.method
        case 'hist'
            O = cell2mat(arrayfun( @(x) sum(reshape(data( find(dtime > (x-beforeT),length(timebins)-1) ), dt, []),1) , trigs2','UniformOutput',false )) ;
            P2 = []; % reshape(O,[],n_reps,size(O,2));
             for sp=1:1:length(stim_unique)
                for rp=unique(rep2)
                    P2(sp,rp,:) =  O(rep2==rp & stimId2==stim_unique(sp),:);
                end
            end
        case 'resample'
            O = cell2mat(arrayfun( @(x) myResample(data( find(dtime >= (x-beforeT),length(timebins)-1) ), sr, newSr) , trigs2','UniformOutput',false )) ;
            P2 = reshape(O,[],n_reps,size(O,2));
%             O = cell2mat(arrayfun( @(x) interp1(timebins, data( find(dtime >= (x-beforeT),length(timebins)) ), lag_t) , trigs2','UniformOutput',false )) ;
%             P2 = reshape(O,[],n_reps,size(O,2));
%             
    end
    
  
%     % correlations
    if plots
        fprintf( strcat([ ' calculating noise ceiling and cross context correlation \n' ] ));
    %     
        o1 = P1; o2 = P2;
        tm = logical(tril(ones(size(o1 ,2)),-1));

        x = []; y = []; z = [];
        for t=1:size(o1  ,3)
            x(t) = mean(corr(squeeze(o1 (:,:,t)),squeeze(o2 (:,:,t))),'all','omitnan');

            r = corr(squeeze(o1 (:,:,t)));
            y(t) = mean(r(tm),'all','omitnan');

            r = corr(squeeze(o2 (:,:,t)));
            z(t) = mean(r(tm),'all','omitnan');
        end

        q = (y+z)/2;
        subplot(length(unique(M.durId)),1,dur)      
        plot(lag_t,x); hold on; plot(lag_t,q); title(num2str(M.stimuli_param.durations(dur)))
    end
    
    random_random_context(dur).order1 = P1;
    random_random_context(dur).order2 = P2;
    random_random_context(dur).durMs = M.stimuli_param.durations(dur)*1000;
    random_random_context(dur).durN = size(P1,1);
    random_random_context(dur).n_total_segs = size(P1,1);
    random_random_context(dur).cat = M.catId(M.durId(:) == dur & M.order(:) == 1 & M.rep(:) == 1 );
    random_random_context(dur).temp_mod = M.stimuli_param.temporal_modulation_rates(M.stimId(M.durId(:) == dur & M.order(:) == 1 & M.rep(:) == 1 ));
    random_random_context(dur).temp_mod_cat = M.stimuli_param.temporal_modulation_rates_cat(M.stimId(M.durId(:) == dur & M.order(:) == 1 & M.rep(:) == 1 ));
    random_random_context(dur).stimId = stimId2;    
    if isfield(M, 'condId')
        random_random_context(dur).cond = M.condId(M.durId(:) == dur & M.order(:) == 1 & M.rep(:) == 1 )';
    end
    if isfield(M, 'rateId')
        random_random_context(dur).rate = M.rateId(M.durId(:) == dur & M.order(:) == 1 & M.rep(:) == 1 );
    end
    
    % pick up onsets and durations of surroundign segments
    % take only second repetition
    trigs= M.triggers_sec(M.durId(:) == dur & M.rep(:) == 2 ) ;
    % pick up indexes for that segment withn +/- X sec
    surrounding_idx = arrayfun(@(x) find(M.triggers_sec >= (x-6) & ...
                                         M.triggers_sec <= (x+6)),trigs , 'UniformOutput', false);
    
    % pickup onsets and durations
    surrounding_durations = cellfun(@(x) M.stimuli_param.durations(M.durId(x))' ,surrounding_idx, 'UniformOutput', false);
    surrounding_onsets = cellfun(@(x,y) M.triggers_sec(x)-y ,surrounding_idx,num2cell(trigs), 'UniformOutput', false);
    
    random_random_context(dur).surrounding_durations = surrounding_durations; % cell2mat(surrounding_durations);
    random_random_context(dur).surrounding_onsets = surrounding_onsets ; %cell2mat(surrounding_onsets);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % natural context
    if args.compute_natural_context 
%         fprintf( strcat([ ' natural context ... (this may take a while) \n' ] ));

        stimuli = unique(stimId1);
        natural_pairs = M.stimuli_param.natural_pairs(stimuli);
        natural_pairs_exploded = cell2mat(natural_pairs');
        natural_pairs_onsets_sec = M.stimuli_param.natural_pairs_onsets(stimuli);
        natural_pairs_onsets_sec_exploded = cell2mat(natural_pairs_onsets_sec');

        %some stimuli won't have natural pairs, get valid stims
        valid_stimuli = find(~cellfun(@isempty,natural_pairs));
        valid_stim_exploded = cell2mat(arrayfun(@(x) repmat(x,1,length(natural_pairs{x})*2),1:length(natural_pairs),'UniformOutput',false));
       
        if isempty(valid_stim_exploded)
            random_natural_context(dur).durMs = M.stimuli_param.durations(dur)*1000;
            random_natural_context(dur).durN = 0;
            random_natural_context(dur).n_total_segs= 0;
            fprintf( strcat([ ' !!! no natural pairs found ... skipping \n' ] ));
            continue
        end
        
        %%
        rep_matrix = reshape(rep1',length(stimuli)/2,[] );

        % Extract the first column of each unique set of rows
        rep_order = unique(rep_matrix, 'rows')';

        % Flatten the matrix to get the desired outcome vector
        rep_order = rep_order(:)';
        [~ , rep_order ] = sort(rep_order);
        myReorder = @(x,y)(x(y));

          %% 

        P31 = P1(valid_stim_exploded,:,:);
        P32 = P2(valid_stim_exploded,:,:);
   
          
        switch args.method
            case 'hist'
                P4 = cell2mat(arrayfun(@(x) permute(reshape(cell2mat(arrayfun(@(y) sum(reshape(data(  ...
                 find(dtime > ((y + natural_pairs_onsets_sec_exploded(x))-beforeT),length(timebins)-1) ), ...
                 dt, []),1),M.triggers_sec(M.stimId == natural_pairs_exploded(x) )', 'UniformOutput',false )),2,n_reps,[]),[1 2 3]) ...
                 , [1:length(natural_pairs_exploded)]', 'UniformOutput',false ));

             
            case 'resample'
                fprintf("Probably doesn't work!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                P4 = cell2mat(arrayfun(@(x) permute(reshape(cell2mat(arrayfun(@(y) myResample(data(  ...
                 find(dtime > ((y + natural_pairs_onsets_sec_exploded(x))-beforeT),length(timebins)-1)), sr, newSr) ...
                 ,M.triggers_sec(M.stimId == natural_pairs_exploded(x) )', 'UniformOutput',false )),2,n_reps,[]),[1 2 3]) ...
                 , [1:length(natural_pairs_exploded)]', 'UniformOutput',false ));
             

            
        end
    

   
        %    % correlations
        if plots
            fprintf( strcat([ ' calculating noise ceiling and cross context correlation \n' ] ));
        %     
            o1 = [P31;P32]; o2 = [P4 ;P4];
            tm = logical(tril(ones(size(o1 ,2)),-1));

            x = []; y = []; z = [];
            for t=1:size(o1  ,3)
                x(t) = mean(corr(squeeze(o1 (:,:,t)),squeeze(o2 (:,:,t))),'all','omitnan');

                r = corr(squeeze(o1 (:,:,t)));
                y(t) = mean(r(tm),'all','omitnan');

                r = corr(squeeze(o2 (:,:,t)));
                z(t) = mean(r(tm),'all','omitnan');
            end

            q =  (y+z)/2;
            subplot(length(unique(M.durId)),1,dur)      
            plot(lag_t,x); hold on; plot(lag_t,q); title(num2str(M.stimuli_param.durations(dur)))
        end

        random_natural_context(dur).order1 = [P31;P32];
        random_natural_context(dur).order2 = [P4 ;P4];
        random_natural_context(dur).durMs = M.stimuli_param.durations(dur)*1000;
        random_natural_context(dur).durN = size(random_natural_context(dur).order1,1);
        random_natural_context(dur).n_total_segs = size(random_natural_context(dur).order1,1);
        random_natural_context(dur).cat = repmat(random_random_context(dur).cat(valid_stim_exploded),1,2);
        random_natural_context(dur).temp_mod = repmat(random_random_context(dur).temp_mod(valid_stim_exploded),1,2);
        random_natural_context(dur).temp_mod_cat = repmat(random_random_context(dur).temp_mod_cat(valid_stim_exploded),1,2);
        if isfield(M, 'condId')
            random_natural_context(dur).cond = repmat(random_random_context(dur).cond(valid_stim_exploded),1,2);
        end
        if isfield(M, 'rateId')
            random_natural_context(dur).rate = repmat(random_random_context(dur).rate(valid_stim_exploded),1,2);
        end
        %TO DECIDE!
%         % pick up onsets and durations of surroundign segments
%         trigs=[trigs1 trigs2];
%         % pick up indexes for that segment withn +/- X sec
%         surrounding_idx = arrayfun(@(x) find(M.triggers_sec >= (x-3*max(M.stimuli_param.durations)) & ...
%                                              M.triggers_sec <= (x+3*max(M.stimuli_param.durations))),trigs , 'UniformOutput', false);
% 
%         % pickup onsets and durations
%         surrounding_durations = cellfun(@(x) M.stimuli_param.durations(M.durId(x))' ,surrounding_idx, 'UniformOutput', false);
%         surrounding_onsets = cellfun(@(x,y) M.triggers_sec(x)-y ,surrounding_idx,num2cell(trigs), 'UniformOutput', false);
% 
%         random_natural_context(dur).surrounding_durations = cell2mat(surrounding_durations);
%         random_natural_context(dur).surrounding_onsets = cell2mat(surrounding_onsets);
    end
    
end


end

% Custom validation function
function mustBeEqualSize(a,b)
    % Test for equal size
    if ~isequal(size(a),size(b))
        eid = 'Size:notEqual';
        msg = 'Size of first input must equal size of second input.';
        error(eid,msg)
    end
end