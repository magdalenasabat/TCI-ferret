function [diffdur,samedur,lag_t,binWidthMs] = mySAR(data, dtime, sr, M, args)
% Newer version of MYSTIMALIGNTCI which only aligned one channel
% here we align all channels at once to creates a SAR matrix for TCI analysis
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
    dtime {mustBeEqualLenght(data,dtime)} 
    sr (1,1) 
    M struct
    args.method string = 'hist'  
    args.binWidthSec (1,1) double = 0.01
    args.compute_natural_context (1,1) logical = true
    args.beforeT double = 0.1;
    args.afterT double = 0.6; 
end
% pick up args
args.binWidthMs = args.binWidthSec *1000;
newSr = 1/args.binWidthSec;

beforeT = args.beforeT;
afterT = args.afterT;

dt_data = mean(diff(dtime));
dt= round(args.binWidthSec / mean(diff(dtime)),0);

% create timebins for SAR
timebins = -beforeT:dt_data:afterT;
lags = (-beforeT:(dt*dt_data):afterT);
lag_t = lags(1:end-1);

%figure out bin size
binWidthMs = dt*dt_data*1000;  

%pick up number of reps
n_reps = max(M.rep);
   
%initialise output variables
diffdur.order1 = [];
diffdur.order2 = [];
diffdur.durMs = [];
diffdur.durN = [];
diffdur.n_total_segs = NaN(length(M.stimuli_param.durations),1);
samedur.order1 = [];
samedur.order2 = [];
samedur.durMs = [];
samedur.durN = [];
samedur.n_total_segs = NaN(length(M.stimuli_param.durations),1);

% for each duration
for dur=fliplr(unique(M.durId))

    fprintf( strcat([ ' *** Duration: '  num2str(M.stimuli_param.durations(dur))  ' sec\n' ] ));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % align responses to each stimuli of given duraiton for both orders
    % order1

    %figue out which stimuli to pick up
    mask = M.durId(:) == dur & M.order(:) == 1  ;
    trigs1 = M.triggers_sec(mask);
    stimId1 = M.stimId(mask);
    stim_unique = unique(stimId1);
    rep1  =M.rep(mask);
    tic
    % now align 
    for ch=fliplr(1:size(data,1)) % flipping so that on the first iteration th entire matrix is created and it doesnt append after
            
        switch args.method
            case 'hist'
                % for each trigger find responses around, it's going to be
                % vectorized so we need to reshape it, we sum the response
                % over to bin it at higher sampling rate
                O = cell2mat(arrayfun( @(x) sum(reshape(data(ch, find(dtime >= (x-beforeT),length(timebins)-1) ), dt, []),1) , trigs1','UniformOutput',false )) ;
                P1 = []; %reshape(O,[],n_reps,size(O,2));
                %this is usually correct anyways but I'm jjuts making
                %double sure the responses are aligned the same way for
                %both orders
                for sp=1:length(stim_unique)
                    for rp=unique(rep1)
                        P1(sp,rp,:) =  O(rep1==rp & stimId1==stim_unique(sp),:);
                    end
                end
            case 'resample'
                O = cell2mat(arrayfun( @(x) myResample(data(ch, find(dtime >= (x-beforeT),length(timebins)-1) ), sr, newSr) , trigs1','UniformOutput',false )) ;
                P1 = reshape(O,[],n_reps,size(O,2));

        end
        diffdur(dur).order1(:,:,:,ch) = P1;
    end
    toc
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
    
    tic
    % this bit of code is nightmare but this we where we align the reponses
    % to onsets
    for ch=fliplr(1:size(data,1))
      
        switch args.method
            case 'hist'
                O = cell2mat(arrayfun( @(x) sum(reshape(data(ch, find(dtime > (x-beforeT),length(timebins)-1) ), dt, []),1) , trigs2','UniformOutput',false )) ;
                P2 = []; % reshape(O,[],n_reps,size(O,2));
                 for sp=1:1:length(stim_unique)
                    for rp=unique(rep2)
                        P2(sp,rp,:) =  O(rep2==rp & stimId2==stim_unique(sp),:);
                    end
                end
            case 'resample'
                O = cell2mat(arrayfun( @(x) myResample(data(ch, find(dtime >= (x-beforeT),length(timebins)-1) ), sr, newSr) , trigs2','UniformOutput',false )) ;
                P2 = reshape(O,[],n_reps,size(O,2));
    %             O = cell2mat(arrayfun( @(x) interp1(timebins, data( find(dtime >= (x-beforeT),length(timebins)-1) ), lag_t) , trigs2','UniformOutput',false )) ;
    %             P2 = reshape(O,[],n_reps,size(O,2));
    %             
        end
        diffdur(dur).order2(:,:,:,ch) = P2;
    end
    
    
    toc
    
    % save stimuli info
    diffdur(dur).n_reps=n_reps;
    diffdur(dur).lag_t=lag_t;
    diffdur(dur).durMs = M.stimuli_param.durations(dur)*1000;
    diffdur(dur).durN = size(P1,1);
    diffdur(dur).n_total_segs = size(P1,1);
    diffdur(dur).cat = M.catId(M.durId(:) == dur & M.order(:) == 1 & M.rep(:) == 1 );
    diffdur(dur).temp_mod = M.stimuli_param.temporal_modulation_rates(M.stimId(M.durId(:) == dur & M.order(:) == 1 & M.rep(:) == 1 ));
    diffdur(dur).temp_mod_cat = M.stimuli_param.temporal_modulation_rates_cat(M.stimId(M.durId(:) == dur & M.order(:) == 1 & M.rep(:) == 1 ));
    diffdur(dur).stimId = stim_unique;    
    if isfield(M, 'condId')
        diffdur(dur).cond = M.condId(M.durId(:) == dur & M.order(:) == 1 & M.rep(:) == 1 )';
    end
    if isfield(M, 'rateId')
        diffdur(dur).rate = M.rateId(M.durId(:) == dur & M.order(:) == 1 & M.rep(:) == 1 );
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % natural context (when the segment we comare to is surrounded by its natural context (from the source))
    if args.compute_natural_context 
        
        % figure out natural pairs for this segment duration
        stimuli = unique(stimId1);
        natural_pairs = M.stimuli_param.natural_pairs(stimuli);
        natural_pairs_exploded = cell2mat(natural_pairs');
        natural_pairs_onsets_sec = M.stimuli_param.natural_pairs_onsets(stimuli);
        natural_pairs_onsets_sec_exploded = cell2mat(natural_pairs_onsets_sec');

        %some stimuli won't have natural pairs, get valid stims
        valid_stimuli = find(~cellfun(@isempty,natural_pairs));
        valid_stim_exploded = cell2mat(arrayfun(@(x) repmat(x,1,length(natural_pairs{x})*2),1:length(natural_pairs),'UniformOutput',false));
        
        %if there is no natural pais skip
        if isempty(valid_stim_exploded)
            samedur(dur).durMs = M.stimuli_param.durations(dur)*1000;
            samedur(dur).durN = 0;
            samedur(dur).n_total_segs= 0;
            fprintf( strcat([ ' !!! no natural pairs found ... skipping \n' ] ));
            continue
        end

        %%
%         rep_matrix = reshape(rep1',length(stimuli)/2,[] );
% 
%         % Extract the first column of each unique set of rows
%         rep_order = unique(rep_matrix, 'rows')';
% 
%         % Flatten the matrix to get the desired outcome vector
%         rep_order = rep_order(:)';
%         [~ , rep_order ] = sort(rep_order);
%         myReorder = @(x,y)(x(y));

        % loop over channels and align, this bit of code is super
        % convoluted sorry
        tic
        for ch=fliplr(1:size(data,1))
            
   
            switch args.method
                case 'hist'
                    P4 = cell2mat(arrayfun(@(x) permute(reshape(cell2mat(arrayfun(@(y) sum(reshape(data(ch,  ...
                     find(dtime > ((y + natural_pairs_onsets_sec_exploded(x))-beforeT),length(timebins)-1) ), ...
                     dt, []),1),M.triggers_sec(M.stimId == natural_pairs_exploded(x) )', 'UniformOutput',false )),2,n_reps,[]),[1 2 3]) ...
                     , [1:length(natural_pairs_exploded)]', 'UniformOutput',false ));


                case 'resample'
                    fprintf("Probably doesn't work!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                    P4 = cell2mat(arrayfun(@(x) permute(reshape(cell2mat(arrayfun(@(y) myResample(data(ch,  ...
                     find(dtime > ((y + natural_pairs_onsets_sec_exploded(x))-beforeT),length(timebins)-1)), sr, newSr) ...
                     ,M.triggers_sec(M.stimId == natural_pairs_exploded(x) )', 'UniformOutput',false )),2,n_reps,[]),[1 2 3]) ...
                     , [1:length(natural_pairs_exploded)]', 'UniformOutput',false ));



            end

            samedur(dur).order2(:,:,:,ch) = [P4 ;P4];

        end
        toc
        % order 1 i just random_random so we pick it up 
        samedur(dur).order1 = [  diffdur(dur).order1(valid_stim_exploded,:,:,:);...
                                                diffdur(dur).order2(valid_stim_exploded,:,:,:)];
      
        % save stimuli info
        samedur(dur).n_reps=n_reps;
        samedur(dur).lag_t=lag_t;
        samedur(dur).durMs = M.stimuli_param.durations(dur)*1000;
        samedur(dur).durN = size(samedur(dur).order1,1);
        samedur(dur).n_total_segs = size(samedur(dur).order1,1);
        samedur(dur).cat = repmat(diffdur(dur).cat(valid_stim_exploded),1,2);
        samedur(dur).temp_mod = repmat(diffdur(dur).temp_mod(valid_stim_exploded),1,2);
        samedur(dur).temp_mod_cat = repmat(diffdur(dur).temp_mod_cat(valid_stim_exploded),1,2);
        samedur(dur).stimId = repmat(diffdur(dur).stimId (valid_stim_exploded),1,2);
        if isfield(M, 'condId')
            samedur(dur).cond = repmat(diffdur(dur).cond(valid_stim_exploded),1,2);
        end
        if isfield(M, 'rateId')
            samedur(dur).rate = repmat(diffdur(dur).rate(valid_stim_exploded),1,2);
        end

        end
    end
end




% Custom validation function
function mustBeEqualLenght(a,b)
    % Test for equal size
    if ~isequal(size(a,2),size(b,2))
        eid = 'Size:notEqual';
        msg = 'Axis 2 of the first input must equal length of axis 2 of the second input.';
        error(eid,msg)
    end
end