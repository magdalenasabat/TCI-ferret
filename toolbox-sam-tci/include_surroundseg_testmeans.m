clc;clear all; close all

sourceFolder = '/home/magdalena/Documents/projects/CST';

        
addpath(strcat([sourceFolder filesep 'analysis' filesep 'code' filesep 'myFunctions']));
addpath(genpath(strcat([sourceFolder filesep 'analysis' filesep 'code' filesep 'toolbox-sam-tci'])));

%%

M={};

% model window params
distr = 'gamma';
widths = [.015,.032,.064,.125,.250];
widths = [.250];
delays=widths/2;
for w=1:length(widths)
    intper_sec = widths(w);
    delay_sec = delays(w);
    shape = 3;
    intervaltype = 'highest';

    internal_sr=100;
    sr=internal_sr;
    rampdur = .015625;

    outer_tsec = -4.5:1/sr:4.5;

    % onsets of surrounding segments
    num_trials = 100;

    time_range = outer_tsec(end)-2.5; % it's gonna be +/- this sec

    overlap_central = zeros(size(outer_tsec));
   
    central_segment_durs=[0.015, 0.032, 0.062, 0.125, 0.250];
    c=0;

    figure ; 
    for central_segment_dur=central_segment_durs
     
        c=c+1;
        %% surrounding segments 
        unique_surrounding_segs_durs = {[0.015, 0.032, 0.062, 0.125, 0.250], [.6,.4,.4,.1,.1]};
        [onsets, durations] = simulate_trials_with_durations(num_trials, central_segment_dur,unique_surrounding_segs_durs, time_range, sr);

        % in L structure we have total_seg per durtion, we should use it
%         % here
        winpow_diff = zeros(size(onsets,2),size(outer_tsec,2));
        overlap_diff = zeros(size(onsets,2),size(outer_tsec,2));
        b=tic;
        for tr=1:size(onsets,2)
            
           
            [winpow_diff(tr,:), ~, overlap_diff(tr,:)] = win_power_ratio_magda(central_segment_dur, ...
                            'gamma', intper_sec, delay_sec, ...
                            'shape', shape, 'tsec', outer_tsec, ...
                            'rampwin', 'hann', 'rampdur', rampdur, ...
                            'intervalmass', .75, ...
                            'intervaltype', intervaltype, ...
                            'delaypoint', 'start',...
                            'surrounding_segments',[onsets{tr};durations{tr}]);
      
        end
        toc(b)
        subplot(length(central_segment_durs),1,c);
        hold on;
        
       
        plot(outer_tsec,mean(winpow_diff,1)); 
        
          tic
        [winpow_diff2, ~, overlap_diff2] = win_power_ratio_magda(central_segment_dur, ...
                            'gamma', intper_sec, delay_sec, ...
                            'shape', shape, 'tsec', outer_tsec, ...
                            'rampwin', 'hann', 'rampdur', rampdur, ...
                            'intervalmass', .75, ...
                            'intervaltype', intervaltype, ...
                            'delaypoint', 'start',...
                            'surrounding_segments',[cell2mat(onsets);cell2mat(durations)],...
                            'n_total_seg',num_trials);
        toc
        plot(outer_tsec,winpow_diff2); 
        
        %%
          %%
        
        a = cell2mat(onsets);
        b = cell2mat(durations);
        
        tmp = arrayfun(@(x) histcounts(a(b==x), outer_tsec),unique(b),'UniformOutput',false);
        tmp1 =cell2mat(tmp')/num_trials;
        tmp1 =tmp1./sum(tmp1);
        [~, argmax] = max(tmp1);
        probabilities = unique_surrounding_segs_durs{1}(argmax);
        
%         figure; plot(outer_tsec(1:end-1),unique_surrounding_segs_durs{1}(argmax))
        
        
%         figure;
%         plot(outer_tsec(1:end-1),tmp1)
%         legend(arrayfun(@(x) num2str(x),unique(b),'UniformOutput',false))
        
        durations_new = [];
        onsets_new = [];
        
        %positive
        ti=0;
        while ti<outer_tsec(end)
            ix = find(outer_tsec(1:end-1)>=ti,1,'first');
            if ~ isempty(ix)
                new_dur = probabilities(ix);
                new_ons = ti;

                durations_new = [durations_new new_dur ];
                onsets_new = [onsets_new new_ons];
                ti=ti+new_dur;
            else
                break
            end
      
        end
        %negative
        ti=-probabilities(find(outer_tsec==0)-1);
        
        durations_new = [ probabilities(find(outer_tsec==0)-1) durations_new];
        onsets_new = [ -probabilities(find(outer_tsec==0)-1) onsets_new];        
        
        while ti>outer_tsec(1)
            ix = find(outer_tsec(1:end-1)>=ti,1,'first');
            if ~ isempty(ix)
                

                new_dur = probabilities(ix);
                ti=ti-new_dur;
                new_ons = ti;

                durations_new = [ new_dur durations_new];
                onsets_new = [ new_ons onsets_new];
            else
                break
            end
            
      
        end
%         figure;plot(onsets_new,durations_new)
             
%%
        b=tic;
        [winpow_diff_avg, ~, overlap_diff_avg] = win_power_ratio_magda(central_segment_dur, ...
                            'gamma', intper_sec, delay_sec, ...
                            'shape', shape, 'tsec', outer_tsec, ...
                            'rampwin', 'hann', 'rampdur', rampdur, ...
                            'intervalmass', .75, ...
                            'intervaltype', intervaltype, ...
                            'delaypoint', 'start',...
                            'surrounding_segments',[onsets_new;durations_new],...
                            'n_total_seg',1);
        toc(b)
        plot(outer_tsec,winpow_diff_avg); 
        
        %%
        unique_surrounding_segs_durs = {[central_segment_dur],[1]};
        [onsets, durations] = simulate_trials_with_durations(num_trials, central_segment_dur,unique_surrounding_segs_durs, time_range, sr);
        

        b=tic;
        [winpow_same, ~, overlap_same] = win_power_ratio_magda(central_segment_dur, ...
                            'gamma', intper_sec, delay_sec, ...
                            'shape', shape, 'tsec', outer_tsec, ...
                            'rampwin', 'hann', 'rampdur', rampdur, ...
                            'intervalmass', .75, ...
                            'intervaltype', intervaltype, ...
                            'delaypoint', 'start');
        toc(b)
        plot(outer_tsec,winpow_same); title(num2str(intper_sec))
        
      
      
    end
    
    legend('mean of trials','treating all trials as single','fake trial','assuming equal durations')
    

end


%%

    % L.indiv='test2';
    % L.diff_context=repmat(ccc_noise2(outer_tsec>=-.5& outer_tsec<=1)',1,2,2);
    % 
    % L.same_context =repmat(noise_cieling(outer_tsec>=-.5& outer_tsec<=1)',1,2,2);
    % 
    % L.same_context_err= repmat(ccc_noise2(outer_tsec>=-.5& outer_tsec<=1).^2',1,2,2);
    % 
    % M = modelfit_cross_context_corr(L, 'overwrite', true, 'plot_figure', true, ...
    %     'shape', 1, 'intper_range', [1/128, 0.50], 'delay_range', [0, 0.25],  ...
    %     'lossfn', 'sqerr');


%% Example usage
% num_trials = 10;
% duration_segments = [0.015, 0.032, 0.062, 0.125, 0.250];
% time_range = 2; % it's gonna be +/- this sec
% sampling_rate = 1000; % Hz
% central_segment_dur=.25;
% [onsets, durations] = simulate_trials_with_durations(num_trials, central_segment_dur,duration_segments, time_range, sampling_rate);



function [onsets, durations] = simulate_trials_with_durations(num_trials,central_segment_dur, duration_segments, time_range, sampling_rate)
    onsets = cell(1, num_trials);
    durations = cell(1, num_trials);
    
    for i = 1:num_trials
        [trial_onsets, trial_durations] = generate_trial_with_durations(central_segment_dur,duration_segments, time_range, sampling_rate);
        onsets{i} = trial_onsets;
        durations{i} = trial_durations;
    end
end

function [onsets, durations] = generate_trial_with_durations(central_segment_dur,duration_segments, time_range, sampling_rate)
    % Generate segment onsets before time 0
    [onsets_before, durations_before] = generate_onsets(duration_segments, time_range, sampling_rate);
    
    % Generate segment onsets after time 0
    [onsets_after, durations_after] = generate_onsets(duration_segments, time_range, sampling_rate);
    
 
    % Combine onsets and durations
    onsets = [-fliplr([onsets_before(2:end) sum(durations_before)]),0,  onsets_after+central_segment_dur];
    durations = [fliplr(durations_before),central_segment_dur, durations_after];
end

function [onsets, durations] = generate_onsets(duration_segments, time_range, sampling_rate)
    onsets = [];
    durations = [];
    current_time = 0;
    while current_time < time_range
        segment_duration = randsample(duration_segments{1}, 1, true, duration_segments{2});
        
        onset_time = current_time;
        onsets = [onsets, onset_time];
        durations = [durations, segment_duration];
        current_time = current_time + segment_duration;
    end
end

