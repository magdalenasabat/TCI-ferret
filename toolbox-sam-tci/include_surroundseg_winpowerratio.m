clc;clear all; close all

sourceFolder = '/home/magdalena/Documents/projects/CST';

        
addpath(strcat([sourceFolder filesep 'analysis' filesep 'code' filesep 'myFunctions']));
addpath(genpath(strcat([sourceFolder filesep 'analysis' filesep 'code' filesep 'toolbox-sam-tci'])));

%%

M={};

% model window params
distr = 'gamma';
widths = [.064,.124,.250];
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

    L={};
    L.param_string =  'boundary-any_lag_win-0-1';
    L.channels = [1,2];
    L.chnames = {'diff dur','same dur'};
    L.boundary = 'any';
    L.rampwin ='hann';
    L.rampdur = 0.0156;
    L.unique_segs=central_segment_durs*1000;
    L.sr=sr;
    L.lag_t= [outer_tsec(outer_tsec>=-.5& outer_tsec<=1)];
    L.n_total_segs =repmat(num_trials,length(central_segment_durs),1);
    L.figure_directory= '/home/magdalena/Documents/projects/CST/test';
    L.output_directory = '/home/magdalena/Documents/projects/CST/test';
    L.indiv=strcat(['test_multipletrials_width-' num2str(intper_sec*1000) 'ms_delay' num2str(delay_sec) 'ms']);
    figure ; 
    for central_segment_dur=central_segment_durs
     
        c=c+1;
        %% surrounding segments 
        unique_surrounding_segs_durs = {[0.015, 0.032, 0.062, 0.125, 0.250], [.6,.4,.4,.1,.1]};
        [onsets, durations] = simulate_trials_with_durations(num_trials, central_segment_dur,unique_surrounding_segs_durs, time_range, sr);

        % in L structure we have total_seg per durtion, we should use it
        % here
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
        plot(outer_tsec,mean(winpow_diff,1)); title(num2str(intper_sec))
    % ASK SAM : if i do a mean across all delayed segment durations
    % (treating all trials as one) can i just divide it by num_trial? 
        %% calculate CCC
        % ccc_noisefree = overlap_central.^2 ./ (overlap_central.^2 + sum(overlap_surround.^2));
        % removing sum & square on overlap surround because already done
        ccc_noisefree = winpow_diff;

        % sim noise ceiling
        % rng(1)
        a = .05;
        b = .5;
        noise_cieling = a.*randn(size(ccc_noisefree)) + b;
        noise_cieling=ones(size(ccc_noisefree));
        % correct by noise
        ccc_noise = noise_cieling .* ccc_noisefree;



        %% if surrounding segments were equal
        unique_surrounding_segs_durs = {[central_segment_dur],[1]};
        [onsets, durations] = simulate_trials_with_durations(num_trials, central_segment_dur,unique_surrounding_segs_durs, time_range, sr);
        
        
        winpow_same = zeros(size(onsets,2),size(outer_tsec,2));
        overlap_same = zeros(size(onsets,2),size(outer_tsec,2));
        b=tic;
        for tr=1:size(onsets,2)
        
            [winpow_same(tr,:), ~, overlap_same(tr,:)] = win_power_ratio_magda(central_segment_dur, ...
                            'gamma', intper_sec, delay_sec, ...
                            'shape', shape, 'tsec', outer_tsec, ...
                            'rampwin', 'hann', 'rampdur', rampdur, ...
                            'intervalmass', .75, ...
                            'intervaltype', intervaltype, ...
                            'delaypoint', 'start',...
                            'surrounding_segments',[onsets{tr};durations{tr}]);
        end
        toc(b)
        plot(outer_tsec,winpow_same); title(num2str(intper_sec))
        %%


        ccc_noisefree2 = winpow_same;


        % correct by noise
        ccc_noise2 = noise_cieling .* ccc_noisefree2;



        %%
%         subplot(length(central_segment_durs),1,c)
%         plot(outer_tsec(outer_tsec>=-.5& outer_tsec<=1) ,ccc_noise(outer_tsec>=-.5& outer_tsec<=1));  hold on;
%         plot(outer_tsec(outer_tsec>=-.5& outer_tsec<=1) ,ccc_noise2(outer_tsec>=-.5& outer_tsec<=1));  
%         plot(outer_tsec(outer_tsec>=-.5& outer_tsec<=1) ,noise_cieling(outer_tsec>=-.5& outer_tsec<=1),'color',[.5 .5 .5]); 
        

        %% fake L struct
        %         param_string: 'boundary-any_lag_win-0-1'
        %             channels: [1×50 double]
        %              chnames: {50×1 cell}
        %             boundary: 'any'
        %              rampwin: 'hann'
        %              rampdur: 0.0156
        %                indiv: 'simulation'
        %                sessN: 0
        %             taskType: 'contstream'
        %            trialType: 'passive'
        %           signalType: 'stimuli/intper_sec-0.0078125-0.5_delay_sec-0.0039063-0.25intervalmass-0.75_delaypoint-median_intervaltype-center_forcecausal-0'
        %          unique_segs: [15.6250 22.0971 31.2500 44.1942 62.5000 88.3883 125 176.7767 250 500]
        %                   sr: 200
        %                lag_t: [1×140 double]
        %         n_total_segs: [10×1 double]
        %     figure_directory: '/home/magdalena/Documents/projects/CST/analysis/figures/sub-simulation/sess-000/contstream/passive/TCI/stim-aligned/stimuli/intper_sec-0.0078125-0.5_delay_sec-0.0039063-0.25intervalmass-0.75_delaypoint-median_intervaltype-center_forcecausal-0/binWidth-5ms/random_random_context'
        %     output_directory: '/home/magdalena/Documents/projects/CST/analysis/data-work/sub-simulation/sess-000/contstream/passive/TCI/stim-aligned/stimuli/intper_sec-0.0078125-0.5_delay_sec-0.0039063-0.25intervalmass-0.75_delaypoint-median_intervaltype-center_forcecausal-0/binWidth-5ms/random_random_context'
        %                  odd: [1×1 struct]
        %                 even: [1×1 struct]
        %         diff_context: [140×10×50 double]
        %         same_context: [140×10×50 double]
        %     same_context_err: [140×10×50 double]




        L.diff_context(:,c,1)=ccc_noise(outer_tsec>=-.5& outer_tsec<=1);
        L.diff_context(:,c,2)=ccc_noise2(outer_tsec>=-.5& outer_tsec<=1);



        L.same_context(:,c,1)=noise_cieling(outer_tsec>=-.5& outer_tsec<=1);
        L.same_context(:,c,2)=noise_cieling(outer_tsec>=-.5& outer_tsec<=1);


        L.same_context_err(:,c,1)=noise_cieling(outer_tsec>=-.5& outer_tsec<=1).^2;
        L.same_context_err(:,c,2)=noise_cieling(outer_tsec>=-.5& outer_tsec<=1).^2;

        % M = modelfit_cross_context_corr(L, 'overwrite', true, 'plot_figure', true, ...
        %     'shape', 1, 'intper_range', [1/128, 0.50], 'delay_range', [0, 0.25],  ...
        %     'lossfn', 'sqerr');
    end
%     title('predicted CCC corrected by noise')
%         xlabel('time')
%         legend('diff dur surr', 'same dur surr','noise ceil')
%         box 
%     M{w} = modelfit_cross_context_corr(L, 'overwrite', true, 'plot_figure', true, ...
%     'shape', [3], 'intper_range', [1/128, 0.25], 'delay_range', [0, 0.25], ...
%     'lossfn', 'sqerr', 'nintper', 20, 'boundstrength', [1]);
%     
%     M{w} = modelfit_cross_context_corr(L, 'overwrite', true, 'plot_figure', true, ...
%         'shape', [1,3,5,7,9], 'intper_range', [1/128, 0.25], 'delay_range', [0, 0.25], ...
%         'lossfn', 'sqerr', 'nintper', 50, 'boundstrength', [0, 0.25, 0.5, 1, 2]);
end

%%
figure;
hold on;
tmp = [];
for i=1:length(M)
    tmp = [tmp M{i}.best_intper_sec]
end
plot([0 .3],[0 .3])
scatter(tmp(1,:),tmp(2,:))
xlabel('diff dur')
ylabel('same dur')
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

