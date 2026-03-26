clc;clear all; close all

sourceFolder = '/home/magdalena/Documents/projects/CST';

        
addpath(strcat([sourceFolder filesep 'analysis' filesep 'code' filesep 'myFunctions']));
addpath(genpath(strcat([sourceFolder filesep 'analysis' filesep 'code' filesep 'toolbox-sam-tci'])));

%%

M={};

% model window params
distr = 'gamma';
widths = [.015,.032,.064,.124,.250];
delays=widths/2;
for w=1:length(widths)
    intper_sec = widths(w);
    delay_sec = delays(w);
    shape = 3;
    intervaltype = 'highest';

    internal_sr=500;
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

    for central_segment_dur=central_segment_durs
        c=c+1

        %% central segment
        segdur_sec =central_segment_dur;
        time_boundary = [-1,1]*(3*segdur_sec + 4*intper_sec + abs(delay_sec));
        % Perform convolution

        % internal time vectory used for convolution
        internal_tsec = time_boundary(1):1/internal_sr:time_boundary(2);

        % adjust vector so there is a zero timepoint
        [~,zero_tp] = min(abs(internal_tsec-0));
        internal_tsec = internal_tsec - internal_tsec(zero_tp);


        % creates a boxcar function timecourse 
        seg = zeros(size(internal_tsec));
        xi1 = find(internal_tsec <= segdur_sec & internal_tsec >= 0, 1, 'first');
        xi2 = find(internal_tsec <= segdur_sec & internal_tsec >= 0, 1, 'last');
        rampdur_smp = round(rampdur*sr);
        halframpdur_smp1 = floor(rampdur_smp/2);
        halframpdur_smp2 = rampdur_smp - halframpdur_smp1;
        xi = xi1-halframpdur_smp1:xi2+halframpdur_smp2;
        seg(xi) = ramp_hann(ones(length(xi),1), rampdur, sr);


        [h,~,causal] = modelwin(distr, intper_sec, delay_sec,...
            'shape', shape, 'delaypoint', 'start', ...
            'tsec', internal_tsec, 'intervalmass', 0.75, ...unique_surrounding_segs_durs
            'intervaltype', intervaltype);

        % convolve that window the segment
        overlap = myconv(seg', h', 'causal', false, 'central_point', zero_tp);

        seg_onset = find(outer_tsec>=0,1,'first');
        oi1 = seg_onset-xi(1);
        oi2 = oi1+length(overlap)-1;
        overlap_central(oi1:oi2) = overlap_central(oi1:oi2)+overlap';



        %% surrounding segments 
        unique_surrounding_segs_durs = {[0.015, 0.032, 0.062, 0.125, 0.250], [.3,.3,.2,.1,.1]};


        [onsets, durations] = simulate_trials_with_durations(num_trials, central_segment_dur,unique_surrounding_segs_durs, time_range, sr);

        overlap_surround_trials = zeros(num_trials,length(outer_tsec));

        % for each segment in 
        for si=1:length(unique_surrounding_segs_durs{1})

            segdur_sec =unique_surrounding_segs_durs{1}(si);
            time_boundary = [-1,1]*(3*segdur_sec + 4*intper_sec + abs(delay_sec));
            % Perform convolution

            % internal time vectory used for convolution
            internal_tsec = time_boundary(1):1/internal_sr:time_boundary(2);

            % adjust vector so there is a zero timepoint
            [~,zero_tp] = min(abs(internal_tsec-0));
            internal_tsec = internal_tsec - internal_tsec(zero_tp);



            % creates a boxcar function timecourse 
            seg = zeros(size(internal_tsec));
            xi1 = find(internal_tsec <= segdur_sec & internal_tsec >= 0, 1, 'first');
            xi2 = find(internal_tsec <= segdur_sec & internal_tsec >= 0, 1, 'last');
            rampdur_smp = round(rampdur*sr);
            halframpdur_smp1 = floor(rampdur_smp/2);
            halframpdur_smp2 = rampdur_smp - halframpdur_smp1;
            xi = xi1-halframpdur_smp1:xi2+halframpdur_smp2;
            seg(xi) = ramp_hann(ones(length(xi),1), rampdur, sr);


            [h,~,causal] = modelwin(distr, intper_sec, delay_sec,...
                'shape', shape, 'delaypoint', 'start', ...
                'tsec', internal_tsec, 'intervalmass', 0.75, ...
                'intervaltype', intervaltype);

            % convolve that window the segment
            overlap = myconv(seg', h', 'causal', false, 'central_point', zero_tp);

        %     figure; hold on;
            surround_segs_onsets = cellfun(@(ix) onsets{ix}(find(durations{ix}==segdur_sec)) , num2cell(1:1:length(durations)), 'UniformOutput' ,false);
            for i=1:num_trials

                if ~isempty(surround_segs_onsets{i})
                    for j=1:length(surround_segs_onsets{i})
                        % now insert the overlaps into right places
                        seg_onset = find(outer_tsec>=surround_segs_onsets{i}(j),1,'first');
                        oi1 = seg_onset-xi(1);
                        oi2 = oi1+length(overlap)-1;
                        overlap_surround_trials(i,oi1:oi2) = overlap_surround_trials(i,oi1:oi2)+overlap'.^2;
                    end
                end
            end

        end

        %%
    %     figure; 
    %     plot(outer_tsec,overlap_surround_trials'); 
    %     xlabel('time')
    %     title('simualted overlap surround per trial')

        overlap_surround = mean(overlap_surround_trials,1);
    %     figure; 
    %     plot(outer_tsec,overlap_surround); hold on; 
    %     plot(outer_tsec,overlap_central)
    %     legend('overlap surround', 'overlap central')
    %     xlabel('time')
    %     box off



        %% calculate CCC
        % ccc_noisefree = overlap_central.^2 ./ (overlap_central.^2 + sum(overlap_surround.^2));
        % removing sum & square on overlap surround because already done
        ccc_noisefree = (overlap_central.^2) ./ (overlap_central.^2 + overlap_surround);

        % sim noise ceiling
        % rng(1)
        a = .05;
        b = .5;
        noise_cieling = a.*randn(size(ccc_noisefree)) + b;
        % correct by noise
        ccc_noise = noise_cieling .* ccc_noisefree;


        % plot



    %     figure ; 
    %     plot(outer_tsec(outer_tsec>=-.5& outer_tsec<=1) ,ccc_noise(outer_tsec>=-.5& outer_tsec<=1)); 
    %     title('predicted CCC corrected by noise')
    %     xlabel('time')
    %     box off

        %% if surrounding segments were equal
        unique_surrounding_segs_durs = {[central_segment_dur],[1]};
        [onsets, durations] = simulate_trials_with_durations(num_trials, central_segment_dur,unique_surrounding_segs_durs, time_range, sr);

        overlap_surround_trials = zeros(num_trials,length(outer_tsec));

        % for each segment in 
        for si=1:length(unique_surrounding_segs_durs{1})

            segdur_sec =unique_surrounding_segs_durs{1}(si);
            time_boundary = [-1,1]*(3*segdur_sec + 4*intper_sec + abs(delay_sec));
            % Perform convolution

            % internal time vectory used for convolution
            internal_tsec = time_boundary(1):1/internal_sr:time_boundary(2);

            % adjust vector so there is a zero timepoint
            [~,zero_tp] = min(abs(internal_tsec-0));
            internal_tsec = internal_tsec - internal_tsec(zero_tp);



            % creates a boxcar function timecourse 
            seg = zeros(size(internal_tsec));
            xi1 = find(internal_tsec <= segdur_sec & internal_tsec >= 0, 1, 'first');
            xi2 = find(internal_tsec <= segdur_sec & internal_tsec >= 0, 1, 'last');
            rampdur_smp = round(rampdur*sr);
            halframpdur_smp1 = floor(rampdur_smp/2);
            halframpdur_smp2 = rampdur_smp - halframpdur_smp1;
            xi = xi1-halframpdur_smp1:xi2+halframpdur_smp2;
            seg(xi) = ramp_hann(ones(length(xi),1), rampdur, sr);


            [h,~,causal] = modelwin(distr, intper_sec, delay_sec,...
                'shape', shape, 'delaypoint', 'start', ...
                'tsec', internal_tsec, 'intervalmass', 0.75, ...
                'intervaltype', intervaltype);

            % convolve that window the segment
            overlap = myconv(seg', h', 'causal', false, 'central_point', zero_tp);

        %     figure; hold on;
            surround_segs_onsets = cellfun(@(ix) onsets{ix}(find(durations{ix}==segdur_sec)) , num2cell(1:1:length(durations)), 'UniformOutput' ,false);
            for i=1:num_trials

                if ~isempty(surround_segs_onsets{i})
                    for j=1:length(surround_segs_onsets{i})
                        % now insert the overlaps into right places
                        seg_onset = find(outer_tsec>=surround_segs_onsets{i}(j),1,'first');
                        oi1 = seg_onset-xi(1);
                        oi2 = oi1+length(overlap)-1;
                        overlap_surround_trials(i,oi1:oi2) = overlap_surround_trials(i,oi1:oi2)+overlap'.^2;
                    end
                end
            end

        end

        %%

        overlap_surround2 = mean(overlap_surround_trials,1);
        figure; 
%         plot(outer_tsec,overlap_surround); hold on; 
        plot(outer_tsec,overlap_surround2); hold on;
        plot(outer_tsec,overlap_central)
        [onsets, durations] = simulate_trials_with_durations(num_trials, segdur_sec,unique_surrounding_segs_durs, time_range, sr);

        
        [winpow_old, ~, overlap_old] = win_power_ratio_magda(segdur_sec, ...
                        'gamma', intper_sec, delay_sec, ...
                        'shape', shape, 'tsec', outer_tsec, ...
                        'rampwin', 'hann', 'rampdur', rampdur, ...
                        'intervalmass', .75, ...
                        'intervaltype', intervaltype, ...
                        'delaypoint', 'start','surrounding_segments',[onsets{1};durations{1}]);
        plot(outer_tsec,winpow_old)
        legend('overlap surround diff','overlap surround same', 'overlap central','overlap_old')
        xlabel('time')
        box off

        ccc_noisefree2 = (overlap_central.^2) ./ (overlap_central.^2 + overlap_surround2);


        % correct by noise
        ccc_noise2 = noise_cieling .* ccc_noisefree2;



        %%
        figure ; 
        plot(outer_tsec(outer_tsec>=-.5& outer_tsec<=1) ,ccc_noise(outer_tsec>=-.5& outer_tsec<=1));  hold on;
        plot(outer_tsec(outer_tsec>=-.5& outer_tsec<=1) ,ccc_noise2(outer_tsec>=-.5& outer_tsec<=1));  
        plot(outer_tsec(outer_tsec>=-.5& outer_tsec<=1) ,noise_cieling(outer_tsec>=-.5& outer_tsec<=1),'color',[.5 .5 .5]); 
        title('predicted CCC corrected by noise')
        xlabel('time')
        legend('diff dur surr', 'same dur surr','noise ceil')
        box 

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
%     M{w} = modelfit_cross_context_corr(L, 'overwrite', true, 'plot_figure', true, ...
%     'shape', [3], 'intper_range', [1/128, 0.25], 'delay_range', [0, 0.25], ...
%     'lossfn', 'sqerr', 'nintper', 20, 'boundstrength', [1]);
    
    M{w} = modelfit_cross_context_corr(L, 'overwrite', true, 'plot_figure', true, ...
        'shape', [1,3,5,7,9], 'intper_range', [1/128, 0.25], 'delay_range', [0, 0.25], ...
        'lossfn', 'sqerr', 'nintper', 50, 'boundstrength', [0, 0.25, 0.5, 1, 2]);
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

