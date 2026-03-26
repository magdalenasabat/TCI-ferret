clc;clear all; close all

sourceFolder = '/home/magdalena/Documents/projects/CST';

        
addpath(strcat([sourceFolder filesep 'analysis' filesep 'code' filesep 'myFunctions']));
addpath(genpath(strcat([sourceFolder filesep 'analysis' filesep 'code' filesep 'toolbox-sam-tci'])));

%%
% model window params
distr = 'gamma';
intper_sec = 0.1;
delay_sec = 0.1;
shape = 3;
intervaltype = 'highest';

internal_sr=100;
sr=100;
rampdur = .015625;

outer_tsec = -3.5:1/sr:3.5;
overlap_surround = zeros(size(outer_tsec));
overlap_central = zeros(size(outer_tsec));

% onsets of surrounding segments
surround_segs_onsets = [-1.1310,-0.8810,-0.8190,-.6940,-.344 -0.0940,-.032,.250,.375,.437, .469,.7190];
surround_seg_durs = [.25,0.062,0.125,0.25,0.25,0.062, 0.032, 0.125, 0.062, 0.032, 0.125,0.25];


%% central segment
segdur_sec =.25;
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
    'shape', shape, 'delaypoint', 'median', ...
    'tsec', internal_tsec, 'intervalmass', 0.75, ...
    'intervaltype', intervaltype);

% convolve that window the segment
overlap = myconv(seg', h', 'causal', false, 'central_point', zero_tp);

seg_onset = find(outer_tsec>=0,1,'first');
oi1 = seg_onset-xi(1);
oi2 = oi1+length(overlap)-1;
overlap_central(oi1:oi2) = overlap_central(oi1:oi2)+overlap';




%% surrounding segments 



% for each segment in 
for si=1:length(surround_seg_durs)
    
    segdur_sec =surround_seg_durs(si);
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
        'shape', shape, 'delaypoint', 'median', ...
        'tsec', internal_tsec, 'intervalmass', 0.75, ...
        'intervaltype', intervaltype);

    % convolve that window the segment
    overlap = myconv(seg', h', 'causal', false, 'central_point', zero_tp);
    
    seg_onset = find(outer_tsec>=surround_segs_onsets(si),1,'first');
    oi1 = seg_onset-xi(1);
    oi2 = oi1+length(overlap)-1;
    overlap_surround(oi1:oi2) = overlap_surround(oi1:oi2)+overlap'.^2;
  
end
figure; 
plot(outer_tsec,overlap_surround); hold on; 
plot(outer_tsec,overlap_central)
legend('overlap surround', 'overlap central')
xlabel('time')
box off


    
%% calculate CCC
% ccc_noisefree = overlap_central.^2 ./ (overlap_central.^2 + sum(overlap_surround.^2));
% removing sum & square on overlap surround because already done
ccc_noisefree = (overlap_central.^2) ./ (overlap_central.^2 + overlap_surround);

% sim noise ceiling
rng(1)
a = .05;
b = .5;
noise_cieling = a.*randn(size(ccc_noisefree)) + b;
% correct by noise
ccc_noise = noise_cieling .* ccc_noisefree;


% plot



figure ; 
plot(outer_tsec(outer_tsec>=-.5& outer_tsec<=1) ,ccc_noise(outer_tsec>=-.5& outer_tsec<=1)); 
title('predicted CCC corrected by noise')
xlabel('time')
box off

%% if surrounding segments were equal

overlap_surround2 = zeros(size(outer_tsec));


% onsets of surrounding segments
surround_segs_onsets = [-1.25,-1,-.75,-.5,-.25,.25,.5,.75,1];
surround_seg_durs = [.25,.25,.25,.25,.25,.25,.25,.25,.25];



% for each segment in 
for si=1:length(surround_seg_durs)
    
    segdur_sec =surround_seg_durs(si);
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
        'shape', shape, 'delaypoint', 'median', ...
        'tsec', internal_tsec, 'intervalmass', 0.75, ...
        'intervaltype', intervaltype);

    % convolve that window the segment
    overlap = myconv(seg', h', 'causal', false, 'central_point', zero_tp);
    
    seg_onset = find(outer_tsec>=surround_segs_onsets(si),1,'first');
    oi1 = seg_onset-xi(1);
    oi2 = oi1+length(overlap)-1;
    overlap_surround2(oi1:oi2) = overlap_surround2(oi1:oi2)+overlap'.^2;
  
end

ccc_noisefree2 = (overlap_central.^2) ./ (overlap_central.^2 + overlap_surround2);


% correct by noise
ccc_noise2 = noise_cieling .* ccc_noisefree2;




figure ; 
plot(outer_tsec(outer_tsec>=-.5& outer_tsec<=1) ,ccc_noise(outer_tsec>=-.5& outer_tsec<=1));  hold on;
plot(outer_tsec(outer_tsec>=-.5& outer_tsec<=1) ,ccc_noise2(outer_tsec>=-.5& outer_tsec<=1));  
title('predicted CCC corrected by noise')
xlabel('time')
legend('diff dur surr', 'same dur surr')
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

L.param_string =  'boundary-any_lag_win-0-1';
L.channels = [1,2];
L.chnames = {'1','2'};
L.boundary = 'any';
L.rampwin ='hann';
L.rampdur = 0.0156;
L.unique_segs=[250,250];
L.sr=sr;
L.lag_t= [outer_tsec(outer_tsec>=-.5& outer_tsec<=1)];
L.n_total_segs =[100;100];
L.figure_directory= '/home/magdalena/Documents/projects/CST/test';
L.output_directory = '/home/magdalena/Documents/projects/CST/test';
L.indiv='test';

L.diff_context = [];
L.diff_context(:,1,1)=ccc_noise(outer_tsec>=-.5& outer_tsec<=1);
L.diff_context(:,2,1)=ccc_noise(outer_tsec>=-.5& outer_tsec<=1);
L.diff_context(:,1,2)=ccc_noise2(outer_tsec>=-.5& outer_tsec<=1);
L.diff_context(:,2,2)=ccc_noise2(outer_tsec>=-.5& outer_tsec<=1);

L.same_context = [];
L.same_context =repmat(noise_cieling(outer_tsec>=-.5& outer_tsec<=1)',1,2,2);


L.same_context_err =[];
L.same_context_err(:,1,1)=ccc_noise(outer_tsec>=-.5& outer_tsec<=1).^2;
L.same_context_err(:,2,1)=ccc_noise(outer_tsec>=-.5& outer_tsec<=1).^2;
L.same_context_err(:,1,2)=ccc_noise2(outer_tsec>=-.5& outer_tsec<=1).^2;
L.same_context_err(:,2,2)=ccc_noise2(outer_tsec>=-.5& outer_tsec<=1).^2;

% M = modelfit_cross_context_corr(L, 'overwrite', true, 'plot_figure', true, ...
%     'shape', 1, 'intper_range', [1/128, 0.50], 'delay_range', [0, 0.25],  ...
%     'lossfn', 'sqerr');

  M = modelfit_cross_context_corr(L, 'overwrite', true, 'plot_figure', true, ...
        'shape', [3], 'intper_range', [1/32, 0.25], 'delay_range', [0, 0.25], ...
        'lossfn', 'sqerr', 'nintper', 50, 'boundstrength', [0, 0.25, 0.5, 1, 2]);


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