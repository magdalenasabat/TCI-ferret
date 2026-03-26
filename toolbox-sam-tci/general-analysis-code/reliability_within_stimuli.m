function R = reliability_within_stimuli(D, varargin)

% Wrapper for reliability.m
% 
% Assumes we have a cell array where each cell
% contains the data matrix for a single neuron
% formatted as [samples, rep, stimuli].
% 
% Reliability is calculated as averaged across stimuli.
%
% if chunksize is set permutationas are shuffling timeseries
% if chunksize=0 permutations shuffle association between repetition and
% stimuli
%
% 2024-03-21: Sam NH
% 2024-09-09: Magdalena S. adapted to loop over neurons

[n_smps, n_reps, n_vars] = size(D);

I.splithalf = true;
I.spearbrown = false;
I.nperms = 0;
I.chunksize = 0;
I.removeNaNs = false;
I.maxoverlap = NaN;
I = parse_optInputs_keyvalue(varargin, I);

n_neuron = length(D);
for i = 1:n_neuron
    fprintf('\n# neuron %d\n', i); drawnow;

    R_single_var{i} = reliability(D{i}, 'splithalf', I.splithalf, 'spearbrown', I.spearbrown, ...
        'nperms', I.nperms, 'chunksize', I.chunksize, 'removeNaNs', I.removeNaNs, ...
        'maxoverlap', I.maxoverlap);
end

% average correlation across stimuli
R.corr = zeros(1,n_neuron);
for i = 1:n_neuron
    R.corr(i) = mean( R_single_var{i}.corr,'omitnan');
end

if I.nperms > 0
    % average null distribution across stimuli
    R.null = zeros(size(R_single_var{i}.null,1),n_neuron);
    for i = 1:n_neuron
        R.null(:,i) = mean( R_single_var{i}.null,2,'omitnan');
    end

    [R.logP_gauss, R.z] = sig_via_null_gaussfit(R.corr, R.null);
    [R.logP_counts] = sig_via_null_counts(R.corr, R.null);
    R.logP_gauss(R.logP_gauss>100) = 100;
end

