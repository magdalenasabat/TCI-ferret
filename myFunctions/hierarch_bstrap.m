function stat_bstrap = hierarch_bstrap(D, subject_indices, statistic, n_bstrap, varargin)

% Computes a sample statistic (e.g., mean) using hierarchical bootstrapping
% first sampling subjects with replacement and then sampling within subject
% datapoints with replacement. Statistic is computed on data concatenated
% across all samples
% 
% -- Arguments -- 
% 
% D: [sample, feature] where the feature axis could be anything (e.g., time, electrode)
% 
% subject_indices: vector with integers or string labels indicating the subject each sample
%
% statistic: 'mean' or 'median' or handle to custom function (assumed to operate over the first dimension)
% 
% n_bstrap: number of bootstrapped samples
% 
% -- Example with 10 data points from 5 unique subjects --
% 
% data = randn(10, 3);
% subject_indices = {'S1', 'S1', 'S2', 'S2', 'S2', 'S4', 'S4', 'S7', 'S7', 'S7'};
% n_bstrap = 100;
% y = hierarch_bstrap(data, subject_indices, 'median', n_bstrap);
% median(y)
% std(y)

if strcmp(class(statistic), 'char')
	switch statistic
		case 'mean'
			fn = @(X)mean(X,1);
		case 'median'
			fn = @(X)median(X,1);
		otherwise
			error('No valid statistic name');
	end
else
	fn = statistic;
end

[n_samples, n_features] = size(D);
assert(n_samples==length(subject_indices))

% turn subject indices into column vector
% sort the indices so that they go from 1 to the number of subjects
% e.g., [2, 2, 2, 4, 4, 5] -> [1, 1, 1, 2, 2, 3]
subject_indices = subject_indices(:);
[~,~,subject_indices] = unique(subject_indices);
n_subjects = max(subject_indices);

% number of electrodes per subject
n_elec_per_subject = nan(1, n_subjects);
for i = 1:n_subjects
	n_elec_per_subject(i) = sum(subject_indices==i);
end

ResetRandStream2(1)

% bootstrapped statistic
stat_bstrap = nan(n_bstrap, n_features);
for i = 1:n_bstrap

	% sample subjects with replacement
	subject_samples = randi(n_subjects, [n_subjects, 1]);

	% sample datapoints within this subject with replacement
	% concatenate samples across subjects
	D_sample = [];
	for j = 1:n_subjects
		s = subject_samples(j);
		D_within_subject = D(subject_indices==s,:);
		within_subject_samples = randi(n_elec_per_subject(s), [n_elec_per_subject(s), 1]);
		D_sample = cat(1, D_sample, D_within_subject(within_subject_samples,:));
	end

	stat_sample = fn(D_sample);
	if i == 1
		stat_bstrap = nan(n_bstrap, length(stat_sample));
	end
	stat_bstrap(i,:) = stat_sample;
end