function [h_relpower, t_sec, h, causal] = win_power_ratio_magda(...
    segdur_sec, distr, intper_sec, delay_sec, varargin)

% Calculates the expected correlation for TCI correlation analyses by
% measuring the variance a given segment should contribute relative to the
% other segments.
% 
% 2019-13-11: Created, Sam NH

% window and sampling rate
% note the window specifies the window of the output
% which can / should be different from the window used for convolution
I.win = [];
I.sr = NaN;

% alternatively you can just specify a vector of timepoints
% again this is the vector of timepoints of the output
% which can / should be different from the time-points used for convolution
I.tsec = [];

% this is the target sampling rate used for internal computations
I.target_sr = 100/(min([intper_sec, segdur_sec]));

% parameters of the window see modelwin.m
I.shape = 1;
I.delaypoint = 'peak';
I.forcecausal = false;
I.intervalmass = 0.75;
I.intervaltype = 'center';

% window applied to beginning / end of segment
% see winoverlap.m
I.rampwin = 'none';
I.rampdur = 0;

% correlation factor
I.corrfac = 0;

% boundary parameters
I.boundstrength = 0;
I.boundexp = 1;

% can optionally plot result
I.plot = false;

% added by Magda, 
I.surrounding_segments = [NaN;segdur_sec];
I.n_total_seg = [1];
% overwrite default parameters with those specified
[I,C] = parse_optInputs_keyvalue(varargin, I);

%%  compute window

% ensure there are an integer number of samples in the segment
internal_sr = ceil(segdur_sec * I.target_sr) / segdur_sec;
segdur_smps = checkint(segdur_sec * internal_sr);

% set time so that we can fit the longest surrounding segments -  Magda
time_boundary = [-1,1]*(3*max(I.surrounding_segments(2,:)) + 4*intper_sec + abs(delay_sec));
internal_tsec = time_boundary(1):1/internal_sr:time_boundary(2);

% measure overlap
[h, t, causal] = winoverlap(segdur_sec, distr, intper_sec, delay_sec, ...
     'shape', I.shape, 'forcecausal', I.forcecausal, 'delaypoint', I.delaypoint, ...
     'plot', false,'internal_sr', internal_sr, 'rampwin', I.rampwin, 'rampdur', I.rampdur, ...
     'intervalmass', I.intervalmass, 'intervaltype', I.intervaltype);
h(h<0) = 0;

% % interpolate h at the same sr as the central segment -  Magda
% h = interp1(t, h, internal_tsec,'nearest',0)';
% t = internal_tsec;

%added by Magda to account for different durations of surrounding segments

if size(I.surrounding_segments,2)==1 &  I.surrounding_segments(2,:)== segdur_sec
    % create delayed copies
    N_delays = ceil(length(h)/segdur_smps);
    shifts = (-N_delays:N_delays)*segdur_smps;
    h_delayed = add_delays(h, shifts);
else
%     mask_time = (I.surrounding_segments(1,:) >= (t(1)-3*max(I.surrounding_segments(2,:)))) & ...
%                 (I.surrounding_segments(1,:) <= (t(end)+ 3*max(I.surrounding_segments(2,:)) ));
    h_delayed=zeros(numel(internal_tsec),size(I.surrounding_segments(2,:),2));
    % for each unique segment
  
    for segdur_sec_this = unique(I.surrounding_segments(2,:))
        mask = I.surrounding_segments(2,:)==segdur_sec_this ;
%         if sum(mask)==0; continue; end
        % create H
        % ensure there are an integer number of samples in the segment
        % internal_sr % reuse

        % measure overlap
        [h_this,t_this, causal_this] = winoverlap(segdur_sec_this, distr, intper_sec, delay_sec, ...
             'shape', I.shape, 'forcecausal', I.forcecausal, 'delaypoint', I.delaypoint, ...
             'plot', false,'internal_sr', internal_sr, 'rampwin', I.rampwin, 'rampdur', I.rampdur, ...
             'tsec',internal_tsec,'intervalmass', I.intervalmass, 'intervaltype', I.intervaltype);
        h_this(h_this<=0) = 0.0;
        
%         % interpolate h at the same sr as the central segment
%         h_this = interp1(t_this, h_this, t,'nearest',0)';
        
        % add delays
        shifts =  floor(I.surrounding_segments(1,mask) * internal_sr);
        
%         h_delayed(:,mask) = add_delays(h_this, shifts);
        h_delayed(:,mask) = add_delays(h_this', shifts);
     
    end

    h_delayed = interp1( internal_tsec,h_delayed,t,'nearest',0);
    
end

%% Boundary effects

% total overlap of adjacent segments
adjacent_sum = h_delayed(:,1:end-1) + h_delayed(:,2:end);

% fraction of that total pairwise overlap
% due to the first segment
split_frac = h_delayed(:,1:end-1) ./ adjacent_sum;
split_frac(adjacent_sum < 1e-5) = 0;

% map split through raised cosine so that 0.5 is 1 and 0/1 are 0
boundary_effect = (0.5-cos(split_frac*2*pi)/2).^I.boundexp;
h_boundary = boundary_effect .* adjacent_sum;

% figure;
% plot([h_boundary, sum(h_boundary,2)]);

% impulse response corresponding to overlap with boundary
% irf = modelwin(distr, intper_sec, delay_sec, ...
%      'shape', I.shape, 'forcecausal', I.forcecausal, 'delaypoint', I.delaypoint, ...
%      'plot', false,'tsec', t, ...
%      'intervalmass', I.intervalmass, 'intervaltype', I.intervaltype);

% irf_delayed = add_delays(irf'/max(irf), shifts);

% delayed copies
% plot(t, [h, irf_h'/max(irf_h)])

%% Calculate predicted correlation


% no correlation
% denom = sum(h_delayed.^2,2) + 11*sum(irf_delayed.^2,2);
denom = sum(h_delayed.^2,2) + I.boundstrength * sum(h_boundary.^2,2);
h_nocorr = bsxfun(@times, h.^2, 1./denom);

% perfect correlation
h_perfcorr = h .* sum(h_delayed,2) ./ sqrt(sum(h_delayed.^2,2) .* sum(h_delayed,2).^2);
h_relpower = I.corrfac * h_perfcorr + (1-I.corrfac) * h_nocorr;
% plot(t, h_relpower)


%% Interpolate

if ~C.sr && ~C.win && ~C.tsec
    t_sec = t;
else
    % set time vectory or window depending upon input
    if C.tsec % if time-vector is specified, adjust window/sampling rate
        t_sec = I.tsec;
        I.win = I.tsec([1 end]);
        I.sr = 1/diff(I.tsec(1:2));
        assert(all(abs(diff(I.tsec)-diff(I.tsec(1:2)))<1e-6));
    else % if time-vector is not specified, determine based on window/sampling rate
        if isempty(I.win)
            I.win = t([1,end]);
        end
        t_sec = (round(I.win(1)*I.sr):round(I.win(2)*I.sr))/I.sr;
    end
    if any(isnan(h_relpower(:))) || any(isnan(h(:))) || any(isnan(h_delayed(:)))
        error('NaN values');
    end
    h_relpower = interp1(t, h_relpower, t_sec(:), 'pchip');
    h = interp1(t, h, t_sec(:), 'pchip');
    xi = t_sec<t(1) | t_sec>t(end);
    h_relpower(xi) = 0;
    h(xi) = 0;
end

%% Plot

if I.plot
    figure;
    subplot(2,1,1);
    ph = plot(t_sec, [h, h_relpower], 'LineWidth', 2);
    xlim(t_sec([1,end]));
    legend(ph, {'Orig', 'Rel Power'}, 'Location', 'Best');
    subplot(2,1,2);
    plot(t, h_delayed, 'LineWidth', 2);
    hold on;
    plot(t, sum(h_delayed,2), 'k--', 'LineWidth', 2);
    xlim(t([1,end]));
    title('All win');
end