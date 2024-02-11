function [msrate, msrateshuffled] = cross_trial_ms_rate(block_ms_present, smoothmethod, onset)%rate_window_length, stepsize)
% cross-trial microsaccade rate
% block_ms_present: ntimepoints x ntrials ms present
if nargin < 2
    smoothmethod = 'alpha'; % default smoothing is alpha weighting.
end
if nargin < 3 % input type not specified (ongoing ms present or ms onset indices)
    onset = 0; % default: ongoing ms present
end

ntps = size(block_ms_present, 1);
ntrials = size(block_ms_present, 2);

if ~onset % derive ms onsets
    block_ms_onsets = nan(size(block_ms_present));
    % extract onsets only
    for trial = 1:ntrials
        block_ms_onsets(:, trial) = [0; diff(block_ms_present(:, trial))];
    end
    block_ms_onsets(block_ms_onsets < 0) = 0; % get rid of offsets
elseif onset
    block_ms_onsets = block_ms_present;
end


%% generate shuffled ms onsets per trial
shuffled_block_ms_onsets = nan(size(block_ms_onsets));
for trial = 1:ntrials
    curr_ms_onsets = block_ms_onsets(:, trial);
    usable_tp_idx = ~isnan(curr_ms_onsets);
    usable_timepoints = curr_ms_onsets(usable_tp_idx);
    shuffled_timepoints = usable_timepoints(randperm(length(usable_timepoints)));
    shuffled_block_ms_onsets(usable_tp_idx, trial) = shuffled_timepoints;
end

%% windowing/weighting is done via causal alpha function.
if ~strcmp(smoothmethod, 'amplitude') % most cases
    msrate = smooth_rate(block_ms_onsets, smoothmethod);
    msrateshuffled = smooth_rate(shuffled_block_ms_onsets, smoothmethod);
else % if using special method for computing ms amplitude across time
    amplitudes = nanmean(block_ms_present, 2);
    window_length = 100; stepsize = 100;
    nwindows = (length(amplitudes) - window_length-1) / stepsize + 1;
    r = nan(nwindows, 1);
    win_idx = 0;
    for win_start = 1:stepsize:(length(amplitudes)-window_length)
        win_idx = win_idx + 1;
        curr_window_idx = [win_start:win_start+(window_length-1)];
        curr_window = amplitudes(curr_window_idx);
        curr_amp_sum = nansum(curr_window);
        curr_amp = curr_amp_sum / sum(~isnan(curr_window)); % divide by # of timepoints used
        r(win_idx) = curr_amp;
    end
end

%% rate function.
    function r = smooth_rate(block_onsets, smoothmethod)
        ntrialspertimepoint = sum(~isnan(block_onsets), 2);

        % combine across trials
        onsets = nansum(block_onsets, 2);

        switch smoothmethod
            case 'alpha'
                % causal alpha kernel
                winsz = 1001; % window length for window function. (essentially inf)
                alph = 1/100;
                taus = 1:winsz;
                t = median(taus);
                taus = (taus - t) * -1;
                taus(taus<0) = 0; % prevent impact of tau > time t
                w = zeros(size(taus));
                w = alph^2 .* taus .* exp(-alph .* taus);
                w = w(:) / sum(w);
                r = convn(onsets, w, 'same');

                % convert to ms/sec
                r = r * 1000 ./ ntrialspertimepoint;

                % account for decay / delay (shift r by 1/alpha ms)
                % and remove start/end transients
                dlay = 1/alph;
                r(1:dlay) = nan;
                r(end-dlay:end) = nan;
                r = [nan(dlay, 1); r(1:end-dlay)];

            case 'savitzky-golay'
                r = sgolayfilt(onsets, 1, 101); % Savitzky-Golay filter
                r = r * (1000 / rate_window_length) ./ ntrialspertimepoint;
            case 'gaussian'
                r = smoothdata(onsets, 'gaussian', 100);
                r = r * (1000 / rate_window_length) ./ ntrialspertimepoint;
            case 'window' % just a sliding window, take the average
                r = smoothdata(onsets, 'movmean', 100);
                r = r * (1000 / rate_window_length) ./ ntrialspertimepoint;
        end
        r(isinf(r)) = nan; % ensure timepoints with no usable trials are set to nan.

    end
end