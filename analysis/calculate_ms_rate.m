function [ms_counts, rate_window_length] = calculate_ms_rate(ms_present, rate_window_length, stepsize)
% calculate microsaccade rate per second
if nargin < 2 % window length not specified
    rate_window_length = 500; % in timepoints (ms at 1000 Hz)
end
if nargin < 3 % step size not specified
    stepsize = 50; % in timepoints (ms at 1000 Hz)
end

%% step = window length
% nwindows = (length(ms_present) - 1) / window_length;
% ntrials = size(ms_present, 2);
% ms_counts_all = zeros(nwindows, ntrials);
% win_idx = 0;
% for win_start = 1:window_length:(length(ms_present)-window_length)
%     win_idx = win_idx + 1;
%     curr_window_idx = [win_start:win_start+window_length];
%     for trial = 1:ntrials
%     curr_window = ms_present(curr_window_idx, trial);
%     if ~isnan(curr_window(1)) && curr_window(1) % ms is ongoing as window begins
%         start_ms = 1;
%     else
%         start_ms = 0;
%     end
%     curr_n_ms = sum(diff(curr_window) == 1) + start_ms;
%     ms_counts_all(win_idx, trial) = curr_n_ms;
%     end
% end
% ms_counts_all = ms_counts_all * (1000 / window_length); % convert to ms/sec
% ms_counts = mean(ms_counts_all, 2);
% end

%% sliding window, variable step size (first try step = 1 ms)
ntrials = size(ms_present, 2);
nwindows = (length(ms_present) - rate_window_length-1) / stepsize + 1;
ms_counts_all = nan(nwindows, ntrials);
win_idx = 0;
for win_start = 1:stepsize:(length(ms_present)-rate_window_length)
    win_idx = win_idx + 1;
    curr_window_idx = [win_start:win_start+rate_window_length];
    %for trial = 1:ntrials
    curr_window = ms_present(curr_window_idx);
    if sum(isnan(curr_window)) < rate_window_length/2 % if less than half of timepoints are nans, use the window
    if ~isnan(curr_window(1)) && curr_window(1) % ms is ongoing as window begins
        start_ms = 1;
    else
        start_ms = 0;
    end
    curr_n_ms = sum(diff(curr_window) == 1) + start_ms;
    else % otherwise if too many nans, skip this window
        curr_n_ms = nan;
    end
        ms_counts_all(win_idx) = curr_n_ms;

end
    ms_counts = ms_counts_all * (1000 / rate_window_length); % convert to ms/sec
%ms_counts = mean(ms_counts_all, 2);
end