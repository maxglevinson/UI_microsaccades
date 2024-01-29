
% store all trial data in .csv table, for use in R.
% run this after doing all the processing in full_eyelink_analysis.m
nonan_idx = ~isnan(rts);
snums = subjs(nonan_idx);
rt_nonan = rts(nonan_idx) ./ 1000;
try % if we computed retinal slip (sometimes not because it takes a long time)
baseline_slip_nonan = baseline_retinal_slip_trials(nonan_idx);
catch
    baseline_retinal_slip_trials = nan(size(min_slip_window));
    baseline_slip_nonan = baseline_retinal_slip_trials(nonan_idx);
end
baseline_nms_nonan = nms_baseline(nonan_idx);
clrs = contrasts(nonan_idx);
szs = eccentricities(nonan_idx);
immobs_nonan = immobilization_durations(nonan_idx) ./ 1000 .* -1; % make positive and in seconds
nblinks_nonan = nblinks_baseline(nonan_idx); % basically 0 or 1 blink always!!
blinkrate_nonan = trial_blinkrate_baseline(nonan_idx);
ms_exist_nonan = ms_exist(nonan_idx);
early_ms_exist_nonan = early_ms_exist(nonan_idx);
late_ms_exist_nonan = late_ms_exist(nonan_idx);
blink_exist_nonan = blink_exist(nonan_idx);
session_num_nonan = session_nums(nonan_idx);
immob_ms_amps_nonan = immob_ms_amplitudes(nonan_idx);
last_ms_times_nonan = last_ms_times(nonan_idx) ./ 1000 .* -1;
last_ms_amps_nonan = last_ms_amplitudes(nonan_idx);
trial_counter_nonan = trial_counter(nonan_idx);
meanmsampl_baseline_nonan = meanmsampl_baseline(nonan_idx);

tabledata = [snums, rt_nonan, baseline_nms_nonan, meanmsampl_baseline_nonan, baseline_slip_nonan, immobs_nonan, immob_ms_amps_nonan, last_ms_times_nonan, last_ms_amps_nonan, clrs, szs, ms_exist_nonan, early_ms_exist_nonan, late_ms_exist_nonan, blink_exist_nonan, nblinks_nonan, session_num_nonan, trial_counter_nonan];
T = array2table(tabledata, 'VariableNames', {'subject', 'RT', 'baselinenms', 'baselinemeanmsampl', 'baselineretinalslip', 'immobduration', 'immobmsamplitude', 'lastmstime', 'lastmsamplitude', 'contrast', 'eccentricity', 'mspresent', 'earlymspresent', 'latemspresent', 'blinkpresent', 'nblinks', 'session', 'trial'});
writetable(T, 'bhv_data/trial_data_withblinks.csv');
