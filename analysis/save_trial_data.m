
% store all min rate data in .csv table.
%nonan_idx = ~isnan(rts) & ~isnan(nms_baseline) & ~isnan(nms_window) & ~isnan(immobilization_durations);
nonan_idx = ~isnan(rts);% & subjs ~= 11 & subjs ~= 12;
%nonan_idx = 1:numel(rts);
snums = subjs(nonan_idx);
rt_nonan = rts(nonan_idx) ./ 1000;
min_nms_nonan = nms_window(nonan_idx) ./ 700 .* 1000;
min_slip_nonan = min_slip_window(nonan_idx);
try
baseline_slip_nonan = baseline_retinal_slip_trials(nonan_idx);
catch
    baseline_retinal_slip_trials = nan(size(min_slip_window));
    baseline_slip_nonan = baseline_retinal_slip_trials(nonan_idx);
end
baseline_nms_nonan = nms_baseline(nonan_idx);% .* 1000;
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

%in_baseline_nms_nonan = in_nms_baseline(nonan_idx);
%out_baseline_nms_nonan = out_nms_baseline(nonan_idx);


%yrt = rt_nonan;
%min_amps = baseline_slip_nonan;

tabledata = [snums, rt_nonan, min_nms_nonan, baseline_nms_nonan, meanmsampl_baseline_nonan, min_slip_nonan, baseline_slip_nonan, immobs_nonan, immob_ms_amps_nonan, last_ms_times_nonan, last_ms_amps_nonan, clrs, szs, ms_exist_nonan, early_ms_exist_nonan, late_ms_exist_nonan, blink_exist_nonan, nblinks_nonan, session_num_nonan, trial_counter_nonan];%, in_baseline_nms_nonan, out_baseline_nms_nonan];
T = array2table(tabledata, 'VariableNames', {'subject', 'RT', 'minnms', 'baselinenms', 'baselinemeanmsampl', 'minretinalslip', 'baselineretinalslip', 'immobduration', 'immobmsamplitude', 'lastmstime', 'lastmsamplitude', 'contrast', 'eccentricity', 'mspresent', 'earlymspresent', 'latemspresent', 'blinkpresent', 'nblinks', 'session', 'trial'});%, 'inbaselinenms', 'outbaselinenms'});
writetable(T, 'bhv_data/trial_data_withblinks.csv');
%writetable(T, 'trial_data_fiximmob.csv');


% normalize by each subj mean
for sj = 1:nsubjects
    subj_idx = snums == sj;
    rt_nonan(subj_idx) = rt_nonan(subj_idx) ./ nanmean(rt_nonan(subj_idx));
    %min_nms_nonan(subj_idx) = baseline_nms_nonan(subj_idx) - min_nms_nonan(subj_idx);
    min_nms_nonan(subj_idx) = min_nms_nonan(subj_idx) ./ nanmean(baseline_nms_nonan(subj_idx));
    min_slip_nonan(subj_idx) = min_slip_nonan(subj_idx) ./ nanmean(baseline_slip_nonan(subj_idx));
    baseline_slip_nonan(subj_idx) = baseline_slip_nonan(subj_idx) ./ nanmean(baseline_slip_nonan(subj_idx));
    baseline_nms_nonan(subj_idx) = baseline_nms_nonan(subj_idx) ./ nanmean(baseline_nms_nonan(subj_idx));
end

% normalize by each block mean (baseline rate & rt)
for sj = 1:nsubjects
    for cl = 1:3
        for sz = 2:2:6
            block_idx = snums == sj & clrs == cl & szs == sz;
            %rt_nonan(block_idx) = rt_nonan(block_idx) ./ nanmean(rt_nonan(block_idx));
            %min_nms_nonan(block_idx) = baseline_nms_nonan(block_idx) - min_nms_nonan(block_idx);
            %baseline_nms_nonan(block_idx) = baseline_nms_nonan(block_idx) ./ nanmean(baseline_nms_nonan(block_idx));
            baseline_nms_nonan(block_idx) = baseline_nms_nonan(block_idx) - nanmean(baseline_nms_nonan(block_idx)); % zero-mean
            min_nms_nonan(block_idx) = min_nms_nonan(block_idx) - nanmean(min_nms_nonan(block_idx)); % zero-mean
            baseline_slip_nonan(block_idx) = baseline_slip_nonan(block_idx) - nanmean(baseline_slip_nonan(block_idx));
        end
    end
end

% remove trials with zero baseline microsaccades?
nozero_idx = baseline_nms_nonan > 0;% & min_nms_nonan > 0;
snums = snums(nozero_idx);
rt_nonan = rt_nonan(nozero_idx);
min_nms_nonan = min_nms_nonan(nozero_idx);
min_slip_nonan = min_slip_nonan(nozero_idx);
baseline_slip_nonan = baseline_slip_nonan(nozero_idx);
baseline_nms_nonan = baseline_nms_nonan(nozero_idx);
clrs = clrs(nozero_idx);
szs = szs(nozero_idx);



figure; hold on;
for i = clr_idx
for sz = sz_idx
scatter(yrt(clrs == i & szs == sz), min_amps(clrs == i & szs == sz))
end
end
figure; hold on;
for i = clr_idx
scatter(yrt(clrs == i), min_amps(clrs == i));
end
figure; hold on;
for sz = sz_idx
scatter(yrt(szs == sz), min_amps(szs == sz));
end


% linear regression
x = baseline_nms_nonan;
y = rt_nonan;
X = [ones(length(x), 1) x];
b = X \ y
ypred = X*b;

figure;
scatter(x,y);
hold on;
plot(x, ypred);


%% raster plots

    ntps_before = 5000; ntps_after = 1000;
    useclrs = {[1 0 0], [0 1 0], [0 0 1]}; % optional colors

which2plot = 'color'; % size, color, or all (black)
shift_subjs = 1; % shift each subject's plot by estimated rt?
subj_shifts = [ -125  -251   -26  -251   -65  -251   249  -251    50  -251    78   249   242  -251   249   -53   137   249  -170     4  -251]; % after 500ms smoothing, bw 300-800ms


figure; hold on; all_subs = 1; plot_iter = 0; % if plotting all subjs together

for subj = 1:nsubjects
subj_idx = snums == subj;
subj_ms_present = ms_present_concat(:, subj_idx);
subj_ntrials = size(subj_ms_present, 2);
subj_clrs = clrs(subj_idx);
subj_szs = szs(subj_idx);
ms_onsets_concat = nan(size(subj_ms_present));
% extract onsets only
for trial = 1:subj_ntrials
    ms_onsets_concat(:, trial) = [0; diff(subj_ms_present(:, trial))];
end
ms_onsets_concat(ms_onsets_concat < 0) = 0; % get rid of offsets
%ms_onsets_concat(isnan(ms_onsets_concat)) = 0; % get rid of nans.
ms_onsets_concat(ms_onsets_concat == 0) = nan; % get rid of 0s
% for plotting size:
    conv_subj_szs = abs(subj_szs / 2 - 4); % 1 = largest, 3 = smallest
    [sort_szs, sort_order] = sort(conv_subj_szs);

if ~all_subs % if plotting each subject separately
    figure; hold on;
end
for trial = 1:subj_ntrials
    switch which2plot
        case 'all'
            curr_color = [0 0 0];
            spikes = ms_onsets_concat(:, trial);
            lgnd = '';
        case 'color'
            curr_color = useclrs{subj_clrs(trial)};
            spikes = ms_onsets_concat(:, trial);
        lgnd = {'low', 'medium', 'high'};
        case 'size'
            curr_color = useclrs{sort_szs(trial)};
            spikes = ms_onsets_concat(:, sort_order == trial);
                    lgnd = {'large', 'medium', 'small'};

    end
    if shift_subjs
        if subj_shifts(subj) > 0 % shift rightward
                        spikes = [nan(subj_shifts(subj), 1); spikes(1:end-subj_shifts(subj))];
                    elseif subj_shifts(subj) < 0 % shift leftward
                        spikes = [spikes((abs(subj_shifts(subj))+1):end); nan(abs(subj_shifts(subj)), 1)];
        end
    end
    if all_subs % plotting all subjects together
        plot_iter = plot_iter + 1;
    else % plotting each subject separately
        plot_iter = trial;
    end
    plot(-ntps_before:ntps_after, plot_iter*spikes, '.', 'Color', curr_color, 'MarkerSize', 5);
end
curr_ylim = get(gca, 'ylim');
hold on; zeroline = plot([0 0], curr_ylim, '-k');
set(get(get(zeroline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
xlabel('time');
%legend(lgnd);
set(gca, 'FontSize', 30);
if ~all_subs
    title(['subj ', num2str(subjects(subj)), ', sorted by ', which2plot]);
else
    title(['all subjects, sorted by ', which2plot]);
end
end



%% raster plot by a certain condition, but shuffling trials randomly across subjects
which2plot = 'size';

% for plotting size:
    conv_szs = abs(szs / 2 - 4); % 1 = largest, 3 = smallest
    [sort_szs, sort_order] = sort(conv_szs);
    
    plot_iter = 0;
figure; hold on;
for cond = 1:3
    switch which2plot
        case 'color'
    cond_idx = Shuffle(find(clrs == cond));
        %cond_idx = find(clrs == cond);

        case 'size'
                    cond_idx = Shuffle(find(sort_szs == cond));
        %cond_idx = find(sort_szs == cond);
    end
cond_ms_present = ms_present_concat(:, cond_idx);
cond_ntrials = size(cond_ms_present, 2);
ms_onsets_concat = nan(size(cond_ms_present));
% get ms rate
nspikes = cross_trial_ms_rate(cond_ms_present, 'alpha');
% extract onsets only
for trial = 1:cond_ntrials
    ms_onsets_concat(:, trial) = [0; diff(cond_ms_present(:, trial))];
end
ms_onsets_concat(ms_onsets_concat < 0) = 0; % get rid of offsets
ms_onsets_concat(ms_onsets_concat == 0) = nan; % get rid of 0s

for trial = 1:cond_ntrials
    curr_color = useclrs{cond};
    spikes = ms_onsets_concat(:, trial);
    if shift_subjs
        subj = snums(cond_idx(trial));
        if subj_shifts(subj) > 0 % shift rightward
                        spikes = [nan(subj_shifts(subj), 1); spikes(1:end-subj_shifts(subj))];
                    elseif subj_shifts(subj) < 0 % shift leftward
                        spikes = [spikes((abs(subj_shifts(subj))+1):end); nan(abs(subj_shifts(subj)), 1)];
        end
    end
        plot_iter = plot_iter + 1;
    plot(-ntps_before:ntps_after, plot_iter*spikes, '.', 'Color', curr_color, 'MarkerSize', 5);
end
% plot nspikes/t within each condition
curr_xlim = get(gca, 'xlim');
%nspikes = cross_trial_ms_rate(ms_onsets_concat, 'alpha', 1);
plot(-ntps_before:ntps_after, nspikes*2000+plot_iter-3000, 'Color', curr_color, 'LineWidth', 3);
%nspikes = nansum(ms_onsets_concat, 2);
%plot(-ntps_before:ntps_after, smoothdata(nspikes, 'movmean', 200)*1000+plot_iter-3000, 'LineWidth', 2);

end
curr_ylim = get(gca, 'ylim');
hold on; zeroline = plot([0 0], curr_ylim, '-k');
set(get(get(zeroline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
xlabel('time');
%legend(lgnd);
set(gca, 'FontSize', 30);
    title(['all subjects, sorted by ', which2plot]);
    