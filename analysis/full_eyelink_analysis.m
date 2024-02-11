%% group analysis for eyelink data
addpath(genpath(pwd));
subjects = [201, 202, 203, 204, 205, 206, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 222, 223, 224];

nsubjects = numel(subjects);
nblocks = 9;
ms_present_all_unsorted = cell(1, nsubjects);
BinocularSaccades_all_unsorted = cell(1, nsubjects);
trial_indices_all_unsorted = cell(1, nsubjects);
posdata_all_unsorted = cell(1, nsubjects);
trial_counter_all_unsorted = cell(1, nsubjects);
extra = cell(nsubjects, 3); % anything extra I want to output and store to check
for subj = 1:nsubjects
    subj_id = subjects(subj);
    for session_id = 1:3
        [curr_ms_present, fps, curr_BinocularSaccades, curr_trial_indices, curr_posdata, curr_trial_counter, extra{subj, session_id}] = analyze_eyelink_fulltrial(subj_id, session_id);
        ms_present_all_unsorted{subj} = [ms_present_all_unsorted{subj}, curr_ms_present];
        BinocularSaccades_all_unsorted{subj} = [BinocularSaccades_all_unsorted{subj}, curr_BinocularSaccades];
        trial_indices_all_unsorted{subj} = [trial_indices_all_unsorted{subj}, curr_trial_indices];
        posdata_all_unsorted{subj} = [posdata_all_unsorted{subj}, curr_posdata]; % t, x pos, y pos
        trial_counter_all_unsorted{subj} = [trial_counter_all_unsorted{subj}, curr_trial_counter];
    end
end


% set order to: sizes  2 4 6 2 4 6 2 4 6
%               colors 1 1 1 2 2 2 3 3 3
% Every subject originally has a randomized block order. Here we align them
% all.
blockorder = [4 8 2 3 1 9 5 7 6; ... % subj 201
    9 1 4 6 8 5 3 7 2; ... % subj 202
    2 5 8 9 4 7 6 1 3; ... % subj 203
    6 8 3 7 5 4 2 1 9; ... % subj 204
    8 7 1 2 4 3 9 5 6; ... % subj 205
    6 8 3 7 5 4 2 1 9; ... %subj 206
    6 8 3 7 5 4 2 1 9; ... % subj 208
    7 6 8 3 1 9 2 5 4; ... % subj 209
    6 8 2 7 5 4 1 3 9; ... % subj 210
    8 6 3 5 7 4 1 2 9; ... % subj 211
    7 8 1 4 2 3 5 9 6; ... % subj 212
    1 6 3 7 9 5 8 2 4; ... % subj 213
    8 6 2 1 7 3 5 4 9; ... % subj 214
    1 6 7 5 3 8 4 9 2; ... % subj 215
    8 4 7 9 3 6 2 5 1; ... % subj 216
    1 3 6 2 5 9 7 4 8; ... % subj 217
    2 8 3 9 4 7 1 6 5; ... % subj 218
    7 1 3 2 9 5 8 4 6; ... % subj 219
    6 3 1 9 4 8 7 2 5; ... % subj 222
    4 1 5 6 7 9 8 2 3; ... % subj 223
    9 1 7 4 2 8 6 5 3]; % subj 224

ntrials_all = zeros(nsubjects, nblocks); % store how many trials are used in each block / subject

% for every millisecond, get the amplitude and direction of any ongoing microsaccade (otherwise NaN)
ms_amplitudes_all_unsorted = repmat({cell(1, nblocks)}, 1, nsubjects);
% directions: 2 = left inwards, 1 = left outwards, -1 = right outwards, -2 = right inwards
ms_directions_all_unsorted = repmat({cell(1, nblocks)}, 1, nsubjects);
for subj = 1:nsubjects
    for block = 1:nblocks
        nblocktrials = size(ms_present_all_unsorted{subj}{block}, 2);
        ms_amplitudes_all_unsorted{subj}{block} = cell(1, nblocktrials);
        ms_directions_all_unsorted{subj}{block} = cell(1, nblocktrials);
        for trial = 1:nblocktrials
            ms_amplitudes_all_unsorted{subj}{block}{trial} = nan(size(ms_present_all_unsorted{subj}{block}{trial}));
            ms_directions_all_unsorted{subj}{block}{trial} = nan(size(ms_present_all_unsorted{subj}{block}{trial}));
            ms_present_indices = find(ms_present_all_unsorted{subj}{block}{trial});
            ntrialms = numel(BinocularSaccades_all_unsorted{subj}{block}{trial}.Start);
            for ms = ms_present_indices
                for nms = 1:ntrialms
                    curr_ms_indices = BinocularSaccades_all_unsorted{subj}{block}{trial}.Start(nms):(BinocularSaccades_all_unsorted{subj}{block}{trial}.End(nms));
                    if ismember(ms, curr_ms_indices) % if current timepoint is within this microsaccade
                        ms_amplitudes_all_unsorted{subj}{block}{trial}(ms) = BinocularSaccades_all_unsorted{subj}{block}{trial}.Amplitude(nms); % assign tp with microsaccade's amplitude
                        ms_directions_all_unsorted{subj}{block}{trial}(ms) = BinocularSaccades_all_unsorted{subj}{block}{trial}.Direction(nms); % assign tp with microsaccade's direction
                    end
                end
            end
        end
    end

    % reorder blocks
    ms_present_all{subj} = ms_present_all_unsorted{subj}(blockorder(subj, :));
    BinocularSaccades_all{subj} = BinocularSaccades_all_unsorted{subj}(blockorder(subj, :));
    trial_indices_all{subj} = trial_indices_all_unsorted{subj}(blockorder(subj, :));
    ms_amplitudes_all{subj} = ms_amplitudes_all_unsorted{subj}(blockorder(subj, :));
    ms_directions_all{subj} = ms_directions_all_unsorted{subj}(blockorder(subj, :));
    posdata_all{subj} = posdata_all_unsorted{subj}(blockorder(subj, :));
    trial_counter_all{subj} = trial_counter_all_unsorted{subj}(blockorder(subj, :));
    for block = 1:nblocks
        ntrials_all(subj, block) = size(ms_present_all{subj}{block}, 2);
    end
end
ntrialstotal_byblock = sum(ntrials_all);

% convert position data to degrees of visual angle?
% calculate deg per pix from screen size (cm), distance, resolution
horzcm = 60.8;
distcm = 60;
pxbycm_x = 2560 / horzcm;
pixperdeg = pi * 2560 / atan(horzcm/distcm/2) / 360; % pixels per degree of visual angle
degperpix = 1 / pixperdeg;

% make blinks (missing data, nan) = 1
blink_present_all = ms_present_all;
for subj = 1:nsubjects
    for block = 1:nblocks
        nblocktrials = size(blink_present_all{subj}{block}, 2);
        for trial = 1:nblocktrials
            ms_present_one = find(ms_present_all{subj}{block}{trial} == 1);
            ms_present_nan = find(isnan(posdata_all{subj}{block}{trial}(2, :)));
            blink_present_all{subj}{block}{trial}(ms_present_one) = 0;
            blink_present_all{subj}{block}{trial}(ms_present_nan) = 1;
        end
    end
end

% get total # excluded trials because of blinks or saccades
badblinktrials = 0;
badsaccadetrials = 0;
for subj = 1:nsubjects
    for sess = 1:3
        badblinktrials = badblinktrials + extra{subj, sess}.excludedblinktrials;
        badsaccadetrials = badsaccadetrials + extra{subj, sess}.excludedsaccadetrials;
    end
end

filename = 'full_eye_data.mat'; % exclude trials with blinks <2s before button press. Use these data for plotting microsaccade rates.
%filename = 'full_eye_data_withblinks.mat'; % for generating trial-wise data to use in linear modeling.
%filename = 'full_eye_data_replay_withblinks.mat'; % same but for replay catch trials.
%filename = 'full_eye_data_sharp_withblinks.mat'; % same but for sharp catch trials.

save(filename, 'subjects', 'nsubjects', 'ntrials_all', ...
    'ntrialstotal_byblock', 'ms_present_all', 'blink_present_all', ...
    'BinocularSaccades_all', 'posdata_all', 'trial_indices_all', ...
    'ms_amplitudes_all', 'ms_directions_all', 'trial_counter_all', ...
    'degperpix', 'fps', 'extra', 'nblocks', 'blockorder', '-v7.3');

%% just load the saved data
datadir = './'; % replace with local directory of experiment code
%filename = [datadir, 'full_eye_data.mat']; % exclude trials with blinks <2s before button press. Use these data for plotting microsaccade rates.
filename = [datadir, 'full_eye_data_withblinks.mat']; % for generating trial-wise data to use in linear modeling.
%filename = [datadir, 'full_eye_data_replay_withblinks.mat']; % same but for replay catch trials.
%filename = [datadir, 'full_eye_data_sharp_withblinks.mat']; % same but for sharp catch trials.
load(filename);

%% epoch trials into block averages
% Run this section, then run the next plot section (for Figure 2b).
% This section also produces .csv file to use in linear modeling and
% making related figures.

rts = []; subjs = []; contrasts = []; eccentricities = []; baseline_retinal_slip_trials = []; ctrst_list=[1 1 1 2 2 2 3 3 3]; ecc_list=[2 4 6 2 4 6 2 4 6];
nms_baseline = []; immobilization_durations = []; nblinks_baseline = []; ms_exist = []; immob_ms_amplitudes = []; trial_counter = [];

% whole trial length: between t = -1000 to 20000

ntps_before = 5000; ntps_after = 1000; % around button press
%ntps_before = 1000; ntps_after = 20000; % from stim onset?
ntimepoints = ntps_before + ntps_after + 1;
reference_point = 'button_press'; % button_press or full_trial (around stim onset) or onset_to_bp (around stim onset but stop at bp)
% We only analyzed button_press, but other options included if of interest.

% baselines: avg for each block / subject
baseline_ms_rate = nan(nsubjects, nblocks);
baseline_Nb = nan(nsubjects, nblocks);
baseline_amplitudes = nan(nsubjects, nblocks);
ntrials_used = zeros(ntimepoints, nsubjects, nblocks);

ms_rate_all = nan(ntimepoints, nsubjects, nblocks);
rate_timepoints = (1:ntimepoints) - ntps_before;
mean_Nb_all = nan((ntimepoints-1)/50, nsubjects, nblocks); % retinal slip box counts
Nb_timepoints = 50/2 + (0:50:(ntimepoints-1-50)) - ntps_before; % stepsize of 50 ms
mean_amplitudes_all = nan(ntimepoints, nsubjects, nblocks);
mean_rts_all = nan(nsubjects, nblocks);

for block = 1:nblocks
    for subj = 1:nsubjects
        % get indices
        nblocktrials = size(ms_present_all{subj}{block}, 2);
        block_ms_present = nan(ntimepoints, nblocktrials);
        fulltrial_block_ms_present = nan(21001, nblocktrials);
        block_ms_amplitudes = nan(ntimepoints, nblocktrials);
        block_ms_rate = nan(numel(rate_timepoints), nblocktrials);
        block_Nb = nan(numel(Nb_timepoints), nblocktrials);
        if nblocktrials > 0 % if not missing whole block
            for trial = 1:nblocktrials
                rt = diff(trial_indices_all{subj}{block}(trial, 2:3));
                if rt > 0%&& rt < 20000 % only analyze trials with certain RT?
                    stim_on = trial_indices_all{subj}{block}(trial, 2);
                    stim_off = trial_indices_all{subj}{block}(trial, 4);
                    button_press = trial_indices_all{subj}{block}(trial, 3);
                    % set full trial timing
                    nan_pad_start_fulltrial = 1000 - (stim_on-1); % if stim onset is a few ms less than ntps_before
                    nan_pad_end_fulltrial = 20000 - (stim_off - stim_on); % if trial didn't last full 20 seconds (i.e. always when button is pressed)
                    if nan_pad_start_fulltrial <= 0
                        first_idx_fulltrial = stim_on - 1000;
                        nan_pad_start_fulltrial = 0;
                    elseif nan_pad_start_fulltrial > 0
                        first_idx_fulltrial = 1;
                    end
                    if nan_pad_end_fulltrial < 0 % if it goes over 20s for some reason
                        last_idx_fulltrial = stim_on + 20000;
                    else
                        last_idx_fulltrial = stim_off;
                    end
                    switch reference_point % epoch around button press or around stim onset?
                        case 'button_press'
                            nan_pad_start = ntps_before - rt; % if RT is faster than ntps_before
                            nan_pad_end = ntps_after-(stim_off - button_press); % if stim offset is a few ms before ntps_after
                            if nan_pad_start >= 0
                                first_idx = stim_on;
                            else
                                first_idx = button_press - ntps_before;
                                nan_pad_start = 0;
                            end
                            if nan_pad_end < 0
                                last_idx = button_press + ntps_after;
                                nan_pad_end = 0;
                            else
                                last_idx = stim_off;
                            end
                        case 'full_trial'
                            nan_pad_start = ntps_before - (stim_on-1); % if stim onset is a few ms less than ntps_before
                            nan_pad_end = ntps_after - (stim_off - stim_on); % if trial didn't last full 20 seconds (i.e. always when button is pressed)
                            if nan_pad_start <= 0
                                first_idx = stim_on - ntps_before + 1;
                                nan_pad_start = 0;
                            elseif nan_pad_start > 0
                                first_idx = 1;
                            end
                            if nan_pad_end < 0 % if it goes over 20s for some reason
                                last_idx = stim_on + ntps_after;
                            else
                                last_idx = stim_off;
                            end
                        case 'onset_to_bp'
                            nan_pad_start = ntps_before - (stim_on-1); % if stim onset is a few ms less than ntps_before
                            nan_pad_end = ntps_after - ((button_press-500) - stim_on); % end after button press
                            if nan_pad_start <= 0
                                first_idx = stim_on - ntps_before + 1;
                                nan_pad_start = 0;
                            elseif nan_pad_start > 0
                                first_idx = 1;
                            end
                            if nan_pad_end < 0 % if it goes over 20s for some reason
                                last_idx = stim_on + ntps_after;
                            else
                                last_idx = (button_press-500);
                            end
                    end

                    rts = [rts; rt]; subjs = [subjs; subj]; contrasts = [contrasts; ctrst_list(block)]; eccentricities = [eccentricities; ecc_list(block)];
                    trial_counter = [trial_counter; trial_counter_all{subj}{block}(trial)];

                    ms_exist = [ms_exist; logical(nansum(ms_present_all{subj}{block}{trial}((stim_on+0):(button_press-300))))];
                    curr_ms_present = ms_present_all{subj}{block}{trial}; % get current trial's timecourse for manipulation
                    curr_posdata = posdata_all{subj}{block}{trial};
                    % remove stim onset transient in ms rate. ONLY for visualization; not removed in the data forwarded on to linear modeling.
                    if strcmp(reference_point, 'button_press')
                        curr_ms_present(stim_on:(stim_on+700)) = nan;
                        curr_posdata(:, stim_on:(stim_on+700)) = nan;
                    end
                    block_ms_present(:, trial) = [nan(1, nan_pad_start), curr_ms_present(first_idx:last_idx), nan(1, nan_pad_end)];
                    curr_blink_present = blink_present_all{subj}{block}{trial};
                    block_blink_present(:, trial) = [nan(1, nan_pad_start), curr_blink_present(first_idx:last_idx), nan(1, nan_pad_end)];
                    fulltrial_block_ms_present(:, trial) = [nan(1, nan_pad_start_fulltrial), ms_present_all{subj}{block}{trial}(first_idx_fulltrial:last_idx_fulltrial), nan(1, nan_pad_end_fulltrial)];
                    % set baseline nms to number of microsaccades from stim onset until button press-300
                    unfilled_ms_onsets = [0, diff(ms_present_all{subj}{block}{trial}(stim_on:button_press-300))];
                    unfilled_ms_onsets(unfilled_ms_onsets < 0) = 0;
                    nms_baseline = [nms_baseline; sum(unfilled_ms_onsets)];
                    % same for blinks
                    unfilled_blink_onsets = [0, diff(blink_present_all{subj}{block}{trial}(stim_on:button_press-300))];
                    unfilled_blink_onsets(unfilled_blink_onsets < 0) = 0;
                    nblinks_baseline = [nblinks_baseline; sum(unfilled_blink_onsets)];
                    % store microsaccade amplitudes
                    block_ms_amplitudes(:, trial) = [nan(1, nan_pad_start), ms_amplitudes_all{subj}{block}{trial}(first_idx:last_idx), nan(1, nan_pad_end)];
                    % get immobilization duration.
                    if strcmp(reference_point, 'button_press')
                        unfilled_ms_offsets = [0, diff(ms_present_all{subj}{block}{trial}(stim_on:button_press))];
                        unfilled_ms_offsets(unfilled_ms_offsets > 0) = 0;
                        ms_offset_idx = find(unfilled_ms_offsets) - numel(unfilled_ms_offsets);
                        last_ms_time = max(ms_offset_idx); % actual last tp of ongoing microsaccade before button press.
                        curr_immobilization_duration = max(ms_offset_idx(ms_offset_idx < -300)); % immobilization time ignores 300ms buffer
                        % and we want the last microsaccade to actually
                        % finish before -300 ms, no microsaccades allowed between -300 and 0.
                        % and we want trials with no blinks (nans) in between the
                        % last microsaccade and button press:
                        if isempty(curr_immobilization_duration) || sum(ms_offset_idx >= -300 & ms_offset_idx <= 0) || sum(isnan(ms_present_all{subj}{block}{trial}((curr_immobilization_duration:0)+button_press)))
                            % no microsaccade OFFSETS between -300 and 0. It's okay if a microsaccade starts at like
                            % t=-10 ms, that's not going to affect anything.
                            curr_immobilization_duration = nan;
                            immob_ms_amplitude = nan;
                        else
                            immob_ms_idx = curr_immobilization_duration + button_press; % restore to timepoint idx (aka in milliseconds)
                            immob_ms_amplitude = ms_amplitudes_all{subj}{block}{trial}(immob_ms_idx-1);
                        end
                        immobilization_durations = [immobilization_durations; curr_immobilization_duration];
                        immob_ms_amplitudes = [immob_ms_amplitudes; immob_ms_amplitude];
                    end
                    trial_xPos = [nan(1, nan_pad_start), nanmean(curr_posdata(2:3, first_idx:last_idx), 1), nan(1, nan_pad_end)];
                    trial_xPos = (trial_xPos - nanmean(trial_xPos)) * degperpix;
                    trial_yPos = [nan(1, nan_pad_start), nanmean(curr_posdata(4:5, first_idx:last_idx), 1), nan(1, nan_pad_end)];
                    trial_yPos = (trial_yPos - nanmean(trial_yPos)) * degperpix;
                    block_Nb(:, trial) = retinal_slip(trial_xPos, trial_yPos, block_ms_present(:, trial), 0); % last arg 1: include microsaccades
                    fulltrial_xPos = nanmean(curr_posdata(2:3, :), 1); fulltrial_yPos = nanmean(curr_posdata(4:5, :), 1);
                    fulltrial_xPos = (fulltrial_xPos - nanmean(fulltrial_xPos)) * degperpix; fulltrial_yPos = (fulltrial_yPos - nanmean(fulltrial_yPos)) * degperpix;
                    fulltrial_block_Nb = retinal_slip(fulltrial_xPos, fulltrial_yPos, ms_present_all{subj}{block}{trial}', 0);
                    baseline_retinal_slip_trials = [baseline_retinal_slip_trials; nanmean(fulltrial_block_Nb(floor(stim_on/50):floor((button_press-300)/50)))];
                    if isnan(baseline_retinal_slip_trials(end))
                        error('tt')
                    end
                end
                ms_rate_all(:, subj, block) = cross_trial_ms_rate(block_ms_present, 'alpha');
                %[~, ms_rate_all(:, subj, block)] = cross_trial_ms_rate(block_ms_present, 'alpha'); % investigate shuffled ms instead.
                mean_amplitudes_all(:, subj, block) = nanmean(block_ms_amplitudes, 2);
                baseline_ms_rate(subj, block) = nanmean(ms_rate_all(1:2000, subj, block)); % first 2000 of epoch = baseline.
                mean_Nb_all(:, subj, block) = nanmean(block_Nb, 2);
                baseline_Nb(subj, block) = nanmedian(mean_Nb_all(1:2000/50, subj, block)); % first 2000ms of epoch = baseline
                baseline_amplitudes(subj, block) = nanmean(mean_amplitudes_all(1:2000, subj, block)); % first 2000 of epoch = baseline.
            end
        end
    end
end

% normalize everything by subject baseline means
norm_ms_rate_all = nan(size(ms_rate_all));
norm_Nb_all = nan(size(mean_Nb_all));
norm_amplitudes_all = nan(size(mean_amplitudes_all));
for subj = 1:nsubjects
    norm_ms_rate_all(:, subj, :) = (ms_rate_all(:, subj, :) ./ nanmean(baseline_ms_rate(subj, :))) .* nanmean(nanmean(baseline_ms_rate));
    norm_Nb_all(:, subj, :) = mean_Nb_all(:, subj, :) ./ nanmean(baseline_Nb(subj, :));
    norm_amplitudes_all(:, subj, :) = mean_amplitudes_all(:, subj, :) ./ nanmean(baseline_amplitudes(subj, :));
end

% save trial data in .csv for linear modeling/plotting in R
fts = rts ./ 1000;
immobs_seconds = immobilization_durations ./ 1000 .* -1; % make positive and in seconds

tabledata = [subjs, fts, nms_baseline, baseline_retinal_slip_trials, immobs_seconds, immob_ms_amplitudes, contrasts, eccentricities, ms_exist, nblinks_baseline, trial_counter];
T = array2table(tabledata, 'VariableNames', {'subject', 'FT', 'baselinenms', 'baselineretinalslip', 'immobduration', 'immobmsamplitude', 'contrast', 'eccentricity', 'mspresent', 'nblinks', 'trial'});
filename = 'trial_data_withblinks.csv'; % main illusion trials
%filename = 'trial_data_replay_withblinks.csv'; % replay catch trials
%filename = 'trial_data_sharp_withblinks.csv'; % sharp catch trials
writetable(T, filename);

%% plot averaged by size or contrast

window_length = 200;
% average contrasts together or averaged sizes together
size_indices = [3 2 1 3 2 1 3 2 1]; % 1=size 6, 2=size 4, 3=size 2
color_indices = [1 1 1 2 2 2 3 3 3]; % 1=low contrast, 2=medium, 3=high contrast
all_indices = [1 1 1 1 1 1 1 1 1];

which2use = 'all'; % size or color or all (average across all blocks) or separate (plot each condition)
which2plot = 'rate'; % rate, amplitude, retinal slip. Only rate is analyzed but other options included if of interest.
subjects2use = 1:nsubjects; % default: all
overlay = 1; % overlay the plots on one axis or use subplots
plotcolors = {[1 0 0], [0 1 0], [0 0 1]};
nobounds = 0; % just plot line, no stderr bounds?

f = figure; hold on;
%for subj = 1:nsubjects % plot each subject separately?
%    subjects2use = subj;
%figure; hold on;
for i = 1:9
    switch which2use
        case 'size'
            if overlay; hold on; else; subplot(1, 3, i); end
            curr_block_idx = find(size_indices == i); % if plotting by size
        case 'color'
            if overlay; hold on; else; subplot(1, 3, i); end
            curr_block_idx = find(color_indices == i); % if plotting by color
        case 'all'
            curr_block_idx = find(all_indices);
            plotcolors = {[0 0 0]};
        case 'separate'
            if overlay; hold on; else; subplot(3, 3, i); end
            curr_block_idx = i;
            plotcolors = {[0 0 .33], [0 .33 0], [.33 0 0], [0 0 .66], [0 .66 0], [.66 0 0], [0 0 1], [0 1 0], [1 0 0]}; % colors = sizes
            %plotcolors = {[.33 0 0], [.66 0 0], [1 0 0], [0 .33 0], [0 .66 0], [0 1 0], [0 0 .33], [0 0 .66], [0 0 1]}; % colors = contrasts
    end

    switch which2plot
        case 'rate'
            norm_ms_rate_all(isinf(norm_ms_rate_all)) = nan;
            ms_rate_corr = norm_ms_rate_all(:, subjects2use, curr_block_idx);
            %ms_rate_corr = ms_rate_all(:, subjects2use, curr_block_idx); % not subject normalized
            ms_rate_curr = nanmean(ms_rate_corr, 3);
            % plot ms rate
            mean_ms_rate = nanmean(ms_rate_curr, 2);
            stderr_ms_rate = nanstd(ms_rate_curr, 0, 2) ./ sqrt(numel(subjects2use));
            if nobounds
                bp = plot(rate_timepoints(250:end-250), mean_ms_rate(250:end-250), 'Color', plotcolors{i});
            else
                bp = boundedline(rate_timepoints(250:end-250), mean_ms_rate(250:end-250), stderr_ms_rate(250:end-250), 'cmap', plotcolors{i}, 'alpha', 'transparency', 0.05);
            end
            ylbl = 'microsaccade rate (Hz)';

        case 'amplitude'
            % plot average ms amplitudes
            %ms_amplitude_corr = mean_amplitudes_all(:, subjects2use, curr_block_idx); % not subject normalized
            ms_amplitude_corr = norm_amplitudes_all(:, subjects2use, curr_block_idx) .* nanmean(nanmean(baseline_amplitudes)); % normalize, then multiply by total mean to interpret it
            ms_amplitude_curr = nanmean(ms_amplitude_corr, 3);
            mean_ms_amplitude = nanmean(ms_amplitude_curr, 2);
            stderr_ms_amplitude = nanstd(ms_amplitude_curr, 0, 2) ./ sqrt(nsubjects);
            smoothed_ampl = smoothdata(mean_ms_amplitude, 'movmean', window_length);
            smoothed_stderr_ampl = smoothdata(stderr_ms_amplitude, 'movmean', window_length);
            bp = boundedline((-ntps_before:ntps_after)' * 1000/fps, smoothed_ampl, smoothed_stderr_ampl, 'cmap', plotcolors{i}, 'alpha', 'transparency', 0.05);
            ylbl = 'm.s. amplitude (deg)';

        case 'retinal slip'
            % plot retinal slip (nboxes of 0.01 length traversed in 50 ms steps)
            mean_Nb_curr = nanmean(norm_Nb_all(:, subjects2use, curr_block_idx), 3) .* nanmean(nanmean(baseline_Nb)); % normalize, then multiply by total mean to interpret it
            %mean_Nb_curr = nanmean(mean_Nb_all(:, subjects2use, curr_block_idx), 3); % raw mean
            mean_Nb = nanmean(mean_Nb_curr, 2);
            stderr_Nb = nanstd(mean_Nb_curr, 0, 2) ./ sqrt(nsubjects);
            Nb_window_length = window_length / 50;
            smoothed_mean_Nb = smoothdata(mean_Nb, 'movmean', Nb_window_length/stepsize);
            smoothed_stderr_Nb = smoothdata(stderr_Nb, 'movmean', Nb_window_length/stepsize);
            % convert Nb to dva: multiply by 0.01 dva
            smoothed_mean_Nb = smoothed_mean_Nb * 0.01; smoothed_stderr_Nb = smoothed_stderr_Nb * 0.01;
            bp = boundedline(Nb_timepoints, smoothed_mean_Nb, smoothed_stderr_Nb, 'cmap', plotcolors{i}, 'alpha', 'transparency', 0.05);
            ylbl = 'retinal slip (deg / 50ms)';

    end
    set(bp, 'LineWidth', 1);

    H = gca;
    H.LineWidth = 2;

    % plot vertical line at t=0 (button press)
    curr_ylim = get(gca, 'ylim');
    hold on; zeroline = plot([0 0], curr_ylim, '--k');
    set(get(get(zeroline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    % set axis limits
    curr_xlim = get(gca, 'xlim');
    if numel(subjects2use) > 1
        switch which2plot
            case 'rate'
                if strcmp(which2use, 'separate')
                    set(gca, 'ylim', [0 1.8]);
                else
                    set(gca, 'ylim', [0.3 1.3]);
                end
            case 'amplitude'
                set(gca, 'ylim', [0.2 0.5]);
            case 'retinal slip'
                set(gca, 'ylim', [0.15 0.2]);
        end
    end


    % add y/x axis labels
    if ~overlay && ~strcmp(which2use, 'all')
        if i == 1
            ylabel('small');
        elseif i == 2
            ylabel('medium');
        elseif i == 3
            ylabel('large');
            xlabel('time (s)');
        end
    else
        ylabel(which2plot);
        xlabel('time (s)');
        if strcmp(which2use, 'all')
            set(bp, 'Color', [0 0 0]);
        end
    end
    if strcmp(which2use, 'all')
        break
    elseif ~strcmp(which2use, 'separate') && i == 3
        break
    end
end
switch which2use
    case 'size'
        if overlay
            legend({'large', 'medium', 'small'});
            sl = title([which2plot ' by boundary eccentricity']);
        else
            sl = suplabel('center radius', 'y');
        end
        set(sl, 'FontSize', 12);
    case 'color'
        if overlay
            legend({'low', 'medium', 'high'});
            sl = title([which2plot ' by color contrast']);
        else
            sl = suplabel('color contrast', 'y');
        end
        set(sl, 'FontSize', 12);
    case 'separate'
        if overlay
            lgnd = legend({'small low', 'medium low', 'large low', 'small medium', 'medium medium', 'large medium', 'small high', 'medium high', 'large high'});
            sl = title([which2plot ' by condition']);
        else
            sly = suplabel('color contrast', 'y');
            slx = suplabel('center size', 'x');
        end
end
set(gca, 'FontSize', 9);

% set subject color if plotting just 1 subject at a time
if numel(subjects2use) == 1
    clrs = colormap;
    set(bp, 'Color', clrs(subjects2use * 3, :));
    set(bp, 'LineWidth', 1); % or make line smaller
end

hold on;
xxlim = get(gca, 'xlim');

%end % uncomment if plotting each subject individually

ylabel(ylbl);
legend boxoff

% Override some plotting parameters for manuscript figures
ylbl = ylabel('Normalized Microsaccade Rate (Hz)');
xlbl = xlabel('Time (ms)');
ttl = title('Microsaccades Time-Locked to Fading Reports', 'FontWeight', 'Normal');
set(gca, 'FontSize', 9, 'LineWidth', 1);
set(bp, 'LineWidth', 1);
set(ylbl, 'FontSize', 10); set(xlbl, 'FontSize', 10); set(ttl, 'FontSize', 12);

switch which2use
    case 'all' % Figure 2b
        legend off
        set(gcf, 'Units','Inches', 'Position', [0, 0, 4, 3]);
        savefig(f, 'fig2b.fig');
        exportgraphics(f, 'fig2b.pdf', 'ContentType', 'vector', 'Resolution', 600);
    case 'separate' % Supplementary Figure 2
        set(gcf, 'Units', 'Inches', 'Position', [0, 0, 6, 3]);
        legend('Location', 'Southwest', 'NumColumns', 3, 'FontSize', 9);
        default_legendlinelength = lgnd.ItemTokenSize;
        lgnd.ItemTokenSize = default_legendlinelength / 2;
        savefig(f, 'figs2.fig');
        exportgraphics(f, 'figs2.pdf', 'ContentType', 'vector', 'Resolution', 600);
end

%% microsaccade characteristics
% plot xy amplitude vs peak xy velocity
all_xy_amplitudes = [];
all_xy_peakvelocities = [];
all_durations = [];
all_subjnums = []; % label each ms with its subject
for subj = 1:nsubjects
    for block = 1:nblocks
        for trial = 1:numel(BinocularSaccades_all{subj}{block})
            if ~isempty(BinocularSaccades_all{subj}{block}{trial})
                all_xy_amplitudes = [all_xy_amplitudes; BinocularSaccades_all{subj}{block}{trial}.Amplitude];
                all_xy_peakvelocities = [all_xy_peakvelocities; BinocularSaccades_all{subj}{block}{trial}.vPeak];
                all_durations = [all_durations; BinocularSaccades_all{subj}{block}{trial}.Duration];
                all_subjnums = [all_subjnums; repmat(subj, numel(BinocularSaccades_all{subj}{block}{trial}.Duration'), 1)];
            end
        end
    end
end

% remove too big or too small
all_xy_amplitudes_sort = all_xy_amplitudes(all_xy_amplitudes < 2 & all_xy_peakvelocities < 300);
all_xy_peakvelocities_sort = all_xy_peakvelocities(all_xy_amplitudes < 2 & all_xy_peakvelocities < 300);
all_subjnums_sort = all_subjnums(all_xy_amplitudes < 2 & all_xy_peakvelocities < 300);
all_durations_sort = all_durations(all_xy_amplitudes < 2 & all_xy_peakvelocities < 300);
f = figure; hold on;
scatter(all_xy_amplitudes_sort(all_durations_sort > 0), all_xy_peakvelocities_sort(all_durations_sort > 0), 2, 'filled', 'k', 'MarkerFaceAlpha', 0.05);

% plotting parameters for Figure 2a
set(gca, 'xscale', 'log'); set(gca, 'yscale', 'log');
set(gca, 'FontSize', 9, 'LineWidth', 1);
set(gca, 'FontName', 'Helvetica')
xlbl = xlabel('Amplitude (dva)'); ylbl = ylabel('Velocity (dva/s)');
ttl = title('Microsaccade Main Sequence', 'FontWeight', 'Normal');
set(gcf, 'Units','Inches', 'Position', [0, 0, 4, 3]);
set(ylbl, 'FontSize', 10); set(xlbl, 'FontSize', 10); set(ttl, 'FontSize', 12);

savefig(f, 'fig2a.fig');
exportgraphics(f, 'fig2a.pdf', 'ContentType', 'vector', 'Resolution', 600);
%% plot each trial gaze timecourses for quality control
% for one subject at a time, plot the gaze paths for each trial, with
% detected microsaccades overlaid. Assess if we should raise or lower
% the noise threshold for microsaccade detection.
subj = 21;
startidx = 2; % which of trial_indices should we use
stopidx = 3;
% 1: fixation onset, 2: stim onset, 3: button press, 4: stim offset
eyecolors = {'r', 'b'};
dircolors = {'r', 'b'};
window_length_tps = 31; % for smoothing

for iBlock = 1:nblocks
    ntrials = ntrials_all(subj, iBlock);
    plotdim1 = floor(sqrt(ntrials));
    plotdim2 = ceil(ntrials/plotdim1);
    fh = figure; suptitle(['s', num2str(subjects(subj)), ' block ', num2str(iBlock)]);
    fh.WindowState = 'maximized';
    for iTrial = 1:ntrials
        starttp = trial_indices_all{subj}{iBlock}(iTrial, startidx);
        stoptp = trial_indices_all{subj}{iBlock}(iTrial, stopidx);
        gazeX = posdata_all{subj}{iBlock}{iTrial}(2:3, starttp:stoptp);
        gazeY = posdata_all{subj}{iBlock}{iTrial}(4:5, starttp:stoptp);
        % smooth data
        gazeX = smoothdata(gazeX, 2, 'movmean', window_length_tps);
        gazeY = smoothdata(gazeY, 2, 'movmean', window_length_tps);
        gazeX = (gazeX - nanmean(gazeX, 2)) * degperpix; % coordinates from center (0) (pixels)
        gazeY = (gazeY - nanmean(gazeY, 2)) * degperpix;
        gazeXms = gazeX; gazeXms(:, ms_present_all{subj}{iBlock}{iTrial}(starttp:stoptp) ~= 1) = NaN;
        gazeYms = gazeY; gazeYms(:, ms_present_all{subj}{iBlock}{iTrial}(starttp:stoptp) ~= 1) = NaN;
        subplot(plotdim1, plotdim2, iTrial); hold on;
        % plot each eye?
        %         for iEye = 1:2
        %             % plot x vs y
        %             %plot(gazeX(iEye, :), gazeY(iEye, :), eyecolors{iEye});
        %             %plot(gazeXms(iEye, :), gazeYms(iEye, :), 'k', 'LineWidth', 5);
        %             % plot x vs time, y vs time
        %             time = 1:size(gazeX(1, :), 2);
        %             plot(time, gazeX(iEye, :), dircolors{1}, time, gazeY(iEye, :), dircolors{2}, 'LineWidth', 1);
        %             plot(time, gazeXms(eye, :), 'k', time, gazeYms(iEye, :), 'k', 'LineWidth', 3);
        %         end
        % average eyes
        gazeYavg = nanmean(gazeY, 1); gazeYmsavg = nanmean(gazeYms, 1);
        gazeXavg = nanmean(gazeX, 1); gazeXmsavg = nanmean(gazeXms, 1);
        % plot x vs y
        %plot(gazeXavg, gazeYavg, eyecolors{1});
        %plot(gazeXmsavg, gazeYmsavg, 'k');
        % plot x vs time, y vs time
        time = 1:size(gazeXavg, 2);
        plot(time, gazeXavg, dircolors{1}, time, gazeYavg, dircolors{2}, 'LineWidth', 1);
        plot(time, gazeXmsavg, 'k', time, gazeYmsavg, 'k', 'LineWidth', 3);
    end
end
