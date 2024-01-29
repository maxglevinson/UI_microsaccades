%% group analysis for eyelink data
addpath(genpath(pwd));
subjects = [201, 202, 203, 204, 205, 206, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 222, 223, 224];

nsubjects = numel(subjects);
nblocks = 9;
ms_present_all_unsorted = cell(1, nsubjects);
pupil_trace_all_unsorted = cell(1, nsubjects);
BinocularSaccades_all_unsorted = cell(1, nsubjects);
trial_indices_all_unsorted = cell(1, nsubjects);
xyVelocity_all_unsorted = cell(1, nsubjects);
posdata_all_unsorted = cell(1, nsubjects);
trial_counter_all_unsorted = cell(1, nsubjects);
extra = cell(nsubjects, 3); % anything extra I want to output and store to check
for subj = 1:nsubjects
    subj_id = subjects(subj);
    for session_id = 1:3
        [curr_ms_present, fps, curr_pupil_trace, curr_BinocularSaccades, curr_trial_indices, curr_xyVelocity, curr_posdata, curr_trial_counter, extra{subj, session_id}] = analyze_eyelink_fulltrial(subj_id, session_id);
        ms_present_all_unsorted{subj} = [ms_present_all_unsorted{subj}, curr_ms_present];
        pupil_trace_all_unsorted{subj} = [pupil_trace_all_unsorted{subj}, curr_pupil_trace];
        BinocularSaccades_all_unsorted{subj} = [BinocularSaccades_all_unsorted{subj}, curr_BinocularSaccades];
        trial_indices_all_unsorted{subj} = [trial_indices_all_unsorted{subj}, curr_trial_indices];
        xyVelocity_all_unsorted{subj} = [xyVelocity_all_unsorted{subj}, curr_xyVelocity];
        posdata_all_unsorted{subj} = [posdata_all_unsorted{subj}, curr_posdata]; % t, x pos, y pos
        trial_counter_all_unsorted{subj} = [trial_counter_all_unsorted{subj}, curr_trial_counter];
    end
end


% set order to: sizes  2 4 6 2 4 6 2 4 6
%               colors 1 1 1 2 2 2 3 3 3
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
% for every millisecond, get the amplitude of any ongoing microsaccade
% (otherwise NaN)
ms_amplitudes_all_unsorted = repmat({cell(1, nblocks)}, 1, nsubjects);
% directions: 2 = left inwards, 1 = left outwards, -1 = right outwards, -2 = right inwards
ms_directions_all_unsorted = repmat({cell(1, nblocks)}, 1, nsubjects);
ms_durations_all_unsorted = repmat({cell(1, nblocks)}, 1, nsubjects);
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
                        ms_durations_all_unsorted{subj}{block}{trial}(ms) = BinocularSaccades_all_unsorted{subj}{block}{trial}.Duration(nms); % assign tp with microsaccade's duration
                        ms_directions_all_unsorted{subj}{block}{trial}(ms) = BinocularSaccades_all_unsorted{subj}{block}{trial}.Direction(nms); % assign tp with microsaccade's direction
                        end
                    end
                end
        end
    end
    
    % reorder blocks
    ms_present_all{subj} = ms_present_all_unsorted{subj}(blockorder(subj, :));
    pupil_trace_all{subj} = pupil_trace_all_unsorted{subj}(blockorder(subj, :));
    BinocularSaccades_all{subj} = BinocularSaccades_all_unsorted{subj}(blockorder(subj, :));
    trial_indices_all{subj} = trial_indices_all_unsorted{subj}(blockorder(subj, :));
    ms_amplitudes_all{subj} = ms_amplitudes_all_unsorted{subj}(blockorder(subj, :));
    ms_durations_all{subj} = ms_durations_all_unsorted{subj}(blockorder(subj, :));
    ms_directions_all{subj} = ms_directions_all_unsorted{subj}(blockorder(subj, :));
    xyVelocity_all{subj} = xyVelocity_all_unsorted{subj}(blockorder(subj, :));
    posdata_all{subj} = posdata_all_unsorted{subj}(blockorder(subj, :));
    trial_counter_all{subj} = trial_counter_all_unsorted{subj}(blockorder(subj, :));
    for block = 1:nblocks
        ntrials_all(subj, block) = size(ms_present_all{subj}{block}, 2);
    end
end
ntrialstotal_byblock = sum(ntrials_all);

% convert position data to degrees of visual angle? actually will need
% to zero mean first
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


% save('current_data.mat', 'subjects', 'nsubjects', 'ntrials_all', ...
%     'ntrialstotal_byblock', 'ms_present_all', 'blink_present_all', 'pupil_trace_all', ...
%     'BinocularSaccades_all', 'posdata_all', 'trial_indices_all', 'ms_amplitudes_all', ...
%     'ms_directions_all', 'ms_durations_all', 'xyVelocity_all', 'trial_counter_all', ...
%     'degperpix', 'fps', 'extra', 'nblocks', 'blockorder', '-v7.3');

%% just load the saved data
datadir = '/export04/data/mlevin/UI_eyelink';
addpath(genpath(datadir));
load('current_data_withblinks.mat');

%% epoch trials into block averages

rts = []; subjs = []; contrasts = []; eccentricities = []; baseline_retinal_slip_trials = []; ctrst_list=[1 1 1 2 2 2 3 3 3]; ecc_list=[2 4 6 2 4 6 2 4 6];
nms_window = []; min_slip_window = []; min_window = -1000:-300; Nb_min_window = -1000/50:-300/50; nms_baseline = []; immobilization_durations = [];
trial_blinkrate_baseline = []; nblinks_baseline = []; ms_exist = []; blink_exist = []; session_nums = []; last_ms_times = []; immob_ms_amplitudes = []; last_ms_amplitudes = [];
early_ms_exist = []; late_ms_exist = []; trial_counter = []; meanmsampl_baseline = [];

subj_shifts = [ -125  -251   -26  -251   -65  -251   249  -251    50  -251    78   249   242  -251   249   -53   137   249  -170     4  -251]; % after 500ms smoothing, bw 300-800ms
subj_rt_window = [-800:-300]; % indiv subj RT estimation can fall within this window.
subj_min_idx = subj_rt_window(250 - subj_shifts);
% whole trial length: between t = -1000 to 20000
baseline_correct = 1; % baseline correct? yes or no
weighsubjects = 0; % weigh subjects by usable trial counts? yes or no
subset_sample = 0; % only keep subset of trials? (keep every XX trial, etc)

sample_counter = 1; % start at 1.
sample_incr = 1; %(0.9 to ignore t=0-750ms)
sample_num = 6; % include trial when we reach this counter number (4 to ignore t=0-750ms)

ntps_before = 5000; ntps_after = 1000; % around button press
%ntps_before = 1000; ntps_after = 20000; % from stim onset
ntimepoints = ntps_before + ntps_after + 1;
reference_point = 'button_press'; % button_press or full_trial (around stim onset) or onset_to_bp (around stim onset but stop at bp)
buffer_baseline_idx = 1400:10000;%200:1400; % this is longer than below to account for transients in ms rate computation.
baseline_sub_idx = 401:8201; % within whole baseline window ^, which do we actually want to use?
ntimepointsbufferbaseline = numel(buffer_baseline_idx);

ntimepointsbaseline = numel(baseline_sub_idx);
baseline_probability = nan(nsubjects, nblocks); % avg for each block / subject
baseline_pupil = nan(nsubjects, nblocks);
%baseline_ms_counts = nan(nsubjects, nblocks);
baseline_Nb = nan(nsubjects, nblocks);
baseline_amplitudes = nan(nsubjects, nblocks);
baseline_xyVelocity = nan(nsubjects, nblocks);
before_press_probability = nan(nsubjects, nblocks);
ntrials_used = zeros(ntimepoints, nsubjects, nblocks);

window_length = 200;
rate_window_length = 1; % for sliding window ms/sec calculation

ms_probability_all = nan(ntimepoints, nsubjects, nblocks);
stepsize = 1;
ms_counts_all = nan(ntimepoints, nsubjects, nblocks);%ms_counts_all = nan((ntimepoints-rate_window_length-1)/stepsize + 1, nsubjects, 9);
count_timepoints = (1:ntimepoints) - ntps_before;%count_timepoints = rate_window_length/2 + (0:stepsize:(ntimepoints-1-rate_window_length)) - ntps_before; % time stamps for each time window in ms counts / second
count_timepoints_baseline = 1:ntimepointsbaseline;%count_timepoints_baseline = rate_window_length/2 + (0:stepsize:(ntimepointsbaseline-1-rate_window_length));
blink_rate_all = nan(ntimepoints, nsubjects, nblocks);
mean_Nb_all = nan((ntimepoints-1)/50, nsubjects, nblocks); % retinal slip box counts
Nb_timepoints = 50/2 + (0:50:(ntimepoints-1-50)) - ntps_before; % stepsize of 50 ms
Nb_timepoints_baseline = 50/2 + (0:50:(ntimepointsbaseline-1-50));
mean_pupil_all = nan(ntimepoints, nsubjects, nblocks);
mean_amplitudes_all = nan(ntimepoints, nsubjects, nblocks);
mean_xyVelocity_all = nan(ntimepoints, nsubjects, nblocks);
% mean_ten50power_all = nan(nsubjects, 9);
% mean_eighty200power_all = nan(nsubjects, 9);
mean_rts_all = nan(nsubjects, nblocks);

trans_ntps_before = 900;
trans_ntps_after = 4000;
transition_probabilities_all = nan(trans_ntps_before+1+trans_ntps_after, nsubjects, nblocks);
baseline_transition_probabilities = transition_probabilities_all; % shuffled model
norm_hazard = 0; % normalize transition probabilities by hazard rate per block.

weighting = nan(nsubjects, nblocks);

ms_present_concat = [];%cell(1, nsubjects); % literally all trials concatenated.

for block = 1:nblocks
    %figure;
    for subj = 1:nsubjects
        subj_blockorder = reshape(blockorder(subj, :), 3, 3); % columns: days
        % get indices
        nblocktrials = size(ms_present_all{subj}{block}, 2);
        baseline_block_ms_present = nan(ntimepointsbufferbaseline, nblocktrials);
        block_ms_present = nan(ntimepoints, nblocktrials);
        fulltrial_block_ms_present = nan(21001, nblocktrials);
        baseline_block_pupil_trace = nan(ntimepointsbaseline, nblocktrials);
        block_pupil_trace = nan(ntimepoints, nblocktrials);
        block_ms_amplitudes = nan(ntimepoints, nblocktrials);
        baseline_block_ms_amplitudes = nan(ntimepointsbaseline, nblocktrials);
        block_xyVelocity = nan(ntimepoints, nblocktrials);
        baseline_block_xyVelocity = nan(ntimepointsbaseline, nblocktrials);
        block_ms_counts = nan(numel(count_timepoints), nblocktrials);
        baseline_block_ms_counts = nan(numel(count_timepoints_baseline), nblocktrials);
        block_blink_present = nan(ntimepoints, nblocktrials);
        block_blink_rate = nan(numel(count_timepoints), nblocktrials);
        baseline_block_blink_rate = nan(numel(count_timepoints_baseline), nblocktrials);
        block_Nb = nan(numel(Nb_timepoints), nblocktrials);
        baseline_block_Nb = nan(numel(Nb_timepoints_baseline), nblocktrials);
        baseline_block_xPos = nan(ntimepointsbaseline, nblocktrials);
        baseline_block_yPos = nan(ntimepointsbaseline, nblocktrials);
        block_rts = nan(1, nblocktrials);
        block_transition_probability = [];
        baseline_block_transition_probability = [];
        curr_subj_reactiontime = subj_min_idx(subj); % estimated rt from experience -> button press for that particular subject
        if nblocktrials > 0 % if not missing whole block
            for trial = 1:nblocktrials
                rt = diff(trial_indices_all{subj}{block}(trial, 2:3));
                if rt > 0 && rt < 19000%20000 % only use trials with certain RT?
                    stim_on = trial_indices_all{subj}{block}(trial, 2);
                    stim_off = trial_indices_all{subj}{block}(trial, 4);
                    button_press = trial_indices_all{subj}{block}(trial, 3);% -randi(2000, 1); % shift by 1000ms for code testing.?
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
                            nan_pad_start = ntps_before - rt;%(button_press - stim_on); % if RT is faster than ntps_before
                            nan_pad_end = ntps_after-(stim_off - button_press);% if stim offset is a few ms before ntps_after
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
                                bp500_idx = (ntps_before-700):(ntps_before-300); % avg for each block / subject
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
                    % ONLY TAKE TRIALS WITH ZERO MICROSACCADES BETWEEN STIM ON+700 AND BUTTON PRESS?
                    % add some buffer for after button press motor command is executed
%                                         if nansum(ms_present_all{subj}{block}{trial}((stim_on):(button_press-300)))
%                                             continue % skip the trial if it has a microsaccade.
%                                         end
                    
                    %                     ONLY TAKE TRIALS WITH ZERO BLINKS / MISSING DATA BETWEEN STIM ON+750 AND BUTTON PRESS?
%                                         if sum(isnan(ms_present_all{subj}{block}{trial}((stim_on):(button_press-300))))
%                                             continue % skip the trial if it has a blink.
%                                         end
                    
                    %                     % ONLY TAKE TRIALS WITH AT LEAST ONE MICROSACCADE?
%                                         if ~nansum(ms_present_all{subj}{block}{trial}((stim_on):(button_press-300)))
%                                             continue % skip the trial if it doesn't have a microsaccade.
%                                         end
                                        
                                        
                                        if subset_sample
                                            if round(sample_counter) < sample_num % if taking subset, only keep every 6 trials
                                                sample_counter = sample_counter + sample_incr;
                                                continue
                                            elseif round(sample_counter) == sample_num
                                                sample_counter = sample_counter - sample_num + 1;
                                            end
                                        end
                    [~, session_num] = find(block == subj_blockorder); % which day was this block on? day 1, 2, 3
                    rts = [rts; rt]; subjs = [subjs; subj]; contrasts = [contrasts; ctrst_list(block)]; eccentricities = [eccentricities; ecc_list(block)];
                    block_rts(trial) = rt; session_nums = [session_nums; session_num]; trial_counter = [trial_counter; trial_counter_all{subj}{block}(trial)];
                    
                    %buffer_baseline_idx = (button_press - 2000):(button_press - 1400);
                    extrananbuffer = buffer_baseline_idx(end) - numel(ms_present_all{subj}{block}{trial}); curr_ms_present_baseline = [ms_present_all{subj}{block}{trial}, nan(1, extrananbuffer)]; curr_ms_present_baseline((button_press-3000):end) = nan;
                    baseline_block_ms_present(:, trial) = curr_ms_present_baseline(buffer_baseline_idx);
                    curr_blink_present_baseline = [blink_present_all{subj}{block}{trial}, nan(1, extrananbuffer)]; curr_blink_present_baseline((button_press-3000):end) = nan;
                    baseline_block_blink_present(:, trial) = curr_blink_present_baseline(buffer_baseline_idx);
                    curr_nms_baseline = [0; diff(baseline_block_ms_present(:, trial))];
                    %nms_baseline = [nms_baseline; sum(curr_nms_baseline(baseline_sub_idx)==1)] / sum(~isnan(curr_nms_baseline(baseline_sub_idx)))];
                
                    trial_nblinks_baseline = [0; diff(baseline_block_blink_present(:, trial))]; nblink_idx = find(trial_nblinks_baseline); if isempty(nblink_idx), nblink_idx = 1; end
                    trial_blinkrate_baseline = [trial_blinkrate_baseline; sum(trial_nblinks_baseline == 1) / nblink_idx(end)]; % last idx = length (end) of baseline window
                    %nblinks_baseline = [nblinks_baseline; sum(trial_nblinks_baseline == 1)];
                    ms_exist = [ms_exist; logical(nansum(ms_present_all{subj}{block}{trial}((stim_on+0):(button_press-300))))];
                    early_ms_exist = [early_ms_exist; logical(nansum(ms_present_all{subj}{block}{trial}(stim_on:stim_on+500)))];
                    late_ms_exist = [late_ms_exist; logical(nansum(ms_present_all{subj}{block}{trial}((stim_on+501):(button_press-300))))];
                    blink_exist = [blink_exist; logical(nansum(blink_present_all{subj}{block}{trial}((stim_on+0):(button_press-300))))];
                    curr_ms_present = ms_present_all{subj}{block}{trial}; % get current trial's timecourse for manipulation
                    curr_blink_present = blink_present_all{subj}{block}{trial};
                    curr_posdata = posdata_all{subj}{block}{trial};
                    if strcmp(reference_point, 'button_press')
                        curr_ms_present(stim_on:(stim_on+700)) = nan; % remove stim onset transient in ms rate.
                        curr_posdata(:, stim_on:(stim_on+700)) = nan;
                    end
                    block_ms_present(:, trial) = [nan(1, nan_pad_start), curr_ms_present(first_idx:last_idx), nan(1, nan_pad_end)];
                    block_blink_present(:, trial) = [nan(1, nan_pad_start), curr_blink_present(first_idx:last_idx), nan(1, nan_pad_end)];
                    curr_nms_window = [0; diff(block_ms_present(ntps_before+1+min_window, trial))];
                    nms_window = [nms_window; sum(curr_nms_window==1)];
                    fulltrial_block_ms_present(:, trial) = [nan(1, nan_pad_start_fulltrial), ms_present_all{subj}{block}{trial}(first_idx_fulltrial:last_idx_fulltrial), nan(1, nan_pad_end_fulltrial)];
                    %%% OPTION regress out global stim-evoked pupil response (see last section in this script)
                    %                     P = pupil_trace_all{subj}{block}{trial}(31:3000);%end-30); % 30 on either side are nans from pupil preprocessing
                    %                     curr_stimresponse = first2000avgpupil(31:numel(P)+30);
                    %                     P0 = P - mean(P);
                    %                     curr_stimresponse = curr_stimresponse - mean(curr_stimresponse);
                    %                     Pres = P - curr_stimresponse'; % just subtract
                    %                     Pres0 = Pres + mean(P); % add back the orig mean
                    %                     Pfull = [nan(1, 30), Pres0, pupil_trace_all{subj}{block}{trial}(3001:end)];
                    %                     %B = P' \ curr_stimresponse;
                    %                     %Pres = P - curr_stimresponse' * B;
                    %                     block_pupil_trace(:, trial) = [nan(1, nan_pad_start), Pfull(first_idx:last_idx), nan(1, nan_pad_end)];
                    %%%
                    % set baseline nms to number of microsaccades from stim onset until button press
                    unfilled_ms_onsets = [0, diff(ms_present_all{subj}{block}{trial}(stim_on:button_press-300))];
                    unfilled_ms_onsets(unfilled_ms_onsets < 0) = 0;
                    nms_baseline = [nms_baseline; sum(unfilled_ms_onsets)];
                    curr_amplitudes = ms_amplitudes_all{subj}{block}{trial}(stim_on:button_press-300);
                    if sum(unfilled_ms_onsets)
                        meanmsampl_baseline = [meanmsampl_baseline; mean(curr_amplitudes(logical(unfilled_ms_onsets)))];
                    else
                        meanmsampl_baseline = [meanmsampl_baseline; nan];
                    end
                    % same for blinks
                    unfilled_blink_onsets = [0, diff(blink_present_all{subj}{block}{trial}(stim_on:button_press-300))];
                    unfilled_blink_onsets(unfilled_blink_onsets < 0) = 0;
                    nblinks_baseline = [nblinks_baseline; sum(unfilled_blink_onsets)];
%                     % get probability of filling report after microsaccade onset
%                     all_ms_onsets = [0, diff(ms_present_all{subj}{block}{trial}(stim_on:button_press+trans_ntps_before))];
%                     shuffled_ms_onsets = all_ms_onsets(randperm(length(all_ms_onsets)));
%                     ms_onset_indices = find(all_ms_onsets);
%                     shuffled_ms_onset_indices = find(shuffled_ms_onsets);
%                     trial_transition_probability = zeros(trans_ntps_before+1+trans_ntps_after, numel(ms_onset_indices));
%                     shuffled_trial_transition_probability = trial_transition_probability;
%                     for ms = 1:numel(ms_onset_indices)
%                         ms_button_latency = numel(all_ms_onsets) - trans_ntps_before - ms_onset_indices(ms);
%                         if ms_button_latency <= trans_ntps_after && ms_button_latency >= -trans_ntps_before % if within window we're plotting
%                             if norm_hazard
%                         % normalize by hazard rate. Accounts for increasing hazard with time, and hazard difference between stimulus conditions.
%                         hazard_times = xouts(:, subj, block); % bin centers for hazard function
%                         [~, hazard_idx] = min(abs(hazard_times - ms_onset_indices(ms))); % which hazard bin is this microsaccade contained in
%                         %hazard_value = lambdas(hazard_idx, subj, block); if isinf(hazard_value) || ~hazard_value, hazard_value = nan; end
%                         hazard_value = cdfs(hazard_idx, subj, block);
%                         trial_transition_probability(ms_button_latency + trans_ntps_before + 1, ms) = 1 * hazard_value;
%                             else
%                                 trial_transition_probability(ms_button_latency + trans_ntps_before + 1, ms) = 1;
%                             end
%                         end
%                         shuffled_ms_button_latency = numel(all_ms_onsets) - trans_ntps_before - shuffled_ms_onset_indices(ms);
%                         if shuffled_ms_button_latency <= trans_ntps_after && shuffled_ms_button_latency >= -trans_ntps_before % if within window we're plotting
%                         if norm_hazard
%                         [~, shuffled_hazard_idx] = min(abs(hazard_times - shuffled_ms_onset_indices(ms))); % which hazard bin is this microsaccade contained in
%                         %shuffled_hazard_value = lambdas(shuffled_hazard_idx, subj, block); if isinf(shuffled_hazard_value) || ~shuffled_hazard_value, shuffled_hazard_value = nan; end
%                         shuffled_hazard_value = cdfs(shuffled_hazard_idx, subj, block);
%                         shuffled_trial_transition_probability(shuffled_ms_button_latency + trans_ntps_before + 1, ms) = 1 * shuffled_hazard_value;
%                         else
%                             shuffled_trial_transition_probability(shuffled_ms_button_latency + trans_ntps_before + 1, ms) = 1;
%                         end
%                         end
%                     end
%                     block_transition_probability = [block_transition_probability, trial_transition_probability];
%                     baseline_block_transition_probability = [baseline_block_transition_probability, shuffled_trial_transition_probability];
                    block_pupil_trace(:, trial) = [nan(1, nan_pad_start), pupil_trace_all{subj}{block}{trial}(first_idx:last_idx), nan(1, nan_pad_end)];
                    block_ms_amplitudes(:, trial) = [nan(1, nan_pad_start), ms_amplitudes_all{subj}{block}{trial}(first_idx:last_idx), nan(1, nan_pad_end)];
                    if strcmp(reference_point, 'button_press') % get immobilization duration and last ms amplitude.
                        %ms_present_idx = find(block_ms_present(:, trial) == 1) - ntps_before; % idx of 1, or nan (NOT 0 (immobile) !), relative to button press (0)
                        unfilled_ms_offsets = [0, diff(ms_present_all{subj}{block}{trial}(stim_on:button_press))];
                        unfilled_ms_offsets(unfilled_ms_offsets > 0) = 0;
                        ms_offset_idx = find(unfilled_ms_offsets) - numel(unfilled_ms_offsets);
                        last_ms_time = max(ms_offset_idx);
                        curr_immobilization_duration = max(ms_offset_idx(ms_offset_idx < -300));
                        %last_ms_time = max(ms_present_idx(ms_present_idx < 1)); % actual last tp of ongoing microsaccade before button press.
                        %curr_immobilization_duration = max(ms_present_idx(ms_present_idx < -300));%curr_subj_reactiontime)); % this just gets the last ms before the subject's reaction time, but aligned to t=0. So need to regress out subject random effects of what t=0 means
                        % if it's from a nan = either a blink or stim onset.
                        % now I set it to == 1, so we are only looking at microsaccades.
                        % and we want the last microsaccade to actually
                        % finish before -300 ms, no microsaccades between -300 and 0.
                        % and we want trials with no blinks (nans) in between the
                        % last microsaccade and button press:
                        %if isempty(curr_immobilization_duration) || sum(ms_present_idx >= -300 & ms_present_idx <= 0) || sum(isnan(block_ms_present((curr_immobilization_duration:0)+ntps_before, trial)))
                        if isempty(curr_immobilization_duration) || sum(ms_offset_idx >= -300 & ms_offset_idx <= 0) || sum(isnan(ms_present_all{subj}{block}{trial}((curr_immobilization_duration:0)+button_press)))
                            curr_immobilization_duration = nan;
                            immob_ms_amplitude = nan;
                        else
                            %immob_ms_idx = curr_immobilization_duration + ntps_before; % restore to timepoint idx (aka in milliseconds)
                            %immob_ms_amplitude = block_ms_amplitudes(immob_ms_idx, trial);
                            immob_ms_idx = curr_immobilization_duration + button_press; % restore to timepoint idx (aka in milliseconds)
                            immob_ms_amplitude = ms_amplitudes_all{subj}{block}{trial}(immob_ms_idx-1);
                        end
                        if isempty(last_ms_time)
                            last_ms_time = nan;
                            last_ms_amplitude = nan;
                        else
                            %last_ms_idx = last_ms_time + ntps_before;
                            %last_ms_amplitude = block_ms_amplitudes(last_ms_idx, trial);
                            last_ms_idx = last_ms_time + button_press;
                            last_ms_amplitude = ms_amplitudes_all{subj}{block}{trial}(last_ms_idx-1);
                        end
                        immobilization_durations = [immobilization_durations; curr_immobilization_duration];
                        immob_ms_amplitudes = [immob_ms_amplitudes; immob_ms_amplitude];
                        last_ms_times = [last_ms_times; last_ms_time];
                        last_ms_amplitudes = [last_ms_amplitudes; last_ms_amplitude];
                    end
                    
                    extranan = buffer_baseline_idx(baseline_sub_idx(end)) - numel(ms_present_all{subj}{block}{trial});
                    curr_pupil_baseline = [pupil_trace_all{subj}{block}{trial}, nan(1, extranan)];
                    baseline_block_pupil_trace(:, trial) = curr_pupil_baseline(buffer_baseline_idx(baseline_sub_idx));
                    curr_amp_baseline = [ms_amplitudes_all{subj}{block}{trial}, nan(1, extranan)]; curr_amp_baseline((button_press-2000):end) = nan;
                    baseline_block_ms_amplitudes(:, trial) = curr_amp_baseline(buffer_baseline_idx(baseline_sub_idx));
                    block_xyVelocity(:, trial) = [nan(1, nan_pad_start), nanmean(xyVelocity_all{subj}{block}{trial}(:, first_idx:last_idx)), nan(1, nan_pad_end)];
                    block_xyVelocity(block_ms_present(:, trial) > 0, trial) = nan; % remove microsaccade timepoints
                    curr_xyVelocity_baseline = [nanmean(xyVelocity_all{subj}{block}{trial}, 1), nan(1, extranan)]; curr_xyVelocity_baseline((button_press-2000):end) = nan;
                    baseline_block_xyVelocity(:, trial) = curr_xyVelocity_baseline(buffer_baseline_idx(baseline_sub_idx));
                    baseline_block_xyVelocity(baseline_block_ms_present(baseline_sub_idx, trial) > 0, trial) = nan; % remove microsaccade timepoints
                    %                    block_ms_counts(:, trial) = calculate_ms_rate(block_ms_present(:, trial), rate_window_length, stepsize);
                    %                    baseline_block_ms_counts(:, trial) = calculate_ms_rate(baseline_block_ms_present(:, trial), rate_window_length, stepsize);
                    trial_xPos = [nan(1, nan_pad_start), nanmean(curr_posdata(2:3, first_idx:last_idx), 1), nan(1, nan_pad_end)];
                    trial_xPos = (trial_xPos - nanmean(trial_xPos)) * degperpix;
                    trial_yPos = [nan(1, nan_pad_start), nanmean(curr_posdata(4:5, first_idx:last_idx), 1), nan(1, nan_pad_end)];
                    trial_yPos = (trial_yPos - nanmean(trial_yPos)) * degperpix;
                    %                     trial_xyPos = hypot(trial_xPos, trial_yPos); trial_xyPos = trial_xyPos(~isnan(trial_xyPos));
                    curr_trial_xPos_baseline = [nanmean(posdata_all{subj}{block}{trial}(2:3, :), 1), nan(1, extranan)]; curr_trial_xPos_baseline((button_press-2000):end) = nan;
                    baseline_block_xPos(:, trial) = curr_trial_xPos_baseline(buffer_baseline_idx(baseline_sub_idx)); baseline_block_xPos(:, trial) = (baseline_block_xPos(:, trial) - nanmean(baseline_block_xPos(:, trial))) * degperpix;
                    curr_trial_yPos_baseline = [nanmean(posdata_all{subj}{block}{trial}(4:5, :), 1), nan(1, extranan)]; curr_trial_yPos_baseline((button_press-2000):end) = nan;
                    baseline_block_yPos(:, trial) = curr_trial_yPos_baseline(buffer_baseline_idx(baseline_sub_idx)); baseline_block_yPos(:, trial) = (baseline_block_yPos(:, trial) - nanmean(baseline_block_yPos(:, trial))) * degperpix;
                    %                     N = ntimepoints-1;
                    %                     trial_fft = (1/N)*fft(trial_xyPos, N)';
                    %                     trial_fft = trial_fft(1:round(N/2)+1); % remove negative frequencies
                    %                     trial_fft = abs(trial_fft); % keep only real amplitudes
                    %                     trial_fft(2:end) = trial_fft(2:end) .* 2; % "make one-sided"
                    %                     frqs = 1000*(0:round(N/2))/N;
                    %                     block_ten50power(trial) = sum(trial_fft(frqs>10 & frqs<50));
                    %                     block_eighty200power(trial) = sum(trial_fft(frqs>80 & frqs<200));
                %                          block_Nb(:, trial) = retinal_slip(trial_xPos, trial_yPos, block_ms_present(:, trial), 0); % last arg 1: include microsaccades
                    min_slip_window = [min_slip_window; nanmean(block_Nb(ntps_before/50+Nb_min_window+1, trial))];
                %                           baseline_block_Nb(:, trial) = retinal_slip(baseline_block_xPos(:, trial), baseline_block_yPos(:, trial), baseline_block_ms_present(baseline_sub_idx, trial), 0);
                    fulltrial_xPos = nanmean(curr_posdata(2:3, :), 1); fulltrial_yPos = nanmean(curr_posdata(4:5, :), 1);
                    fulltrial_xPos = (fulltrial_xPos - nanmean(fulltrial_xPos)) * degperpix; fulltrial_yPos = (fulltrial_yPos - nanmean(fulltrial_yPos)) * degperpix;
                %    fulltrial_block_Nb = retinal_slip(fulltrial_xPos, fulltrial_yPos, ms_present_all{subj}{block}{trial}', 0);
                    %baseline_retinal_slip_trials = [baseline_retinal_slip_trials; nanmean(baseline_block_Nb(:, trial))];
                %    baseline_retinal_slip_trials = [baseline_retinal_slip_trials; nanmean(fulltrial_block_Nb(floor(stim_on/50):floor((button_press-300)/50)))];
                %    if isnan(baseline_retinal_slip_trials(end))
                %        error('tt')
                %    end
                    %                     fulltrial_xPos = [nan(1, nan_pad_start_fulltrial), nanmean(posdata_all{subj}{block}{trial}(2:3, first_idx_fulltrial:last_idx_fulltrial), 1), nan(1, nan_pad_end_fulltrial)];
                    %                     fulltrial_xPos = (fulltrial_xPos - nanmean(fulltrial_xPos)) * degperpix;
                    %                     fulltrial_yPos = [nan(1, nan_pad_start_fulltrial), nanmean(posdata_all{subj}{block}{trial}(4:5, first_idx_fulltrial:last_idx_fulltrial), 1), nan(1, nan_pad_end_fulltrial)];
                    %                     fulltrial_yPos = (fulltrial_yPos - nanmean(fulltrial_yPos)) * degperpix;
                    %                     fulltrial_xyPos = hypot(fulltrial_xPos, fulltrial_yPos);
                    
                    
                    if baseline_correct
                        % pupil
                        T = block_pupil_trace(:, trial);
                        A = pupil_trace_all{subj}{block}{trial}(first_idx_fulltrial:last_idx_fulltrial); % normalize against whole trial
                        curr_pupil_baseline = [pupil_trace_all{subj}{block}{trial}, nan(1, extranan)];
                        %A = curr_pupil_baseline(buffer_baseline_idx(baseline_sub_idx));% normalize against baseline
                        % z score ignoring nans
                        %block_pupil_trace(:, trial) = (T - nanmean(A)) ./ nanstd(A);
                        % robust z score (now done in analyze_eyelink_fulltrial already)
                        %block_pupil_trace(:, trial) = (T - nanmedian(A)) ./ nanmedian(abs(T - nanmedian(A)));
                        % percent change from baseline idx
                        %block_pupil_trace(:, trial) = (T - nanmean(A)) ./ nanmean(A) .* 100;
                        % just zero mean within epoch
                        %block_pupil_trace(:, trial) = T - nanmean(T);
                        % zero mean to surrounding button press
                        %block_pupil_trace(:, trial) = T - nanmean(T(end-3000:end-1000)); % -500 ms to 0 ms
                        % subtract fixation baseline
                        %A = pupil_trace_all{subj}{block}{trial}(800:1000); % last 200 ms of fixation
                        %block_pupil_trace(:, trial) = T - nanmean(A);
                        % subtract baseline = first 50ms following stimulus onset
                        %A = pupil_trace_all{subj}{block}{trial}(1001:1050);
                        %block_pupil_trace(:, trial) = T - nanmean(A);
                        block_pupil_trace(:, trial) = (T - nanmean(A)) ./ nanstd(A); % z
                        
                        % amplitude, z score ignoring nans
                        %A = ms_amplitudes_all{subj}{block}{trial}(first_idx_fulltrial:last_idx_fulltrial);
                        %block_ms_amplitudes(:, trial) = (block_ms_amplitudes(:, trial) - nanmean(A)) ./ nanstd(A);
                        
                        % xy velocity, z score ignoring nans
                        %A = xyVelocity_all{subj}{block}{trial}(first_idx_fulltrial:last_idx_fulltrial);
                        %block_xyVelocity(:, trial) = (block_xyVelocity(:, trial) - nanmean(A)) ./ nanstd(A);
                        
                        % rate, z score ignoring nans
                        %                        A = calculate_ms_rate(fulltrial_block_ms_present(:, trial), rate_window_length, stepsize);
                        %                        block_ms_counts(:, trial) = (block_ms_counts(:, trial) - nanmean(A)) ./ nanstd(A);
                        
                        % Nb, z score ignoring nans
                        %                        A = retinal_slip(fulltrial_xPos, fulltrial_yPos, fulltrial_block_ms_present(:, trial));
                        %                        block_Nb(:, trial) = (block_Nb(:, trial) - nanmean(A)) ./ nanstd(A);
                    end
                end
            end
            baseline_probability(subj, block) = nanmean(nanmean(baseline_block_ms_present, 2)); % avg over second 600 ms of fixation
            ms_probability_all(:, subj, block) = nanmean(block_ms_present, 2);
            %ms_counts_all(:, subj, block) = nanmean(block_ms_counts, 2);%calculate_ms_rate(block_ms_present, window_length, stepsize);
            %ms_counts_all(:, subj, block) = calculate_ms_rate(nansum(block_ms_present, 2), rate_window_length, stepsize);
            ms_counts_all(:, subj, block) = cross_trial_ms_rate(block_ms_present, 'alpha');
            blink_rate_all(:, subj, block) = cross_trial_ms_rate(block_blink_present, 'alpha');
            %[~, ms_counts_all(:, subj, block)] = cross_trial_ms_rate(block_ms_present, 'alpha'); % investigate shuffled ms instead.
            %baseline_ms_counts(subj, block) = nanmean(nanmean(baseline_block_ms_counts, 2));
            %curr_baseline_r = cross_trial_ms_rate(baseline_block_ms_present, 'alpha');
            %curr_baseline_r = curr_baseline_r(baseline_sub_idx);
            %[~, curr_baseline_r] = cross_trial_ms_rate(fulltrial_block_ms_present, 'alpha'); % shuffled
            mean_amplitudes_all(:, subj, block) = nanmean(block_ms_amplitudes, 2);
            mean_xyVelocity_all(:, subj, block) = nanmean(block_xyVelocity, 2);
            %baseline_ms_counts(subj, block) = nanmean(curr_baseline_r);
            baseline_ms_counts(subj, block) = nanmean(ms_counts_all(1:2000, subj, block)); % first 2000 of epoch = baseline.
            curr_baseline_blinks = cross_trial_ms_rate(baseline_block_blink_present, 'alpha');
            curr_baseline_blinks = curr_baseline_blinks(baseline_sub_idx);
            baseline_blink_rate(subj, block) = nanmean(curr_baseline_blinks);
            %baseline_blink_rate(subj, block) = nanmedian(blink_rate_all(1:2000, subj, block));
            mean_pupil_all(:, subj, block) = nanmean(block_pupil_trace, 2);
            %mean_pupil_all(:, subj, block) = nanmean((block_pupil_trace - nanmedian(mean_pupil_all(:, subj, block))) ./ nanmedian(abs(block_pupil_trace - nanmedian(mean_pupil_all(:, subj, block)))), 2); % robust z scored from block mean
            mean_Nb_all(:, subj, block) = nanmean(block_Nb, 2);
            %baseline_Nb(subj, block) = nanmean(nanmedian(baseline_block_Nb)); % median across time, mean across trials
            baseline_Nb(subj, block) = nanmedian(mean_Nb_all(1:2000/50, subj, block)); % first 2000ms of epoch = baseline
            %baseline_amplitudes(subj, block) = nanmean(nanmean(baseline_block_ms_amplitudes));
            baseline_amplitudes(subj, block) = nanmean(mean_amplitudes_all(1:2000, subj, block)); % first 2000 of epoch = baseline.
            %baseline_xyVelocity(subj, block) = nanmean(nanmean(baseline_block_xyVelocity));
            baseline_xyVelocity(subj, block) = nanmean(mean_xyVelocity_all(1:3000, subj, block));
            %ms_present_concat = [ms_present_concat, block_ms_present];
            try transition_probabilities_all(:, subj, block) = nanmean(block_transition_probability, 2); end
            try baseline_transition_probabilities(:, subj, block) = nanmean(baseline_block_transition_probability, 2); end
            if baseline_correct
                % probability, z score ignoring nans
                %A = nanmean(fulltrial_block_ms_present, 2);
                %ms_probability_all(:, subj, block) = (ms_probability_all(:, subj, block) - nanmean(A)) ./ nanstd(A);
                
                % rate, z score ignoring nans
                %A = cross_trial_ms_rate(fulltrial_block_ms_present, 'alpha');
                %ms_counts_all(:, subj, block) = (ms_counts_all(:, subj, block) - nanmean(A)) ./ nanstd(A);
                
                % rate, robust z-score
                %A = cross_trial_ms_rate(fulltrial_block_ms_present, 'alpha');
                %ms_counts_all(:, subj, block) = (ms_counts_all(:, subj, block) - nanmedian(A)) ./ nanmedian(abs(ms_counts_all(:, subj, block) - nanmedian(A)));
                
                % rate, normalize by baseline
                %A = cross_trial_ms_rate(baseline_block_ms_present, 'alpha');
                %A = curr_baseline_r;
                %ms_counts_all(:, subj, block) = ms_counts_all(:, subj, block) ./ nanmean(A);
                %ms_counts_all(:, subj, block) = (ms_counts_all(:, subj, block) - nanmean(A)) ./ nanstd(A);
                %ms_counts_all(:, subj, block) = (ms_counts_all(:, subj, block) - nanmedian(A)) ./ nanmedian(abs(ms_counts_all(:, subj, block) - nanmedian(A)));
            end
            
            
            %baseline_pupil(subj, block) = nanmean(nanmean(baseline_block_pupil_trace));
            baseline_pupil(subj, block) = nanmean(mean_pupil_all(1:3000, subj, block)); % first 2000 of epoch = baseline.
            %baseline_pupil(subj, block) = nanmean(mean_pupil_all(:, subj, block)); % baseline = mean of whole epoch.
            %plotpupil = block_pupil_trace;
            %plot(plotpupil); hold on;
            %             mean_ten50power_all(subj, block) = nanmean(block_ten50power);
            %             mean_eighty200power_all(subj, block) = nanmean(block_eighty200power);
            %ms_probability_all(:, subj, block) = (nanmean(ms_present_all{subj}{block}, 2));% - baseline_probability(subj,block)) ./ baseline_probability(subj, block) .* 100;
            mean_rts_all(subj, block) = nanmean(block_rts);
            % weigh by trial counts?
            weighting(subj, block) = ntrials_all(subj, block) ./ ntrialstotal_byblock(block) .* nsubjects;
            if weighsubjects
                ms_probability_all(:, subj, block) = ms_probability_all(:, subj, block) .* weighting(subj, block);
                ms_counts_all(:, subj, block) = ms_counts_all(:, subj, block) .* weighting(subj, block);
                mean_pupil_all(:, subj, block) = mean_pupil_all(:, subj, block) .* weighting(subj, block);
                mean_amplitudes_all(:, subj, block) = mean_amplitudes_all(:, subj, block) .* weighting(subj, block);
                mean_xyVelocity_all(:, subj, block) = mean_xyVelocity_all(:, subj, block) .* weighting(subj, block);
                baseline_probability(subj, block) = baseline_probability(subj, block) .* weighting(subj, block);
                baseline_pupil(subj, block) = baseline_pupil(subj, block) .* weighting(subj, block);
                baseline_ms_counts(subj, block) = baseline_ms_counts(subj, block) .* weighting(subj, block);
                mean_Nb_all(:, subj, block) = mean_Nb_all(:, subj, block) .* weighting(subj, block);
            end
            ntrials_used(:, subj, block) = sum(~isnan(block_ms_present), 2);
        else
        end
    end
end

% normalize everything by subject baseline means
norm_ms_counts_all = nan(size(ms_counts_all));
norm_blink_rate_all = nan(size(blink_rate_all));
norm_Nb_all = nan(size(mean_Nb_all));
norm_amplitudes_all = nan(size(mean_amplitudes_all));
norm_xyVelocity_all = nan(size(mean_xyVelocity_all));
norm_baseline_ms_counts = nan(size(baseline_ms_counts));
norm_baseline_Nb = nan(size(baseline_Nb));
percchange_baseline_ms_counts = nan(size(baseline_ms_counts));
norm_rts_all = nan(size(mean_rts_all));
norm_pupil_all = nan(size(mean_pupil_all));
for subj = 1:nsubjects
    norm_ms_counts_all(:, subj, :) = (ms_counts_all(:, subj, :) ./ nanmean(baseline_ms_counts(subj, :))) .* nanmean(nanmean(baseline_ms_counts));
    %norm_ms_counts_all(:, subj, :) = (ms_counts_all(:, subj, :) - nanmean(baseline_ms_counts(subj, :))) ./ nanmean(baseline_ms_counts(subj, :)) .* 100;
    %norm_ms_counts_all(:, subj, :) = ms_counts_all(:, subj, :);
    %norm_ms_counts_all(:, subj, :) = ms_counts_all(:, subj, :) - nanmean(baseline_ms_counts(subj, :)); % difference
    %norm_blink_rate_all(:, subj, :); = blink_rate_all(:, subj, :);
    norm_blink_rate_all(:, subj, :) = (blink_rate_all(:, subj, :) ./ nanmean(baseline_blink_rate(subj, :))) .* nanmean(nanmean(baseline_blink_rate));
    norm_Nb_all(:, subj, :) = mean_Nb_all(:, subj, :) ./ nanmean(baseline_Nb(subj, :));
    %norm_Nb_all(:, subj, :) = mean_Nb_all(:, subj, :);
    norm_amplitudes_all(:, subj, :) = mean_amplitudes_all(:, subj, :) ./ nanmean(baseline_amplitudes(subj, :));
    %norm_amplitudes_all(:, subj, :) = mean_amplitudes_all(:, subj, :);
    norm_xyVelocity_all(:, subj, :) = mean_xyVelocity_all(:, subj, :) ./ nanmean(baseline_xyVelocity(subj, :));
    %norm_xyVelocity_all(:, subj, :) = mean_xyVelocity_all(:, subj, :);
    %norm_baseline_ms_counts(subj, :) = baseline_ms_counts(subj, :);
    norm_baseline_ms_counts(subj, :) = (baseline_ms_counts(subj, :) ./ nanmean(baseline_ms_counts(subj, :))) .* nanmean(nanmean(baseline_ms_counts));
    %norm_baseline_ms_counts(subj, :) = baseline_ms_counts(subj, :) - nanmean(baseline_ms_counts(subj, :)); % difference
    norm_baseline_amplitudes(subj, :) = baseline_amplitudes(subj, :) ./ nanmean(baseline_amplitudes(subj, :));
    norm_baseline_xyVelocity(subj, :) = baseline_xyVelocity(subj, :) ./ nanmean(baseline_xyVelocity(subj, :));
    norm_baseline_Nb(subj, :) = baseline_Nb(subj, :) ./ nanmean(baseline_Nb(subj, :));
    percchange_baseline_ms_counts(subj, :) = (baseline_ms_counts(subj, :) - nanmean(baseline_ms_counts(subj, :))) ./ nanmean(baseline_ms_counts(subj, :));
    norm_rts_all(subj, :) = mean_rts_all(subj, :) ./ nanmean(mean_rts_all(subj, :));
    norm_pupil_all(:, subj, :) = mean_pupil_all(:, subj, :);
    %norm_pupil_all(:, subj, :) = mean_pupil_all(:, subj, :) ./ nanmean(baseline_pupil(subj, :)); % div normalization by subject mean
    %norm_pupil_all(:, subj, :) = mean_pupil_all(:, subj, :) - reshape(baseline_pupil(subj, :), 1, 1, 9); % zero mean per block
    %norm_pupil_all(:, subj, :) = mean_pupil_all(:, subj, :) ./ reshape(baseline_pupil(subj, :), 1, 1, 9); % div normalization per block
end

all_ntrials_used = nansum(nansum(ntrials_used, 3), 2);

% maybe include this again at some point (need to rewrite a little I think)
%if strcmp(reference_point, 'button_press') % if centering around button press, get -500ms value
%                before_press_probability(subj, block) = nanmean(nanmean(block_ms_present(bp500_idx, :), 2)); % avg for -700 ms : -300 ms
%            end

%% plot averaged by size or contrast
% average contrasts together or averaged sizes together
size_indices = [3 2 1 3 2 1 3 2 1]; % 1=size 6, 2=size 4, 3=size 2
color_indices = [1 1 1 2 2 2 3 3 3]; % 1=low contrast, 2=medium, 3=high contrast
all_indices = [1 1 1 1 1 1 1 1 1];

which2use = 'separate'; % size or color or all (average across all blocks) or separate (plot each condition)
which2plot = 'rate'; % probability, rate, amplitude, pupil, velocity, retinal slip
subjects2use = 1:nsubjects;%[1:10, 13:nsubjects];%[1:nsubjects]; % normally all 
baseline_correct = 0; % baseline correct each block? yes or no
shift_subjs = 0; % shift each subjects timeseries by their avg minimum (to somewhat account for variable rts)
shift_subjs_block = 0; % same but each block's shift calculated on the other 8
overlay = 1; % overlay the plots on one axis or use subplots
plotcolors = {[1 0 0], [0 1 0], [0 0 1]};
nobounds = 0; % just plot line, no stderr bounds?
weighsubjs = 0; % weigh subjects by trial counts?

% only plot subjects that have highbaseline > lowbaseline
if strcmp(which2use, 'color')
    baseline_sequence_subjs = nanmean(baseline_ms_counts(:, color_indices==3), 2) > nanmean(baseline_ms_counts(:, color_indices==1), 2);
elseif strcmp(which2use, 'size')
    baseline_sequence_subjs = nanmean(baseline_ms_counts(:, size_indices==3), 2) > nanmean(baseline_ms_counts(:, size_indices==1), 2);
end
%subjects2use = find(~baseline_sequence_subjs);

figure; hold on;
%for subj = 1:nsubjects
%    subjects2use = subj;
%figure; hold on;
for i = 1:9
    switch which2use
        case 'size'
            if overlay; hold on; else; subplot(1, 3, i); end
            if size(ms_counts_all, 3) == 9 % if all blocks epoched separately
                curr_block_idx = find(size_indices == i); % if plotting by size
            elseif size(ms_counts_all, 3) == 3 % if jointly combined blocks within size
                curr_block_idx = i;
            end
        case 'color'
            if overlay; hold on; else; subplot(1, 3, i); end
            if size(ms_counts_all, 3) == 9 % if all blocks epoched separately
                curr_block_idx = find(color_indices == i); % if plotting by color
            elseif size(ms_counts_all, 3) == 3 % if jointly combined blocks within color
                curr_block_idx = i;
            end
        case 'all'
            if size(ms_counts_all, 3) == nblocks % if all blocks epoched separately
                curr_block_idx = find(all_indices);
            elseif size(ms_counts_all, 3) == 3 % if jointly combined across 1 dimension
                curr_block_idx = [1 2 3];
            end
            plotcolors = {[0 0 0]};
        case 'separate'
            if overlay; hold on; else; subplot(3, 3, i); end
            curr_block_idx = i;
            plotcolors = {[0 0 .33], [0 .33 0], [.33 0 0], [0 0 .66], [0 .66 0], [.66 0 0], [0 0 1], [0 1 0], [1 0 0]}; % colors = sizes
            %plotcolors = {[.33 0 0], [.66 0 0], [1 0 0], [0 .33 0], [0 .66 0], [0 1 0], [0 0 .33], [0 0 .66], [0 0 1]}; % colors = contrasts
    end
    
    switch which2plot
        case 'probability'
            % plot ms probability
            if baseline_correct
                ms_probability_uncorr = ms_probability_all(:, subjects2use, curr_block_idx);
                curr_prob_baseline = baseline_probability(subjects2use, curr_block_idx);
                ms_probability_corr = nan(size(ms_probability_uncorr));
                for ss = 1:size(ms_probability_uncorr, 2) % per subject
                    for bb = 1:size(ms_probability_uncorr, 3) % per block
                        curr_uncorr = ms_probability_uncorr(:, ss, bb);
                        % baseline correct
                        curr_base = curr_prob_baseline(ss, bb); % from computed baseline
                        %curr_base = nanmean(curr_uncorr(1:1000)); % from this actual epoch
                        %ms_probability_corr(:, ss, bb) = (curr_uncorr - curr_base) ./ curr_base .* 100; % percent change
                        ms_probability_corr(:, ss, bb) = curr_uncorr - curr_base; % just subtract baseline
                        % z score
                        %ms_probability_corr(:, ss, bb) = (curr_uncorr - nanmean(curr_uncorr)) ./ nanstd(curr_uncorr);
                        % and weigh by trial counts
                        % ms_probability_corr(:, ss, bb) = ms_probability_corr(:, ss, bb) .* ntrials_all(ss, bb) ./ ntrialstotal_byblock(bb) .* nsubjects;
                    end
                end
                ms_probability_corr(isinf(ms_probability_corr)) = nan;
                ms_probability_curr = nanmean(ms_probability_corr, 3);
            else
                ms_probability_curr = nanmean(ms_probability_all(:, subjects2use, curr_block_idx), 3);
            end
            mean_ms_probability = nanmean(ms_probability_curr, 2); % average by subject
            stderr_ms_probability = nanstd(ms_probability_curr, 0, 2) ./ sqrt(nsubjects);
            smoothed_prob = smoothdata(mean_ms_probability, 'movmean', window_length);
            smoothed_stderr = smoothdata(stderr_ms_probability, 'movmean', window_length);
            %plot((-ntps_before:ntps_after)' * 1000/fps, smoothed_prob);
            bp = boundedline((-ntps_before:ntps_after)' * 1000/fps, smoothed_prob, smoothed_stderr, 'cmap', plotcolors{i}, 'alpha', 'transparency', 0.05);
            %ylim([0.0, 0.06]);
            %ylim([-50 100]);
            %ylim([-.5 .5]);
            
        case 'rate'
            if baseline_correct
                ms_counts_uncorr = norm_ms_counts_all(:, subjects2use, curr_block_idx);
                curr_ms_counts_baseline = norm_baseline_ms_counts(subjects2use, curr_block_idx);
                %curr_weights = weighting(subjects2use, curr_block_idx);
                ms_counts_corr = nan(size(ms_counts_uncorr));
                for bb = 1:size(ms_counts_uncorr, 3) % per block
                    for ss = 1:size(ms_counts_uncorr, 2) % per subject
                        curr_uncorr = ms_counts_uncorr(:, ss, bb);
                        % baseline correct
                        curr_base = curr_ms_counts_baseline(ss, bb); % from computed baseline
                        %curr_base = nanmean(curr_uncorr(1:(1000/stepsize))); % from this actual epoch
                        %ms_counts_corr(:, ss, bb) = (curr_uncorr - curr_base) ./ curr_base .* 100; % percent change
                        ms_counts_corr(:, ss, bb) = curr_uncorr - curr_base; % just subtract baseline
                        %ms_counts_corr(:, ss, bb) = curr_uncorr ./ curr_base; % normalize by baseline (base=1)
                        % z score
                        %ms_counts_corr(:, ss, bb) = (curr_uncorr - nanmean(curr_uncorr)) ./ nanstd(curr_uncorr);
                    end
                    % normalize by subject mean of pre-stim baseline (per block).
                    %curr_uncorr = ms_counts_uncorr(:, :, bb);
                    %curr_base = nanmean(curr_ms_counts_baseline(:, bb)); % from computed baseline
                    %ms_counts_corr(:, :, bb) = curr_uncorr ./ curr_base;
                end
                ms_counts_corr(isinf(ms_counts_corr)) = nan;
                % make baseline = the total avg baseline, instead of 0
                ms_counts_corr = ms_counts_corr + nanmean(nanmean(norm_baseline_ms_counts));
            else
                norm_ms_counts_all(isinf(norm_ms_counts_all)) = nan;
                ms_counts_corr = norm_ms_counts_all(:, subjects2use, curr_block_idx);
                %ms_counts_corr = ms_counts_all(:, subjects2use, curr_block_idx); % not subject normalized
            end
            if weighsubjs
                curr_weights = weighting(subjects2use, curr_block_idx);
                for ss = 1:numel(subjects2use)
                    for bb = 1:numel(curr_block_idx)
                        ms_counts_corr(:, ss, bb) = ms_counts_corr(:, ss, bb) * curr_weights(ss, bb);
                    end
                end
            end
            if shift_subjs_block
                for ss = 1:numel(subjects2use)
                    sj = subjects2use(ss);
                    for bb = 1:numel(curr_block_idx)
                        bl = curr_block_idx(bb);
                        if all_subj_shifts(sj, bl) > 0 % shift rightward
                            ms_counts_corr(:, ss, bb) = [nan(all_subj_shifts(sj, bl), 1, 1); ms_counts_corr(1:end-all_subj_shifts(sj, bl), ss, bb)];
                        elseif all_subj_shifts(sj, bl) < 0 % shift leftward
                            ms_counts_corr(:, ss, bb) = [ms_counts_corr((abs(all_subj_shifts(sj, bl))+1):end, ss, bb); nan(abs(all_subj_shifts(sj, bl)), 1, 1)];
                        end
                    end
                end
            end
            ms_counts_curr = nanmean(ms_counts_corr, 3);
            %ms_counts_curr = smoothdata(ms_counts_curr, 'movmean', 500);
            % plot ms rate
            if shift_subjs
                % shift by data-driven subject rt estimate?
                % this is done by aligning minima of indiv subject mean ms rate
                % (between -800 and -300 ms) to -550 ms.
                %minimum_window = ntps_before - (800:-1:300);
                %[~, subj_min_idx] = min(ms_counts_curr(minimum_window, :)); % calculate using this
                %subj_shifts = 250 - subj_min_idx;
                %subj_shifts = [-30 91 -83 -251 47 -98 221 -251 -28 -250 189 249 249 -250 249 47 76 240 -48 136]; bw 250-750ms, no extra smoothing
                %subj_shifts = [ -125  -251   -26  -251   -65  -251   249  -251    50  -251    78   249   242  -251   249   -53   137   249  -170     4  -251]; % after 500ms smoothing, bw 300-800ms
                for ss = 1:numel(subjects2use)
                    if subj_shifts(ss) > 0 % shift rightward
                        ms_counts_curr(:, ss) = [nan(subj_shifts(ss), 1); ms_counts_curr(1:end-subj_shifts(ss), ss)];
                    elseif subj_shifts(ss) < 0 % shift leftward
                        ms_counts_curr(:, ss) = [ms_counts_curr((abs(subj_shifts(ss))+1):end, ss); nan(abs(subj_shifts(ss)), 1)];
                    end
                end
            end
            mean_ms_counts = nanmean(ms_counts_curr, 2);
            stderr_ms_counts = nanstd(ms_counts_curr, 0, 2) ./ sqrt(numel(subjects2use));
            if nobounds
                bp = plot(count_timepoints(1:end-250), mean_ms_counts(1:end-250), 'Color', plotcolors{i});
            else
                bp = boundedline(count_timepoints(1:end-250), mean_ms_counts(1:end-250), stderr_ms_counts(1:end-250), 'cmap', plotcolors{i}, 'alpha', 'transparency', 0.05);
            end
            %ylim([0, 3]);
            %ylim([-1.25 1]);
            ylbl = 'microsaccade rate (Hz)';%'normalized microsaccade rate (Hz)';
            
        case 'amplitude'
            % plot average ms amplitudes
            if baseline_correct
                ms_amplitude_uncorr = norm_amplitudes_all(:, subjects2use, curr_block_idx);
                curr_amp_baseline = norm_baseline_amplitudes(subjects2use, curr_block_idx);
                ms_probability_corr = nan(size(ms_amplitude_uncorr));
                for ss = 1:size(ms_amplitude_uncorr, 2) % per subject
                    for bb = 1:size(ms_amplitude_uncorr, 3) % per block
                        curr_uncorr = ms_amplitude_uncorr(:, ss, bb);
                        % baseline correct
                        curr_base = curr_amp_baseline(ss, bb); % from computed baseline
                        %curr_base = nanmean(curr_uncorr(1:1000)); % from this actual epoch
                        %ms_amplitude_corr(:, ss, bb) = (curr_uncorr - curr_base) ./ curr_base .* 100; % percent change
                        %ms_amplitude_corr(:, ss, bb) = curr_uncorr - curr_base; % just subtract baseline
                        ms_amplitude_corr(:, ss, bb) = curr_uncorr ./ curr_base; % just divide
                        % z score
                        %ms_amplitude_corr(:, ss, bb) = (curr_uncorr - nanmean(curr_uncorr)) ./ nanstd(curr_uncorr);
                    end
                end
                ms_amplitude_corr(isinf(ms_amplitude_corr)) = nan;
            else
                %ms_amplitude_corr = mean_amplitudes_all(:, subjects2use, curr_block_idx);
                ms_amplitude_corr = norm_amplitudes_all(:, subjects2use, curr_block_idx) .* nanmean(nanmean(baseline_amplitudes)); % normalize, then multiply by total mean to interpret it
            end
            if shift_subjs_block
                for ss = 1:numel(subjects2use)
                    sj = subjects2use(ss);
                    for bb = 1:numel(curr_block_idx)
                        bl = curr_block_idx(bb);
                        if all_subj_shifts(sj, bl) > 0 % shift rightward
                            ms_amplitude_corr(:, ss, bb) = [nan(all_subj_shifts(sj, bl), 1, 1); ms_amplitude_corr(1:end-all_subj_shifts(sj, bl), ss, bb)];
                        elseif all_subj_shifts(sj, bl) < 0 % shift leftward
                            ms_amplitude_corr(:, ss, bb) = [ms_amplitude_corr((abs(all_subj_shifts(sj, bl))+1):end, ss, bb); nan(abs(all_subj_shifts(sj, bl)), 1, 1)];
                        end
                    end
                end
            end
            ms_amplitude_curr = nanmean(ms_amplitude_corr, 3);
            mean_ms_amplitude = nanmean(ms_amplitude_curr, 2);
            stderr_ms_amplitude = nanstd(ms_amplitude_curr, 0, 2) ./ sqrt(nsubjects);
            smoothed_ampl = smoothdata(mean_ms_amplitude, 'movmean', window_length);
            smoothed_stderr_ampl = smoothdata(stderr_ms_amplitude, 'movmean', window_length);
            %plot((-ntps_before:ntps_after)' * 1000/fps, smoothed_ampl);
            bp = boundedline((-ntps_before:ntps_after)' * 1000/fps, smoothed_ampl, smoothed_stderr_ampl, 'cmap', plotcolors{i}, 'alpha', 'transparency', 0.05);
            %ylim([0.0, 0.05]);
            %ylim([-50 100]);
            %ylim([-1 1]);
            ylbl = 'm.s. amplitude (deg)';
            
        case 'pupil'
            % plot pupil area (both eyes)
            % baseline correction is already done per trial, during epoching.
            if baseline_correct == 1 % if I want to do per-block correction instead or also.
                mean_pupil_uncorr = mean_pupil_all(:, subjects2use, curr_block_idx);
                curr_pupil_baseline = baseline_pupil(subjects2use, curr_block_idx);
                mean_pupil_corr = nan(size(mean_pupil_uncorr));
                for ss = 1:size(mean_pupil_uncorr, 2) % per subject
                    for bb = 1:size(mean_pupil_uncorr, 3) % per block
                        curr_uncorr = mean_pupil_uncorr(:, ss, bb);
                        % baseline correct
                        %curr_base = curr_pupil_baseline(ss, bb); % from computed baseline
                        curr_base = nanmean(curr_uncorr(1:1000)); % from this actual epoch
                        %mean_pupil_corr(:, ss, bb) = (curr_uncorr - curr_base) ./ curr_base .* 100; % percent change
                        mean_pupil_corr(:, ss, bb) = curr_uncorr - curr_base; % just subtract baseline
                        % z score
                        %mean_pupil_corr(:, ss, bb) = (curr_uncorr - nanmean(curr_uncorr)) ./ nanstd(curr_uncorr);
                    end
                end
                mean_pupil_curr = nanmean(mean_pupil_corr, 3); % avg across blocks
            else
                %mean_pupil_curr = nanmean(mean_pupil_all(:, subjects2use, curr_block_idx), 3);
                mean_pupil_curr = nanmean(norm_pupil_all(:, subjects2use, curr_block_idx), 3);
            end
            mean_pupil_area = nanmean(mean_pupil_curr, 2);
            stderr_pupil_area = nanstd(mean_pupil_curr, 0, 2) ./ sqrt(nsubjects);
            smoothed_area = smoothdata(mean_pupil_area, 'movmean', window_length);
            smoothed_stderr_area = smoothdata(stderr_pupil_area, 'movmean', window_length);
            %bp = boundedline((-ntps_before:ntps_after)' * 1000/fps, smoothed_area, smoothed_stderr_area, plotcolors(i), 'cmap', plotcolors{i}, 'alpha', 'transparency', 0.05);
            bp = boundedline((-ntps_before:ntps_after-100)' * 1000/fps, mean_pupil_area(1:end-100), stderr_pupil_area(1:end-100), 'cmap', plotcolors{i}, 'alpha', 'transparency', 0.05);
            %ylim([1000 3000]);
            %ylim([-50 300]);
            ylbl = 'pupil size (a.u.)';
            
        case 'velocity'
            % plot xy velocity (both eyes)
            mean_xyVelocity_curr = nanmean(norm_xyVelocity_all(:, subjects2use, curr_block_idx), 3);
            mean_xyVelocity = nanmean(mean_xyVelocity_curr, 2);
            stderr_xyVelocity = nanstd(mean_xyVelocity_curr, 0, 2) ./ sqrt(nsubjects);
            smoothed_xyVelocity = smoothdata(mean_xyVelocity, 'movmean', window_length);
            smoothed_stderr_xyVelocity = smoothdata(stderr_xyVelocity, 'movmean', window_length);
            bp = boundedline((-ntps_before:ntps_after)' * 1000/fps, smoothed_xyVelocity, smoothed_stderr_xyVelocity, 'cmap', plotcolors{i}, 'alpha', 'transparency', 0.05);
            %ylim([-1 1]);
            
        case 'retinal slip'
            % plot retinal slip (nboxes of 0.01 length traversed in 50 ms steps)
            if baseline_correct
                mean_Nb_uncorr = norm_Nb_all(:, subjects2use, curr_block_idx);
                mean_Nb_corr = nan(size(mean_Nb_uncorr));
                mean_baselines = nan(numel(subjects2use, numel(curr_block_idx)));
                for ss = 1:size(mean_Nb_uncorr, 2) % per subject
                    for bb = 1:size(mean_Nb_uncorr, 3) % per block
                        curr_uncorr = mean_Nb_uncorr(:, ss, bb);
                        % baseline correct
                        %curr_base = nanmean(curr_uncorr(1:(2000/50))); % from this actual epoch
                        curr_base = norm_baseline_Nb(ss, bb);
                        mean_baselines(ss, bb) = curr_base;
                        %mean_Nb_corr(:, ss, bb) = (curr_uncorr - curr_base) ./ curr_base .* 100; % percent change
                        mean_Nb_corr(:, ss, bb) = curr_uncorr - curr_base; % just subtract baseline
                        %mean_Nb_corr(:, ss, bb) = curr_uncorr ./ curr_base; % just divide
                        %mean_Nb_corr(:, ss, bb) = curr_uncorr; % don't correct yet
                        % z score
                        %mean_Nb_corr(:, ss, bb) = (curr_uncorr - nanmean(curr_uncorr)) ./ nanstd(curr_uncorr);
                    end
                end
                mean_Nb_corr(isinf(mean_Nb_corr)) = nan;
                mean_Nb_curr = nanmean(mean_Nb_corr, 3);
                %mean_Nb_curr = mean_Nb_curr - nanmean(mean_baselines, 2)';
                mean_Nb_curr = mean_Nb_curr .* nanmean(nanmean(baseline_Nb));
            else
                mean_Nb_curr = nanmean(norm_Nb_all(:, subjects2use, curr_block_idx), 3) .* nanmean(nanmean(baseline_Nb)); % normalize, then multiply by total mean to interpret it
                %mean_Nb_curr = nanmean(mean_Nb_all(:, subjects2use, curr_block_idx), 3); % raw mean
            end
            mean_Nb = nanmean(mean_Nb_curr, 2);
            stderr_Nb = nanstd(mean_Nb_curr, 0, 2) ./ sqrt(nsubjects);
            Nb_window_length = window_length / 50;
            smoothed_mean_Nb = smoothdata(mean_Nb, 'movmean', Nb_window_length/stepsize);
            smoothed_stderr_Nb = smoothdata(stderr_Nb, 'movmean', Nb_window_length/stepsize);
            % convert Nb to dva? (multiply by 0.01 dva)
            smoothed_mean_Nb = smoothed_mean_Nb * 0.01; smoothed_stderr_Nb = smoothed_stderr_Nb * 0.01;
            bp = boundedline(Nb_timepoints, smoothed_mean_Nb, smoothed_stderr_Nb, 'cmap', plotcolors{i}, 'alpha', 'transparency', 0.05);
            ylbl = 'retinal slip (deg / 50ms)';
            
    end
    set(bp, 'LineWidth', 4);
    
    H = gca;
    H.LineWidth = 2;
    
    % plot vertical line at t=0
    curr_ylim = get(gca, 'ylim');
    hold on; zeroline = plot([0 0], curr_ylim, '--k');
    set(get(get(zeroline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    % and for baseline correction, horizontal line at t=0
    curr_xlim = get(gca, 'xlim');
    if baseline_correct
        if strcmp(which2plot, 'rate')
            %set(gca, 'ylim', [-.8 0.4]) % subtract baseline
            set(gca, 'ylim', [.3 1.4]); % make baseline = actual baseline
            %set(gca, 'ylim', [-60, 20]) % percent change
        end
        %baseline_line = 0;
        baseline_line = nanmean(nanmean(norm_baseline_ms_counts));
    elseif ~baseline_correct & numel(subjects2use) > 1
        switch which2plot
            case 'rate'
                set(gca, 'ylim', [0.3 1.3])
                baseline_line = nanmean(nanmean(norm_baseline_ms_counts));
            case 'amplitude'
                baseline_line = nanmean(nanmean(baseline_amplitudes));
                set(gca, 'ylim', [0.2 0.5]);
            case 'retinal slip'
                baseline_line = nanmean(nanmean(baseline_Nb)) * 0.01; % in dva
                set(gca, 'ylim', [0.15 0.2]);
            case 'pupil'
                baseline_line = 0;
                set(gca, 'ylim', [-50, 50]);%[-0.02, 0.02]);
            otherwise
                baseline_line = 1;
        end
    elseif numel(subjects2use) == 1
        switch which2plot
            case 'rate'
                baseline_line = nanmean(baseline_ms_counts(subj, :));
        end
    end
    hold on; horzline = plot(curr_xlim, [baseline_line baseline_line], '--k');
    
    set(get(get(horzline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    
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
        if ~baseline_correct
            ylabel([which2plot]);
        elseif baseline_correct
            ylabel([which2plot, ' change from baseline']);
        end
        xlabel('time (s)');
        if strcmp(which2use, 'all')
            set(bp, 'Color', [0 0 0]);
            %break % stop for loop
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
            sl = suplabel('center size', 'y');
        end
        set(sl, 'FontSize', 26);
    case 'color'
        if overlay
            legend({'low', 'medium', 'high'});
            sl = title([which2plot ' by color contrast']);
        else
            sl = suplabel('color contrast', 'y');
        end
        set(sl, 'FontSize', 26);
    case 'separate'
        if overlay
            legend({'small low', 'medium low', 'large low', 'small medium', 'medium medium', 'large medium', 'small high', 'medium high', 'large high'});
            sl = title([which2plot ' by condition']);
        else
            sly = suplabel('color contrast', 'y');
            slx = suplabel('center size', 'x');
        end
end
set(gca, 'FontSize', 26);

% set subject color if plotting just 1 subject at a time
if numel(subjects2use) == 1
    clrs = colormap;
    set(bp, 'Color', clrs(subjects2use * 3, :));
    set(bp, 'LineWidth', 1); % or make line smaller
    %set(bp, 'Color', [0.5 0.5 0.5]);
end
%set(gca, 'ylim', [-1000 1000]);

hold on;
xxlim = get(gca, 'xlim');
%plot(xxlim, [nanmean(baseline_xyVelocity(subj, :)), nanmean(baseline_xyVelocity(subj, :))])
%plot(xxlim, [nanmean(nanmean(baseline_Nb)), nanmean(nanmean(baseline_Nb))])

%end

switch reference_point
    case 'button_press'
%xlim([-3500, 750]);
end
ylabel(ylbl);
legend boxoff

%set(gcf, 'WindowState', 'maximized')
set(gcf, 'Position', [440 314 563 433]);
%yticks([-.6, -.3, 0, .3]);
%ylabel('rate change (Hz)');

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
figure; hold on;
scatter(all_xy_amplitudes_sort(all_durations_sort > 0), all_xy_peakvelocities_sort(all_durations_sort > 0), '.k', 'MarkerFaceAlpha', .05, 'MarkerEdgeAlpha', .05);

set(gca, 'xscale', 'log'); set(gca, 'yscale', 'log');
set(gca, 'FontSize', 26);
   H = gca;
    H.LineWidth = 2;
set(gca, 'FontName', 'Helvetica')
xlabel('amplitude (dva)', 'FontName', 'Helvetica'); ylabel('velocity (dva/s)', 'FontName', 'Helvetica');

%% plot each trial gaze timecourses for quality control
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
        velocity = xyVelocity_all{subj}{iBlock}{iTrial}(:, starttp:stoptp);
        velocity = smoothdata(velocity, 2, 'movmean', window_length_tps);
        velocityms = velocity; velocityms(:, ms_present_all{subj}{iBlock}{iTrial}(starttp:stoptp) ~= 1) = NaN;
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
        %             % plot xy velocity vs time
        %             %plot(time, velocity(iEye, :), eyecolors{iEye}, 'LineWidth', 1);
        %             %plot(time, velocityms(iEye, :), 'k', 'LineWidth', 3);
        %         end
        % average eyes
        gazeYavg = nanmean(gazeY, 1); gazeYmsavg = nanmean(gazeYms, 1);
        gazeXavg = nanmean(gazeX, 1); gazeXmsavg = nanmean(gazeXms, 1);
        velocityavg = nanmean(velocity, 1); velocitymsavg = nanmean(velocityms, 1);
        % plot x vs y
        %plot(gazeXavg, gazeYavg, eyecolors{1});
        %plot(gazeXmsavg, gazeYmsavg, 'k');
        % plot x vs time, y vs time
        time = 1:size(gazeXavg, 2);
        %plot(time, gazeXavg, dircolors{1}, time, gazeYavg, dircolors{2}, 'LineWidth', 1);
        %plot(time, gazeXmsavg, 'k', time, gazeYmsavg, 'k', 'LineWidth', 3);
        % plot xy velocity vs time
        plot(time, velocityavg, eyecolors{1}, 'LineWidth', 1);
        plot(time, velocitymsavg, 'k', 'LineWidth', 3);
    end
end
