
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
            % martinez-conde
            if isfield(BinocularSaccades_all_unsorted{subj}{block}{trial}, 'start')
                ntrialms = numel(BinocularSaccades_all_unsorted{subj}{block}{trial}.start);
                for ms = ms_present_indices
                    for nms = 1:ntrialms
                        curr_ms_indices = BinocularSaccades_all_unsorted{subj}{block}{trial}.start(nms):(BinocularSaccades_all_unsorted{subj}{block}{trial}.start(nms)+BinocularSaccades_all_unsorted{subj}{block}{trial}.duration(nms));
                        if ismember(ms, curr_ms_indices) % if current timepoint is within this microsaccade
                            ms_amplitudes_all_unsorted{subj}{block}{trial}(ms) = BinocularSaccades_all_unsorted{subj}{block}{trial}.XYAmplitude(nms); % assign tp with microsaccade's amplitude
                            ms_durations_all_unsorted{subj}{block}{trial}(ms) = BinocularSaccades_all_unsorted{subj}{block}{trial}.duration(nms); % assign tp with microsaccade's duration
                        end
                    end
                end
                % engbert
            elseif isfield(BinocularSaccades_all_unsorted{subj}{block}{trial}, 'msStarts')
                ntrialms = numel(BinocularSaccades_all_unsorted{subj}{block}{trial}.msStarts);
                for ms = ms_present_indices
                    for nms = 1:ntrialms
                        curr_ms_indices = BinocularSaccades_all_unsorted{subj}{block}{trial}.msStarts(nms):BinocularSaccades_all_unsorted{subj}{block}{trial}.msEnds(nms);
                        if ismember(ms, curr_ms_indices) % if current timepoint is within this microsaccade
                            ms_amplitudes_all_unsorted{subj}{block}{trial}(ms) = BinocularSaccades_all_unsorted{subj}{block}{trial}.msXYAmplitude(nms); % assign tp with microsaccade's amplitude
                        end
                    end
                end
                % edfImport edfExtractMicrosaccades
            elseif isfield(BinocularSaccades_all_unsorted{subj}{block}{trial}, 'Start')
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
% for subj = 1:nsubjects
%     for block = 1:nblocks
%         ntr = size(posdata_all{subj}{block}, 2);
%         for tr = 1:ntr
%             posdata_all{subj}{block}{tr}(2:5, :) = posdata_all{subj}{block}{tr}(2:5, :) * degperpix;
%         end
%     end
% end
%

% % make saccades = 0 and make blinks (missing data, nan) = 1? (current_data_blinks.mat)
% for subj = 1:nsubjects
%     for block = 1:nblocks
%         nblocktrials = size(ms_present_all{subj}{block}, 2);
%         for trial = 1:nblocktrials
%             ms_present_one = find(ms_present_all{subj}{block}{trial} == 1);
%             ms_present_nan = find(isnan(posdata_all{subj}{block}{trial}(2, :)));
%             ms_present_all{subj}{block}{trial}(ms_present_one) = 0;
%             ms_present_all{subj}{block}{trial}(ms_present_nan) = 1;
%         end
%     end
% end

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

%
% get total # excluded trials because of blinks or saccades
badblinktrials = 0;
badsaccadetrials = 0;
for subj = 1:nsubjects
    for sess = 1:3
        badblinktrials = badblinktrials + extra{subj, sess}.excludedblinktrials;
        badsaccadetrials = badsaccadetrials + extra{subj, sess}.excludedsaccadetrials;
    end
end
%
%
% save('current_data.mat', 'subjects', 'nsubjects', 'ntrials_all', ...
%     'ntrialstotal_byblock', 'ms_present_all', 'blink_present_all', 'pupil_trace_all', ...
%     'BinocularSaccades_all', 'posdata_all', 'trial_indices_all', 'ms_amplitudes_all', ...
%     'ms_directions_all', 'ms_durations_all', 'xyVelocity_all', 'trial_counter_all', ...
%     'degperpix', 'fps', 'extra', 'nblocks', '-v7.3');

%% just load the saved data
datadir = '/Users/maxlevinson/Documents/McGill/OneDrive - McGill University/neurospeed/UI_eyelink';
%datadir = '/export04/data/mlevin/UI_eyelink';
addpath(genpath(datadir));
load('current_data_withblinks.mat');

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

%% get hazard rates for each block
% I don't think this is necessary or useful at all anymore.

% hazard rate code from Qi Li
% we normalize every subject's RTs, combine across subjects for each
% block, generate hazard rate, then un-normalize back to subject values.
subj_meanrts = zeros(1, nsubjects);
for subj = 1:nsubjects
    subj_rts = [];
    for block = 1:nblocks
        subj_rts = [subj_rts; diff(trial_indices_all{subj}{block}(:, 2:3), 1, 2)];
    end
    subj_meanrts(subj) = nanmean(subj_rts);
end
norm_subjmeans = subj_meanrts ./ mean(subj_meanrts);

block_reactiontimes = cell(1, nblocks);
nbins = 21;
xout = nan(nbins, nblocks); xouts = nan(nbins, nsubjects, nblocks);
pdf = nan(nbins, nblocks); pdfs = nan(nbins, nsubjects, nblocks);
cdf = nan(nbins, nblocks); cdfs = nan(nbins, nsubjects, nblocks);
lambda = nan(nbins, nblocks); lambdas = nan(nbins-1, nsubjects, nblocks); % hazard rate
for block = 1:nblocks
    for subj = 1:nsubjects
        block_reactiontimes{block} = [block_reactiontimes{block}; diff(trial_indices_all{subj}{block}(:, 2:3), 1, 2) ./ norm_subjmeans(subj)];
    end
    [n, xout(:, block)] = hist(block_reactiontimes{block}, nbins);
        pdf(:, block) = n / sum(n);
        cdf(:, block) = cumsum(pdf(:, block));
        lambda(:, block) = pdf(:, block) ./ (1 - cdf(:, block)); % hazard rate
end
lambda = lambda(1:end-1, :);
%xout = xout(1:end-1, :);

% restore back to subject-specific value ranges
for block = 1:nblocks
for subj = 1:nsubjects
    lambdas(:, subj, block) = lambda(:, block);
    pdfs(:, subj, block) = pdf(:, block);
    cdfs(:, subj, block) = cdf(:, block);
    xouts(:, subj, block) = xout(:, block) .* norm_subjmeans(subj);
end
end



% define hazard rate per subject, per block
% reactiontimes = cell(nsubjects, nblocks);
% nbins = 11;
% xouts = nan(nbins, nsubjects, nblocks);
% pdfs = nan(nbins, nsubjects, nblocks);
% cdfs = nan(nbins, nsubjects, nblocks);
% lambdas = nan(nbins, nsubjects, nblocks); % hazard rate
% for subj = 1:nsubjects
%     for block = 1:nblocks
%         reactiontimes{subj, block} = diff(trial_indices_all{subj}{block}(:, 2:3), 1, 2);
%         [n, xouts(:, subj, block)] = hist(reactiontimes{subj, block}, nbins);
%         pdfs(:, subj, block) = n / sum(n);
%         cdfs(:, subj, block) = cumsum(pdfs(:, subj, block));
%         lambdas(:, subj, block) = pdfs(:, subj, block) ./ (1 - cdfs(:, subj, block)); % hazard rate
%     end
% end
% lambdas = lambdas(1:end-1, :, :); % remove last one (inf)
% xouts = xouts(1:end-1, :, :);

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

%% epoch rates by colors / sizes (collapsing across the other dimension)
% the point is to maximize signal - include more trials in rate estimate.

%weighsubjects = 1; % weigh subjects by usable trial counts? yes or no

ntps_before = 5000; ntps_after = 1000; % around button press
%ntps_before = 1000; ntps_after = 2000;%20000; % from stim onset
ntimepoints = ntps_before + ntps_after + 1;
reference_point = 'button_press'; % button_press or full_trial (around stim onset) or onset_to_bp (around stim onset but stop at bp)

rate_window_length = 1; % for sliding window ms/sec calculation
stepsize = 1;
%ms_counts_all = nan((ntimepoints-rate_window_length-1)/stepsize + 1, nsubjects, 3);
%count_timepoints = rate_window_length/2 + (0:stepsize:(ntimepoints-1-rate_window_length)) - ntps_before; % time stamps for each time window in ms counts / second
ms_counts_all = nan(ntimepoints, nsubjects, 3);
baseline_ms_counts = nan(nsubjects, 3);
count_timepoints = (1:ntimepoints) - ntps_before;
block_indices = {[1:3], [4:6], [7:9]}; % blocks per color (1, 2, 3)
%block_indices = {[3 6 9], [2 5 8], [1 4 7],}; % blocks per size (6, 4, 2)
for clr = 1:3
    for subj = 1:nsubjects
        %threeblocks_ms_present = cell(3, 1);
        clr_ms_present = [];
        baseline_clr_ms_present = [];
        for block = block_indices{clr}
            % get indices
            nblocktrials = size(ms_present_all{subj}{block}, 2);
            block_ms_present = nan(ntimepoints, nblocktrials);
            baseline_block_ms_present = nan(ntimepointsbaseline, nblocktrials);
            if nblocktrials > 0 % if not missing whole block
                for trial = 1:nblocktrials
                    rt = diff(trial_indices_all{subj}{block}(trial, 2:3));
                    if rt > 0 && rt < 20000 % only use trials with certain RT?
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
                                nan_pad_start = ntps_before - (button_press - stim_on); % if RT is faster than ntps_before
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
                                nan_pad_end = ntps_after - (button_press - stim_on); % end after button press
                                if nan_pad_start <= 0
                                    first_idx = stim_on - ntps_before + 1;
                                    nan_pad_start = 0;
                                elseif nan_pad_start > 0
                                    first_idx = 1;
                                end
                                if nan_pad_end < 0 % if it goes over 20s for some reason
                                    last_idx = stim_on + ntps_after;
                                else
                                    last_idx = button_press;
                                end
                        end
                        block_ms_present(:, trial) = [nan(1, nan_pad_start), ms_present_all{subj}{block}{trial}(first_idx:last_idx), nan(1, nan_pad_end)];
                        %extranan = buffer_baseline_idx(end) - numel(ms_present_all{subj}{block}{trial}); curr_ms_present_baseline = [ms_present_all{subj}{block}{trial}, nan(1, extranan)]; curr_ms_present_baseline((button_press-2000):end) = nan;
                        %baseline_block_ms_present(:, trial) = curr_ms_present_baseline(buffer_baseline_idx);
                    end
                end
            end
            clr_ms_present = [clr_ms_present, block_ms_present];
            %baseline_clr_ms_present = [baseline_clr_ms_present, baseline_block_ms_present];
        end
        ms_counts_all(:, subj, clr) = cross_trial_ms_rate(clr_ms_present, 'alpha');
        %curr_baseline_r = cross_trial_ms_rate(baseline_clr_ms_present, 'alpha');
        %curr_baseline_r = curr_baseline_r(baseline_sub_idx);
        %baseline_ms_counts(subj, clr) = nanmean(curr_baseline_r);
        baseline_ms_counts(subj, clr) = nanmean(ms_counts_all(1:2000, subj, clr)); % first 2000 of epoch = baseline.
    end
end

% normalize everything by subject baseline means
norm_ms_counts_all = nan(size(ms_counts_all));
norm_baseline_ms_counts = nan(size(baseline_ms_counts));
for subj = 1:nsubjects
    %norm_ms_counts_all(:, subj, :) = ms_counts_all(:, subj, :);
    norm_ms_counts_all(:, subj, :) = ms_counts_all(:, subj, :) ./ nanmean(baseline_ms_counts(subj, :));
    %norm_ms_counts_all(:, subj, :) = (ms_counts_all(:, subj, :) - nanmean(baseline_ms_counts(subj, :))) ./ nanmean(baseline_ms_counts(subj, :)) .* 100;
    norm_baseline_ms_counts(subj, :) = baseline_ms_counts(subj, :) ./ nanmean(baseline_ms_counts(subj, :));
end

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
                % engbert
                %                 all_xy_amplitudes = [all_xy_amplitudes, BinocularSaccades_all{subj}{block}{trial}.msXYAmplitude];
                %                 all_xy_peakvelocities = [all_xy_peakvelocities, BinocularSaccades_all{subj}{block}{trial}.msMaxXYVel];
                %                 all_durations = [all_durations, BinocularSaccades_all{subj}{block}{trial}.msLengths];
                % martinez-conde
                if isfield(BinocularSaccades_all{subj}{block}{trial}, 'XYAmplitude')
                all_xy_amplitudes = [all_xy_amplitudes; BinocularSaccades_all{subj}{block}{trial}.XYAmplitude];
                all_xy_peakvelocities = [all_xy_peakvelocities; BinocularSaccades_all{subj}{block}{trial}.peakXYVelocity];
                all_durations = [all_durations; BinocularSaccades_all{subj}{block}{trial}.duration];
                all_subjnums = [all_subjnums; repmat(subj, numel(BinocularSaccades_all{subj}{block}{trial}.duration), 1)];
            % edfImport
                elseif isfield(BinocularSaccades_all{subj}{block}{trial}, 'Amplitude')
                all_xy_amplitudes = [all_xy_amplitudes; BinocularSaccades_all{subj}{block}{trial}.Amplitude];
                all_xy_peakvelocities = [all_xy_peakvelocities; BinocularSaccades_all{subj}{block}{trial}.vPeak];
                all_durations = [all_durations; BinocularSaccades_all{subj}{block}{trial}.Duration];
                all_subjnums = [all_subjnums; repmat(subj, numel(BinocularSaccades_all{subj}{block}{trial}.Duration'), 1)];
                end
            end
        end
    end
end

% % convert to dva if needed
% all_xy_amplitudes = all_xy_amplitudes * degperpix;
% all_xy_peakvelocities = all_xy_peakvelocities * degperpix;

% remove too big or too small
all_xy_amplitudes_sort = all_xy_amplitudes(all_xy_amplitudes < 2 & all_xy_peakvelocities < 300);
all_xy_peakvelocities_sort = all_xy_peakvelocities(all_xy_amplitudes < 2 & all_xy_peakvelocities < 300);
all_subjnums_sort = all_subjnums(all_xy_amplitudes < 2 & all_xy_peakvelocities < 300);
all_durations_sort = all_durations(all_xy_amplitudes < 2 & all_xy_peakvelocities < 300);
figure; hold on;
scatter(all_xy_amplitudes_sort(all_durations_sort > 0), all_xy_peakvelocities_sort(all_durations_sort > 0), '.k', 'MarkerFaceAlpha', .05, 'MarkerEdgeAlpha', .05);
% clrs = colormap;
% for subj = 1:nsubjects
%     scatter(all_xy_amplitudes_sort(all_subjnums_sort == subj), all_xy_peakvelocities_sort(all_subjnums_sort == subj), '.', 'MarkerFaceColor', clrs(subj * 3, :));
% end
%ylim([3 40]); xlim([0.1 2]); % martinez-conde scales
set(gca, 'xscale', 'log'); set(gca, 'yscale', 'log');
set(gca, 'FontSize', 26);
   H = gca;
    H.LineWidth = 2;
set(gca, 'FontName', 'Helvetica')
%H.XAxis.TickLength = [0 0];
xlabel('amplitude (dva)', 'FontName', 'Helvetica'); ylabel('velocity (dva/s)', 'FontName', 'Helvetica');

%% split by direction
% 2 = left inwards, 1 = left outwards, -1 = right outwards, -2 = right inwards
backup_ms_present_all = ms_present_all;
% for analysis: then set ms_present_all equal to one of the below
out_ms_present_all = ms_present_all;
in_ms_present_all = ms_present_all;
left_ms_present_all = ms_present_all;
right_ms_present_all = ms_present_all;

for subj = 1:nsubjects
    for block = 1:nblocks
        nblocktrials = size(ms_present_all{subj}{block}, 2);
        for trial = 1:nblocktrials
            out_ms_present_all{subj}{block}{trial}(abs(ms_directions_all{subj}{block}{trial}) == 2) = 0;
            in_ms_present_all{subj}{block}{trial}(abs(ms_directions_all{subj}{block}{trial}) == 1) = 0;
            left_ms_present_all{subj}{block}{trial}(ms_directions_all{subj}{block}{trial} < 0) = 0;
            right_ms_present_all{subj}{block}{trial}(ms_directions_all{subj}{block}{trial} > 0) = 0;
        end
    end
end

%% comparing baseline to -500ms period before button press
mean_baseline_rate = nanmean(baseline_ms_counts, 1);
stderr_baseline_rate = nanstd(baseline_ms_counts, 0, 1) ./ sqrt(nsubjects);
mean_bp_rate = nanmean(min_rates, 1);
stderr_bp_rate = nanstd(min_rates, 0, 1) ./ sqrt(nsubjects);

perc_change_prob = (before_press_probability - baseline_probability) ./ baseline_probability .* 100;
perc_change_avg_prob = (mean_bp_prob - mean_baseline_prob) ./ mean_baseline_prob .* 100;
figure;
barwitherr([stderr_baseline_rate', stderr_bp_rate'], [mean_baseline_rate', mean_bp_rate']);
figure;
bar(perc_change_avg_prob);

%% look at minimum rate per condition
min_rates = nan(nsubjects, nblocks); % take minimum/min mean of timeseries per subject, per block.
post_rates = nan(nsubjects, nblocks); % take post-button press window
min_amps = nan(nsubjects, nblocks);
min_slip = nan(nsubjects, nblocks);
decr_start = nan(nsubjects, nblocks); % 700ms time window surrounding the "start" of the microsaccade decrease
window_use = (ntps_before-800):(ntps_before-300); % -800ms to -300ms
post_window_use = (ntps_before+150):(ntps_before+550); % +150ms to +550ms
start_window_use = (ntps_before-2300):(ntps_before-1800);
Nb_window_use = (ntps_before/50-800/50):(ntps_before/50-300/50);
Nb_post_window_use = (ntps_before/50+150/50):(ntps_before/50+550/50);
for ss = 1:nsubjects
    for bb = 1:nblocks
        %min_rates(ss, bb) = nanmin(ms_counts_all(window_use, ss, bb)); % raw MIN
        min_rates(ss, bb) = nanmean(ms_counts_all(window_use, ss, bb)); % raw.
        %min_rates(ss, bb) = nanmean(ms_counts_all(window_use+all_subj_shifts(ss, bb), ss, bb)); % raw with subj shift
        %min_rates(ss, bb) = nanmean(norm_ms_counts_all(window_use, ss, bb)); % normalized to subject mean.
        %min_rates(ss, bb) = nanmean(norm_ms_counts_all(window_use+subj_shifts(ss), ss, bb)); % normalized to subject mean with subj shift
        %min_rates(ss, bb) = (nanmean(ms_counts_all(window_use, ss, bb)) - baseline_ms_counts(ss, bb)) ./ baseline_ms_counts(ss, bb);% .* 100; % percent change from block baseline mean.
        %min_rates(ss, bb) = (nanmean(norm_ms_counts_all(window_use, ss, bb)) - norm_baseline_ms_counts(ss, bb)) ./ norm_baseline_ms_counts(ss, bb);% .* 100; % normalized percent change from norm block baseline mean.
        %min_rates(ss, bb) = (nanmean(ms_counts_all(window_use, ss, bb)) - nanmean(baseline_ms_counts(ss, :))) ./ nanmean(baseline_ms_counts(ss, :));% .* 100; % percent change from subject baseline mean.
        %min_rates(ss, bb) = nanmean(norm_ms_counts_all(window_use, ss, bb)) - norm_baseline_ms_counts(ss, bb); % norm minus normalized block baseline.
        %min_rates(ss, bb) = nanmean(norm_ms_counts_all(window_use, ss, bb)) ./ norm_baseline_ms_counts(ss, bb); % norm divided by normalized block baseline.
        %min_rates(ss, bb) = nanmean(norm_ms_counts_all(window_use+subj_shifts(ss), ss, bb)) - norm_baseline_ms_counts(ss, bb); % norm minus normalized block baseline with subj shift
        post_rates(ss, bb) = nanmean(ms_counts_all(post_window_use, ss, bb)); % raw.
        %post_rates(ss, bb) = nanmean(norm_ms_counts_all(post_window_use, ss, bb)); % normalized to subject mean
        min_amps(ss, bb) = nanmean(mean_amplitudes_all(window_use, ss, bb)); % raw
        %min_amps(ss, bb) = nanmean(mean_amplitudes_all(window_use+all_subj_shifts(ss, bb), ss, bb)); % raw with subj shift
        post_amps(ss, bb) = nanmean(mean_amplitudes_all(post_window_use, ss, bb)); % raw
        min_slip(ss, bb) = nanmean(mean_Nb_all(Nb_window_use, ss, bb)); % raw
        post_slip(ss, bb) = nanmean(mean_Nb_all(Nb_post_window_use, ss, bb));
        %min_rates(ss, bb) = nanmean(ms_counts_all(window_use, ss, bb)) - baseline_ms_counts(ss, bb); % raw difference
        %min_rates(ss, bb) = nanmean(ms_counts_all(window_use, ss, bb)) ./ baseline_ms_counts(ss, bb); % raw divide by block baseline (gives some inf)
        decr_start(ss, bb) = nanmean(ms_counts_all(start_window_use+subj_shifts(ss), ss, bb)); % raw
        decr_start(ss, bb) = nanmean(norm_ms_counts_all(start_window_use+subj_shifts(ss), ss, bb)); % subj normalized
        %min_rates(ss, bb) = min_rates(ss, bb) - decr_start(ss, bb); % subtract from decrease start
    end
end
mean_min_rates = nanmean(min_rates, 1);
stderr_min_rates = nanstd(min_rates, 0, 1) ./ sqrt(nsubjects);

% group by color
color_min_rates = [mean(min_rates(:, 1:3), 2), mean(min_rates(:, 4:6), 2), mean(min_rates(:, 7:9), 2)];
size_min_rates = [mean(min_rates(:, [3 6 9]), 2), mean(min_rates(:, [2 5 8]), 2), mean(min_rates(:, [1 4 7]), 2)];
color_baseline_rates = [mean(baseline_ms_counts(:, 1:3), 2), mean(baseline_ms_counts(:, 4:6), 2), mean(baseline_ms_counts(:, 7:9), 2)];
size_baseline_rates = [mean(baseline_ms_counts(:, [3 6 9]), 2), mean(baseline_ms_counts(:, [2 5 8]), 2), mean(baseline_ms_counts(:, [1 4 7]), 2)];

% check prevalence of an effect
sum(color_min_rates(:, 1) < color_min_rates(:, 2) & color_min_rates(:, 2) < color_min_rates(:, 3))
sum(size_min_rates(:, 1) < size_min_rates(:, 2) & size_min_rates(:, 2) < size_min_rates(:, 3))
sum(size_min_rates(:, 1) > size_min_rates(:, 3))
sum(color_baseline_rates(:, 1) < color_baseline_rates(:, 2) & color_baseline_rates(:, 2) < color_baseline_rates(:, 3))
sum(color_baseline_rates(:, 1) < color_baseline_rates(:, 3))



% normalize min_rates by subject mean of min_rates?
%min_rates(ss, :) = min_rates(ss, :) ./ nanmean(min_rates(ss, :)); % normalize min_rates by subject mean of min_rates.


% % averaged across subjects first
% mean_min_rates = nan(1, nblocks);
% stderr_min_rates = nan(1, nblocks);
% for bb = 1:nblocks
%     %mean_min_rates(bb) = nanmin(nanmean(ms_counts_all(window_use, :, bb), 2));
%     mean_min_rates(bb) = nanmean(nanmean(norm_ms_counts_all(window_use, :, bb), 2));
%     stderr_min_rates(bb) = nanmean(nanstd(norm_ms_counts_all(window_use, :, bb), 0, 2)) ./ sqrt(nsubjects);
% end
%
% % then avged by color (1, 2, 3) and size (6, 4, 2)
% color_mean_min_rates = [mean(mean_min_rates(1:3)), mean(mean_min_rates(4:6)), mean(mean_min_rates(7:9))]
% size_mean_min_rates = [mean(mean_min_rates([3 6 9])), mean(mean_min_rates([2 5 8])), mean(mean_min_rates([1 4 7]))]
%
% % scatter plot min rate by mean rt (per block/subject)
% %first run analyze_group_data.m
% %permuted_rts = squeeze(permute(all_mean_rts(1, :, :), [3 2 1]));
% permuted_rts = squeeze(permute(norm_all_mean_rts(1, :, :), [3 2 1]));
% %plotrates = min_rates(:, [1:9]);
% %plotrates = norm_baseline_ms_counts(:);
% %plotrates = nanmean(norm_baseline_ms_counts, 2);
% plotrates = post_rates(:);
% %plotrts = permuted_rts(:, [1:9]);
% %plotrts = nanmean(permuted_rts, 2);
% %plotrts = permuted_rts(:);
% %plotrts = norm_rts_all(:);
% plotrts = post_slip(:);
% x = plotrates(~isnan(plotrates) & ~isnan(plotrts));
% y = plotrts(~isnan(plotrates) & ~isnan(plotrts));
% X = [ones(length(x), 1) x];
% b = X \ y;
% ypred = X*b;
%
% figure;
% scatter(x,y);
% hold on;
% plot(x, ypred);
% xlabel('minimum ms rate'); ylabel('mean rt');
%
%
% % average within conditions then scatter
% rates = baseline_ms_counts;%min_rates;
% rts = permuted_rts;
%
% colormeanrates = [nanmean(rates(:, 1:3), 2), nanmean(rates(:, 4:6), 2), nanmean(rates(:, 7:9), 2)];
% sizemeanrates = [nanmean(rates(:, [3 6 9]), 2), nanmean(rates(:, [2 5 8]), 2), nanmean(rates(:, [1 4 7]), 2)];
%
% colormeanrts = [nanmean(rts(:, 1:3), 2), nanmean(rts(:, 4:6), 2), nanmean(rts(:, 7:9), 2)];
% sizemeanrts = [nanmean(rts(:, [3 6 9]), 2), nanmean(rts(:, [2 5 8]), 2), nanmean(rts(:, [1 4 7]), 2)];
%
% x = colormeanrates(:); y = colormeanrts(:);
% figure; scatter(colormeanrates(:), colormeanrts(:));
%
% x = sizemeanrates(:); y = sizemeanrts(:);
% figure; scatter(sizemeanrates(:), sizemeanrts(:));
%
% X = [ones(length(x), 1) x];
% b = X \ y;
% ypred = X*b;
%
% figure;
% scatter(x,y);
% hold on;
% plot(x, ypred);
% xlabel('minimum ms rate'); ylabel('mean rt');
%
%
% % anova on min rates
% %snums = repmat((1:18)', 9, 1);
% %y = min_rates(:);
% %clrs = [ones(18*3, 1); 2*ones(18*3, 1); 3*ones(18*3, 1)];
% %szs = repmat([ones(18, 1); 2*ones(18, 1); 3*ones(18, 1)], 3, 1);
% within = table(categorical([1 1 1 2 2 2 3 3 3])', categorical([3 2 1 3 2 1 3 2 1]'), 'VariableNames', {'color', 'size'});
% blocknames = {'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9'};
% trates = array2table(min_rates, 'VariableNames', blocknames);
% rm = fitrm(trates,'b1-b9~1', 'WithinDesign', within);
% ranovatbl = ranova(rm, 'WithinModel', 'color*size')
%
% % post hoc
%
%
% % anova on baseline rates
% norm_baseline_ms_counts = nan(size(baseline_ms_counts));
% for subj = 1:nsubjects
%     norm_baseline_ms_counts(subj, :) = baseline_ms_counts(subj, :) ./ nanmean(baseline_ms_counts(subj, :));
% end
% trates = array2table(baseline_ms_counts, 'VariableNames', blocknames);
% rm = fitrm(trates,'b1-b9~1', 'WithinDesign', within);
% ranovatbl = ranova(rm, 'WithinModel', 'color*size')


% store all min rate data in .csv table.
snums = repmat((1:nsubjects)', 9, 1);
y = min_rates(:);
basel = baseline_ms_counts(:);
baseamp = baseline_amplitudes(:);
baseslip = baseline_Nb(:);
baseblink = baseline_blink_rate(:);
nonan_idx = find(~isnan(y) & ~isinf(y) & ~isnan(basel));
rtimes = mean_rts_all(:) ./ 1000;
%nonan_idx = find(~isnan(rtimes)); % if analyzing only trials w/o microsaccades
dstarts = decr_start(:);
clrs = [ones(nsubjects*3, 1); 2*ones(nsubjects*3, 1); 3*ones(nsubjects*3, 1)];
szs = repmat([2*ones(nsubjects, 1); 4*ones(nsubjects, 1); 6*ones(nsubjects, 1)], 3, 1);
% scale szs by cortical magnification
A = 0.063; % from Engel et al
l0 = 36.54 * A;
%M_szs = 1 ./ (A .* szs); % dl
M_szs = (log(szs) - l0) ./ A; % l
% extra measures
rate_post = post_rates(:);
ampl_min = min_amps(:);
ampl_post = post_amps(:);
slip_min = min_slip(:);
slip_post = post_slip(:);

tabledata = [snums(nonan_idx), y(nonan_idx), basel(nonan_idx), baseblink(nonan_idx), dstarts(nonan_idx), rtimes(nonan_idx), clrs(nonan_idx), szs(nonan_idx), M_szs(nonan_idx), rate_post(nonan_idx), ampl_min(nonan_idx), ampl_post(nonan_idx), slip_min(nonan_idx), slip_post(nonan_idx), baseamp(nonan_idx), baseslip(nonan_idx)];
T = array2table(tabledata, 'VariableNames', {'subject', 'msratemin', 'baselinerate', 'baselineblinks', 'msratestart', 'RT', 'contrast', 'eccentricity', 'scaled_eccentricity', 'msratepost', 'amplmin', 'amplpost', 'slipmin', 'slippost', 'baselineampl', 'baselineslip'});
writetable(T, 'bhv_data/minrates.csv');
%writetable(T, 'bhv_data/noms_minrates.csv');
%writetable(T, 'bhv_data/onlyms_subset_minrates.csv');

%% regress out stim-evoked pupil response.
% first compute & plot full trial epoch
% then extract the grand avg timeseries up to 2000s post-stim
% should capture the whole transient.
% (start at 51 ms to avoid artifact at very beginning)
first2000avgpupil = [zeros(50, 1); mean_pupil_area(51:3000); zeros(18001, 1)];

% then compute button_press epoch

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

%% plot transition probability following each microsaccade
smooth_transition_probabilities_all = smoothdata(transition_probabilities_all, 'movmean', 100);
%smooth_transition_probabilities_all = smoothdata(baseline_transition_probabilities, 'movmean', 100);
size_indices = [3 2 1 3 2 1 3 2 1]; % 1=size 6, 2=size 4, 3=size 2
color_indices = [1 1 1 2 2 2 3 3 3]; % 1=low contrast, 2=medium, 3=high contrast
all_indices = [1 1 1 1 1 1 1 1 1];

which2use = 'size'; % size or color or all (average across all blocks)
subjects2use = 1:nsubjects;%[1:10, 13:nsubjects];%[1:nsubjects]; % normally all 
overlay = 1; % overlay the plots on one axis or use subplots
plotcolors = {[1 0 0], [0 1 0], [0 0 1]};
nobounds = 0; % just plot line, no stderr bounds?
baseline_correct = 1;

figure; hold on;
%for subj = 1:nsubjects
%    subjects2use = subj;
%figure; hold on;
for i = 1:3
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
    end
    
    if baseline_correct
        transition_probability_uncorr = smooth_transition_probabilities_all(:, subjects2use, curr_block_idx);
        curr_probability_baseline = squeeze(nanmean(baseline_transition_probabilities(trans_ntps_before:trans_ntps_before+2000, subjects2use, curr_block_idx), 1));
                transition_probability_corr = nan(size(transition_probability_uncorr));
                for bb = 1:size(transition_probability_uncorr, 3) % per block
                    for ss = 1:size(transition_probability_uncorr, 2) % per subject
                        curr_uncorr = transition_probability_uncorr(:, ss, bb);
                        % baseline correct
                        curr_base = curr_probability_baseline(ss, bb); % from computed baseline
                        %curr_base = nanmean(curr_uncorr(1:(1000/stepsize))); % from this actual epoch
                        transition_probability_corr(:, ss, bb) = (curr_uncorr - curr_base) ./ curr_base .* 100; % percent change
                        %transition_probability_corr(:, ss, bb) = curr_uncorr - curr_base; % just subtract baseline
                        %transition_probability_corr(:, ss, bb) = curr_uncorr ./ curr_base; % normalize by baseline (base=1)
                        % z score
                        %transition_probability_corr(:, ss, bb) = (curr_uncorr - nanmean(curr_uncorr)) ./ nanstd(curr_uncorr);
                    end
                end
                transition_probability_corr(isinf(transition_probability_corr)) = nan;
                transition_probability_curr = nanmean(transition_probability_corr, 3);
    else
        transition_probability_curr = nanmean(smooth_transition_probabilities_all(:, subjects2use, curr_block_idx), 3);
    end
    
            mean_transition_probability = nanmean(transition_probability_curr, 2); % average by subject
            stderr_transition_probability = nanstd(transition_probability_curr, 0, 2) ./ sqrt(nsubjects);
            %plot((-ntps_before:ntps_after)' * 1000/fps, smoothed_prob);
            bp = boundedline([-trans_ntps_before:trans_ntps_after], mean_transition_probability, stderr_transition_probability, 'cmap', plotcolors{i}, 'alpha', 'transparency', 0.05);

            set(bp, 'LineWidth', 4);
            
            
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
            ylabel('probability');
        elseif baseline_correct
            ylabel(['probability change from baseline']);
        end
        xlabel('time (s)');
        if strcmp(which2use, 'all')
            set(bp, 'Color', [0 0 0]);
            %break % stop for loop
        end
    end
    if strcmp(which2use, 'all')
        break
    end
end
switch which2use
    case 'size'
        if overlay
            legend({'large', 'medium', 'small'});
            sl = title(['by boundary eccentricity']);
        else
            sl = suplabel('center size', 'y');
        end
        set(sl, 'FontSize', 26);
    case 'color'
        if overlay
            legend({'low', 'medium', 'high'});
            sl = title(['by color contrast']);
        else
            sl = suplabel('color contrast', 'y');
        end
        set(sl, 'FontSize', 26);
end
set(gca, 'FontSize', 26);

% set subject color if plotting just 1 subject at a time
if numel(subjects2use) == 1
    clrs = colormap;
    set(bp, 'Color', clrs(subjects2use * 3, :));
    set(bp, 'LineWidth', 1); % or make line smaller
    %set(bp, 'Color', [0.5 0.5 0.5]);
end

%end

legend boxoff

%set(gcf, 'WindowState', 'maximized')
set(gcf, 'Position', [440 314 563 433]);

%% plot immobilization duration scatter, colored by size or contrast
plotcolors = {[1 0 0], [0 1 0], [0 0 1]};
% by size
lowidx = eccentricities == 6;
medidx = eccentricities == 4;
highidx = eccentricities == 2;
% by color
lowidx = contrasts == 1;
medidx = contrasts == 2;
highidx = contrasts == 3;

figure; hold on;
scatter(ones(size(immobilization_durations(lowidx))), immobilization_durations(lowidx), 'MarkerFaceColor', plotcolors{1}, 'MarkerFaceAlpha', 0.05, 'MarkerEdgeColor', plotcolors{1}, 'MarkerEdgeAlpha', 0.05, 'jitter', 'on');
scatter(ones(size(immobilization_durations(medidx))), immobilization_durations(medidx), 'MarkerFaceColor', plotcolors{2}, 'MarkerFaceAlpha', 0.05, 'MarkerEdgeColor', plotcolors{2}, 'MarkerEdgeAlpha', 0.05, 'jitter', 'on')
scatter(ones(size(immobilization_durations(highidx))), immobilization_durations(highidx), 'MarkerFaceColor', plotcolors{3}, 'MarkerFaceAlpha', 0.05, 'MarkerEdgeColor', plotcolors{3}, 'MarkerEdgeAlpha', 0.05, 'jitter', 'on');
