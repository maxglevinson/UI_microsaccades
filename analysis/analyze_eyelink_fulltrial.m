% get eyelink data, but save entire trial and indices for trial events
% (ONLY SAVE FROM FIX_ONSET to STIM_OFFSET)
function [ms_present, fps, pupil_trace, all_BinocularSaccades, trial_indices, all_xyVelocity, posdata, block_trial_counter, subj_extra] = analyze_eyelink_fulltrial(subj_id, session_id)

% load data
nblocks = 3; % blocks per session
ntrialsperblock = 50;
subj_id = num2str(subj_id);
session_id = num2str(session_id);
trials_file = ['bhv_data/trials_sub' subj_id 's' session_id '.mat'];

% if processing first from raw
%edffile = ['bhv_data/sub' subj_id 's' session_id '.edf'];
%Trials = edfImport(edffile, [1 1 1]); % last 1 parameter is to include all samples
%save(trials_file, 'Trials'); % save Trials struct for future use

% note for subject 210, I manually combined sessions 3-4 and changed the block
% numbers of session 4 from [1,2] to [2,3]. Original edfs are unchanged;
% just the new trials_file is changed. Same general thing done for subj 223.

% if already processed before from raw, load processed version
load(trials_file); % load Trials struct instead of full edf

fps = Trials(1).Header.rec.sample_rate;

% get indices of good trials (not aborted due to bad fixation or blink)
good_trials = ones(1, numel(Trials));
for iTrial = 1:numel(Trials)
    if sum(contains(Trials(iTrial).Events.message, 'aborted'))
        good_trials(iTrial) = 0;
    end
end
% index them!
good_trials = find(good_trials);
nblocks_actual = numel(good_trials)/ntrialsperblock;
trial_counter = nan(1, numel(Trials));
trial_counter(good_trials) = repmat(1:ntrialsperblock, 1, nblocks_actual);

% also select only main illusion trials, not control / catch trials
bhv_data = load(['bhv_data/color_' subj_id '_allsessions_data.mat']);
block_idx = (str2num(session_id)-1) * 3 + (1:3);
%bhv_data = load(['bhv_data/color_' subj_id '_session_4_data.mat']); block_idx=1:50; % to look at subj 201 block4 control
session_trials = find(ismember(bhv_data.data(:, 3), block_idx));
illusion_trials = bhv_data.data(session_trials, 2) == 1; %1 = illusion, 2 = replay, 3 = catch
% also exclude trials with median fixation displacement > 1 dva
illusion_trials = illusion_trials & bhv_data.data(session_trials, 8) < 1;
illusion_trials = illusion_trials(1:numel(good_trials)); % in case there are more data entries than actual recorded trials (if experiment ended earlier than planned)
good_trials = good_trials(illusion_trials);

% extract interesting events: fixations, saccades, blinks, button presses
Trials = edfExtractInterestingEvents(Trials, 'STIM_ONSET'); % starts trial at STIM_ONSET message; default is "!MODE RECORD" message

% extract trial variables: subject, trial num, block, luminance, radius
Trials = edfExtractVariables(Trials);

% extract key events
Trials = edfExtractKeyEventsTiming(Trials);

% get indices of good trials with button press
button_trials = zeros(1, numel(Trials));
for iTrial = 1:numel(Trials)
    if ismember(iTrial, good_trials)
        if isfield(Trials(iTrial).KeyEvents, 'BUTTON_PRESS')
            % BUT ONLY TAKE TRIALS WITH RT > 2.
            rt = Trials(iTrial).KeyEvents.BUTTON_PRESS - Trials(iTrial).KeyEvents.STIM_ONSET;
            if rt > 2000
                button_trials(iTrial) = 1;
            end
        end
    end
end

% get rid of trials in which blink or large saccade occurs within some cutoff time before button
% press (e.g. 2 seconds: within 1.5 seconds before subjective filling-in, presumably)
exclude_blinks = 1; % default yes, but set to 0 if we want to analyze even trials with blinks in this time window.

subj_extra.excludedsaccadetrials = 0;
subj_extra.excludedblinktrials = 0;
if exclude_blinks
%blink_cutoff = -2000; % for cleaner ms rate curve visual - avoid blinks within 2 seconds before button press
blink_cutoff = -300; % just during motor reaction time. Better for linear modeling.
for iTrial = 1:numel(Trials)
    if button_trials(iTrial)
        button_idx = Trials(iTrial).KeyEvents.BUTTON_PRESS - Trials(iTrial).KeyEvents.STIM_ONSET; % indexing starts at stim onset
        blink_indices = get_event_indices(Trials(iTrial), 150, 150, 'blink');
        saccade_indices = get_event_indices(Trials(iTrial), 150, 150, 'saccade');
        exclude_indices = button_idx + (blink_cutoff:0); % cuts off at button (t=0), can extend later too
        if sum(ismember(exclude_indices, blink_indices))
            button_trials(iTrial) = 0; % remove trial from analysis.
            subj_extra.excludedblinktrials = subj_extra.excludedblinktrials+1;
        elseif sum(ismember(exclude_indices, saccade_indices))
            button_trials(iTrial) = 0;
            subj_extra.excludedsaccadetrials = subj_extra.excludedsaccadetrials+1;
        end
    end
end
clear blink_indices
end

trial_counter(~button_trials) = nan;

% sort trials by block
block_button_trials = cell(1, nblocks);
block_trial_counter = cell(1, nblocks);
for bl = 1:nblocks
    curr_button_trials = button_trials;
    for iTrial = 1:numel(Trials)
        if button_trials(iTrial) && Trials(iTrial).Variables.block ~= bl
            curr_button_trials(iTrial) = 0;
        end
    end
    block_button_trials{bl} = find(curr_button_trials);
    block_trial_counter{bl} = trial_counter(logical(curr_button_trials));
end


% prepare for microsaccade algorithm

% extract microsaccades
Trials(1).Microsaccades = [];
if sum(strcmp(subj_id, {'219', '223', '224'}))
    eta = 6; % higher threshold
else
    eta = 5; % default
end
Trials(logical(button_trials)) = UI_extractMicrosaccades(Trials(logical(button_trials)), eta);

% initialize data structures
ms_present = cell(1, nblocks);
all_BinocularSaccades = cell(1, nblocks);
trial_indices = cell(1, nblocks);
posdata = cell(1, nblocks);
for bl = 1:nblocks
    ms_present{bl} = cell(1, numel(block_button_trials{bl}));
    pupil_trace{bl} = cell(1, numel(block_button_trials{bl}));
    all_BinocularSaccades{bl} = cell(1, numel(block_button_trials{bl}));
    trial_indices{bl} = nan(numel(block_button_trials{bl}), 4);
    posdata{bl} = cell(1, numel(block_button_trials{bl}));
end

for bl = 1:nblocks
    button_trials_use = block_button_trials{bl};
    for iTrial = 1:numel(button_trials_use)
        trial_idx = button_trials_use(iTrial);
        trial_start = Trials(trial_idx).KeyEvents.FIX_ONSET; % tracker time, ms
        stim_start = Trials(trial_idx).KeyEvents.STIM_ONSET; % tracker time, ms
        button_press = Trials(trial_idx).KeyEvents.BUTTON_PRESS; % tracker time, ms
        
        % replace button press for shift onset, if we want to look at data
        % surrounding that timepoint
%                 if isfield(Trials(trial_idx).KeyEvents, 'SHIFT_ONSET')
%                     button_press = Trials(trial_idx).KeyEvents.SHIFT_ONSET;
%                 end
        
        trial_end = Trials(trial_idx).KeyEvents.STIM_OFFSET; % tracker time, ms
        full_time = Trials(trial_idx).Samples.time; % tracker time
        trial_start_idx = find(full_time == trial_start); stim_start_idx = find(full_time == stim_start); trial_end_idx = find(full_time == trial_end); button_press_idx = find(full_time == button_press);
        trial_time_idx = full_time >= trial_start & full_time <= trial_end; % extract full trial indices
        trial_data = zeros(7, sum(trial_time_idx)); % 1st row = time, 2nd-3rd = x pos (L, R), 4rd-5th = y pos (L, R), 6th-7th = pupil area (L, R)
        trial_data(1, :) = [full_time(trial_time_idx)]; % tracker time, ms
        trial_data(1, :) = trial_data(1, :) - stim_start; % set stim onset to t=0 (so fixation time starts in negative)
        trial_data([2:3], :) = [Trials(trial_idx).Samples.gx(:, trial_time_idx)]; % x pos for both eyes
        trial_data([4:5], :) = [Trials(trial_idx).Samples.gy(:, trial_time_idx)]; % y pos for both eyes
        trial_data([6:7], :) = [Trials(trial_idx).Samples.pa(:, trial_time_idx)]; % pupil area for both eyes
        
        % remove timepoints during and 150ms before/after blinks
        % Blinks measured normally for stimulus period
        blink_indices = [];
        if isfield(Trials(trial_idx).Blinks, 'sttime')
            blink_indices = get_event_indices(Trials(trial_idx), 150, 150, 'blink');
            blink_indices = blink_indices + stim_start_idx; % align to recording onset
            blink_indices = blink_indices(blink_indices <= (trial_end_idx - trial_start_idx)); % keep within trial. indices are already aligned to stim onset..
            blink_indices = blink_indices - trial_start_idx; % align to fixation onset (t = ~-1000) and trial data matrix
        end
        % for fixation period, have to manually find them
        for eyepos = 2:5 % both eyes, x and y pos
            large_value_indices = find(trial_data(eyepos, :) > 10000); % outrageous eye position values
            for v = large_value_indices
                blink_indices = unique([blink_indices, v-150:v+150]);
            end
        end
        blink_indices = blink_indices(blink_indices > 0); % keep within trial
        
        trial_data([2:7], blink_indices) = NaN;
        
        trial_data = double(trial_data);
        posdata{bl}{iTrial} = trial_data(1:5, :);
        
        % preprocess pupil area over time
        curr_pupil = nanmean(trial_data([6:7], :), 1); % average both eyes
        % interpolate missing timepoints
        if sum(isnan(curr_pupil))
            nanstarts = []; nanends = []; currently_in_nan = 0;
            for t = 1:numel(curr_pupil)
                if isnan(curr_pupil(t))
                    if ~currently_in_nan
                        nanstarts = [nanstarts; t]; % mark start of nan sequence
                        currently_in_nan = 1;
                    end
                else
                    if currently_in_nan
                        nanends = [nanends; t-1]; % mark end of nan sequence
                        currently_in_nan = 0;
                    end
                end
            end
            if isnan(curr_pupil(1)); nanstarts(1) = []; nanends(1) = []; end % if trial begins with nan, don't interpolate the first section
            if isnan(curr_pupil(end)); nanstarts(end) = []; end % if trial ends with nan, don't interpolate the last section
            %linearly interpolate
            for nanseq = 1:numel(nanstarts)
                curr_pupil(nanstarts(nanseq)-1:nanends(nanseq)+1) = linspace(curr_pupil(nanstarts(nanseq)-1), curr_pupil(nanends(nanseq)+1), nanends(nanseq)-nanstarts(nanseq)+3);
            end
        end
        
        % filter
        nonnanidx = ~isnan(curr_pupil); % account for trial starting/ending with nan sequences
        curr_pupil_highfilt = highpass(curr_pupil(nonnanidx), 0.01, fps);% 0.01 Hz to remove low drift & zero mean
        curr_pupil_filt = lowpass(curr_pupil_highfilt, 6, fps); % 6 hz a la Brascamp et al 2021
        curr_pupil(nonnanidx) = [nan(1, 30), curr_pupil_filt(31:end-30), nan(1, 30)]; % replace transients with nans
        % resample?
        %curr_pupil_resample = resample(curr_pupil_filt, 10, fps); % 10 Hz a la Brascamp et al 2021
        % z score?, ignoring nans, relative to t=2000 onwards to avoid stim onset transient.
        %curr_pupil = (curr_pupil - nanmean(curr_pupil(3000:end))) ./ nanstd(curr_pupil(3000:end));
        % or by median (robust z score)
        %curr_pupil = (curr_pupil - nanmedian(curr_pupil)) ./ nanmedian(abs(curr_pupil - nanmedian(curr_pupil)));
        
        pupil_trace{bl}{iTrial} = curr_pupil;
        
        
        % assign each timepoint with whether or not microsaccade is ongoing
        ms_present{bl}{iTrial} = zeros(1, size(trial_data, 2));
        % run actual algorithm
                [~, all_xyVelocity{bl}{iTrial}] = UI_extractMicrosaccades(Trials(trial_idx));
                if ~isempty(Trials(trial_idx).Microsaccades)
                for ms = 1:numel(Trials(trial_idx).Microsaccades.Start)
                    msStart = Trials(trial_idx).Microsaccades.Start(ms);
                    msEnd = Trials(trial_idx).Microsaccades.End(ms);
                    all_ms_indices = (msStart:msEnd);% - trial_start_idx + 1;
                    all_ms_indices = all_ms_indices(all_ms_indices > 0 & all_ms_indices <= trial_end_idx);
                    if ~isempty(all_ms_indices)
                        ms_present{bl}{iTrial}(all_ms_indices) = 1;
                    end
                end
                end
                all_BinocularSaccades{bl}{iTrial} = Trials(trial_idx).Microsaccades;
        % get event indices, where trial start (fixation onset) = 1
        start_idx_trial = trial_start_idx - trial_start_idx + 1;
        stim_on_idx_trial = stim_start_idx - trial_start_idx + 1;
        button_idx_trial = button_press_idx - trial_start_idx + 1;
        end_idx_trial = trial_end_idx - trial_start_idx + 1; % should be last timepoint
        trial_indices{bl}(iTrial, :) = [start_idx_trial, stim_on_idx_trial, button_idx_trial, end_idx_trial];
    end
    if isempty(button_trials_use)
        all_xyVelocity{bl} = cell(1, 0); % account for a fully missing block
    end
end