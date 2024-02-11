function Trials = UI_extractMicrosaccades(Trials, eta)
% custom microsaccade extraction from eyelink data
% Max Levinson, 2023
% adapted from edfImport library (A. Pastukhov),
% Martinez-Conde et al. 2006, and Engbert & Kliegl 2003.

if ~isfield(Trials, 'KeyEvents')
    Trials = edfExtractKeyEventsTiming(Trials); % needed for resolving a bug with blink detection (see get_event_indices)
end

if nargin < 2 % eta not specified
    eta = 5; % default velocity threshold multiplier
end

for iTrial = 1:length(Trials)
    badtrial = 0;
    try trial_start = Trials(iTrial).KeyEvents.FIX_ONSET; % tracker time, ms
    catch badtrial = 1;
    end
    if badtrial
        continue
    end

    fps = Trials(iTrial).Header.rec.sample_rate;
    dt = (1/fps) * 1000; % in ms
    xPos = Trials(iTrial).Samples.gx;
    yPos = Trials(iTrial).Samples.gy;

    % remove timepoints during and 150ms before/after blinks
    trial_start = Trials(iTrial).KeyEvents.FIX_ONSET; % tracker time, ms
    stim_start = Trials(iTrial).KeyEvents.STIM_ONSET; % tracker time, ms
    trial_end = Trials(iTrial).KeyEvents.STIM_OFFSET; % tracker time, ms
    full_time = Trials(iTrial).Samples.time; % tracker time
    trial_start_idx = find(full_time == trial_start); stim_start_idx = find(full_time == stim_start); trial_end_idx = find(full_time == trial_end);
    trial_start_idx = trial_start_idx(1);
    % Blinks measured normally for stimulus period
    blink_indices = [];
    if isfield(Trials(iTrial).Blinks, 'sttime')
        blink_indices = get_event_indices(Trials(iTrial), 150, 150, 'blink');
        blink_indices = blink_indices + stim_start_idx; % align to recording onset
        blink_indices = blink_indices(blink_indices <= (trial_end_idx - trial_start_idx)); % keep within trial. indices are already aligned to stim onset..
        blink_indices = blink_indices - trial_start_idx; % align to fixation onset (t = ~-1000) and trial data matrix
    end
    % for fixation period, have to manually find them
    for eye = 1:2 % both eyes, x and y pos
        large_value_indicesX = find(xPos(eye, :) > 10000); % outrageous eye position values
        large_value_indicesY = find(yPos(eye, :) > 10000);
        for v = [large_value_indicesX, large_value_indicesY]
            blink_indices = unique([blink_indices, v-150:v+150]);
        end
    end
    blink_indices = blink_indices(blink_indices > 0); % keep within trial

    xPosClean = xPos;
    yPosClean = yPos;
    xPosClean(:, blink_indices) = NaN;
    yPosClean(:, blink_indices) = NaN;

    % only keep timepoints within trial (trial start to trial end)
    trial_time_idx = full_time >= trial_start & full_time <= trial_end; % extract full trial indices
    xPos = xPos(:, trial_time_idx);
    yPos = yPos(:, trial_time_idx);
    xPosClean = xPosClean(:, trial_time_idx);
    yPosClean = yPosClean(:, trial_time_idx);

    xPos = xPos - nanmean(xPosClean, 2); % coordinates from center (0) (pixels)
    yPos = yPos - nanmean(yPosClean, 2);

    % calculate deg per pix from screen size (cm), distance, resolution
    % MNI eyelink room
    horzcm = 60.8;
    distcm = 60;
    pxbycm_x = 2560 / horzcm;
    pixperdeg = pi * 2560 / atan(horzcm/distcm/2) / 360; % pixels per degree of visual angle
    degperpix = 1 / pixperdeg;
    xPos = xPos * degperpix; % convert to dva
    yPos = yPos * degperpix;

    % computing velocity (adapted from edfExtractMicrosaccades.m, in edfImport)
    xVelocity = nan(size(xPos));
    yVelocity = nan(size(yPos));
    xyAngle = nan(size(xVelocity));

    window_length_ms = 31; % in ms, surrounding timepoint (should be odd number); 31 taken from Martinez-Conde 2006
    window_length_tps = window_length_ms / dt;

    dx = zeros(size(xPos));
    dy = zeros(size(yPos));
    xVelocity = dx;
    yVelocity = dy;
    xyVelocity = xVelocity;
    xyAngleRad = xyVelocity;
    xyAngleDeg = xyAngleRad;
    for iEye = 1:2
        dx(iEye, :) = [0, diff(xPos(iEye, :))];
        xVelocity(iEye, :) = smoothdata(dx(iEye, :), 'movmean', window_length_tps) * 1000;
        dy(iEye, :) = [0, diff(yPos(iEye, :))];
        yVelocity(iEye, :) = smoothdata(dy(iEye, :), 'movmean', window_length_tps) * 1000;
        xyVelocity(iEye, :) = hypot(xVelocity(iEye, :), yVelocity(iEye, :));
        xyAngleRad(iEye, :) = atan(yVelocity(iEye, :) ./ xVelocity(iEye, :));
        xyAngleDeg(iEye, :) = xyAngleRad(iEye, :) * 180 / pi;
    end

    % remove blink periods from velocity traces
    blink_indices = blink_indices(blink_indices < trial_end_idx); % keep within trial
    xVelocity(:, blink_indices) = NaN;
    yVelocity(:, blink_indices) = NaN;
    xyVelocity(:, blink_indices) = NaN;

    % detect possible microsaccades
    part_of_microsaccade = ones(size(xyVelocity)); % true or false for each timepoint and eye
    part_of_microsaccade(:, isnan(xyVelocity(1, :))) = 0;
    velInThreshold = zeros(size(xyVelocity));
    angle_change_threshold = 15; % if change in xyAngleDeg is > threshold, mark the timepoint as NOT a microsaccade
    microsaccade_details = [];
    for iEye = 1:2
        angle_change_deg = abs(diff(xyAngleDeg(iEye, :)));
        angle_change_deg(isnan(angle_change_deg)) = 0;
        angle_change_deg = [0, angle_change_deg]; % first timepoint change = 0
        part_of_microsaccade(iEye, angle_change_deg > angle_change_threshold) = 0;

        % set velocity threshold - eye is considered 'stopped' when velocity is
        % too slow (calculation taken from Engbert & Kliegl 2003)
        % multiplier eta is applied to noise level (sigma below)
        xsigma = sqrt(nanmedian(xVelocity .^ 2, 2) - (nanmedian(xVelocity, 2) .^ 2));
        ysigma = sqrt(nanmedian(yVelocity .^ 2, 2) - (nanmedian(yVelocity, 2) .^ 2));
        xthresh = eta * xsigma;
        ythresh = eta * ysigma;
        % which velocities are in threshold? AKA if output here is >1.
        velInThreshold(iEye, :) = hypot(xVelocity(iEye, :) ./ xthresh(iEye), yVelocity(iEye, :) ./ ythresh(iEye));
        part_of_microsaccade(iEye, velInThreshold(iEye, :) < 1) = 0;

        microsaccade_start = [0, diff(part_of_microsaccade(iEye, :))];
        microsaccade_start(microsaccade_start < 1) = 0;
        start_idx = find(microsaccade_start);
        microsaccade_details(iEye).Start = nan(numel(start_idx), 1);
        microsaccade_details(iEye).Duration = nan(numel(start_idx), 1);
        microsaccade_details(iEye).End = nan(numel(start_idx), 1);
        microsaccade_details(iEye).Amplitude = nan(numel(start_idx), 1);
        microsaccade_details(iEye).Phi = nan(numel(start_idx), 1);
        microsaccade_details(iEye).DeltaX = nan(numel(start_idx), 1);
        microsaccade_details(iEye).DeltaY = nan(numel(start_idx), 1);
        microsaccade_details(iEye).peakXVelocity = nan(numel(start_idx), 1);
        microsaccade_details(iEye).peakYVelocity = nan(numel(start_idx), 1);
        microsaccade_details(iEye).vPeak = nan(numel(start_idx), 1);
        for ms_idx = 1:numel(start_idx)
            xamp = 0;
            yamp = 0;
            start_timepoint = start_idx(ms_idx);
            microsaccade_details(iEye).Start(ms_idx) = start_timepoint;
            tp = start_timepoint;
            while part_of_microsaccade(iEye, tp) && tp < size(xVelocity, 2)
                xamp = xamp + xVelocity(iEye, tp) / fps;
                yamp = yamp + yVelocity(iEye, tp) / fps;
                tp = tp + 1;
            end
            amplitude = hypot(xamp, yamp);
            angle_rad = atan(yamp / xamp);
            angle_deg = angle_rad * 180 / pi;
            microsaccade_details(iEye).Duration(ms_idx) = tp - start_timepoint; % minimum = 1 timepoint
            microsaccade_details(iEye).End(ms_idx) = tp - 1;
            microsaccade_details(iEye).Amplitude(ms_idx) = amplitude;
            microsaccade_details(iEye).Phi(ms_idx) = angle_deg;
            microsaccade_details(iEye).DeltaX(ms_idx) = xamp;
            microsaccade_details(iEye).DeltaY(ms_idx) = yamp;
            microsaccade_details(iEye).peakXVelocity(ms_idx) = max(xVelocity(iEye, (start_timepoint):(start_timepoint+microsaccade_details(iEye).Duration(ms_idx))));
            microsaccade_details(iEye).peakYVelocity(ms_idx) = max(yVelocity(iEye, (start_timepoint):(start_timepoint+microsaccade_details(iEye).Duration(ms_idx))));
            microsaccade_details(iEye).vPeak(ms_idx) = max(xyVelocity(iEye, (start_timepoint):(start_timepoint+microsaccade_details(iEye).Duration(ms_idx))));
        end

        % merge microsaccades that occur too close to each other
        min_separation = 12; % at least 12 ms apart (tho Martinez Conde uses 20, to merge ms with overshoot)
        n_ms = numel(microsaccade_details(iEye));
        merged_list = false(1, numel(n_ms)); % keep track of which microsaccades we are merging
        for ms_idx = 2:numel(n_ms) % start at 2nd microsaccade
            if microsaccade_details(iEye).Start(ms_idx) - microsaccade_details(iEye).End(ms_idx-1) < min_separation
                merged_list(ms_idx) = true; % note that we are merging the current microsaccade.
                curr_merged_idx = 1; % default to merging with the previous one
                if merged_list(ms_idx-1) % if we've already merged the previous microsaccade, add this one too
                    previous_merged = 1;
                    while previous_merged % if there are many short ones in a row, keep going back until we hit the one that wasn't merged
                        curr_merged_idx = curr_merged_idx + 1; % have to go one more microsaccade back
                        if ~merged_list(ms_idx-curr_merged_idx)
                            previous_merged = 0;
                        end
                    end
                end
                microsaccade_details(iEye).End(ms_idx-curr_merged_idx) = microsaccade_details(iEye).End(ms_idx);
                microsaccade_details(iEye).Duration(ms_idx-curr_merged_idx) = microsaccade_details(iEye).End(ms_idx-curr_merged_idx) - microsaccade_details(iEye).Start(ms_idx-curr_merged_idx) + 1;
                microsaccade_details(iEye).DeltaX(ms_idx-curr_merged_idx) = microsaccade_details(iEye).DeltaX(ms_idx-curr_merged_idx) + microsaccade_details(iEye).DeltaX(ms_idx);
                microsaccade_details(iEye).DeltaY(ms_idx-curr_merged_idx) = microsaccade_details(iEye).DeltaY(ms_idx-curr_merged_idx) + microsaccade_details(iEye).DeltaY(ms_idx);
                microsaccade_details(iEye).Amplitude(ms_idx-curr_merged_idx) = hypot(microsaccade_details(iEye).DeltaX(ms_idx-curr_merged_idx), microsaccade_details(iEye).DeltaY(ms_idx-curr_merged_idx));
                angle_rad = atan(microsaccade_details(iEye).DeltaY(ms_idx-curr_merged_idx) / microsaccade_details(iEye).DeltaX(ms_idx-curr_merged_idx));
                angle_deg = angle_rad * 180 / pi;
                microsaccade_details(iEye).Phi(ms_idx-curr_merged_idx) = angle_deg;
                microsaccade_details(iEye).peakXVelocity(ms_idx-curr_merged_idx) = max(xVelocity(iEye, (microsaccade_details(iEye).Start(ms_idx-curr_merged_idx)):(microsaccade_details(iEye).End(ms_idx-curr_merged_idx))));
                microsaccade_details(iEye).peakYVelocity(ms_idx-curr_merged_idx) = max(yVelocity(iEye, (microsaccade_details(iEye).Start(ms_idx-curr_merged_idx)):(microsaccade_details(iEye).End(ms_idx-curr_merged_idx))));
                microsaccade_details(iEye).vPeak(ms_idx-curr_merged_idx) = max(xyVelocity(iEye, (microsaccade_details(iEye).Start(ms_idx-curr_merged_idx)):(microsaccade_details(iEye).End(ms_idx-curr_merged_idx))));
            end
        end
        allfields = fields(microsaccade_details(iEye));
        for f = 1:numel(allfields)
            microsaccade_details(iEye).(allfields{f})(merged_list) = [];
        end

        % remove eye movements too small / large, or too short in duration
        % important so that we don't get a "spurious binocular" microsaccade where
        % it's really monocular, but there's like 1 or 2 noisy timepoints in the other eye.
        % in Otero-Millan et al., 2008, avg magnitude = 0.5-0.6 deg
        % Otero-Millan et al 2008 has mean ~12 ms, range 5 - 30 ms
        % use essentially retinal slip for amplitude, so we remove eye movements
        % that traverse too much. Accounts for large merged 'square' microsaccades that
        % go in and out, so summed amplitude is closer to 0.
        min_duration = 3; % ms % this gets duration of monocular microsaccades. Then combined later to binocular, with final min duration
        min_duration_timepoints = round(min_duration * fps / 1000);
        max_amplitude = 2; % degrees
        min_amplitude = 3/60; % 3 min of arc, 0.05 deg
        all_amplitudes = microsaccade_details(iEye).Amplitude;
        bad_amplitude = (all_amplitudes < min_amplitude) | (all_amplitudes > max_amplitude);
        for i_ms = 1:numel(microsaccade_details(iEye).Start)
            curr_start = microsaccade_details(iEye).Start(i_ms);
            curr_end = microsaccade_details(iEye).End(i_ms);
            curr_slip = sum(abs(xyVelocity(iEye, curr_start:curr_end))) ./ fps; % total velocity / total time = total distance covered
            if curr_slip > max_amplitude
                bad_amplitude(i_ms) = 1;
            end
        end
        bad_amplitude = false(numel(microsaccade_details(iEye).Start), 1);
        all_durations = microsaccade_details(iEye).Duration;
        bad_duration = (all_durations < min_duration_timepoints);
        microsaccade_details(iEye).Start(bad_amplitude | bad_duration) = [];
        microsaccade_details(iEye).Duration(bad_amplitude | bad_duration) = [];
        microsaccade_details(iEye).End(bad_amplitude | bad_duration) = [];
        microsaccade_details(iEye).Amplitude(bad_amplitude | bad_duration) = [];
        microsaccade_details(iEye).Phi(bad_amplitude | bad_duration) = [];
        microsaccade_details(iEye).DeltaX(bad_amplitude | bad_duration) = [];
        microsaccade_details(iEye).DeltaY(bad_amplitude | bad_duration) = [];
        microsaccade_details(iEye).peakXVelocity(bad_amplitude | bad_duration) = [];
        microsaccade_details(iEye).peakYVelocity(bad_amplitude | bad_duration) = [];
        microsaccade_details(iEye).vPeak(bad_amplitude | bad_duration) = [];
    end


    %% check binocular correspondence
    % adapted from edfImport edfExtractMicrosaccades.m
    % we set the "left" eye to the eye with fewest microsaccades.
    % I guess assuming that all the binocular microsaccades we want are
    % dominated by "right" eye, and we just take those that also have some
    % timepoints in the "left."

    % so first make sure that "left" eye has fewer microsaccades
    if (length(microsaccade_details(2).Start)<length(microsaccade_details(1).Start))
        EyeTemp= microsaccade_details(2);
        microsaccade_details(2)= microsaccade_details(1);
        microsaccade_details(1)= EyeTemp;
    end

    % now find the binocular microsaccades
    bi_microsaccade_details = microsaccade_details(1); % start with left eye, then look for where left microsaccades are included in right ones.
    iBadLeftEyeSaccades= [];
    for iLS= 1:length(bi_microsaccade_details.Start)
        % find right eye microsaccade that starts before/when current left one ends, and ends after/on current left one starts.
        iOverlap= find(microsaccade_details(2).Start<=bi_microsaccade_details.End(iLS) & microsaccade_details(2).End>=bi_microsaccade_details.Start(iLS));
        if (~isempty(iOverlap)) % there can be more than one iOverlap, if left saccade is long and two short right saccades are included.
            bi_microsaccade_details.Start(iLS)= min([bi_microsaccade_details.Start(iLS); microsaccade_details(2).Start(iOverlap)]);
            bi_microsaccade_details.End(iLS)= max([bi_microsaccade_details.End(iLS); microsaccade_details(2).End(iOverlap)]);
            bi_microsaccade_details.Duration(iLS) = bi_microsaccade_details.End(iLS) - bi_microsaccade_details.Start(iLS) + 1;
            bi_microsaccade_details.vPeak(iLS)= mean([bi_microsaccade_details.vPeak(iLS); microsaccade_details(2).vPeak(iOverlap)]);
            bi_microsaccade_details.Amplitude(iLS)= mean([bi_microsaccade_details.Amplitude(iLS); microsaccade_details(2).Amplitude(iOverlap)]);
            bi_microsaccade_details.DeltaX(iLS)= mean([bi_microsaccade_details.DeltaX(iLS); microsaccade_details(2).DeltaX(iOverlap)]);
            bi_microsaccade_details.DeltaY(iLS)= mean([bi_microsaccade_details.DeltaY(iLS); microsaccade_details(2).DeltaY(iOverlap)]);
            angle_rad = atan(bi_microsaccade_details.DeltaY(iLS) / bi_microsaccade_details.DeltaX(iLS));
            angle_deg = angle_rad * 180 / pi;
            bi_microsaccade_details.Phi(iLS) = angle_deg; % pos = to the right, neg = to the left
        else
            iBadLeftEyeSaccades= [iBadLeftEyeSaccades iLS];
        end
    end

    bi_microsaccade_details.Start(iBadLeftEyeSaccades)= [];
    bi_microsaccade_details.Duration(iBadLeftEyeSaccades)= [];
    bi_microsaccade_details.End(iBadLeftEyeSaccades)= [];
    bi_microsaccade_details.peakXVelocity(iBadLeftEyeSaccades)= [];
    bi_microsaccade_details.peakYVelocity(iBadLeftEyeSaccades)= [];
    bi_microsaccade_details.vPeak(iBadLeftEyeSaccades)= [];
    bi_microsaccade_details.DeltaX(iBadLeftEyeSaccades)= [];
    bi_microsaccade_details.DeltaY(iBadLeftEyeSaccades)= [];
    bi_microsaccade_details.Amplitude(iBadLeftEyeSaccades)= [];
    bi_microsaccade_details.Phi(iBadLeftEyeSaccades)= [];


    %% merge microsaccades close together again
    % merge microsaccades that occur too close to each other
    min_separation = 12; % at least 12 ms apart (tho Martinez Conde uses 20, to merge ms with overshoot)
    merged_list = false(1, numel(bi_microsaccade_details.Start)); % keep track of which microsaccades we are merging
    for ms_idx = 2:numel(bi_microsaccade_details.Start) % start at 2nd microsaccade
        if bi_microsaccade_details.Start(ms_idx) - bi_microsaccade_details.End(ms_idx-1) < min_separation
            merged_list(ms_idx) = true; % note that we are merging the current microsaccade.
            curr_merged_idx = 1; % default to merging with the previous one
            if merged_list(ms_idx-1) % if we've already merged the previous microsaccade, add this one too
                previous_merged = 1;
                while previous_merged % if there are many short ones in a row, keep going back until we hit the one that wasn't merged
                    curr_merged_idx = curr_merged_idx + 1; % have to go one more microsaccade back
                    if ~merged_list(ms_idx-curr_merged_idx)
                        previous_merged = 0;
                    end
                end
            end
            bi_microsaccade_details.End(ms_idx-curr_merged_idx) = bi_microsaccade_details.End(ms_idx);
            bi_microsaccade_details.Duration(ms_idx-curr_merged_idx) = bi_microsaccade_details.End(ms_idx-curr_merged_idx) - bi_microsaccade_details.Start(ms_idx-curr_merged_idx) + 1;
            bi_microsaccade_details.DeltaX(ms_idx-curr_merged_idx) = bi_microsaccade_details.DeltaX(ms_idx-curr_merged_idx) + bi_microsaccade_details.DeltaX(ms_idx);
            bi_microsaccade_details.DeltaY(ms_idx-curr_merged_idx) = bi_microsaccade_details.DeltaY(ms_idx-curr_merged_idx) + bi_microsaccade_details.DeltaY(ms_idx);
            bi_microsaccade_details.Amplitude(ms_idx-curr_merged_idx) = hypot(bi_microsaccade_details.DeltaX(ms_idx-curr_merged_idx), bi_microsaccade_details.DeltaY(ms_idx-curr_merged_idx));
            angle_rad = atan(bi_microsaccade_details.DeltaY(ms_idx-curr_merged_idx) / bi_microsaccade_details.DeltaX(ms_idx-curr_merged_idx));
            angle_deg = angle_rad * 180 / pi;
            bi_microsaccade_details.Phi(ms_idx-curr_merged_idx) = angle_deg;
            bi_microsaccade_details.peakXVelocity(ms_idx-curr_merged_idx) = max(xVelocity(iEye, (bi_microsaccade_details.Start(ms_idx-curr_merged_idx)):(bi_microsaccade_details.End(ms_idx-curr_merged_idx))));
            bi_microsaccade_details.peakYVelocity(ms_idx-curr_merged_idx) = max(yVelocity(iEye, (bi_microsaccade_details.Start(ms_idx-curr_merged_idx)):(bi_microsaccade_details.End(ms_idx-curr_merged_idx))));
            bi_microsaccade_details.vPeak(ms_idx-curr_merged_idx) = max(xyVelocity(iEye, (bi_microsaccade_details.Start(ms_idx-curr_merged_idx)):(bi_microsaccade_details.End(ms_idx-curr_merged_idx))));
        end
    end
    allfields = fields(bi_microsaccade_details);
    for f = 1:numel(allfields)
        bi_microsaccade_details.(allfields{f})(merged_list) = [];
    end

    if isempty(bi_microsaccade_details.Start)
        continue
    end

    %% remove eye movements too large too small, or too short in duration
    % in Otero-Millan et al., 2008, avg magnitude = 0.5-0.6 deg
    % Otero-Millan et al 2008 has mean ~12 ms, range 5 - 30 ms
    min_duration = 8; % ms
    min_duration_timepoints = round(min_duration * fps / 1000);
    max_amplitude = 2; % degrees
    min_amplitude = 3/60; % 3 min of arc, 0.05 deg
    all_amplitudes = bi_microsaccade_details.Amplitude;
    bad_amplitude = (all_amplitudes < min_amplitude) | (all_amplitudes > max_amplitude);
    for i_ms = 1:numel(bi_microsaccade_details.Start)
        curr_start = bi_microsaccade_details.Start(i_ms);
        curr_end = bi_microsaccade_details.End(i_ms);
        curr_slip = nanmax(sum(abs(xyVelocity(:, curr_start:curr_end)), 2)) ./ fps; % total velocity / total time = total distance covered
        bi_microsaccade_details.Slip(i_ms) = curr_slip;
        if curr_slip > max_amplitude
            bad_amplitude(i_ms) = 1;
        end
    end
    all_durations = bi_microsaccade_details.Duration;
    bad_duration = (all_durations < min_duration_timepoints);
    bi_microsaccade_details.Start(bad_amplitude | bad_duration) = [];
    bi_microsaccade_details.Duration(bad_amplitude | bad_duration) = [];
    bi_microsaccade_details.End(bad_amplitude | bad_duration) = [];
    bi_microsaccade_details.Amplitude(bad_amplitude | bad_duration) = [];
    bi_microsaccade_details.Phi(bad_amplitude | bad_duration) = [];
    bi_microsaccade_details.DeltaX(bad_amplitude | bad_duration) = [];
    bi_microsaccade_details.DeltaY(bad_amplitude | bad_duration) = [];
    bi_microsaccade_details.peakXVelocity(bad_amplitude | bad_duration) = [];
    bi_microsaccade_details.peakYVelocity(bad_amplitude | bad_duration) = [];
    bi_microsaccade_details.vPeak(bad_amplitude | bad_duration) = [];
    bi_microsaccade_details.Slip(bad_amplitude | bad_duration) = [];

    for iS = 1:numel(bi_microsaccade_details.Start)
        bi_microsaccade_details.StartTime(iS) = Trials(iTrial).Samples.time(bi_microsaccade_details.Start(iS) + trial_start_idx - 1);
        bi_microsaccade_details.EndTime(iS) = Trials(iTrial).Samples.time(bi_microsaccade_details.End(iS) + trial_start_idx - 1);
    end

    % set microsaccade direction
    % positive = leftwards, negative = rightwards
    % 1 = outwards to periphery, 2 = inwards to center fixation
    lrDirections = sign(bi_microsaccade_details.DeltaX); % 1 = right, -1 = left
    ioDirections = nan(numel(bi_microsaccade_details.Start), 1);
    for iS = 1:numel(bi_microsaccade_details.Start)
        startTime = bi_microsaccade_details.Start(iS);
        endTime = bi_microsaccade_details.End(iS);
        startXY = hypot(nanmean(xPos(:, startTime)), nanmean(yPos(:, startTime)));
        endXY = hypot(nanmean(xPos(:, endTime)), nanmean(yPos(:, endTime)));
        if endXY > startXY % outwards microsaccade
            ioDirections(iS) = 1;
        else
            ioDirections(iS) = 2;
        end
    end
    bi_microsaccade_details.Direction = lrDirections .* ioDirections;

    if isempty(bi_microsaccade_details.Start)
        continue
    end

    %% save in Trials struct
    Trials(iTrial).Microsaccades = bi_microsaccade_details;
end


%% plot microsaccades for quality control

% plot each binocular microsaccade from center
%figure;plot([zeros(1, numel(bi_microsaccade_details.DeltaX)); bi_microsaccade_details.DeltaX'], [zeros(1, numel(bi_microsaccade_details.DeltaY)); bi_microsaccade_details.DeltaY']); % normalized to (0,0)

% plot binocular microsaccade amplitude against peak velocity
%figure;scatter(bi_microsaccade_details.Amplitude, bi_microsaccade_details.vPeak);
