function [data, timing, exitNow, bad_fix] = UI_runTrial(winMain, params, data, timing, stimulus, exitNow)

% To stop and restart the experiment if needed
[keyDown, keyTime, keyCode] = KbCheck(-1);
if keyDown
    if keyCode(params.keys.qkey)
        disp('Experiment on pause. Press any key to continue!');
        pause;
    end
end

uniform_frame_start = timing(4); % start of control shift (if control trial)
uniform_frame_end = timing(4) + floor(params.physical_shift_period ./ params.IFI);
if isnan(uniform_frame_start) % if illusion trial
    uniform_frame_start = Inf;
    uniform_frame_end = Inf;
end
shift_start = NaN; % default for illusion trial

% Start trial
nframes=round(timing(end)/params.IFI);

% display fixation cross at beginning of trial
if params.shrinkDisplay
    Screen('FillRect', winMain, [0;0;0]); % black underlying everything
    Screen('FillRect', winMain, params.colBackground, params.totalRect);
else
    Screen('FillRect', winMain, params.colBackground);
end
if params.fixation
    UI_drawFixation(winMain, params);
end
t_display=Screen('Flip', winMain); % Start of the display
if params.eyeLink
    Eyelink('Message', 'KEY_EVENT FIX_ONSET');
end


%% Get stimulus
colorCenter = data(12);
        colorPeriphery = data(13);
        [stimulus.fillColors, stimulus.frameColors] = UI_color_create(colorCenter, colorPeriphery, params);

nframeRects = size(stimulus.frameRects, 2);

%% additional params

% eyelink fixation tracking parameters
bad_fix = 0; % set to 1 if bad fixation occurs, which restarts trial
if params.eyeLink
    % if miss_fix counter becomes params.frametolerance + 1, restart trial.
    miss_fix = 0; % counter of how many bad frames you've had
    gaze_distances = zeros(1, nframes); % store distance from fixation point across frames
    gaze_distances(1:30) = nan; % starts measuring at frame 31
end

illusion_perceived_before=0;
start_uniformity = NaN; % assume no illusion perceived until button press
start_uniformity_frame = NaN;
lastframe = nframes + 1;

% CREATE MOVIE FROM STIMULUS
%stimMovie = Screen('CreateMovie', winMain, ['movie_', params.stimtype, '_', num2str(data(trial_idx, 2)), '.mov'], params.DisplaySizeX*2, params.DisplaySizeY*2, [], ':CodecSettings=EncodingQuality=0.5');

%% display stimulus
[keyDown, keyTime, keyCode] = KbCheck(-1); % set empty KeyDown

for frame=1:nframes
    if frame == lastframe; break; end % if we've reached 1 second past button press
    % Timestamp for this frame (also check whether a key has been pressed)
    if frame==1
        previous_frame_onset = t_display+params.fixationTime - params.IFI; % wait 1 second before frame starts
    end
    if frame > 30
            [keyDown, keyTime, keyCode] = KbCheck(-1);
        if params.eyeLink && ~illusion_perceived_before % check fixation if before button press
            % Eye tracking and gaze-contingent correction
            error=Eyelink('CheckRecording');
            if error ~= 0
                break;
            end
            if ~params.dummymode
                if Eyelink( 'NewFloatSampleAvailable') > 0
                    % get the sample in the form of an event structure
                    evt = Eyelink( 'NewestFloatSample');
                    if params.eye_used ~= -1 % do we know which eye to use yet?
                        % if we do, get current gaze position from sample
                        if params.eye_used < 2
                            xgaze = evt.gx(params.eye_used+1); % +1 as we're accessing MATLAB array
                            xgaze = evt.gy(params.eye_used+1);
                        else % binocular
                            xgaze = mean(evt.gx);
                            ygaze = mean(evt.gy);
                        end
                        % do we have valid data and is the pupil visible?
                        if xgaze~=params.MISSING_DATA && ygaze~=params.MISSING_DATA && max(evt.pa)>0 % if valid data
                            gaze_distances(frame) = sqrt((xgaze-params.xmid)^2 + (ygaze-params.ymid)^2);
                            if gaze_distances(frame) >= params.shiftolerance
                                miss_fix = miss_fix + 1; % count 1 more frame of bad fixation (saccade)
                            else
                                miss_fix = 0; % reset counter back to 0 if fixation is good
                            end
                        else % if invalid data (blink)
                            %miss_fix = miss_fix + 1; % count 1 more frame of bad fixation (blink)
                            miss_fix = miss_fix-10; % allow blinking; we don't want to be so evil
                        end
                    end
                end
            else % dummy mode, get mouse position
                [mx, my] = GetMouse(0); xgaze = mx; ygaze = my; % arg 1 with 2 monitors, arg 0 with just macbook
                %mx = params.xmid; my = params.ymid; % force good mouse fixation
                gaze_distances(frame) = sqrt((xgaze-params.xmid)^2 + (ygaze-params.ymid)^2);
                if gaze_distances(frame) >= params.shiftolerance
                    miss_fix = miss_fix + 1;
                else
                    miss_fix = 0;
                end
            end
            
            if miss_fix > params.frametolerance % 3 frames straight of bad fixation
                bad_fix = 1;
                Eyelink('Message', 'trial aborted');
                break
            end
        end
    end
    
    if params.shrinkDisplay
        Screen('FillRect', winMain, [0;0;0]); % black underlying everything
        Screen('FillRect', winMain, params.colBackground, params.totalRect);
    else
        Screen('FillRect', winMain, params.colBackground);
    end
    
    curr_fillColors = stimulus.fillColors;
    curr_frameColors = stimulus.frameColors;
    if data(2) == 2 % control shift
        target_color = mean(curr_fillColors, 2); % average of the two
        %target_color = curr_fillColors(:, 2); % center shade
        %target_color = curr_fillColors(:, 1); % periphery shade
        frameColor_distances = curr_frameColors(:, :) - target_color; % difference bw frame shade and center shade
        if frame<uniform_frame_end && frame>uniform_frame_start-1 % if control and transitioning to uniform
            frameScale = (uniform_frame_end - frame)/(uniform_frame_end-uniform_frame_start); % slowly scale down luminance of frame rects and periphery
            curr_frameColors(:, :) = curr_frameColors(:, :) - (frameColor_distances .* (1-frameScale));
            peripheryColor_distance = (curr_fillColors(:, 1) - target_color) * (1-frameScale); % slowly morph peripheral shade into target shade
            curr_fillColors(:, 1) = curr_fillColors(:, 1) - peripheryColor_distance; % if below 0, gets set to 0
            centerColor_distance = (curr_fillColors(:, 2) - target_color) * (1-frameScale); % slowly morph center shade into target shade
            curr_fillColors(:, 2) = curr_fillColors(:, 2) - centerColor_distance; % if below 0, gets set to 0
        elseif frame>uniform_frame_end-1 % if control and matching illusory frame now
            curr_fillColors(:, 1) = target_color; % only target shade across display
            curr_fillColors(:, 2) = target_color;
            curr_frameColors(:, :) = curr_frameColors(:, :) - frameColor_distances;
        end
    end
    
        curr_fillRectColors = curr_fillColors(:, 1);
        curr_fillOvalColors = curr_fillColors(:, 2);
        curr_fillRectRects = stimulus.fillRects(:, 1);
        curr_fillOvalRects = stimulus.fillRects(:, 2);
        curr_frameOvalColors = curr_frameColors;
    
    Screen('FillRect', winMain, [0;0;0]); % black underlying everything
    Screen('FillRect', winMain, curr_fillRectColors, curr_fillRectRects);
    Screen('FillOval', winMain, curr_fillOvalColors, curr_fillOvalRects);
    if data(2) ~= 3 % if not sharp catch trial
        Screen('FrameOval', winMain, curr_frameOvalColors, stimulus.frameRects, stimulus.thickness);
    end
    
    if params.fixation
        UI_drawFixation(winMain, params, 0);
    end
    
    
    %     %%%%% DEBUGGING: DRAW FIXATION WINDOW FOR REFERENCE
    %    Screen('FrameRect', winMain, [255;0;0], CenterRectOnPoint([-shiftolerance, -shiftolerance, shiftolerance, shiftolerance], params.xmid, params.ymid));
    %     %%%%%
    %     %%%%% DRAW EYE POSITION
    %     if frame > 60 && params.eyeLink
    %         Screen('FillOval', winMain, [255;255;255], CenterRectOnPoint(2*[-1,-1,1,1], xgaze, ygaze));
    %     end
    
    previous_frame_onset=Screen('Flip', winMain, previous_frame_onset+(params.IFI*(0.5)));
    % SAVE FRAME TO MOVIE
    %Screen('AddFrameToMovie', winMain, [], [], stimMovie);
    
    if frame == 1
        stim_start = previous_frame_onset;
        if params.eyeLink
            Eyelink('Message', 'KEY_EVENT STIM_ONSET');
            Eyelink('Message', num2str(stim_start));
        end
    elseif frame == uniform_frame_start
        shift_start = previous_frame_onset;
        if params.eyeLink
            Eyelink('Message', 'KEY_EVENT SHIFT_ONSET');
            Eyelink('Message', num2str(shift_start));
        end
    end
    
    if frame > 30 && (keyDown) % avoid incidental button presses at beginning of trial
        if keyCode(params.keys.escape)
            exitNow=1; bad_fix = 0; break
        elseif keyCode(params.keys.button2use)
            if illusion_perceived_before==0  % Only if this is the first time the illusion is perceived
                data(6)=1;  % Illusion was perceived at some point during this trial
                    start_uniformity = keyTime;
                    if params.eyeLink
                        Eyelink('Message', 'KEY_EVENT BUTTON_PRESS');
                        Eyelink('Message', num2str(start_uniformity));
                    end
                start_uniformity_frame = frame-1; % frame when button was pressed
                params.post_uniformity_frames = params.post_uniformity / params.IFI; % how many more frames to show before trial end
                lastframe = ceil(frame - 1 + params.post_uniformity_frames);
            end
            illusion_perceived_before=1;
        end
    end
end

if params.shrinkDisplay
    Screen('FillRect', winMain, [0;0;0]); % black underlying everything
    Screen('FillRect', winMain, params.colBackground, params.totalRect);
else
    Screen('FillRect', winMain, params.colBackground);
end

nNoiseFrames = round(lastframe / 2); % show noise for half the trial length
if ~bad_fix & ~exitNow
    [stimulus_offset] = UI_wipeScreen(winMain, params, nNoiseFrames, redCenter, previous_frame_onset, data(12), data(13));
end

%% save data
%Screen('FinalizeMovie', stimMovie); % SAVE MOVIE
if exitNow | bad_fix
    return
end
if params.shrinkDisplay
    Screen('FillRect', winMain, [0;0;0]); % black underlying everything
    Screen('FillRect', winMain, params.colBackground, params.totalRect);
else
    Screen('FillRect', winMain, params.colBackground);
end
% gather timing information
timing(7) = stimulus_offset;
timing(1) = t_display;
timing(2) = stim_start;
timing(3) = shift_start;
timing(5) = start_uniformity;
timing(6) = start_uniformity_frame;

% gather gaze statistics (in dva)
if params.eyeLink
    data(8) = nanmedian(gaze_distances ./ params.ppd);
    data(9) = nanmean(gaze_distances ./ params.ppd);
    data(10) = max(gaze_distances ./ params.ppd);
end

if data(6)==1  % subject has reported experiencing the illusion
    data(7) = timing(5) - timing(2); % reaction time
else
    data(6) = 0;
end
disp(data(7)); % show experimenter the reaction time
end


% -------- %
% FUNCTION
% -------- %
function [stimulus_offset] = UI_wipeScreen(winMain, params, nNoiseFrames, redCenter, previous_frame_onset, centerdata, peripherydata)
% displays noise across entire screen
% may help get rid of after-images and serial dependence
% nframes: how many frames to show noise for
% centerdata and peripherydata = used for color stimuli; dkl azimuths
% otherwise ignored

% choose overlay color - should be relatively opponent to stimulus color
if centerdata < 180
            center_opponent_deg = centerdata + 180;
        else
            center_opponent_deg = centerdata - 180;
        end
        if peripherydata < 180
            periphery_opponent_deg = peripherydata + 180;
        else
            periphery_opponent_deg = peripherydata - 180;
        end
        center_overlay = UI_dkl2rgb(params, center_opponent_deg);
        center_overlay = [center_overlay; .7 * params.white]; % alpha of 0.7
        periphery_overlay = UI_dkl2rgb(params, periphery_opponent_deg);
        periphery_overlay = [periphery_overlay; .7 * params.white]; % alpha of 0.7
    
for f = 1:nNoiseFrames
    noiseimg = 15*randn(floor(params.DisplaySizeX/8), floor(params.DisplaySizeY/8)) + 60; % mean 60, stddev 15
    noisetex = Screen('MakeTexture', winMain, noiseimg, [], 4);
    Screen('DrawTexture', winMain, noisetex, [], params.totalRect, [], 0);
        Screen('FillRect', winMain, periphery_overlay, params.totalRect + [0 0 -1 -1]);
        % (need to remove 1 pixel from each dimension otherwise alpha blending doesn't work)
        if redCenter
            Screen('FillOval', winMain, center_overlay, params.Central_rect);
        end
    % After drawing, we can discard the noise texture.
    Screen('Close', noisetex);
    previous_frame_onset=Screen('Flip', winMain, previous_frame_onset+(params.IFI*(0.5)));
    if f == 1
        stimulus_offset = previous_frame_onset;
        if params.eyeLink
            Eyelink('Message', 'KEY_EVENT STIM_OFFSET');
            Eyelink('Message', num2str(stimulus_offset));
        end
    end
end
end