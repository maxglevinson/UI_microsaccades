%% GENERAL PROPERTIES
clearvars
close all
clear all
params = struct;
params.isDemo = 0; % set to 1 for demo, just to show what filling-in should look like
params.isPractice = 0; % set to 1 for practice (run after isoluminant generation)
params.debug = 0; % debug mode?
params.fixation = 1; % include fixation point?
params.shrinkDisplay = 1; % make stimulus take up only 75% of screen?
params.eyeLink = 1; % leave this on unless doing some sort of testing/debugging
params.useCalibratedIsoluminance = 1; % use isoluminant RGB values generated per subject (1), or online dkl conversion (0)?
params.dummymode = 0; % change to 1 if testing eyelink code without using actual tracker
params.showCursor = 0; % show mouse cursor?
if params.dummymode
    params.showCursor = 1; % always show cursor in dummy mode so we can see it
end
params.writeData = 1; % 1 for experiments, 0 for testing/practice
params.button2use = 'space'; % which button does participant use? KbName code ('space' or '1!')
Priority(1);

% response and experimenter keys
KbName('UnifyKeyNames');
KbCheck;
params.keys.escape = KbName('escape');
params.keys.qkey = KbName('q');
params.keys.button2use = KbName(params.button2use);

params.stimtype = 'color';
params.subj_id = input('Subject number: ', 's');
params.session_id = input('Session number: ', 's');

params.scripts_dir='./'; % replace with local directory of experiment code
addpath(genpath(params.scripts_dir));

% trial timing details
params.fixationTime = 1; % fixation time before stimulus appears
params.shiftstart_range = [4 6]; % range of possible starts of physical uniformity shift, in seconds
params.physical_shift_period = 2; % how many seconds it takes to shift to physical uniformity
%params.noiseTime = 3; % how much time in seconds to show dynamic noise at the end of each trial
% ^^ not used, changed instead to half of trial length
params.trial_length_max = 20; % maximum trial length w/o button press
params.post_uniformity = 1; % how much time to keep stimulus up after button press

% block details
params.nBlocks = 3; % 3 per session
params.nTrials_per_block = 50; % trials per block
params.nIllusion_trials_per_block = 40;
params.nControlShift_trials_per_block = 5;
params.nControlSharp_trials_per_block = 5; % no soft transition bw center & periphery. should at least slightly increase RT
% set block types. 1 = smallest luminance difference, 3 largest luminance difference
% show all periphery values for each center size value
params.blockOrder = [1 2 3]; % for periphery color; randomize offline and set here manually
params.blockSizeOrder = [6 4 2]; % randomize offline and set here manually
% (run randomize_blocks.m once per subject to generate list of random block
% orders)

% set demo example parameters
% shift trials so I can show the stimulus for a while and then it
% demonstrates what the experience should be like.
if params.isDemo
    params.eyeLink = 0;
    params.writeData = 0;
    params.useCalibratedIsoluminance = 0;
    params.nBlocks = 1;
    params.nTrials_per_block = 5;
    params.nIllusion_trials_per_block = 0;
    params.nControlShift_trials_per_block = 5;
    params.nControlSharp_trials_per_block = 0;
    params.blockOrder = 2; % middle color
    params.blockSizeOrder = 4; % middle size
    params.shiftstart_range = [8 10];
end

% set practice parameters. run this after isoluminant generation
if params.isPractice
    params.eyeLink = 1;
    params.writeData = 1;
    params.session_id = '0'; % session 0 means practice
    params.useCalibratedIsoluminance = 1;
    params.nBlocks = 1;
    params.nTrials_per_block = 10;
    params.nIllusion_trials_per_block = 10;
    params.nControlShift_trials_per_block = 0;
    params.nControlSharp_trials_per_block = 0;
    params.blockOrder = 2; % middle color
    params.blockSizeOrder = 4; % middle size
end

% eyelink fixation tracking parameters
if params.eyeLink
    params.fixlimit_dva = 2; % how many dva the gaze can be vertically / horizontally off from center
    params.frametolerance = 2; % how many bad frames you're allowed to have before trial ends (restarts)
end

% set stimulus parameters
switch params.stimtype
    case 'color'
        load('phosphors.mat'); % monitor RGB phosphors (template, not calibrated)
        params.phosphors = phosphors;
        load('fundamentals_ss.mat'); % cone fundamentals, Stockman & Sharpe
        params.cone_fundamentals = fundamentals;
        params.dkl_radius = 0.07;
        params.colorCenter = 270; % dkl azimuth, in degrees
        params.colorPeriphery = [280;290;300];
        % set below for if we want to skip the online dkl calculation.
        % but make sure above azimuth values are still accurate!!
        if params.useCalibratedIsoluminance
            s_colors = load(['rgb_colors_', params.subj_id, '.mat']);
            params.centerRGBforced = s_colors.colors1.center; % 1 x n
            params.peripheryRGBforced = [s_colors.colors1.periphery, s_colors.colors2.periphery, s_colors.colors3.periphery]; % 3 x n
        end
end

if params.writeData
    params.data_dir=[params.scripts_dir, 'bhv_data/'];
    mkdir(params.data_dir);
    params.filestart = [params.data_dir params.stimtype '_' params.subj_id '_session_' params.session_id];
    params.datafile = [params.filestart, '_data.mat'];
end

%% SET UP WINDOW AND CENTER/PERIPHERY

[params, winMain] = UI_openScreen(params);

% System tests
AssertOpenGL;
rng('default'); rng shuffle;
if ~params.showCursor
    HideCursor;
end

%% SET UP EYELINK IF USING
if params.eyeLink
    el=EyelinkInitDefaults(winMain);
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    if ~EyelinkInit(params.dummymode, 1)
        fprintf('Eyelink Init failed.\n');
        Eyelink('Shutdown');
        sca;
        ListenChar(0);
        Priority(0);
        return;
    end
    
    % set params
    el.backgroundcolour = params.colBackground;
    el.foregroundcolour = params.colForeground;
    el.calibrationtargetcolour = params.colForeground;
    EyelinkUpdateDefaults(el);
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, params.DisplaySizeX-1, params.DisplaySizeY-1);
    
    
    % make sure that we get gaze data from the Eyelink (from Yair then amended from Eyelink youtube tutorial)
    Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
    Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,BLINK,SACCADE,BUTTON,FIXUPDATE,INPUT');
    Eyelink('command', 'file_sample_filter = LEFT,RIGHT,GAZE,SACCADE,BLINK,MESSAGE,AREA,GAZERES,STATUS');
    Eyelink('command', 'file_sample_data = LEFT,RIGHT,GAZE,HREF,RAW,AREA,GAZERES,BUTTON,STATUS,INPUT');
    Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
    
    Eyelink('command', 'enable_automatic_calibration = YES');
    Eyelink('command', 'calibration_type = HV9'); % or HV5/9/13
    % shrinked calibration window for sub 210 session 4 (last 2 blocks of session 2)
    % and all following subjects used shrinked area, plus HV9 instead of HV5.
    Eyelink('command', 'calibration_area_proportion = 0.3 0.3');
    Eyelink('command', 'validation_area_proportion = 0.3 0.3');
    
    % open file to record data to
    edfFile=['sub', params.subj_id, 's', params.session_id, '.edf'];
    Eyelink('Openfile', edfFile);
    
    params.eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
    %     if params.eye_used == el.BINOCULAR % if both eyes are tracked
    %         params.eye_used = el.LEFT_EYE; % use left eye
    %     end
    
    params.MISSING_DATA = el.MISSING_DATA;
    
    % calculate allowed gaze offset from fixation, in pixels
    params.shiftolerance = params.fixlimit_dva * params.ppd; % i.e., the dimensions of a square box surrounding the fixation point
    
end

%% ARRAY OF TRIALS

% Create trial array

all_trials=params.nBlocks*params.nTrials_per_block;

data = NaN(all_trials,11);
% C1 = subject number
% C2 = trial type (1=illusion, 2=control shift, 3=no soft border)
% C3 = block number in the session (1-3)
% C4 = block color contrast level (1, 2, 3 = low, medium, high)
% C5 = trial number within block
% C6 = Y/N did subject press the button (uniformity perceived)
% C7 = reaction time from stimulus onset
% C8=median gaze distance to the fixation point (in dva)
% C9=mean distance to the fixation point (in dva)
% C10=max distance to the fixation point (in dva)
% C11=center radius (dva)
% C12=center stimulus color (deg)
% C13=periphery stimulus color (deg)

timing = NaN(all_trials, 8);
% C1 = trial start (fixation)
% C2 = stimulus onset
% C3 = start of physical uniformity shift (control trial only)
% C4 = frame index of start of physical uniformity shift (control trial only)
% C5 = button press time (signaling perceived uniformity)
% C6 = frame index of button press
% C7 = stimulus offset (trial end)
% C8 = maximum display time if no button press

% fill in general characteristics
data(:, 1) = str2num(params.subj_id); % subject number
timing(:, 8) = params.trial_length_max;

% Set Experimental Blocks parameters
trialtypes = [ones(params.nIllusion_trials_per_block, 1); 2 * ones(params.nControlShift_trials_per_block, 1); 3 * ones(params.nControlSharp_trials_per_block, 1)];
%trialtypes = 1 * ones(params.nTrials_per_block, 1); % for testing, make all trials the same type
for block=1:params.nBlocks
    
    % BLOCK, BLOCK TYPE AND TRIAL
    
    blocktrial_1=(block-1)*params.nTrials_per_block+1;
    blocktrials=[blocktrial_1:block*params.nTrials_per_block];
    data(blocktrials, 2) = Shuffle(trialtypes); % randomize trial order within each block
    data(blocktrials, 3) = block;
    data(blocktrials, 4) = params.blockOrder(block);
    controlshiftblocktrials = data(blocktrials, 2) == 2;
    controlcatchblocktrials = data(blocktrials, 2) == 3;
    
    data(blocktrials,5)=[1:params.nTrials_per_block]';
    
    % STIMULUS CHARACTERISTICS
    data(blocktrials, 12) = params.colorCenter;
    data(blocktrials, 13) = params.colorPeriphery(params.blockOrder(block));
    
    % get random start of physical shift, between range given in params.shiftstart_range
    shiftstart_range_frames = (floor(params.shiftstart_range(1) ./ params.IFI):floor(params.shiftstart_range(2) ./ params.IFI));
    timing(blocktrials(controlshiftblocktrials), 4) = datasample(shiftstart_range_frames, sum(controlshiftblocktrials));
    
    % set center size for this block
    data(blocktrials, 11) = params.blockSizeOrder(block);
end

%% INITIALIZE STIMULUS
stimulus = struct;
[stimulus.fillRects, stimulus.frameRects, stimulus.thickness] = UI_color_init(params);

%% DISPLAY INSTRUCTIONS
exitNow = 0;

if params.eyeLink
    EyelinkDoTrackerSetup(el); % can calibrate now before the block starts
    WaitSecs(0.1);
end

start = 0;
while start == 0
    [keyDown, keyTime, keyCode] = KbCheck(-1);
    if params.shrinkDisplay
        Screen('FillRect', winMain, [0;0;0]); % black underlying everything
        Screen('FillRect', winMain, params.colBackground, params.totalRect);
    else
        Screen('FillRect', winMain, params.colBackground);
    end
    UI_drawInstructions('maintask_instructions_color', winMain, params);
    Screen('Flip', winMain);
    if keyDown
        if keyCode(params.keys.escape)
            exitNow = 1;
            break
        end
        start=1;
    end
end

if params.shrinkDisplay
    Screen('FillRect', winMain, [0;0;0]); % black underlying everything
    Screen('FillRect', winMain, params.colBackground, params.totalRect);
else
    Screen('FillRect', winMain, params.colBackground);
end
UI_drawInstructions('maintask_instructions_color', winMain, params);
Screen('Flip', winMain);
WaitSecs(.2);

%% RUN EXPERIMENT
for block = 1:params.nBlocks
    if exitNow == 1
        break
    end
    % To stop and restart the experiment if needed
    [keyDown, keyTime, keyCode] = KbCheck(-1);
    if keyDown
        if keyCode(params.keys.qkey)
            disp('Experiment on pause. Press any key to continue!');
            pause;
        elseif keyCode(params.keys.escape)
            exitNow = 1;
            break
        end
    end
    
    blocktrials_total=params.nTrials_per_block;
    blocktrial_1=(block-1)*params.nTrials_per_block+1;
    blocktrials=[blocktrial_1:block*params.nTrials_per_block];
    
    % because changing center size every block, re-init stimulus
    params.centerdva = params.blockSizeOrder(block);
    params = UI_generate_center_rects(params);
    [stimulus.fillRects, stimulus.frameRects, stimulus.thickness] = UI_color_init(params);
    
    if params.eyeLink && block > 1 % don't re-calibrate for first block
        disp('ready to calibrate'); % let experimenter know they are ready
        EyelinkDoTrackerSetup(el); % can calibrate again
    end
    WaitSecs(0.1);
    
    start = 0;
    WaitSecs(.2);
    while start == 0
        [keyDown, keyTime, keyCode] = KbCheck(-1);
        if params.shrinkDisplay
            Screen('FillRect', winMain, [0;0;0]); % black underlying everything
            Screen('FillRect', winMain, params.colBackground, params.totalRect);
        else
            Screen('FillRect', winMain, params.colBackground);
        end
        UI_drawInstructions('blockstart', winMain, params);
        Screen('Flip', winMain);
        if keyDown
            if keyCode(params.keys.escape)
                exitNow = 1;
                break
            end
            start=1;
            disp(['block ' num2str(block) ' started'])
        end
    end
    
    
    %%%%%%%% Trials %%%%%%%%%
    
    for trial = 1:blocktrials_total
        disp(['trial ' num2str(trial)]);
        bad_fix = 1; % FOR EYELINK ONLY: when a trial is successfully completed with good fixation, this is changed to 0
        % otherwise repeat trial until it is successful
        % when not using eyelink, this is just changed to 0 so the trial
        % is always successful.
        while bad_fix & ~exitNow
            trial_idx = blocktrials(trial);
            if params.shrinkDisplay
                Screen('FillRect', winMain, [0;0;0]); % black underlying everything
                Screen('FillRect', winMain, params.colBackground, params.totalRect);
            else
                Screen('FillRect', winMain, params.colBackground);
            end
            UI_drawInstructions('trialstart', winMain, params);
            Screen('Flip', winMain);
            
            pause(0.3);
            start = 0;
            while start == 0
                [keyDown, keyTime, keyCode] = KbCheck(-1);
                if params.shrinkDisplay
                    Screen('FillRect', winMain, [0;0;0]); % black underlying everything
                    Screen('FillRect', winMain, params.colBackground, params.totalRect);
                else
                    Screen('FillRect', winMain, params.colBackground);
                end
                UI_drawInstructions('trialstart', winMain, params);
                Screen('Flip', winMain);
                if keyDown
                    if keyCode(params.keys.escape)
                        exitNow = 1;
                        break
                    elseif keyCode(params.keys.button2use)
                        start = 1;
                    end
                end
            end
            % eye tracking
            if params.eyeLink
                Eyelink('Message', 'TRIALID %d', trial); % in case we want to use the Eyelink Dataviewer trial integration
                Eyelink('SetOfflineMode'); % first set offline mode before starting
                Eyelink('StartRecording');
                WaitSecs(0.1); % let it accumulate some samples first
                Eyelink('Message', 'TRIAL_VAR subject %s', params.subj_id);
                Eyelink('Message', 'TRIAL_VAR block %d', block);
                Eyelink('Message', 'TRIAL_VAR trial %d', trial);
                Eyelink('Message', 'TRIAL_VAR radius %d', params.blockSizeOrder(block));
                Eyelink('Message', 'TRIAL_VAR periphery_value %s', num2str(data(trial_idx, 13)));
                params.eye_used = Eyelink('EyeAvailable');
            end
            
            [data(trial_idx, :), timing(trial_idx, :), exitNow, bad_fix] = UI_runTrial(winMain, params, data(trial_idx, :), timing(trial_idx, :), stimulus, exitNow);
            
            if params.eyeLink
                Eyelink('Message', 'TRIAL_RESULT 0'); % if using dataviewer thing
                Eyelink('StopRecording');
            end
            
            if bad_fix
                disp('trial restarted');
                if params.shrinkDisplay
                    Screen('FillRect', winMain, [0;0;0]); % black underlying everything
                    Screen('FillRect', winMain, params.colBackground, params.totalRect);
                else
                    Screen('FillRect', winMain, params.colBackground);
                end
                UI_drawInstructions('bad_fix', winMain, params);
                Screen('Flip', winMain);
                pause(2);
            end
            
            if exitNow
                break
            end
        end
    end
    
    %%%%%%%%        %%%%%%%%%
    disp('block over');
    if params.writeData
        save(params.datafile, 'params', 'data', 'timing', 'stimulus'); % save after each block
    end
    if block<params.nBlocks
        if params.shrinkDisplay
            Screen('FillRect', winMain, [0;0;0]); % black underlying everything
            Screen('FillRect', winMain, params.colBackground, params.totalRect);
        else
            Screen('FillRect', winMain, params.colBackground);
        end
        UI_drawInstructions('blockrest', winMain, params);
        Screen('Flip', winMain);
        WaitSecs(10);
        
        start = 0;
        while start == 0
            [keyDown, keyTime, keyCode] = KbCheck(-1);
            if params.shrinkDisplay
                Screen('FillRect', winMain, [0;0;0]); % black underlying everything
                Screen('FillRect', winMain, params.colBackground, params.totalRect);
            else
                Screen('FillRect', winMain, params.colBackground);
            end
            UI_drawInstructions('blockrest', winMain, params);
            Screen('Flip', winMain);
            if keyDown
                if keyCode(params.keys.escape)
                    exitNow = 1;
                    break
                end
                start=1;
            end
        end
        
    elseif block==params.nBlocks
        if params.shrinkDisplay
            Screen('FillRect', winMain, [0;0;0]); % black underlying everything
            Screen('FillRect', winMain, params.colBackground, params.totalRect);
        else
            Screen('FillRect', winMain, params.colBackground);
        end
        UI_drawInstructions('experiment_end', winMain, params);
        Screen('Flip', winMain);
    end
    if params.writeData
        save(params.datafile, 'params', 'data', 'timing', 'stimulus'); % save final data
    end
    if exitNow
        break
    end
end
if params.eyeLink
    Eyelink('CloseFile');
    % download data file
    try
        fprintf('Receiving data file ''%s''\n', edfFile );
        status=Eyelink('ReceiveFile');
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2==exist(edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
        end
    catch rdf
        fprintf('Problem receiving data file ''%s''\n', edfFile );
        rdf;
    end
    Eyelink('Shutdown');
end

WaitSecs(3);
ListenChar(0);
Priority(0);
ShowCursor;
Screen('CloseAll');