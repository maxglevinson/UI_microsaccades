function [params, winMain] = UI_openScreen(params)
% opens screen for UI experiments and sets basic display parameters

scres=[1920 1080];%[1024 768]; % screen resolution we want
FR=120;%100; % frame rate we want

if ismac % max's macbook pro
    params.cm_to_screen=61;
    if params.eyeLink
        params.cm_to_screen=57;
    end
elseif ispc & contains(params.scripts_dir, 'C:') % if eyelink Windows machine
    params.cm_to_screen=60; % fixed with chin rest
elseif ispc & contains(params.scripts_dir, 'D:') % if MEG Windows machine
    params.cm_to_screen=52;
else
    params.cm_to_screen=56; % Participant's distance to the point of fixation.
end

% set up ProPixx projector
if ispc & contains(params.scripts_dir, 'D:') % if Windows machine (MEG stimulus PC)
    Datapixx('Open');
    Datapixx('RegWrRd');
    if params.rapidTag % rapid frequency tagging
        if params.real_FR == 480
            % initialize to use 480 Hz method (4 quadrants of screen)
            Datapixx('SetPropixxDlpSequenceProgram', 2); % 2 for 480, 5 for 1440 Hz, 0 for normal (120)
        elseif params.real_FR == 1440
            Datapixx('SetPropixxDlpSequenceProgram', 5); % each RGB channel becomes 1 grayscale frame
        end
    else
        Datapixx('SetPropixxDlpSequenceProgram', 6); % normal usage: 120 Hz
    end
    Datapixx('RegWrRd');
    % these 2 are taken from propixx 480 Hz instructions
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible'); % will use 32-bit precision if possible in combo with alpha blending
    % but if that's too taxing, it will switch to 16-bit precision
end

% Open main window, double-buffered fullscreen window
if params.debug
    PsychDebugWindowConfiguration;
else
    clear Screen
end
if ispc & contains(params.scripts_dir, 'D:') % MEG computer
    Screen('Preference', 'SkipSyncTests', 0); %  Apparently it is better to set this to 0 when you are running the experiment for real (ie not params.debugging or playing around).
elseif ~ismac & isunix % BIC workstation
    Screen('Preference', 'SkipSyncTests', 1);
elseif params.debug
    Screen('Preference', 'SkipSyncTests', 1);
else
    Screen('Preference', 'SkipSyncTests', 1);%0);
end
%Screen('Preference', 'VisualDebugLevel', 1); % makes startup screen black instead of white
screenNumber=max(Screen('Screens'));
if numel(Screen('Screens')) > 2
    screenNumber = 2;%1; % force a screen in MEG lab or eyelink room
end
if params.debug == 2
    [winMain, rectMain] = Screen('OpenWindow', screenNumber, 0, [0, 0, 1000, 1000], 32, 2, [], [], kPsychNeed32BPCFloat);
else
    %[winMain, rectMain] = Screen('OpenWindow', screenNumber);%, [], [], [], 2, [], 4);%,  kPsychNeed32BPCFloat);
%     if ismac & params.eyeLink && params.dummymode == 1
%          PsychImaging('PrepareConfiguration');
%          PsychImaging('AddTask', 'General', 'UseRetinaResolution');
%     end
    [winMain, rectMain] = PsychImaging('OpenWindow', screenNumber);
    %[winMain, rectMain] = Screen('OpenWindow', screenNumber, 0,[],32,2);
end
if params.useC48 % if using C48 color mode
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', 2); % mode 2
    % mode 2 averages adjacent pixels together to make 1 wide pixel.
end

fullX = rectMain(3);%1920;
fullY = rectMain(4);%1080;
if params.rapidTag
    rectMain = [0 0 fullX/2 fullY/2]; % first quadrant
end
% default alpha blending
Screen('BlendFunction', winMain, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');  % Set up alpha-blending for smooth (anti-aliased) lines
%Screen('BlendFunction', winMain, 'GL_ONE', 'GL_ZERO'); % maybe no blending?
%Screen('BlendFunction', winMain, 'GL_SRC_ALPHA', 'GL_ZERO');

%switch params.stimtype
%    case 'retinotopy'
%Screen('BlendFunction', winMain, 'GL_SRC_ALPHA', 'GL_ZERO');
%case 'shade'
%Screen('BlendFunction', winMain, 'GL_SRC_ALPHA', 'GL_ZERO');
%end
%Screen('LoadNormalizedGammaTable', winMain, gammaTable2*[1 1 1]);
ListenChar(2);

params.IFI = Screen('GetFlipInterval', winMain); % get length of frame
params.frame_rate = round(Screen('FrameRate', winMain));% get refresh rate
if params.frame_rate == 0 % sometimes happens on mac
    params.frame_rate = round(1 / params.IFI);
end

if ismac % force 30 Hz frame rate...
    params.IFI = 1/30;
    params.frame_rate = round(1/params.IFI);
end

params.white = WhiteIndex(winMain);
% if params.useL48
%     params.white = 65535;
% end
params.black = BlackIndex(winMain);

% variables related to display size
params.DisplaySizeX = rectMain(3);
params.DisplaySizeY = rectMain(4);
params.xmid = params.DisplaySizeX / 2;
params.ymid = params.DisplaySizeY / 2;
if params.shrinkDisplay % make stimulus display smaller for eye tracking experiment
    % but the other parameters like displaysizeX, etc refer to entire display
    params.totalRect = CenterRect(rectMain .* .75, rectMain);
else
    params.totalRect = [0 0 fullX fullY]; % i.e. rectMain
end
params.xstart_pos = round(params.DisplaySizeX/15) + rectMain(1);

% Screen resolution (for converting pixels to cm, etc).
set(0,'units','centimeters');
params.CM_SS = get(0,'screensize');
% if params.rapidTag
%     params.CM_SS = params.CM_SS ./ 2; % just get size of one quadrant?
% end
if ispc & contains(params.scripts_dir, 'D:') % MEG computer
    params.CM_SS = [0 0 29*(1920/1080) 29]; % set custom screen size by literally measuring the projected display.
elseif ismac & params.eyeLink % running eyelink experiment
    params.CM_SS = params.CM_SS;
elseif ismac & ~params.eyeLink % testing on mac, avoid retina issues
    params.CM_SS = [0 0 29 29*900/1440];
    %params.CM_SS = [0 0 50.8 31.75];
elseif ispc & contains(params.scripts_dir, 'C:') % eyelink pc (or home pc)
    [widthMM, heightMM] = Screen('DisplaySize', 2);
    params.CM_SS = [0 0 widthMM/10 heightMM/10];
end
params.Res = rectMain./params.CM_SS;
params.pxbycm_x=params.Res(3);
params.pxbycm_y=params.Res(4);
params.ppd = pi * params.DisplaySizeX / atan(params.CM_SS(3)/params.cm_to_screen/2) / 360; % pixels per degree of visual angle

% if params.frame_rate~=FR || isequal( [params.DisplaySizeX params.DisplaySizeY], scres) ==0
%     disp(['PLEASE CHANGE FRAME RATE TO ' num2str(FR) ' AND SCREEN SIZE TO ' num2str(scres(1)) 'X' num2str(scres(2))]);
%     sca;
% end

% Text
textsize = round(params.DisplaySizeX/60);
if params.shrinkDisplay
    textsize = round(textsize * 0.75);
end
params.wrep = round(params.DisplaySizeX);
Screen('TextFont', winMain, 'Arial');
Screen('TextSize', winMain, textsize);
Screen('TextStyle', winMain, 1); % 0=normal, 1=bold, 2=italic, 4=underlined, 8=outline

% Colors
params.colForeground = [160,160,160] .* params.white/255;
params.colFixation=[200,0,0] .* params.white/255;%[0, 0, 200];
switch params.stimtype
    case {'dotblur', 'dotsize'}
        params.colBackground = [.7 .7 0] .* params.white; % yellow
    case 'retinotopy'
        params.colBackground = [50 50 50] .* params.white/255;
    case 'color'
        if ispc & contains(params.scripts_dir, 'D:') % if Windows machine (MEG stimulus PC)
            % account for linear propixx, no gamma applied
            params.colBackground = .1 * [1 1 1] .* params.white;
            params.colForeground = .01 * [1 1 1] .* params.white;
        else
            params.colBackground = 0.3 * [1 1 1] .* params.white;
            params.colForeground = 0.05 * [1 1 1] .* params.white;
        end
    otherwise
        params.colBackground = [20 20 20] .* params.white/255;%[10 10 10];
end


% if params.useL48 % add colors to clut?
%     params.clut(4, :) = params.colBackground;
%     params.clut(5, :) = params.colForeground;
%     params.clut(6, :) = params.colFixation;
% end

%% initialize center area
if ~isfield(params, 'centerdva') % if not already specified
    params.centerdva = 5; % radius
end
params = UI_generate_center_rects(params);

%% set up MEG button triggers
if ispc & contains(params.scripts_dir, 'D:') % MEG lab stim PC
    SetButtonsAndPixelMode('Both');
    Datapixx('Open');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
end

% pixel mode color triggers
params.triggers = struct;
params.triggers.pixel1 = [0 1 0]; % start trial (fixation).
params.triggers.pixel2 = [0 2 0]; % frequency tagging: start flickering (stimulus on)
params.triggers.pixel3 = [0 4 0]; % frequency tagging: end flickering (stimulus off)
params.triggers.pixel4 = [0 8 0]; % end trial.
params.triggers.nopixeltrigger = [0 0 0]; % show this color to prevent a pixel mode trigger
params.triggers.pixeltriggerrect = [0; 0; 2; 2;]; % where to show pixel trigger on screen
if params.rapidTag
    newtriggerrect = [];
    for quad = 1:4
        newtriggerrect = [newtriggerrect, params.triggers.pixeltriggerrect + params.quadshifts(:, quad)];
    end
    params.triggers.pixeltriggerrect = newtriggerrect;
end
% if params.useL48 % add colors to clut?
%     params.clut(7, :) = params.triggers.pixel1;
%     params.clut(8, :) = params.triggers.pixel2;
%     params.clut(9, :) = params.triggers.pixel3;
%     params.clut(10, :) = params.triggers.pixel4;
%     params.clut(11, :) = params.triggers.nopixeltrigger;
% end