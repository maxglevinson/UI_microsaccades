function [params, winMain] = UI_openScreen(params)
% opens screen for UI experiment and sets basic display parameters

params.cm_to_screen = 60; % fixed with chin rest

% Open main window, double-buffered fullscreen window
if params.debug
    PsychDebugWindowConfiguration;
else
    clear Screen
end
if params.debug
    Screen('Preference', 'SkipSyncTests', 1);
else
    Screen('Preference', 'SkipSyncTests', 0);
end
screenNumber=max(Screen('Screens'));
if numel(Screen('Screens')) > 2
    screenNumber = 2; % force second screen
end
if params.debug == 2
    [winMain, rectMain] = Screen('OpenWindow', screenNumber, 0, [0, 0, 1000, 1000], 32, 2, [], [], kPsychNeed32BPCFloat);
else
    [winMain, rectMain] = PsychImaging('OpenWindow', screenNumber);
end

fullX = rectMain(3);
fullY = rectMain(4);
% default alpha blending
Screen('BlendFunction', winMain, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');  % Set up alpha-blending for smooth (anti-aliased) lines
ListenChar(2);

params.IFI = Screen('GetFlipInterval', winMain); % get length of frame
params.frame_rate = round(Screen('FrameRate', winMain));% get refresh rate
if params.frame_rate == 0 % sometimes happens on mac
    params.frame_rate = round(1 / params.IFI);
end

params.white = WhiteIndex(winMain); % should be 255
params.black = BlackIndex(winMain); % should be 0

% variables related to display size
params.DisplaySizeX = rectMain(3);
params.DisplaySizeY = rectMain(4);
params.xmid = params.DisplaySizeX / 2;
params.ymid = params.DisplaySizeY / 2;
if params.shrinkDisplay % make stimulus display smaller
    % but the other parameters like displaysizeX, etc refer to entire display
    params.totalRect = CenterRect(rectMain .* .75, rectMain);
else
    params.totalRect = [0 0 fullX fullY]; % i.e. rectMain
end
params.xstart_pos = round(params.DisplaySizeX/15) + rectMain(1);

% Screen resolution (for converting pixels to cm, etc).
set(0,'units','centimeters');
params.CM_SS = get(0,'screensize');
if ismac
    params.CM_SS = params.CM_SS;
elseif ispc
    [widthMM, heightMM] = Screen('DisplaySize', 2);
    params.CM_SS = [0 0 widthMM/10 heightMM/10];
end
params.Res = rectMain./params.CM_SS;
params.pxbycm_x=params.Res(3);
params.pxbycm_y=params.Res(4);
params.ppd = pi * params.DisplaySizeX / atan(params.CM_SS(3)/params.cm_to_screen/2) / 360; % pixels per degree of visual angle

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
params.colFixation=[200,0,0] .* params.white/255;
params.colBackground = 0.3 * [1 1 1] .* params.white;
params.colForeground = 0.05 * [1 1 1] .* params.white;


%% initialize center area
if ~isfield(params, 'centerdva') % if not already specified
    params.centerdva = 5; % radius default (but it should already be specified in params)
end
params = UI_generate_center_rects(params);