% script to generate isoluminant rgb values
% first computes rgb from dkl, then allows individual to slightly alter
% one rgb to be perceptually isoluminant
% (using heterochromatic flicker technique)

% parameters
useDebug = 0; % debug mode?
%color_idx = 3; % 1, 2 or 3: which periphery color we are using. now input when calling function
center_azimuth = 270;
all_periphery_azimuths = [280; 290; 300];
dkl_radius = 0.07; % 0.07 for eyelink; possibly different for propixx?
bg_lum = 0.3; % background grey luminance to match, between 0 - 1
    % set this to 0.3 for eyelink experiment, probably 0.1 for propixx



subj_id = input('Subject number: ', 's');
color_idx = str2num(input('Which color? 1 2 or 3', 's'));
if ismac % max's macbook pro
    scripts_dir='~/Documents/McGill/neurospeed/Uniformity_Illusion/';
elseif isunix % rosaline BIC workstation
    scripts_dir='/export03/data/mlevin/Uniformity_Illusion/';
else % eyelink room or home pc
    scripts_dir='C:\Users\tneuro3\Documents\MATLAB\Max\Uniformity_Illusion\';
    %params.scripts_dir='C:\Users\max\Documents\neurospeed\Uniformity_Illusion\';
end
addpath(genpath(scripts_dir));
data_dir=[scripts_dir, 'bhv_data/'];
mkdir(data_dir);
datafile = [data_dir 'rgb_colors_' subj_id, '.mat'];
try
    load(datafile);
catch
    save(datafile, 'subj_id'); % init empty data file, will append with results later
end

% load example monitor phosphors and cone fundamentals from toolbox
% note on this circle: 0 deg = "red", 90 = "yellow", 180 = "green", 270 = "blue"
    % 0, 180 = maximum L-M contrast, no S-(L+M) contrast
        % aka S cone activation is equal to the "sum" of L+M
        % red (0): greatest L-M difference. But S = 0 so why is there no S-(L+M) contrast?
            % maybe because there is no M either?
        % green (180): L+M = S
    % 90, 270 = maximum S-(L+M) contrast, no L-M contrast
        % blue (270): S peak, no L or M activation (so L-M = 0)   pos max
        % yellow (90): no S activation, L=M (so L-M = 0)          neg max
load phosphors.mat;
load fundamentals_ss.mat;
bg_rgb = bg_lum * [1; 1; 1];
color1_azimuth = center_azimuth;
color2_azimuth = all_periphery_azimuths(color_idx);

% compute colors
bg_lms = rgb2lms(phosphors, fundamentals, bg_rgb);
color1_rgb = lms2rgb(phosphors, fundamentals, dkl2lms(bg_lms, dkl_sph2cart([dkl_radius, color1_azimuth, 0]))) .* 255;
color2_rgb = lms2rgb(phosphors, fundamentals, dkl2lms(bg_lms, dkl_sph2cart([dkl_radius, color2_azimuth, 0]))) .* 255;

% do minimum motion flicker
% just: take the single target RGB color, and increment or decrement all 3
    % RGB channels by 1 to move a luminance step.
color1_rgb = round(color1_rgb);
color2_rgb = round(color2_rgb);
start_colors = [color1_rgb, color2_rgb];

% open screen
KbName('UnifyKeyNames');
escape = KbName('ESCAPE');
left = KbName('LeftArrow');
right = KbName('RightArrow');
screenid = max(Screen('Screens'));
if useDebug
    PsychDebugWindowConfiguration;
else
    clear Screen
end
Screen('Preference', 'SkipSyncTests', 1);
win = PsychImaging('OpenWindow', screenid);
flipinterval = Screen('GetFlipInterval',win);
Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');  % Set up alpha-blending for smooth (anti-aliased) lines

% set parameters
flicker_rate_hz = 15;
params.real_IFI = flipinterval;
params.trial_length_max = 5*60; % 5 minutes as needed
period_inframes = (1/flicker_rate_hz) / params.real_IFI;
[alpha_reference, alpha_mod] = generateAlphaSinusoids(flicker_rate_hz, flicker_rate_hz, 1, params);
alpha_mod = [alpha_mod((1+period_inframes/2):end), alpha_mod(1:(period_inframes/2+1))]; % alternate the two flickers
alpha_reference = alpha_reference ./ 2 .* 255; % set between 0 and 255
alpha_mod = alpha_mod ./ 2 .* 255;
reference_color_idx = 1; % which color stays constant (the other gets luminance modulated)
Screen('FillRect', win, bg_rgb .* 255);
Screen('Flip', win');
circle_radius = 200;
[xmid, ymid] = WindowCenter(win);
circ_rec = CenterRectOnPoint([0; 0; circle_radius*2; circle_radius*2], xmid, ymid);

% display
reference_rgb = start_colors(:, reference_color_idx); % which stays constant. default = color1 (center)
mod_rgb = start_colors(:, mod(reference_color_idx, 2) + 1); % which gets modulated. default = periphery
orig_mod_rgb = mod_rgb;
for frame = 1:numel(alpha_reference)
    % check keyboard only every 10th frame to save computing power
    if mod(frame, 10) == 0
        [keyDown, keyTime, keyCode] = KbCheck(-1); % check for button press
        if keyDown
            if keyCode(escape)
                %exitNow = 1;
                break
            elseif keyCode(right)
                mod_rgb = mod_rgb + 1;
            elseif keyCode(left)
                mod_rgb = mod_rgb - 1;
            end
        end
    end
    % don't let rgb values go out of range
    mod_rgb(mod_rgb > 255) = 255;
    mod_rgb(mod_rgb < 0) = 0;
    Screen('FillOval', win, [[reference_rgb; alpha_reference(frame)], [mod_rgb; alpha_mod(frame)]], [circ_rec', circ_rec']);
    Screen('Flip', win);
end
Screen('CloseAll');

final_colors = struct;
final_colors.center = reference_rgb;
final_colors.periphery = mod_rgb;
if color_idx == 1
    colors1 = final_colors;
    save(datafile, 'colors1', '-append');
elseif color_idx == 2
    colors2 = final_colors;
    save(datafile, 'colors2', '-append');
elseif color_idx == 3
    colors3 = final_colors;
    save(datafile, 'colors3', '-append');
end