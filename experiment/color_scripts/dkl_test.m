% DKL TEST

%% first attempts - not used
% % load example monitor phosphors and cone fundamentals from toolbox
% load phosphors.mat;
% load fundamentals_ss.mat;

% % ORIGINAL JUST TESTING STUFF
% % set desired colors
% % dkl values come from Shim paper
% % or, dkl values come from textbook circle (so different from Shim's
% % circle? why?)
% % the colors aren't quite what I expect them to be, but relatively close.
%     % see plotDKLspace
% bg_rgb = [0.5; 0.5; 0.5]; % mid grey
% %dkl_radius = 0.8; % constant used in Shim paper, saturation
% dkl_radius = 0.3; % use lower bc output rgb are too high
% dkl_elevation = 0; % isoluminant from the grey background
% blue_dkl = [dkl_radius; 270; dkl_elevation];
% magenta_dkl = [dkl_radius; 315; dkl_elevation];
% yellow_dkl = [dkl_radius; 90; dkl_elevation];
% cyan_dkl = [dkl_radius; 225; dkl_elevation];
% bluegreen_dkl = [dkl_radius; 260; dkl_elevation];
% 
% % compute colors
% bg_lms = rgb2lms(phosphors, fundamentals, bg_rgb); % background from rgb to lms
% % convert dkl from spherical to cartesian, then to lms
% blue_lms = dkl2lms(bg_lms, dkl_sph2cart(blue_dkl));
% magenta_lms = dkl2lms(bg_lms, dkl_sph2cart(magenta_dkl));
% yellow_lms = dkl2lms(bg_lms, dkl_sph2cart(yellow_dkl));
% cyan_lms = dkl2lms(bg_lms, dkl_sph2cart(cyan_dkl));
% bluegreen_lms = dkl2lms(bg_lms, dkl_sph2cart(bluegreen_dkl));
% % convert targets from lms to rgb
% blue_rgb = lms2rgb(phosphors, fundamentals, blue_lms) .* 255;
% magenta_rgb = lms2rgb(phosphors, fundamentals, magenta_lms) .* 255;
% yellow_rgb = lms2rgb(phosphors, fundamentals, yellow_lms) .* 255;
% cyan_rgb = lms2rgb(phosphors, fundamentals, cyan_lms) .* 255;
% bluegreen_rgb = lms2rgb(phosphors, fundamentals, bluegreen_lms) .* 255;


% %% compute colors surrounding a target color
% load phosphors.mat;
% load fundamentals_ss.mat;
% dkl_radius = 0.3;
% dkl_elevation = 0; % isoluminant
% % set desired colors
% bg_rgb = [0.25; 0.25; 0.25];
% basecolor_azimuth = 260; % blue
% small_diff = 5; % small difference in azimuth
% large_diff = 10;
% % color matrices: 3x5
% % rows: radius, azimuth, elevation
% % columns: base color, small plus, small minus, large plus, large minus
% all_dkl = repmat([dkl_radius; nan; dkl_elevation], 1, 5);
% all_dkl(2, :) = basecolor_azimuth + [0, small_diff, -small_diff, large_diff, -large_diff];
% 
% % compute colors
% bg_lms = rgb2lms(phosphors, fundamentals, bg_rgb);
% all_lms = zeros(size(all_dkl));
% all_rgb = zeros(size(all_dkl));
% for target = 1:size(all_dkl, 2)
%     all_lms(:, target) = dkl2lms(bg_lms, dkl_sph2cart(all_dkl(:, target)));
%     all_rgb(:, target) = lms2rgb(phosphors, fundamentals, all_lms(:, target)) .* 255;
% end
% 
% % plot them to just take a look
% figure; hold on;
% for target = 1:size(all_rgb, 2)
%     plot(target, 1, 'o', 'Color', all_rgb(:, target) ./ 255, 'MarkerFace', all_rgb(:, target) ./ 255, 'MarkerSize', 30)
% end
% xlim([0 target+1]);

%% Two target colors: generate RGB for flicker fusion / minimum motion isoluminance
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
bg_rgb = .4 * [1; 1; 1];
color1_azimuth = 0;%60;
color2_azimuth = 45;%150;
dkl_radius1 = 0.05;%0.07;
dkl_radius2 = 0.05;%dkl_radius1; % use different saturation?
% bg_rgb .4, radius 0.09, azimuths 60 & 110 replicate p well the
% uniformillusion.com red/green. But the red is too striking for filling?
% but maybe need 0.05 to facilitate filling.

% compute colors
bg_lms = rgb2lms(phosphors, fundamentals, bg_rgb);
color1_rgb = lms2rgb(phosphors, fundamentals, dkl2lms(bg_lms, dkl_sph2cart([dkl_radius1, color1_azimuth, 0]))) .* 255;
color2_rgb = lms2rgb(phosphors, fundamentals, dkl2lms(bg_lms, dkl_sph2cart([dkl_radius2, color2_azimuth, 0]))) .* 255;

% plot both colors at larger saturation
color1_fullsat = lms2rgb(phosphors, fundamentals, dkl2lms(bg_lms, dkl_sph2cart([0.15, color1_azimuth, 0]))) .* 255;
color2_fullsat = lms2rgb(phosphors, fundamentals, dkl2lms(bg_lms, dkl_sph2cart([0.15, color2_azimuth, 0]))) .* 255;
which_plot = 'fullsat'; % fullsat or rgb (the one we'll show)
switch which_plot
    case 'fullsat'
        color1_plot = color1_fullsat;
        color2_plot = color2_fullsat;
    case 'rgb'
        color1_plot = color1_rgb;
        color2_plot = color2_rgb;
end
figure; hold on;
plot(1, 1, 'o', 'Color', color1_plot ./ 255, 'MarkerFace', color1_plot ./ 255, 'MarkerSize', 30);
plot(2, 1, 'o', 'Color', color2_plot ./ 255, 'MarkerFace', color2_plot ./ 255, 'MarkerSize', 30);
xlim([0 3]);

% 
% % compute dkl given a rgb (TESTING/TROUBLESHOOTING)
% t_rgb = [127,100,99] ./ 255;
% t_lms = rgb2lms(phosphors, fundamentals,t_rgb);
% t_dkl_rad = lms2dkl(bg_lms, t_lms);
% t_dkl_deg = dkl_cart2sph(t_dkl_rad)

% do minimum motion flicker
% just: take the single target RGB color, and increment or decrement all 3
% RGB channels by 1 to move a luminance step.
color1_rgb = round(color1_rgb);
color2_rgb = round(color2_rgb);
start_colors = [color1_rgb, color2_rgb];
%start_colors = [[200; 0; 0], [0; 200; 0]]; 

% open screen
KbName('UnifyKeyNames');
escape = KbName('ESCAPE');
left = KbName('LeftArrow');
right = KbName('RightArrow');
screenid = max(Screen('Screens'));
%PsychDebugWindowConfiguration;
clear Screen
win = PsychImaging('OpenWindow', screenid);
flipinterval = Screen('GetFlipInterval',win);
Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');  % Set up alpha-blending for smooth (anti-aliased) lines

% set parameters
flicker_rate_hz = 15;
params.real_IFI = flipinterval; params.trial_length_max = 60;
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
reference_rgb = start_colors(:, reference_color_idx); % which stays constant
mod_rgb = start_colors(:, mod(reference_color_idx, 2) + 1); % which gets modulated
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
