function UI_drawFixation(winMain, params, colorOrder)
% draws fixation target from Thaler et al., 2013, Vision Research
% 3rd input colorOrder should be 0 or 1: flips colors

whiteColor = [params.white params.white params.white];
blackColor = [params.black params.black params.black];
if ispc & contains(params.scripts_dir, 'D:') % if Windows machine (MEG stimulus PC)
    % make white dimmer because propixx is linear, no gamma applied
    whiteColor = whiteColor .* (150/255);
else % otherwise do make it a bit dimmer but not by as much
    whiteColor = whiteColor .* (200/255);
end

% default settings
if nargin < 3 % no colorOrder specified
    colorOrder = 0;
end
if colorOrder == 0
    colorOval = whiteColor; % black circles
    colorCross = blackColor; % white cross
elseif colorOrder == 1
    colorOval = blackColor;
    colorCross = whiteColor;
end

% particular conditions to change fixation colors
switch params.stimtype
    case 'dotblur'
        if colorOrder == 0
            colorOval = blackColor; % black circles
            colorCross = whiteColor; % white cross
        elseif colorOrder == 1
            colorOval = whiteColor;
            colorCross = blackColor;
        end
    case 'retinotopy'
        colorOval = blackColor;
        if colorOrder == 0
            colorCross = [200 0 0]; % red
            if params.rapidTag && params.real_FR == 1440
                colorCross = [100 100 100]; % grey
            end
        elseif colorOrder == 1
            colorCross = [0 0 200]; % blue
            if params.rapidTag && params.real_FR == 1440
                colorCross = [10 10 10]; % black
            end
        end
end

d1 = 0.6; % diameter of outer circle in degrees
d1_scalar = d1/2*params.ppd;
d2 = 0.2; % diameter of inner circle in degrees
d2_scalar = d2/2*params.ppd;

if params.rapidTag % need to draw fixation in all 4 quadrants
    xmid = zeros(1, 4);
    ymid = zeros(1, 4);
    xmid_lines = zeros(1, 8);
    ymid_lines = zeros(1, 8);
    for quad = 1:4
        [curr_x, curr_y] = RectCenter(params.quadrants.rect(quad, :));
        xmid(quad) = curr_x;
        ymid(quad) = curr_y;
        xmid_lines(quad*2-1:quad*2) = curr_x;
        ymid_lines(quad*2-1:quad*2) = curr_y;
    end
else
    xmid = params.xmid;
    ymid = params.ymid;
end
Screen('FillOval', winMain, colorOval, [xmid-d1_scalar; ymid-d1_scalar; xmid+d1_scalar; ymid+d1_scalar], d1*params.ppd);
if params.rapidTag
    Screen('DrawLines', winMain, [xmid_lines + d1_scalar.*[-1 1 -1 1 -1 1 -1 1]; ymid_lines], d2*params.ppd, colorCross, [0 0]);%[1920/2 1080/2]);
    Screen('DrawLines', winMain, [xmid_lines; ymid_lines + d1_scalar.*[-1 1 -1 1 -1 1 -1 1]], d2*params.ppd, colorCross, [0 0]);%[1920/2 1080/2]);
else
    %Screen('DrawLine', winMain, colorCross, xmid-d1_scalar, ymid, xmid+d1_scalar, ymid, d2*params.ppd);
   % Screen('DrawLine', winMain, colorCross, xmid, ymid-d1_scalar, xmid, ymid+d1_scalar, d2*params.ppd);
    Screen('DrawLines', winMain, [xmid-d1_scalar, xmid+d1_scalar; ymid, ymid], d2*params.ppd, colorCross, [0 0], 0, 1);
    Screen('DrawLines', winMain, [xmid, xmid; ymid-d1_scalar, ymid+d1_scalar], d2*params.ppd, colorCross, [0 0], 0, 1);
end
Screen('FillOval', winMain, colorOval, [xmid-d2_scalar; ymid-d2_scalar; xmid+d2_scalar; ymid+d2_scalar], d2*params.ppd);