function UI_drawFixation(winMain, params, colorOrder)
% draws fixation target from Thaler et al., 2013, Vision Research
% 3rd input colorOrder should be 0 or 1: flips colors

whiteColor = [params.white params.white params.white];
blackColor = [params.black params.black params.black];
% make white a bit dimmer
whiteColor = whiteColor .* (200/255);


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

d1 = 0.6; % diameter of outer circle in degrees
d1_scalar = d1/2*params.ppd;
d2 = 0.2; % diameter of inner circle in degrees
d2_scalar = d2/2*params.ppd;

xmid = params.xmid;
ymid = params.ymid;

Screen('FillOval', winMain, colorOval, [xmid-d1_scalar; ymid-d1_scalar; xmid+d1_scalar; ymid+d1_scalar], d1*params.ppd);
Screen('DrawLines', winMain, [xmid-d1_scalar, xmid+d1_scalar; ymid, ymid], d2*params.ppd, colorCross, [0 0], 0, 1);
Screen('DrawLines', winMain, [xmid, xmid; ymid-d1_scalar, ymid+d1_scalar], d2*params.ppd, colorCross, [0 0], 0, 1);
Screen('FillOval', winMain, colorOval, [xmid-d2_scalar; ymid-d2_scalar; xmid+d2_scalar; ymid+d2_scalar], d2*params.ppd);