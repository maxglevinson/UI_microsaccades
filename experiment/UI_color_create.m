function [fillColors, frameColors] = UI_color_create(colorCenter, colorPeriphery, params)
% generates rgb colors for stimulus
% based on DKL azimuth (colorCenter and colorPeriphery)
    % and DKL radius (params.dkl_radius)
% requires monitor phosphors in params.phosphors (here we just use an example)
% and cone fundamentals in params.cone_fundamentals (we use Stockman & Sharpes)

% first check if we instead just manually set the RGB values
if isfield(params, 'centerRGBforced') % we set centerRGBforced and peripheryRGBforced in this case
    % if so, ignore colorCenter and colorPeriphery arguments
    if colorCenter == params.colorCenter % a normal block
        center_rgb = params.centerRGBforced;
        periphery_idx = find(colorPeriphery == params.colorPeriphery);
        periphery_rgb = params.peripheryRGBforced(:, periphery_idx);
    else % all uniform block; take average color
        center_rgb = params.centerRGBforced;
        periphery_rgb = params.peripheryRGBforced(:, 2);
        center_rgb = mean([center_rgb, periphery_rgb], 2); % average them
        periphery_rgb = center_rgb; % set center and periphery equal
    end
else % we didn't set them, have to calculate them from DKL values
    center_rgb = UI_dkl2rgb(params, colorCenter);
    periphery_rgb = UI_dkl2rgb(params, colorPeriphery);
end

fillColors = [periphery_rgb, center_rgb]; % PERIPHERY FIRST!
% % gamma correct to see what it looks like
% for i = 1:numel(fillColors)
%     fillColors(i) = ((fillColors(i) / 255) ^ (1/2.2)) * 255;
% end

%% make the soft border
% scale from current center & periphery colors
nframeRects = 10;
frameColors = zeros(3, nframeRects);
for rgb = 1:3
    curr_scalars = linspace(fillColors(rgb, 2), fillColors(rgb, 1), nframeRects+2);
    frameColors(rgb, :) = curr_scalars(:, 2:end-1);
end