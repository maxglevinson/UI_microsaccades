function [fillRects, frameRects, thickness] = UI_color_init(params)
% generates 2 rectangles: full display and center
% generates X framed rects, corresponding to fade from center -> periphery
% output order of framed rects: center -> periphery
nframeRects = 10;

fillRects = NaN(4, 2);
fillRects(:, 1) = params.totalRect; % periphery
fillRects(:, 2) = params.Central_rect; % center

% create larger 'central' areas similar to how we created original center
% these will be the framed rectangles
fade_area_deg = params.centerdva * 2 * 0.01 * 6; % dva of total fading area; comes out to 12% the radius
thickness = round((tand(fade_area_deg)*params.cm_to_screen*params.pxbycm_y) / nframeRects); % in pixels
frameRects = NaN(4, nframeRects);
c = params.Central_rect;
for f = 1:nframeRects
    frameRects(1, f) = c(1) - thickness*f;
    frameRects(2, f) = c(2) - thickness*f;
    frameRects(3, f) = c(3) + thickness*f;
    frameRects(4, f) = c(4) + thickness*f;
end
