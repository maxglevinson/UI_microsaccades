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
fade_area_deg = params.centerdva * 2 * 0.01 * 6; % same as online study (circle radius = 4); dva of total fading area
thickness = round((tand(fade_area_deg)*params.cm_to_screen*params.pxbycm_y) / nframeRects); % in pixels
frameRects = NaN(4, nframeRects);
c = params.Central_rect;
for f = 1:nframeRects
    frameRects(1, f) = c(1) - thickness*f;
    frameRects(2, f) = c(2) - thickness*f;
    frameRects(3, f) = c(3) + thickness*f;
    frameRects(4, f) = c(4) + thickness*f;
end

%% if rapid frequency tagging, add more quadrants
if params.rapidTag
    fillRects = squeeze([params.quadrants.rect(1, :)', params.quadrants.rect(2, :)', params.quadrants.rect(3, :)', params.quadrants.rect(4, :)']);
    for quad = 1:4
        fillRects = [fillRects, params.quadrants.center(quad, :)']; % peripheries first, then centers
        if quad ~= 1
        for f = 1:nframeRects
            frameRects = [frameRects, CenterRect(frameRects(:, f)', params.quadrants.rect(quad, :))'];
        end
        end
    end
end

%% IMPORTANT NOTE: MIGHT NEED TO TRY LOWERING # OF FRAMED RECTS IF DISPLAY CANNOT FLIP IN TIME %%