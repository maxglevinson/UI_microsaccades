function params = UI_generate_center_rects(params)

%% CREATE PERIPHERAL/CENTRAL FIELDS
% We imitate Yair's peripheral overflow examples and create a rectangular
% area corresponding to the 'central' visual field.

Limit_central=params.centerdva; % degrees from center to top or bottom of center rect
Half_ht=tand(Limit_central)*params.cm_to_screen*params.pxbycm_y;
%Half_wt=Half_ht*params.DisplaySizeX/params.DisplaySizeY;
Half_wt = Half_ht; % same visual angle for horizontal and vertical
params.Central_rect=[params.xmid-Half_wt params.ymid-Half_ht params.xmid+Half_wt params.ymid+Half_ht];

%% create quadrants if using rapid frequency tagging
if params.rapidTag
    fullX = params.DisplaySizeX * 2;
    fullY = params.DisplaySizeY * 2;
    params.quadrants = struct;
    params.quadrants.rect = nan(4, 4); % quadrants x corners
    params.quadrants.center = nan(4, 4); % quadrants x corners
    params.quadrants.fourrects = nan(4, 4, 4); % corners x rects x quadrants
    params.quadrants.rect(1, :) = [0 0 fullX/2 fullY/2]; % top left total rect
    params.quadrants.center(1, :) = CenterRect(params.Central_rect, params.quadrants.rect(1, :)); % top left central rect
    params.quadrants.rect(2, :) = [fullX/2 0 fullX fullY/2]; % top right
    params.quadrants.center(2, :) = CenterRect(params.Central_rect, params.quadrants.rect(2, :));
    params.quadrants.rect(3, :) = [0 fullY/2 fullX/2 fullY]; % bottom left
    params.quadrants.center(3, :) = CenterRect(params.Central_rect, params.quadrants.rect(3, :));
    params.quadrants.rect(4, :) = [fullX/2 fullY/2 fullX fullY]; % bottom right
    params.quadrants.center(4, :) = CenterRect(params.Central_rect, params.quadrants.rect(4, :));
    params.quadshifts = [[0; 0; 0; 0], [1; 0; 1; 0], [0; 1; 0; 1], [1; 1; 1; 1]];
    for quad = 1:4
        params.quadshifts(:, quad) = params.quadshifts(:, quad) .* [params.DisplaySizeX; params.DisplaySizeY; params.DisplaySizeX; params.DisplaySizeY];
    end
end

