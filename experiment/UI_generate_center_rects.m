function params = UI_generate_center_rects(params)

%% CREATE PERIPHERAL/CENTRAL FIELDS
% adapted from Suarez-Pinilla et al., 2018, i-Perception
% "We imitate Yair's peripheral overflow examples and create a rectangular
% area corresponding to the 'central' visual field."

Limit_central=params.centerdva; % degrees from center to top or bottom of center rect
Half_ht=tand(Limit_central)*params.cm_to_screen*params.pxbycm_y;
Half_wt = Half_ht; % same visual angle for horizontal and vertical
params.Central_rect=[params.xmid-Half_wt params.ymid-Half_ht params.xmid+Half_wt params.ymid+Half_ht];