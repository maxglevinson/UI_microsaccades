function target_rgb = UI_dkl2rgb(params, azm)
% computes rgb values for an input dkl azimuth, given other parameters set in params
% target and background are isoluminant

% set background rgb and lms
bg_rgb = params.colBackground' ./ params.white; % grey DKL 'background' to generate contrast from
bg_lms = rgb2lms(params.phosphors, params.cone_fundamentals, bg_rgb);

% set target dkl
target_dkl = [params.dkl_radius, azm, 0]; % zero elevation = isoluminant to grey background

% convert to rgb (0-255 normally)
target_rgb = lms2rgb(params.phosphors, params.cone_fundamentals, dkl2lms(bg_lms, dkl_sph2cart(target_dkl))) .* params.white;