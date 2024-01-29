function Nb = retinal_slip(xPos, yPos, ms_present, incl_ms)
% retinal slip calculation from Engbert & Mergenthaler 2006
% creates boxes of length 0.01 dva (epsilon, ~diameter of a cone receptive field)
% For every 50ms time window, number of boxes needed to fully cover the
% trajectory (Nb) is used as a quantitative measure of retinal image slip.
%
% assumes xPos and yPos input are in dva, mean-centered, averaged across both eyes.
% ignore any timepoint containing a microsaccade in-flight (unless incl_ms = 1)

if nargin < 4
    incl_ms = 0; % default: exclude microsaccades.
end

epsilon = 0.01; % epsilon in dva
winsize = 50; % time window length in ms

% 5 x 5 dva grid of box centers (points) separated by 0.01 dva
dva_side = 20;
nboxes_side = dva_side ./ epsilon;
boxgrid = zeros(nboxes_side, nboxes_side);
for b = 1:numel(boxgrid)
    boxgrid(b) = b;
end

% 2 x N array listing the edges in dva of each box.
% first row: left/bottom, second row: right/top.
boxsides = [0:epsilon:dva_side-epsilon; epsilon:epsilon:dva_side] - dva_side/2;


% set microsaccades to NaN, and 10ms surrounding each one
if ~incl_ms
ms_changes = [0; diff(ms_present)];
ms_onsets = find(ms_changes == 1);
ms_offsets = find(ms_changes == -1);
xPos(ms_present == 1) = NaN;
for iOns = 1:numel(ms_onsets)
    ons = ms_onsets(iOns);
    if ons > 10
        xPos((ons-10):ons) = NaN;
    elseif ons <= 10
        xPos(1:ons) = NaN;
    end
end
for iOfs = 1:numel(ms_offsets)
    ofs = ms_offsets(iOfs);
    if ofs < numel(xPos) - 10
        xPos(ofs:(ofs+10)) = NaN;
    else
        xPos(ofs:end) = NaN;
    end
end
end

% for each time window: for each timepoint, get the current box. Add
% up total number of boxes traversed.
nwindows = floor(numel(xPos) / winsize);
Nb = nan(1, nwindows);
for win = 1:nwindows
    curr_window = [1:winsize] + winsize*(win-1); % indices into trial_data (timepoints)
    %if sum(isnan(xPos(curr_window))) < winsize/4 % if less than a quarter of timepoints are nans, use the window
    if sum(isnan(xPos(curr_window))) == 0 % if any nans, skip the window
        boxes = nan(1, numel(curr_window)); % add each traversed box number to this vector
        i=0; % increment boxes
        for t = curr_window
            i = i+1;
            if ~isnan(ms_present(t)) && (incl_ms || ~(ms_present(t))) % if usable tp && no microsaccade currently ongoing (unless incl_ms=1)
                curr_box_x = xPos(t) >= boxsides(1, :) & xPos(t) < boxsides(2, :);
                curr_box_y = yPos(t) >= boxsides(1, :) & yPos(t) < boxsides(2, :);
                if sum(curr_box_x) && sum(curr_box_y) % eye position is within the full grid
                    boxes(i) = boxgrid(curr_box_x, curr_box_y); % current box number
                end % otherwise just leave it at nan, if eye position is wild (outside the grid)
            end
        end
        Nb(win) = numel(~isnan(unique(boxes))); % how many unique boxes traversed.
    else
        Nb(win) = nan;
    end
end