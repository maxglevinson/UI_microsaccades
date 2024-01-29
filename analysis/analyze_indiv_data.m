function [mean_rts_by_block, percentage_clicked_by_block, reorder, median_dists, illusion_data, subj_data, mean_shiftstarts] = analyze_indiv_data(subj, sess, ploton)
% code to analyze individual subject main task data.
% reorder changes the order of blocks to match:
%   sizes   2 4 6 2 4 6 2 4 6
%   colors  1 1 1 2 2 2 3 3 3

scripts_dir = '/Users/maxlevinson/Documents/McGill/OneDrive - McGill University/neurospeed/UI_eyelink';
%scripts_dir = '/export04/data/mlevin/UI_eyelink';
data_dir = [scripts_dir, '/bhv_data'];

%subj = '207'; % enter subject code here
subj = num2str(subj);
%sess = 'allsessions'; % enter session number here, or 'allsessions'
sess = num2str(sess);
datafile = [data_dir, '/color_', subj, '_', sess, '_data.mat'];
load(datafile);
%ploton = 0; % 1 to plot, 0 to not

blocks_to_use = 1:max(data(:, 3)); % default = all
blocks = data(:, 3);
use_idx = ismember(blocks, blocks_to_use);
blocks = blocks(use_idx);
rts = data(use_idx, 7);
% generate rolling average of rt density across trials, across time
rts_ms = rts * 1000;
window_step = 250; % ms
window_length = 500; % ms
nwindows = 15000 / window_step; % 15000 = maximum trial length
rt_density = zeros(1, nwindows);
rt_time = zeros(1, nwindows);
for win = 1:nwindows
    curr_window_start = (win-1)*window_step;
    curr_window = [curr_window_start, curr_window_start + window_length] ./ 1000;
    rt_time(win) = mean(curr_window);
    rts_in_window = rts > curr_window(1) & rts < curr_window(2);
    rt_density(win) = sum(rts_in_window);
end
if ploton; figure; plot(rt_time, rt_density); end
    
subj_data = data(use_idx, :); % all data to analyze
trialtypes = data(use_idx, 2);
clicked = data(use_idx, 6);
blockcolor = data(use_idx, 4);
blocksize = data(use_idx, 11);

% get data, split by trial type (normal / controls)
mean_rts_by_block = zeros(3, numel(unique(blocks)));
std_rts_by_block = zeros(3, numel(unique(blocks)));
percentage_clicked_by_block = zeros(3, numel(unique(blocks)));
illusion_data = zeros(sum(trialtypes==1), size(data, 2));%3); % illusion rts, contrasts, center sizes
mean_shiftstarts = nan(1, numel(unique(blocks)));
for block = unique(blocks)'
    for tt = 1:3
        curr_rts = rts(trialtypes == tt & blocks == block);
        mean_rts_by_block(tt, block) = nanmean(curr_rts(curr_rts > 1)); % only take rt > 1 s...
        std_rts_by_block(tt, block) = nanstd(curr_rts(curr_rts > 1)); % only take rt > 1 s...
        percentage_clicked_by_block(tt, block) = sum(~isnan(curr_rts)) / numel(curr_rts);
        if tt == 1 % get illusion trial data matrix
            illusion_data((1:40)+40*(block-1), 1) = curr_rts;
            illusion_data((1:40)+40*(block-1), 2) = blockcolor(trialtypes == tt & blocks == block);
            illusion_data((1:40)+40*(block-1), 3) = blocksize(trialtypes == tt & blocks == block);
        end
        if tt == 2 % get shiftstarts
            curr_shiftstarts = timing(trialtypes == tt & blocks == block, 3) - timing(trialtypes == tt & blocks == block, 2);
            mean_shiftstarts(block) = nanmean(curr_shiftstarts);
        end
    end
end

% get mean distances from fixation for illusion trials
median_dists = data(use_idx, 8);
median_dists = median_dists(trialtypes == 1);

% split by color
mean_rts_by_contrast = zeros(3, 2);
std_rts_by_contrast = zeros(3, 2);
percentage_clicked_by_contrast = zeros(3, 2);
for contrast = 1:3 % columns, color contrast (block types 1, 2, or 3)
    for tt = 1:3 % rows, trial types
        curr_rts = rts(trialtypes == tt & blockcolor == contrast);
        mean_rts_by_contrast(tt, contrast) = nanmean(curr_rts);
        std_rts_by_contrast(tt, contrast) = nanstd(curr_rts);
        percentage_clicked_by_contrast(tt, contrast) = sum(~isnan(curr_rts)) / numel(curr_rts);
    end
end

% split by center size
mean_rts_by_size = zeros(3, 3);
std_rts_by_size = zeros(3, 3);
percentage_clicked_by_size = zeros(3, 3);
sizelist = [2 4 6];
for csize = 1:3 % columns, center size (2, 4, or 6)
    curr_size = sizelist(csize);
    for tt = 1:3 % rows, trial type
        curr_rts = rts(trialtypes == tt & blocksize == curr_size);
        mean_rts_by_size(tt, csize) = nanmean(curr_rts);
        std_rts_by_size(tt, csize) = nanstd(curr_rts);
        percentage_clicked_by_size(tt, csize) = sum(~isnan(curr_rts)) / numel(curr_rts);
    end
end


% plot data by color
if ploton
figure; hold on;
barwitherr(std_rts_by_contrast, [1:3], mean_rts_by_contrast);
xticklabels({'illusion', '', 'replay', '', 'catch'});
legend({'low contrast', 'medium contrast', 'high contrast'});
ylabel('reaction time (s)');
h = gca;
h.XAxis.TickLength = [0 0];
set(gca, 'FontSize', 20);
legend box off
box off
end

% plot data by size
if ploton
figure; hold on;
barwitherr(std_rts_by_size, [1:3], mean_rts_by_size);
xticklabels({'illusion', '', 'replay', '', 'catch'});
legend({'small', 'medium', 'big'});
ylabel('reaction time (s)');
h = gca;
h.XAxis.TickLength = [0 0];
set(gca, 'FontSize', 20);
legend box off
box off
end

% plot illusion trials only, by color / size
% reorder data
colororder = blockcolor(1:50:401);
sizeorder = blocksize(1:50:401);
colorordertarget = [1 1 1 2 2 2 3 3 3];
sizeordertarget = [2 4 6 2 4 6 2 4 6];
nblocks = numel(colororder);
reorder = zeros(1, nblocks);
for b = 1:nblocks
    reorder(b) = find(colororder == colorordertarget(b) & sizeorder == sizeordertarget(b));
end

illusion_mean_rts_by_block = mean_rts_by_block(1, reorder);
illusion_mean_rts_by_block = reshape(illusion_mean_rts_by_block, 3, 3);
illusion_std_rts_by_block = std_rts_by_block(1, reorder);
illusion_std_rts_by_block = reshape(illusion_std_rts_by_block, 3, 3);
if ploton
figure; hold on;
barwitherr(illusion_std_rts_by_block, [1:3], illusion_mean_rts_by_block);
legend({'low contrast', 'medium contrast', 'high contrast'});
xticks(1:3); xticklabels({'small', 'medium', 'large'});
set(gca, 'FontSize', 20);
ylabel('RT (s)'); xlabel('center size');
end

%%

% % consolidate behavior to quickly look at (during brainstorm analyses)
% blocks = 1:6;
% databyblock = nan(30, 13, 6);
% for b = blocks
%     databyblock(:,:,b) = data(data(:, 3) == b, :);
% end
% 
% trialtypebyblock = squeeze(databyblock(:, 2, :));
% curr_b = 6;
% illusions = find(trialtypebyblock(:, curr_b) == 1);
% controls = find(trialtypebyblock(:, curr_b) == 2);
% catches = find(trialtypebyblock(:, curr_b) == 3);