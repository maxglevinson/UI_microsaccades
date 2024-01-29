% analyze group data

subjects = [201 202 203 204 205 206 208 209 210 211 212 213 214 215 216 217 218 219 222 223 224];
sess = 'allsessions'; % use all sessions
ploton = 0; % don't plot when retrieving individual data

nsubjects = numel(subjects);
all_mean_rts = nan(3, 9, nsubjects); % rows: illusion trials, replay trials, catch (harder) trials
all_perc_clicked = nan(3, 9, nsubjects);
subj_means = nan(3, nsubjects); % mean across all trial types, to normalize individual differences (maybe)
subj_means_perc_clicked = nan(3, nsubjects);
all_illusion_data = [];
all_data = [];
all_mean_shiftstarts = [];
for s = 1:nsubjects
    subj = subjects(s);
    [mean_rts_by_block, percentage_clicked_by_block, reorder, median_dists, illusion_data, subj_data, mean_shiftstarts] = analyze_indiv_data(subj, sess, ploton);
    all_mean_rts(:, :, s) = mean_rts_by_block(:, reorder);
    all_perc_clicked(:, :, s) = percentage_clicked_by_block(:, reorder);
    all_median_dists(:, s) = median_dists;
    subj_means(:, s) = nanmean(all_mean_rts(:, :, s), 2);
    subj_means_perc_clicked(:, s) = nanmean(all_perc_clicked(:, :, s), 2);
    all_illusion_data = [all_illusion_data; illusion_data];
    all_data = [all_data; subj_data];
    all_mean_shiftstarts = [all_mean_shiftstarts; mean_shiftstarts];
end
full_data_means = mean(subj_means, 2);
full_data_stderrs = std(subj_means, 0, 2) ./ sqrt(nsubjects);
full_data_means_perc_clicked = mean(subj_means_perc_clicked, 2);
full_data_stds_perc_clicked = std(subj_means_perc_clicked, 0, 2);
% it's reordered to:
%   sizes   2 4 6 2 4 6 2 4 6
%   colors  1 1 1 2 2 2 3 3 3


% store all illusion data in .csv table
full_behavior = all_data;
full_behavior(:, end) = (full_behavior(:, end) / 10) - 27; % convert contrast to 1 (low), 2 (med), 3 (high)
T = array2table(full_behavior, 'VariableNames', {'subject', 'trialtype', 'block', 'ignore', 'trial', 'pressed', 'FT', 'medianfixdist', 'meanfixdist', 'maxfixdist', 'eccentricity', 'naignore', 'contrast'});
writetable(T, 'bhv_data/full_behavior.csv');

% normalize individual differences?
norm_all_mean_rts = nan(size(all_mean_rts));
for s = 1:nsubjects
    curr_scalars = full_data_means ./ subj_means(:, s);
    %norm_all_mean_rts(:, :, s) = all_mean_rts(:, :, s) .* curr_scalars;
    norm_all_mean_rts(:, :, s) = all_mean_rts(:, :, s); % don't normalize
end

% plot mean rts for each individual subject in their own figure
for subj = 1:nsubjects
    figure; hold on;
    bar(reshape(all_mean_rts(1, :, subj), 3, 3));
legend({'low contrast', 'medium contrast', 'high contrast'});
xticks(1:3); xticklabels({'small', 'medium', 'large'});
set(gca, 'FontSize', 20);
ylabel('RT (s)'); xlabel('center size');
title(['subject ', num2str(subjects(subj))]);
set(gcf, 'Position', [1 29 1920 943]);
end

% indiv subjs, avgd across sizes
for subj = 1:nsubjects
    figure; hold on;
    bar(nanmean(reshape(all_mean_rts(1, :, subj), 3, 3), 1));
    xticks(1:3); xticklabels({'low', 'medium', 'high'});
    xlabel('contrast'); ylabel('RT (s)');
    set(gca, 'FontSize', 20);
    title(['subject ', num2str(subjects(subj))]);
end
% these subjects don't show low < medium < high: 
% 202, 205, 211, 212, 213, 214, 217, 223

% indiv subjs, avgd across colors
for subj = 1:nsubjects
    figure; hold on;
    bar(nanmean(reshape(all_mean_rts(1, :, subj), 3, 3), 2));
    xticks(1:3); xticklabels({'small', 'medium', 'large'});
    xlabel('size'); ylabel('RT (s)');
    set(gca, 'FontSize', 20);
    title(['subject ', num2str(subjects(subj))]);
end
% these subjects don't show large < medium < small: 
% 204, 206, 208, 211, 217, 218, 219, 222, 223

% these subjects don't show expected contrast OR size effect:
% 211, 217, 223.


mean_mean_rts = nanmean(norm_all_mean_rts, 3);
stderr_mean_rts = nanstd(norm_all_mean_rts, 0, 3) ./ sqrt(nsubjects);
mean_perc_clicked = mean(all_perc_clicked, 3);
stderr_perc_clicked = std(all_perc_clicked, 0, 3) ./ sqrt(nsubjects);


% plot mean rts across subjects
plotcolors = {[1 0 0], [0 1 0], [0 0 1]}; % red = weakest, green = mid, blue = strongest
figure; hold on;
br = barwitherr(reshape(stderr_mean_rts(1, :), 3, 3), [1:3], reshape(mean_mean_rts(1, :), 3, 3));
legend({'low contrast', 'medium contrast', 'high contrast'});
xticks(1:3); xticklabels({'small', 'medium', 'large'});
set(gca, 'FontSize', 20);
ylabel('RT (s)'); xlabel('center size');
for b = 1:3   
    br(b).FaceColor = 'flat';
    br(b).CData(1:3, :) = repmat(plotcolors{b}, 3, 1);
end

figure;
illusion_rts = squeeze(permute(norm_all_mean_rts(1, :, :), [3 2 1]));
grouping = categorical([[2 1]; [4 1]; [6 1]; [2 2]; [4 2]; [6 2]; [2 3]; [4 3]; [6 3]]); % first col: contrast, 2nd col: size
boxplot(illusion_rts, grouping, 'factorgap', [10 1], 'factorseparator', [1], 'colorgroup', [1 1 1 2 2 2 3 3 3], 'colors', [);

% mean rts switched order
figure; hold on;
barwitherr(reshape(stderr_mean_rts(1, :), 3, 3)', [1:3], reshape(mean_mean_rts(1, :), 3, 3)');
legend({'small size', 'medium size', 'large size'});
xticks(1:3); xticklabels({'low', 'medium', 'high'});
set(gca, 'FontSize', 20);
ylabel('RT (s)'); xlabel('color contrast');


% plot perc clicked across subjects
figure; hold on;
barwitherr(reshape(stderr_perc_clicked(1, :), 3, 3), [1:3], reshape(mean_perc_clicked(1, :), 3, 3));
legend({'low contrast', 'medium contrast', 'high contrast'});
xticks(1:3); xticklabels({'small', 'medium', 'large'});
set(gca, 'FontSize', 20);
ylabel('percentage of trials filled-in'); xlabel('center size');

% plot 3 trial types for each subject, averaged across blocks
subj_mean_rts = squeeze(nanmean(all_mean_rts, 2));
%subj_stderr_rts = nanstd(all_mean_rts, 0, 2) ./ sqrt(nsubjects);
figure;
bar(subj_mean_rts');
legend({'illusion', 'replay', 'catch'});
set(gca, 'FontSize', 20);
xlabel('subject');
ylabel('RT (s)');

% plot 3 trial types for each block, averaged across subjects
figure;
barwitherr(stderr_mean_rts', [1:9], mean_mean_rts');
legend({'illusion', 'replay', 'catch'});
xticks(1:9);
set(gca, 'FontSize', 20);
ylabel('RT (s)'); xlabel('block');

% plot 3 trial types averaged across both blocks and subjects.
figure;
bar(mean(subj_mean_rts, 2)');
xticks(1:3); xticklabels({'illusion', 'replay', 'catch'});
set(gca, 'FontSize', 20);
ylabel('RT (s)');