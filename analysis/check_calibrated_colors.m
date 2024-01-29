% load all calibrated colors for each subject and get summary statistics.

subjects = string([201:206 208:219, 222:224]);
nsubjects = numel(subjects);

data_dir = '/Users/maxlevinson/Documents/McGill/OneDrive - McGill University/neurospeed/UI_eyelink/bhv_data/';

peripheries_low = nan(3, nsubjects);
peripheries_med = nan(3, nsubjects);
peripheries_high = nan(3, nsubjects);

for s = 1:nsubjects
    filename = [data_dir, 'rgb_colors_', subjects{s}, '.mat'];
    load(filename);
    peripheries_low(:, s) = colors1.periphery;
    peripheries_med(:, s) = colors2.periphery;
    peripheries_high(:, s) = colors3.periphery;
end

median_colors = [median(peripheries_low, 2), median(peripheries_med, 2), median(peripheries_high, 2)];

% convert subject colors to offset from median (assuming median = the
% default start color
offsets_low = peripheries_low(1, :) - median_colors(1, 1);
offsets_med = peripheries_med(1, :) - median_colors(1, 2);
offsets_high = peripheries_high(1, :) - median_colors(1, 3);

% median offsets, mean ± standard deviation
median_offsets = [median(offsets_low), median(offsets_med), median(offsets_high)];
mean_offsets = [mean(offsets_low), mean(offsets_med), mean(offsets_high)];
std_offsets = [std(offsets_low), std(offsets_med), std(offsets_high)];
% mean ± std of those that aren't default; absolute-valued
mean_offsets_not0 = [mean(abs(offsets_low(offsets_low ~= 0))), mean(abs(offsets_med(offsets_med ~= 0))), mean(abs(offsets_high(offsets_high ~= 0)))];
std_offsets_not0 = [std(abs(offsets_low(offsets_low ~= 0))), std(abs(offsets_med(offsets_med ~= 0))), std(abs(offsets_high(offsets_high ~= 0)))];

% percentage of subjects with 0 offsets
perc_default_low = sum(offsets_low == 0) / nsubjects * 100;
perc_default_med = sum(offsets_med == 0) / nsubjects * 100;
perc_default_high = sum(offsets_high == 0) / nsubjects * 100;