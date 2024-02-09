% analyze group data

subjects = [201 202 203 204 205 206 208 209 210 211 212 213 214 215 216 217 218 219 222 223 224];
sess = 'allsessions'; % use allsessions, or input 1, 2, or 3

nsubjects = numel(subjects);
all_data = [];
sess = num2str(sess); % in case we input a session number
for s = 1:nsubjects
    subj = num2str(subjects(s));
datafile = [data_dir, '/color_', subj, '_', sess, '_data.mat'];
load(datafile);
blocks_to_use = 1:max(data(:, 3)); % default = all
blocks = data(:, 3);
use_idx = ismember(blocks, blocks_to_use);
subj_data = data(use_idx, :); % all data to analyze
    all_data = [all_data; subj_data];
end

% store all illusion data in .csv table
full_behavior = all_data;
full_behavior(:, end) = (full_behavior(:, end) / 10) - 27; % convert contrast to 1 (low), 2 (med), 3 (high)
T = array2table(full_behavior, 'VariableNames', {'subject', 'trialtype', 'block', 'ignore', 'trial', 'pressed', 'FT', 'medianfixdist', 'meanfixdist', 'maxfixdist', 'eccentricity', 'naignore', 'contrast'});
writetable(T, 'bhv_data/full_behavior.csv');