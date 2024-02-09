% merge several behavioral data files together.

subjects = string([201:206 208:219, 222:224]);
nsubjects = numel(subjects);
scripts_dir = '/export04/data/mlevin/Uniformity_Illusion';
data_dir = [scripts_dir, '/bhv_data'];

for s = 1:nsubjects
    subj = subjects{s};
% subjects 210 and 223 had different session schedules
if strcmp(subj, '210')
    blocks_to_use = {[1:3], [1:3], [1], [1:4]}; % subject 210
elseif strcmp(subj, '223')
    blocks_to_use = {[1:2], [1], [1:3], [1], [1:2]}; % subject 223
else
    blocks_to_use = {[1:3], [1:3], [1:3]}; % default: 3 blocks for 3 sessions
end
rawdata = {};
rawblocks = {};
use_idx = {};
nsessions = numel(blocks_to_use);
data = [];
timing = [];
last_block_idx = 0;
if strcmp(subj, '223')
    session_order = [1 4 2 5 3];
else
    session_order = 1:nsessions; % default
end
for sess = session_order
    curr_datafile = [data_dir, '/color_', subj, '_session_', num2str(sess), '_data.mat'];
    rawdata{sess} = load(curr_datafile);
    rawblocks{sess} = rawdata{sess}.data(:, 3);
    use_idx{sess} = ismember(rawblocks{sess}, blocks_to_use{sess});
    curr_data = rawdata{sess}.data(use_idx{sess}, :);
    curr_data(:, 3) = curr_data(:, 3) + last_block_idx; % make block numbers 1-9
    data = [data; curr_data];
    timing = [timing; rawdata{sess}.timing(use_idx{sess}, :)];
    last_block_idx = last_block_idx + numel(blocks_to_use{sess});
end

final_datafile = [data_dir, '/color_', subj, '_allsessions_data.mat'];
save(final_datafile, 'data', 'timing');