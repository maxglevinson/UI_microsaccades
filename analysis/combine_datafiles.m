% code to merge several data files together.
% currently written for eyetracking data: 3 sessions of 3 blocks each.

scripts_dir = '/Users/maxlevinson/Documents/McGill/neurospeed/Uniformity_Illusion';
%scripts_dir = '/export03/data/mlevin/Uniformity_Illusion';
data_dir = [scripts_dir, '/bhv_data'];

subj = '224';
blocks_to_use = {[1:3], [1:3], [1:3]}; % default: 3 blocks for 3 sessions
%blocks_to_use = {[1:3], [1:3], [1], [1:4]}; % subject 210
%blocks_to_use = {[1:2], [1], [1:3], [1], [1:2]}; % subject 223
rawdata = {};
rawblocks = {};
use_idx = {};
nsessions = numel(blocks_to_use);
data = [];
timing = [];
last_block_idx = 0;
for sess = 1:nsessions
%for sess = [1 4 2 5 3] % subject 223
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