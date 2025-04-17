% get statistics about trial counts / FT in trials with or without
% microsaccades (& blinks)
data = readtable('trial_data_withblinks.csv');

subjects = unique(data.subject);
n_subjects = numel(subjects);

subject_col = data.subject;
ft = data.FT;
mspresent = data.mspresent;
blinkpresent = data.nblinks > 0;
baselinenms = data.baselinenms;
baseline_msrate = baselinenms ./ ft;

[avg_ft_withms, avg_ft_noms, avg_msrate, n_trials_noms] = deal(zeros(1, n_subjects));
all_ms_counts = cell(1, n_subjects);

for i_subj = 1:n_subjects
    withms_idx = subject_col == subjects(i_subj) & (mspresent == 1 | blinkpresent == 1);
    noms_idx = subject_col == subjects(i_subj) & (mspresent == 0 & blinkpresent == 0);
    avg_ft_withms(i_subj) = mean(ft(withms_idx));
    avg_ft_noms(i_subj) = mean(ft(noms_idx));
    avg_msrate(i_subj) = mean(baseline_msrate(subject_col == subjects(i_subj)));
    n_trials_noms(i_subj) = sum(noms_idx);
    all_ms_counts{i_subj} = baselinenms(subject_col == subjects(i_subj));
end

% mean + stderr
perc_trials_noms = n_trials_noms / 360; % 360 main task trials
mean_perc = mean(perc_trials_noms);
stderr_perc = std(perc_trials_noms) ./ sqrt(n_subjects);

mean_msrate = mean(avg_msrate);
stderr_msrate = std(avg_msrate) ./ sqrt(n_subjects);

% plot histograms of ms counts
n_rows = 3;
n_cols = 7;
binedges_tmp = [0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 100];
binedges_plot = binedges_tmp; binedges_plot(end) = 11;
f = figure; hold on;
for i_subj = 1:n_subjects
    subplot(n_rows, n_cols, i_subj);
    bincounts = histcounts(all_ms_counts{i_subj}, binedges_tmp);
    h_plot = histogram('BinEdges', binedges_plot, 'BinCounts', bincounts, 'FaceColor', [0.2 0.2 0.2]);
    xticks([0.5:1:9.5, 11]);
    xticklabels({'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10+'});
    xtickangle(0);
    title(['participant ', num2str(i_subj)]);
end
yl = suplabel('trial count', 'y', [0.12 .08 .84 .84]);
yl.FontSize = 14;
xl = suplabel('number of microsaccades per trial', 'x', [0.08 .09 .84 .84]);
xl.FontSize = 14;
set(f, 'Position', [1 58 1440 739]);
exportgraphics(f, 'ms_count_hists.pdf', 'ContentType', 'vector');