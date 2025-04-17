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

% plot histograms of ms counts