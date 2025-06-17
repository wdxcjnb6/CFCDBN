
clc; clear variables; restoredefaultpath;

%% === Toolbox Initialization ===
matlab_path = matlabroot;
project_path = 'E:\PAC_network';

% Toolbox paths
addpath(fullfile(matlab_path, 'toolbox', 'eeglab2021.0'));
addpath(fullfile(matlab_path, 'toolbox', 'fieldtrip-master'));
ft_defaults;

% Add FieldTrip externals
ext_subfolders = {'mne','brewermap','matplotlib','cmocean','colorcet','gifti','spm12'};
cellfun(@(f) addpath(fullfile(matlab_path, 'toolbox', 'fieldtrip-master', 'external', f)), ext_subfolders);

% Add custom tool paths
addpath(genpath(fullfile(matlab_path, 'toolbox', 'fooof_mat-main')));
addpath(fullfile(project_path, 'Code'));
addpath(fullfile(project_path, 'Data'));
addpath(fullfile(project_path, 'Tools'));
addpath(fullfile(project_path, 'Result'));

%% === Define Paths ===
group_list = {'HC', 'Anxiety'};
input_root = fullfile(project_path, 'Result', 'Preprocess_result', 'EEG_pure');
fooof_output_path = fullfile(project_path, 'Result', 'Fooof_result', 'fooof_data_group_plot');
if ~exist(fooof_output_path, 'dir'), mkdir(fooof_output_path); end

% ROI and channel config
chanloc_file = fullfile(project_path, 'Tools', 'sensor_relevent', 'HBN.sfp');
band_ranges = {'Delta', [1 4]; 'Theta', [4 7]; 'Alpha', [7 13]; 'Beta', [13 30]; 'Gamma', [30 45]};
beta_band = [13 30];

% Read subject list from Excel
file = readtable(fullfile(project_path, 'Result', 'Demo_result', 'result_table_refine.xlsx'));
subject_names = string(file{:, 1});

%% === Run FOOOF on Each Group ===
for g = 1:numel(group_list)
    group = group_list{g};
    group_dir = fullfile(input_root, group);
    file_list = dir(fullfile(group_dir, '*set'));
    fooof_group_data = {};

    for i = 1:numel(file_list)
        sub_file = file_list(i).name;
        sub_id = sub_file(1:end-9);
        if ~ismember(sub_id, subject_names), continue; end

        EEG = f_safe_loadset(sub_file, group_dir);
        EEG = pop_select(EEG, 'nochannel', {'E1','E8','E14','E17','E21','E25','E32','E48','E49',...
            'E56','E63','E68','E73','E81','E88','E94','E99','E107','E113','E119','E125','E126','E127','E128'});
        if isempty(EEG), continue; end

        % Convert to FieldTrip
        cfg = []; cfg.channel = 'all'; cfg.latency = 'all'; cfg.method = 'trial';
        ft_data = eeglab2fieldtrip(EEG, 'preprocessing', cfg);

        % Concatenate trials to continuous and segment
        cfg = []; cfg.trials = 'all'; cfg.continuous = 'yes';
        ft_data = ft_redefinetrial(cfg, ft_data);
        cfg = []; cfg.length = 4; cfg.overlap = 0.8;
        ft_data = ft_redefinetrial(cfg, ft_data);

        % Compute PSD
        cfg = []; cfg.method = 'mtmfft'; cfg.output = 'pow'; cfg.tapsmofrq = 0.5; cfg.pad = 'nextpow2';
        freq_data = ft_freqanalysis(cfg, ft_data);
        psd_avg = mean(freq_data.powspctrm, 1);

        % Run FOOOF
        settings = struct('aperiodic_mode','fixed','max_n_peaks',6,'peak_width_limits',[2 8],'verbose',true);
        fooof_result = fooof(freq_data.freq, psd_avg, [0 45], settings, true);
        fooof_group_data{end+1} = fooof_result;
    end
    save(fullfile(fooof_output_path, [group '_fooof_group.mat']), 'fooof_group_data');
end

%% === Plot FOOOF Fit Across Groups ===
colors = {[0, 0.447, 0.741], [0.85, 0.325, 0.098]};
shadow_color = [0.7, 0.7, 0.7];
lw = 2.5;
band_to_highlight = {'Delta',[1 4]; 'Beta',[13 30]};

all_model_data = cell(1,numel(group_list));
all_model_mean = cell(1,numel(group_list));
all_model_std = cell(1,numel(group_list));
all_freq_x = cell(1,numel(group_list));

for g = 1:numel(group_list)
    group = group_list{g};
    fooof_data = load(fullfile(fooof_output_path, [group '_fooof_group.mat']));
    model_data = load(fullfile(fooof_output_path, ['model_data_' group '.mat'])).model_data;

    freq_x = fooof_data.fooof_group_data{1}.freqs;
    all_model_data{g} = model_data;
    all_model_mean{g} = mean(model_data,1);
    all_model_std{g} = std(model_data,[],1);
    all_freq_x{g} = freq_x;
end

figure; hold on;
for g = 1:numel(group_list)
    freq_x = all_freq_x{g};
    m = all_model_mean{g}; s = all_model_std{g};
    fill([freq_x, fliplr(freq_x)], [m-s, fliplr(m+s)], shadow_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(freq_x, m, 'Color', colors{g}, 'LineWidth', lw);
end
for b = 1:size(band_to_highlight,1)
    r = band_to_highlight{b,2}; yl = ylim;
    fill([r(1),r(2),r(2),r(1)], [yl(1),yl(1),yl(2),yl(2)], [0.9,0.9,0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    text(mean(r), yl(2)*0.9, band_to_highlight{b,1}, 'HorizontalAlignment', 'center', 'FontWeight','bold');
end
legend({'HC','Anxiety'}); xlabel('Frequency (Hz)'); ylabel('Power (a.u.)');
title('FOOOF Model Fits (Mean ± SD)');
grid on; hold off;

%% === Save Group-Level Features ===
save_dir = fullfile(project_path, 'Figure', 'fooof');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
cd(save_dir);

freq_x = all_freq_x{1};
save('fooof_results.mat', 'freq_x', 'all_model_data', 'all_model_mean', 'all_model_std', 'group_list', 'band_to_highlight');

%% === Statistical Comparison of Beta Peaks ===
all_labels = {}; all_peak_freqs = []; all_peak_powers = []; all_subjects = {}; group_stats = [];

for g = 1:numel(group_list)
    group = group_list{g}; model_data = all_model_data{g}; freq_x = all_freq_x{g};
    power_list = []; freq_list = [];

    for s = 1:size(model_data,1)
        idx = freq_x >= beta_band(1) & freq_x <= beta_band(2);
        beta_vals = model_data(s,idx); beta_freqs = freq_x(idx);
        [p, idx_peak] = max(beta_vals);
        f = beta_freqs(idx_peak);

        all_labels{end+1} = group;
        all_peak_powers(end+1) = p;
        all_peak_freqs(end+1) = f;
        all_subjects{end+1} = sprintf('%s_sub%02d', group, s);

        power_list(end+1) = p;
        freq_list(end+1) = f;
    end

    group_stats(end+1).Group = group;
    group_stats(end).Mean_Power = mean(power_list);
    group_stats(end).Mean_Freq = mean(freq_list);
    group_stats(end).Std_Power = std(power_list);
    group_stats(end).Std_Freq = std(freq_list);
    group_stats(end).Power_List = power_list(:);
    group_stats(end).Freq_List = freq_list(:);
end

% Two-sample t-tests
[~, p1, ~, stat1] = ttest2(group_stats(1).Power_List, group_stats(2).Power_List);
[~, p2, ~, stat2] = ttest2(group_stats(1).Freq_List, group_stats(2).Freq_List);
group_stats(1).Tstat_Power = stat1.tstat; group_stats(2).Tstat_Power = stat1.tstat;
group_stats(1).Tstat_Freq = stat2.tstat; group_stats(2).Tstat_Freq = stat2.tstat;
group_stats(1).Pval_Power = p1; group_stats(2).Pval_Power = p1;
group_stats(1).Pval_Freq = p2; group_stats(2).Pval_Freq = p2;

% Save results
individual_table = table(all_subjects', all_labels', all_peak_powers', all_peak_freqs', ...
    'VariableNames', {'Subject', 'Group', 'Beta_Peak_Power', 'Beta_Peak_Frequency'});
writetable(individual_table, 'beta_individual_peak_data.csv');

stat_table = struct2table(group_stats);
stat_table.Power_List = []; stat_table.Freq_List = [];
writetable(stat_table, 'beta_group_statistics.csv');

% Boxplots
figure;
subplot(1,2,1); boxplot(all_peak_powers, all_labels); ylabel('Beta Peak Power'); title('Power Comparison');
subplot(1,2,2); boxplot(all_peak_freqs, all_labels); ylabel('Beta Peak Frequency (Hz)'); title('Frequency Comparison');
saveas(gcf, 'beta_peak_boxplots.png');

disp('✅ All results saved (MAT + CSV + PNG)');
