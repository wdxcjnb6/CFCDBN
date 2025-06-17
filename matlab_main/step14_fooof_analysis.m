
clc; clear variables;
restoredefaultpath; clear RESTOREDEFAULTPATH_EXECUTED;

%% === Set Project and Toolbox Paths ===
project_root = '...';  % Replace with root directory
code_path = fullfile(project_root, 'Code');
data_path = fullfile(project_root, 'Data');
result_path = fullfile(project_root, 'Result');
toolbox_path = fullfile(project_root, 'Tools');

% Add toolboxes
addpath(fullfile(toolbox_path, 'eeglab2021.0'));
addpath(fullfile(toolbox_path, 'fieldtrip-master'));
external_ft = {'mne','brewermap','matplotlib','cmocean','colorcet','gifti','spm12'};
cellfun(@(f)addpath(fullfile(toolbox_path, 'fieldtrip-master','external',f)), external_ft);
addpath(genpath(fullfile(toolbox_path, 'fooof_mat-main')));
ft_defaults;
eeglab('nogui');

% Add project paths
addpath(code_path); addpath(data_path); addpath(result_path); addpath(toolbox_path);

%% === Configuration Parameters ===
group_list = {'HC','Anxiety'};
channel_file = fullfile(toolbox_path, 'sensor', 'HBN.sfp');
band_ranges = {'Delta', [1 4]; 'Theta', [4 7]; 'Alpha', [7 13]; 'Beta', [13 30]; 'Gamma', [30 45]};
threshold = 0.4;
peak_band = [13 30];

% Load subject inclusion table
T = readtable(fullfile(result_path, 'Demo_result', 'result_table_refine.xlsx'));
subject_names = string(T{:,1});

% Set input/output
eeg_input_root = fullfile(result_path, 'Preprocess_result', 'EEG_pure');
fooof_output_path = fullfile(result_path, 'Fooof_result', 'fooof_data_group_plot');
if ~exist(fooof_output_path, 'dir'), mkdir(fooof_output_path); end

%% === Group-wise FOOOF Feature Extraction ===
for g = 1:numel(group_list)
    group = group_list{g};
    group_path = fullfile(eeg_input_root, group);
    eeg_files = dir(fullfile(group_path, '*set'));
    fooof_group_data = {};
    k = 1;

    for i = 1:numel(eeg_files)
        sub_file = eeg_files(i).name;
        sub_id = sub_file(1:end-9);  % Remove _eeg.set or similar

        if ~ismember(sub_id, subject_names), continue; end
        EEG = f_safe_loadset(sub_file, group_path);
        if isempty(EEG), continue; end

        % Channel exclusion
        EEG = pop_select(EEG, 'nochannel', {'E1','E8','E14','E17','E21','E25','E32','E48','E49','E56','E63',...
            'E68','E73','E81','E88','E94','E99','E107','E113','E119','E125','E126','E127','E128'});

        % Convert to FieldTrip format
        cfg = []; cfg.method = 'trial'; cfg.channel = 'all';
        ft_data = eeglab2fieldtrip(EEG, 'preprocessing', cfg);
        cfg = []; cfg.trials = 'all'; cfg.continuous = 'yes';
        ft_data_cont = ft_redefinetrial(cfg, ft_data);
        cfg = []; cfg.length = 4; cfg.overlap = 0.8;
        ft_data_seg = ft_redefinetrial(cfg, ft_data_cont);

        % PSD estimation
        cfg = []; cfg.method = 'mtmfft'; cfg.output = 'pow'; cfg.tapsmofrq = 0.5; cfg.pad = 'nextpow2';
        freq = ft_freqanalysis(cfg, ft_data_seg);

        psd_avg = mean(freq.powspctrm, 1);
        freq_x = freq.freq;

        % Run FOOOF
        settings = struct('aperiodic_mode', 'fixed', 'verbose', true, ...
                          'max_n_peaks', 6, 'peak_width_limits', [2, 8]);
        fooof_result = fooof(freq_x, psd_avg, [0 45], settings, true);
        fooof_group_data{k} = fooof_result;
        k = k + 1;
    end

    save(fullfile(fooof_output_path, [group '_fooof_group.mat']), 'fooof_group_data');
end

%% === Extract and Compare Beta Peak Features ===
group_data = struct();
for g = 1:numel(group_list)
    group = group_list{g};
    load(fullfile(fooof_output_path, [group '_fooof_group.mat']), 'fooof_group_data');
    freqs = []; powers = [];

    for s = 1:numel(fooof_group_data)
        if isempty(fooof_group_data{s}), continue; end
        peaks = fooof_group_data{s}.peak_params;
        beta_idx = peaks(:,1) >= peak_band(1) & peaks(:,1) <= peak_band(2);
        beta_peaks = peaks(beta_idx, :);
        if isempty(beta_peaks), continue; end
        [~, best_idx] = max(beta_peaks(:,3));
        freqs(end+1) = beta_peaks(best_idx,1);
        powers(end+1) = beta_peaks(best_idx,3);
    end

    group_data.(group).freqs = freqs;
    group_data.(group).powers = powers;

    fprintf('[%s] Beta peak freq: %.2f ± %.2f Hz\n', group, mean(freqs), std(freqs));
    fprintf('[%s] Beta peak power: %.2f ± %.2f\n', group, mean(powers), std(powers));
end

% === Independent Sample t-Test ===
if all(isfield(group_data, {'Anxiety','HC'}))
    [~, p_freq, ~, s_freq] = ttest2(group_data.Anxiety.freqs, group_data.HC.freqs);
    [~, p_pow, ~, s_pow] = ttest2(group_data.Anxiety.powers, group_data.HC.powers);
    fprintf('\n[Freq] t(%.0f)=%.2f, p=%.4f\n', s_freq.df, s_freq.tstat, p_freq);
    fprintf('[Power] t(%.0f)=%.2f, p=%.4f\n', s_pow.df, s_pow.tstat, p_pow);
end
