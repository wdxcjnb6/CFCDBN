%% CLEAN AND SET PATHS
clc; clear variables; restoredefaultpath; clear RESTOREDEFAULTPATH_EXECUTED;

% === Set base paths ===
matlab_path = matlabroot;
current_project_path = fullfile(matlab_path, 'project_folder');

% === Add toolbox paths ===
eeglab_path = fullfile(matlab_path, 'toolbox', 'eeglab2021.0');
fieldtrip_path = fullfile(matlab_path, 'toolbox', 'fieldtrip-master');
addpath(eeglab_path);
addpath(fieldtrip_path);
ft_defaults;

% Add FieldTrip externals
external_subfolders = {'mne','brewermap','matplotlib','cmocean','colorcet','gifti','spm12'};
for i = 1:length(external_subfolders)
    addpath(fullfile(fieldtrip_path, 'external', external_subfolders{i}));
end

% Add project paths
addpath(fullfile(current_project_path, 'Code'));
addpath(fullfile(current_project_path, 'Data'));
addpath(fullfile(current_project_path, 'Tools'));
addpath(fullfile(current_project_path, 'Result'));

% Add FOOOF package
addpath(genpath(fullfile(matlab_path, 'toolbox', 'fooof_mat-main')));

%% SET I/O PATHS AND PARAMETERS
group = {'HC','Anxiety'};
input_root = fullfile(current_project_path, 'Result');
input_paths = {'Preprocess_result/EEG_pure/'};
output_paths = {'Fooof_result/fooof_data/'};
input_tool_root = fullfile(current_project_path, 'Tools');
input_tool_paths = {'souce_relevent/','excel_relevent/','sensor_relevent/'};

% Create output folders
f_createOutputFolders(input_root, output_paths);

% Define paths
EEG_pure_path = fullfile(input_root, input_paths{1});
fooof_data_path = fullfile(input_root, output_paths{1});
Channel_location_tool_path = fullfile(input_tool_root, input_tool_paths{3});
Channel_location = [Channel_location_tool_path 'HBN.sfp'];

band_ranges = {
    'Delta', [1 4], 'c';
    'Theta', [4 7], 'b';
    'Alpha', [7 13], 'm';
    'Beta', [13 30], 'g';
    'Gamma', [30 45], 'r'
};

threshold = 0.4;
eeglab nogui;

% Read subject list
subject_table = readtable(fullfile(current_project_path, 'Result', 'Demo_result', 'result_table_refine.xlsx'));
subname = subject_table{:,1};

%% FOOOF ANALYSIS LOOP
for group_n = 1:numel(group)
    fooof_group = {}; k = 1;
    group_name = group{group_n};
    sub_group_path = fullfile(EEG_pure_path, group_name);
    subject_files = dir(fullfile(sub_group_path, '*.set'));
    fprintf('Subject count in %s: %d\n', group_name, numel(subject_files));

    for sub_n = 1:numel(subject_files)
        sub_file = subject_files(sub_n).name;
        sub_id = sub_file(1:end-9);
        if ~ismember(sub_id, subname), continue; end

        EEG = f_safe_loadset(sub_file, sub_group_path);
        EEG = pop_select(EEG, 'nochannel', {'E1','E8','E14','E17','E21','E25','E32','E48','E49','E56','E63','E68',...
                                            'E73','E81','E88','E94','E99','E107','E113','E119','E125','E126','E127','E128'});
        if isempty(EEG), warning('Skipping file due to load error.'); continue; end

        % Convert to FieldTrip format
        cfg = []; cfg.channel = 'all'; cfg.latency = 'all'; cfg.method = 'trial';
        ft_data = eeglab2fieldtrip(EEG, 'preprocessing', cfg);

        % Concatenate trials into continuous data
        cfg = []; cfg.trials = 'all'; cfg.continuous = 'yes';
        ft_data_cont = ft_redefinetrial(cfg, ft_data);

        % Redefine trials with overlap
        cfg = []; cfg.length = 4; cfg.overlap = 0.8;
        ft_data_redef = ft_redefinetrial(cfg, ft_data_cont);

        % Power spectrum estimation
        cfg = [];
        cfg.pad = 'nextpow2';
        cfg.tapsmofrq = 0.5;
        cfg.method = 'mtmfft';
        cfg.output = 'pow';
        fre_data = ft_freqanalysis(cfg, ft_data_redef);

        freq_x = fre_data.freq;
        psd_y = fre_data.powspctrm;
        psd_avg = mean(psd_y, 1);

        % FOOOF parameter settings
        settings = struct();
        settings.aperiodic_mode = 'fixed';
        settings.verbose = true;
        settings.max_n_peaks = 6;
        settings.peak_width_limits = [2, 8];
        fooof_result = fooof(freq_x, psd_avg, [1 45], settings, true);

        fooof_group{k} = fooof_result;
        k = k + 1;
    end

    save(fullfile(fooof_data_path, [group_name '_fooof_group.mat']), 'fooof_group');
end

%% PLOT FOOOF MODEL FIT ACROSS GROUPS
group_name = {'HC', 'Anxiety'};
colors = {[1, 0, 0], [0, 0, 1]}; % red for HC, blue for Anxiety
figure; hold on;

for group_n = 1:numel(group_name)
    group = group_name{group_n};
    fooof_group = load(fullfile(fooof_data_path, [group '_fooof_group.mat'])).fooof_group;
    model_data = [];

    for sub_n = 1:numel(fooof_group)
        if isempty(fooof_group{sub_n}), continue; end
        freq_x = fooof_group{sub_n}.freqs;
        model_y = fooof_group{sub_n}.fooofed_spectrum;
        model_data(sub_n, :) = model_y;
    end

    model_mean = mean(model_data, 1);
    model_std = std(model_data, [], 1);

    fill([freq_x, fliplr(freq_x)], ...
         [model_mean - model_std, fliplr(model_mean + model_std)], ...
         colors{group_n}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(freq_x, model_mean, 'Color', colors{group_n}, 'LineWidth', 2.5);
end

legend({'HC ± STD', 'HC Mean', 'Anxiety ± STD', 'Anxiety Mean'}, 'Location', 'northeast');
xlabel('Frequency (Hz)');
ylabel('Power (a.u.)');
title('FOOOF Model Fit Across Groups');
grid on;
hold off;
