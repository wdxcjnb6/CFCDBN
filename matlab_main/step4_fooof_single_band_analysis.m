%% CLEAN AND SET PATHS
clc; clear variables; restoredefaultpath; clear RESTOREDEFAULTPATH_EXECUTED;

% === Base path setup ===
matlab_path = matlabroot;
project_root = fullfile(matlab_path, 'project_folder');

% === Toolbox paths ===
eeglab_path = fullfile(matlab_path, 'toolbox', 'eeglab2021.0');
fieldtrip_path = fullfile(matlab_path, 'toolbox', 'fieldtrip-master');
addpath(eeglab_path);
addpath(fieldtrip_path);
ft_defaults;

external_subfolders = {'mne','brewermap','matplotlib','cmocean','colorcet','gifti','spm12'};
for i = 1:numel(external_subfolders)
    addpath(fullfile(fieldtrip_path, 'external', external_subfolders{i}));
end

addpath(fullfile(project_root, 'Code'));
addpath(fullfile(project_root, 'Data'));
addpath(fullfile(project_root, 'Tools'));
addpath(fullfile(project_root, 'Result'));
addpath(genpath(fullfile(matlab_path, 'toolbox', 'fooof_mat-main')));

%% DEFINE I/O PATHS
group_list = {'HC', 'Anxiety'};
input_root = fullfile(project_root, 'Result');
input_paths = {'Preprocess_result/EEG_pure/'};
tool_root = fullfile(project_root, 'Tools');
tool_subdirs = {'souce_relevent/', 'excel_relevent/', 'sensor_relevent/'};
output_paths = {
    'Fooof_result/fooof_data/', ...
    'Fooof_result/fooof_fig/', ...
    'Fooof_result/fooof_excel/', ...
    'Fooof_result/fooof_sub_select_periodic/', ...
    'Fooof_result/fooof_sub_select_aperiodic/' ...
};

f_createOutputFolders(input_root, output_paths);

% Paths
EEG_pure_path = fullfile(input_root, input_paths{1});
fooof_data_path = fullfile(input_root, output_paths{1});
fooof_fig_path = fullfile(input_root, output_paths{2});
fooof_excel_path = fullfile(input_root, output_paths{3});
fooof_sub_select_periodic_path = fullfile(input_root, output_paths{4});
fooof_sub_select_aperiodic_path = fullfile(input_root, output_paths{5});
channel_location_file = fullfile(tool_root, tool_subdirs{3}, 'HBN.sfp');

%% PARAMETERS
band_ranges = {
    'Delta', [1 4], 'c';
    'Theta', [4 7], 'b';
    'Alpha', [7 13], 'm';
    'Beta', [13 30], 'g';
    'Gamma', [30 45], 'r'
};
threshold = 0.4;
excel_headers = {'Subject', 'Group', 'Freq Peak Fr', 'Freq Peak Pw', 'Freq Peak Bw', 'Freq Ratio'};
excel_data = {};
excel_filename = fullfile(fooof_sub_select_periodic_path, 'Single_Band_FOOOF_Results.xlsx');
freq_pairs = nchoosek(band_ranges(:,1), 1);

%% BEGIN FOOOF LOOP
eeglab nogui;
for group_n = 1:numel(group_list)
    group = group_list{group_n};
    sub_group_path = fullfile(EEG_pure_path, group);
    file_list = dir(fullfile(sub_group_path, '*.set'));
    file_list(1:2) = [];
    sub_names = {file_list.name}';

    fprintf('Group %s: %d subjects\n', group, numel(file_list));

    for s = 1:numel(sub_names)
        sub_name = sub_names{s};
        EEG = f_safe_loadset(sub_name, sub_group_path);
        EEG = pop_select(EEG, 'nochannel', {'E1','E8','E14','E17','E21','E25','E32','E48','E49','E56','E63','E68',...
                                             'E73','E81','E88','E94','E99','E107','E113','E119','E125','E126','E127','E128'});
        if isempty(EEG), warning('Skipped empty EEG'); continue; end

        cfg = []; cfg.channel = 'all'; cfg.latency = 'all'; cfg.method = 'trial';
        ft_data = eeglab2fieldtrip(EEG, 'preprocessing', cfg);
        channel_label_all = ft_data.label;

        cfg = []; cfg.trials = 'all'; cfg.continuous = 'yes';
        ft_data_cont = ft_redefinetrial(cfg, ft_data);

        cfg = []; cfg.length = 4; cfg.overlap = 0.8;
        ft_data_seg = ft_redefinetrial(cfg, ft_data_cont);

        cfg = []; cfg.method = 'mtmfft'; cfg.output = 'pow'; cfg.pad = 'nextpow2'; cfg.tapsmofrq = 0.5;
        fre_data = ft_freqanalysis(cfg, ft_data_seg);

        freq_x = fre_data.freq;
        psd_y = fre_data.powspctrm;
        psd_avg = mean(psd_y, 1);
        if any(isnan(psd_avg)) || any(isinf(psd_avg)) || any(psd_avg <= 0), continue; end

        settings = struct('aperiodic_mode','fixed','verbose',true,'max_n_peaks',6,'peak_width_limits',[2,8]);
        fooof_results_all = fooof(freq_x, psd_avg, [1 40], settings, true);

        fooof_results_aperiodic.exp = fooof_results_all.aperiodic_params(2);
        fooof_results_aperiodic.psd = fooof_results_all.power_spectrum;
        fooof_results_aperiodic.ap_fit = fooof_results_all.ap_fit;

        aperiodic_save_dir = fullfile(fooof_sub_select_aperiodic_path, group);
        f_create_directories({aperiodic_save_dir});
        save(fullfile(aperiodic_save_dir, sprintf('%s.mat', sub_name(1:min(12,end)))), 'fooof_results_aperiodic');

        fig_save_dir = fullfile(fooof_fig_path, group);
        f_create_directories({fig_save_dir});
        fooof_plot(fooof_results_all, false); drawnow;
        saveas(gcf, fullfile(fig_save_dir, sprintf('FOOOF_%s.png', sub_name(1:12))));
        close(gcf);

        fig2 = figure('Visible','off');
        plot(freq_x, smooth(fooof_results_all.fooofed_spectrum - fooof_results_all.ap_fit), 'k', 'LineWidth', 2);
        hold on; f_plot_fooof_peaks(fooof_results_all, band_ranges);
        title(['Oscillatory Spectrum for ', sub_name(1:12)]);
        xlabel('Frequency (Hz)'); ylabel('Power'); grid on; hold off;
        saveas(fig2, fullfile(fig_save_dir, sprintf('Spectrum_%s.png', sub_name(1:12))));
        close(fig2);

        for pair_n = 1:size(freq_pairs,1)
            low_band = freq_pairs{pair_n,1};
            high_band = freq_pairs{pair_n,1};

            fooof_results = f_process_fooof_and_calculate_ratios(freq_x, psd_y, [1 40], ...
                settings, band_ranges, channel_label_all, low_band, high_band);

            if fooof_results.low_freq.ratio >= threshold && fooof_results.high_freq.ratio >= threshold
                periodic_save_dir = fullfile(fooof_sub_select_periodic_path, group);
                f_create_directories({periodic_save_dir});
                save(fullfile(periodic_save_dir, sprintf('%s_%s.mat', sub_name(1:min(12,end)), low_band)), 'fooof_results');

                excel_data(end+1,:) = {
                    sub_name(1:min(12,end)), group, ...
                    num2str(fooof_results.low_freq.avg(1),'%.2f'), ...
                    num2str(fooof_results.low_freq.avg(2),'%.2f'), ...
                    num2str(fooof_results.low_freq.avg(3),'%.2f'), ...
                    num2str(fooof_results.low_freq.ratio,'%.2f')};

                sheet_name = ['_' high_band '_'];
                writetable(cell2table(excel_data(end,:), 'VariableNames', excel_headers), ...
                    excel_filename, 'Sheet', sheet_name, 'WriteMode', 'append');
            end
        end
    end
end
