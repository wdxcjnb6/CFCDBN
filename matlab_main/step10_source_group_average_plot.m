
clc; clear variables;
restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED;

% === Set MATLAB and project paths ===
matlab_path = matlabroot;
current_project_path = fullfile(matlab_path, 'your_project_root');  % Replace with actual relative or absolute path

% === Toolbox paths ===
toolbox_paths = {
    fullfile(matlab_path, 'toolbox', 'eeglab2021.0');
    fullfile(matlab_path, 'toolbox', 'fieldtrip-20240110');
};
external_subfolders = {'mne','brewermap','matplotlib','cmocean','colorcet','gifti','spm12'};
fieldtrip_external_path = fullfile(matlab_path, 'toolbox', 'fieldtrip-master', 'external');
cellfun(@addpath, toolbox_paths);
cellfun(@(f) addpath(fullfile(fieldtrip_external_path,f)), external_subfolders);
addpath(genpath(fullfile(matlab_path, 'toolbox', 'fooof_mat-main')));
ft_defaults;

% === Project folder structure ===
addpath(fullfile(current_project_path, 'Code'));
addpath(fullfile(current_project_path, 'Data'));
addpath(fullfile(current_project_path, 'Tools'));
addpath(genpath(fullfile(current_project_path, 'Result')));
cd(fullfile(current_project_path, 'Code'));

% === I/O paths ===
input_root = fullfile(current_project_path, 'Result');
output_root = fullfile(current_project_path, 'Result');
source_location_path = fullfile(input_root, 'Source_result', 'source_location_data');
source_location_ana_save_path = fullfile(output_root, 'Source_result', 'source_location_group_plot');
f_createOutputFolders(output_root, {'Source_result/source_location_group_plot'});

% === Tool paths ===
source_tool_path = fullfile(current_project_path, 'Tools', 'souce_relevent');
channel_location_path = fullfile(current_project_path, 'Tools', 'sensor_relevent', 'HBN.sfp');
atlas = ft_read_atlas(fullfile(source_tool_path, 'ROI_MNI_V4.nii'));
atlas = ft_convert_units(atlas, 'mm');
elec = ft_read_sens(channel_location_path);
elec = ft_convert_units(elec, 'mm');

% === Group and band definitions ===
group = {'HC', 'Anxiety'};
band_roi = {'Delta', 'Beta'};
band_roi_data = containers.Map({'Delta', 'Beta'}, {[1 4], [13 30]});
pow_grop_summary = struct();

% === Subject list from external Excel ===
excel_table_path = fullfile(output_root, 'Demo_result', 'result_table_refine.xlsx');
file = readtable(excel_table_path);
subname = string(file{:, 1});

% === Load, clean and average subject-level data ===
for band_n = 1:numel(band_roi)
    band_name = band_roi{band_n};
    for group_n = 1:numel(group)
        group_name = group{group_n};
        save_file_path = fullfile(source_location_ana_save_path, [group_name, '_', band_name, '.mat']);
        if exist(save_file_path, 'file'), fprintf('File exists: %s\n', save_file_path); continue; end
        source_group_path = fullfile(source_location_path, group_name);
        file_list = dir(fullfile(source_group_path, ['*' band_name '*']));
        file_names = {file_list.name};
        valid_files = file_names(contains(file_names, subname));
        if isempty(valid_files), continue; end
        temp = load(fullfile(source_group_path, valid_files{1}));
        data_size = size(temp.source.avg.pow);
        source_pow_all = NaN([data_size, numel(valid_files)]);
        for i = 1:numel(valid_files)
            try
                d = load(fullfile(source_group_path, valid_files{i}));
                source_pow_all(:, :, i) = d.source.avg.pow;
            catch
                fprintf('Warning: failed to load %s\n', valid_files{i});
            end
        end
        pow_grop_summary.(group_name).(band_name) = source_pow_all;
        save(save_file_path, 'source_pow_all');
    end
end
save(fullfile(source_location_ana_save_path, 'all_group.mat'), 'pow_grop_summary');

% === Construct source structure for plotting ===
source_struct = struct();
source_stru = load(fullfile(source_tool_path, 'source_stru_HBN')).sub_source_stru;
source_stru = ft_convert_units(source_stru, 'mm');

for g = 1:numel(group)
    for b = 1:numel(band_roi)
        group_name = group{g}; band_name = band_roi{b};
        mat_file = fullfile(source_location_ana_save_path, [group_name, '_', band_name, '.mat']);
        if ~exist(mat_file, 'file'), continue; end
        dat = load(mat_file).source_pow_all;
        avg_pow = mean(dat, 3);
        source_struct.(group_name).(band_name) = source_stru;
        source_struct.(group_name).(band_name).avg.pow = avg_pow;
    end
end

% === Plot settings ===
save_path = fullfile(source_location_ana_save_path, 'figures');
if ~exist(save_path, 'dir'), mkdir(save_path); end
view_angles = {[90 0], [-90 0]};
view_names = {'rightside', 'leftside'};
hemis = {'left', 'right'};

cfg = [];
cfg.method = 'surface';
cfg.funparameter = 'pow';
cfg.projmethod = 'sphere_avg';
cfg.sphereradius = 10;
cfg.maskparameter = 'mask';

for g = 1:numel(group)
    for b = 1:numel(band_roi)
        group_name = group{g}; band_name = band_roi{b};
        if ~isfield(source_struct, group_name) || ~isfield(source_struct.(group_name), band_name)
            continue;
        end
        data = source_struct.(group_name).(band_name);
        for h = 1:numel(hemis)
            hemi = hemis{h};
            cfg.surffile = sprintf('surface_white_%s.mat', hemi);
            cfg.surfinflated = sprintf('surface_inflated_%s_caret.mat', hemi);
            for v = 1:numel(view_angles)
                fig = figure('Visible', 'off'); set(fig, 'Renderer', 'painters');
                ft_sourceplot(cfg, data);
                view(view_angles{v}); lighting gouraud; camlight('headlight');
                material dull; shading interp; axis vis3d off;
                if ft_hastoolbox('brewermap', 1)
                    colormap(flipud(brewermap(64, 'RdBu')));
                else
                    colormap(parula);
                end
                title(sprintf('%s | %s | %s | %s hemi', group_name, band_name, view_names{v}, hemi), ...
                    'FontSize', 14, 'FontWeight', 'bold');
                exportgraphics(fig, fullfile(save_path, ...
                    sprintf('source_%s_%s_%s_%s.png', group_name, band_name, view_names{v}, hemi)), 'Resolution', 600);
                close(fig);
            end
        end
    end
end