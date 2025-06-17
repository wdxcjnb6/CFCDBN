

clc; clear variables;
restoredefaultpath;

% === Path Initialization ===
matlab_path = matlabroot;
project_root = fullfile('..', 'PAC_network');  % Replace with actual relative or absolute path
addpath(project_root);

% Toolbox paths
toolbox_paths = {
    fullfile(matlab_path, 'toolbox', 'eeglab2021.0'),...
    fullfile(matlab_path, 'toolbox', 'fieldtrip-master')};
external_subfolders = {'mne', 'brewermap', 'matplotlib', 'cmocean', 'colorcet', 'gifti', 'spm12'};
fieldtrip_ext = fullfile(matlab_path, 'toolbox', 'fieldtrip-master', 'external');
cellfun(@addpath, toolbox_paths);
cellfun(@(f) addpath(fullfile(fieldtrip_ext, f)), external_subfolders);

% Initialize FieldTrip
ft_defaults;

% Project subfolders
addpath(fullfile(project_root, 'Code'));
addpath(fullfile(project_root, 'Data'));
addpath(fullfile(project_root, 'Tools'));
addpath(genpath(fullfile(project_root, 'Result')));

% Output paths
output_root = fullfile(project_root, 'Result');
source_fig_path = fullfile(output_root, 'Direct_net_result', 'Direct_net_figure', 'Figure');
if ~exist(source_fig_path, 'dir'), mkdir(source_fig_path); end

% SourceMesh toolbox
addpath(genpath(fullfile(matlab_path, 'toolbox', 'SourceMesh-master')));

% Load AAL label list
AAL_labels = load('AAL_labels.mat').labels;

%% === Define ROI Networks ===
DMN_roi = {'Hippocampus_L', 'Hippocampus_R', 'Frontal_Sup_Medial_L', 'Frontal_Sup_Medial_R', ...
           'Precuneus_L', 'Precuneus_R', 'Angular_L', 'Angular_R', ...
           'Cingulum_Post_L', 'Cingulum_Post_R'};

SN_roi  = {'Insula_L', 'Insula_R', 'Cingulum_Ant_L', 'Cingulum_Ant_R', ...
           'Amygdala_L', 'Amygdala_R', 'Thalamus_L', 'Thalamus_R'};

FPN_roi = {'Frontal_Mid_L', 'Frontal_Mid_R', 'Parietal_Inf_L', 'Parietal_Inf_R', ...
           'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R'};

DAN_roi = {'Frontal_Sup_L', 'Frontal_Sup_R', 'Parietal_Sup_L', 'Parietal_Sup_R'};

SMN_roi = {'Postcentral_L', 'Postcentral_R', 'Precentral_L', 'Precentral_R', ...
           'Supp_Motor_Area_L', 'Supp_Motor_Area_R'};

roi_list = {DMN_roi, SN_roi, FPN_roi, DAN_roi, SMN_roi};
network_names = {'DMN', 'SN', 'FPN', 'DAN', 'SMN'};

%% === Visualize Networks with AAL Template ===
for net_idx = 1:length(roi_list)
    overlay = zeros(numel(AAL_labels), 1);
    roi = roi_list{net_idx};
    roi = sort(roi(:));
    roi_indices = zeros(length(roi), 1);
    color = 1;

    for i = 1:length(roi)
        region_name = roi{i};
        idx = find(strcmp(AAL_labels, region_name));
        if isempty(idx)
            fprintf('Warning: ROI not found in AAL_labels -> %s\n', region_name);
        else
            roi_indices(i) = idx;
            overlay(idx) = color;
        end
    end

    figure('Position', [480 160 900 800]);
    atemplate('mesh', 'def4', 'overlay', overlay, 'method', {'aal_super', []});

    % Create color map: gray (default) + distinct parula tone
    gray_map = repmat([1 1 1], 1, 1);  % gray
    parula_map = parula(256);
    cmap_custom = [gray_map; parula_map(50*net_idx,:)];
    colormap(cmap_custom);

    % Color axis + transparency
    caxis([min(overlay) max(overlay)]);
    alpha 0.7;

    sgtitle(['Network: ' network_names{net_idx}], 'FontSize', 16);

    % Save as .svg
    filename = fullfile(source_fig_path, [network_names{net_idx} '.svg']);
    saveas(gcf, filename);
    fprintf('✅ Saved: %s\n', filename);
    close all;
end
