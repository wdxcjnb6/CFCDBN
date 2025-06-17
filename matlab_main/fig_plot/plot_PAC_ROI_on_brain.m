%% CLEAN AND SET PATHS
clc; clear; restoredefaultpath;

% === Basic project paths ===
matlab_path = matlabroot;
project_root = 'E:\';  % Replace with your base path
project_path = fullfile(project_root, 'PAC_network');
addpath(project_path);

% === Toolbox paths ===
toolbox_paths = {
    fullfile(matlab_path, 'toolbox', 'eeglab2021.0'), ...
    fullfile(matlab_path, 'toolbox', 'fieldtrip-master')
};
external_folders = {'mne', 'brewermap', 'matplotlib', 'cmocean', 'colorcet', 'gifti', 'spm12'};
fieldtrip_ext = fullfile(matlab_path, 'toolbox', 'fieldtrip-master', 'external');
cellfun(@addpath, toolbox_paths);
cellfun(@(f) addpath(fullfile(fieldtrip_ext, f)), external_folders);
ft_defaults;

% === Add project subpaths ===
addpath(fullfile(project_path, 'Code'));
addpath(fullfile(project_path, 'Data'));
addpath(fullfile(project_path, 'Tools'));
addpath(genpath(fullfile(project_path, 'Result')));

% === Define input/output paths ===
data_method = 'TE';
input_root = fullfile(project_path, 'Result');
input_dir = fullfile(input_root, 'Direct_net_result', 'Direct_net_data', data_method);
output_dir = fullfile(input_root, 'PAC_result');
pac_data_path = fullfile(output_dir, 'Data');
pac_fig_path = fullfile(output_dir, 'Fig');

% === Create output folders ===
f_createOutputFolders(input_root, {'PAC_result\Data\', 'PAC_result\Fig\'});

% === Load tools ===
source_tool = fullfile(project_path, 'Tools', 'souce_relevent');
excel_tool = fullfile(project_path, 'Tools', 'excel_relevent');
sensor_tool = fullfile(project_path, 'Tools', 'sensor_relevent');
addpath(genpath(fullfile(matlab_path, 'toolbox', 'SourceMesh-master')));
addpath(genpath(fullfile(matlab_path, 'toolbox', 'DrosteEffect-BrewerMap-3.2.5.0')));

% === Load AAL labels ===
AAL_Labels = load('AAL_labels.mat').labels;

%% ROI definition and reordering
channels_labels = load(fullfile(source_tool, 'virtual_labels.mat')).virtual_label;

% Define region groups from AAL90
regions_AAL90_part = {
    {'S_Hippocampus_L', 'S_Hippocampus_R'}, {'S_Frontal_Sup_Medial_L', 'S_Frontal_Sup_Medial_R'}, ...
    {'S_Precuneus_L', 'S_Precuneus_R'}, {'S_Angular_L', 'S_Angular_R'}, {'S_Cingulum_Post_L', 'S_Cingulum_Post_R'}, ...
    {'S_Insula_L', 'S_Insula_R'}, {'S_Cingulum_Ant_L', 'S_Cingulum_Ant_R'}, ...
    {'S_Amygdala_L', 'S_Amygdala_R'}, {'S_Thalamus_L', 'S_Thalamus_R'}, ...
    {'S_Frontal_Mid_L', 'S_Frontal_Mid_R'}, {'S_Parietal_Inf_L', 'S_Parietal_Inf_R'}, ...
    {'S_Frontal_Inf_Tri_L', 'S_Frontal_Inf_Tri_R'}, ...
    {'S_Frontal_Sup_L', 'S_Frontal_Sup_R'}, {'S_Parietal_Sup_L', 'S_Parietal_Sup_R'}, ...
    {'S_Postcentral_L', 'S_Postcentral_R'}, {'S_Precentral_L', 'S_Precentral_R'}, ...
    {'S_Supp_Motor_Area_L', 'S_Supp_Motor_Area_R'}
};

% Flatten region list
all_regions = horzcat(regions_AAL90_part{:});
roi_indices = find(ismember(channels_labels, all_regions));
roi_names = channels_labels(roi_indices);

% Load AAL abbreviations
aal90_abbr = readtable('Node_AAL90.txt', 'Delimiter', '\t', 'ReadVariableNames', false);
roi_names_abbr = table2cell(aal90_abbr(roi_indices, 6));

% Sort ROI names
[~, sorted_idx] = ismember(all_regions, roi_names);
sorted_idx = sorted_idx(sorted_idx > 0);
sorted_roi_names = roi_names(sorted_idx);
sorted_roi_names_abbr = roi_names_abbr(sorted_idx);

%% Load PAC data from Excel files
group = {'Anxiety', 'HC'};
input_dir_excel = fullfile(output_dir, 'Fig');
group_pac_data = cell(1, numel(group));

for g = 1:numel(group)
    T = readtable(fullfile(input_dir_excel, [group{g} '_pac_data.xlsx']), 'ReadRowNames', true);
    group_pac_data{g} = table2array(T);
    fprintf('%s group data loaded: %d × %d\n', group{g}, size(group_pac_data{g},1), size(group_pac_data{g},2));
end

%% Compute mean PAC and normalize
group_avg_pac = cellfun(@(x) mean(x, 2, 'omitnan'), group_pac_data, 'UniformOutput', false);
all_data = [group_avg_pac{1}; group_avg_pac{2}];
rescaled = rescale(all_data, 0.1, 1);
n1 = length(group_avg_pac{1});
group_avg_pac_rescaled = {rescaled(1:n1), rescaled(n1+1:end)};

%% Project PAC to source mesh (example: Anxiety group)
group_index = 2;  % 1=Anxiety, 2=HC
pac_vector = group_avg_pac_rescaled{group_index};

% Map ROI names to AAL indices
aal_indices = zeros(length(sorted_roi_names), 1);
for i = 1:length(sorted_roi_names)
    idx = find(strcmp(AAL_Labels, sorted_roi_names{i}));
    if ~isempty(idx), aal_indices(i) = idx; end
end

% Fill PAC vector across AAL space
whole_brain_map = zeros(numel(AAL_Labels), 1);
whole_brain_map(aal_indices) = pac_vector;

% Load pastel colormap
mp = load('pastelmap'); cmap = flipud(mp.map);

% === Left hemisphere (lateral) ===
figure('Position', [100 100 1000 800]);
atemplate('hemi','l','mesh','def4','overlay',whole_brain_map,'method',{'aal_super', []});
colormap(cmap); caxis([-max(whole_brain_map), max(whole_brain_map)]);
view([90 0]); title('Left Hemisphere - Lateral View');

% === Left hemisphere (medial) ===
figure('Position', [100 100 1000 800]);
atemplate('hemi','l','mesh','def4','overlay',whole_brain_map,'method',{'aal_super', []});
colormap(cmap); caxis([-max(whole_brain_map), max(whole_brain_map)]);
view([-90 0]); title('Left Hemisphere - Medial View');

% === Right hemisphere (lateral) ===
figure('Position', [100 100 1000 800]);
atemplate('hemi','r','mesh','def4','overlay',whole_brain_map,'method',{'aal_super', []});
colormap(cmap); caxis([-max(whole_brain_map), max(whole_brain_map)]);
view([-90 0]); title('Right Hemisphere - Lateral View');

% === Right hemisphere (medial) ===
figure('Position', [100 100 1000 800]);
atemplate('hemi','r','mesh','def4','overlay',whole_brain_map,'method',{'aal_super', []});
colormap(cmap); caxis([-max(whole_brain_map), max(whole_brain_map)]);
view([90 0]); title('Right Hemisphere - Medial View');
