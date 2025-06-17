
clc; clear variables; restoredefaultpath;

% === PATH INITIALIZATION ===
matlab_project_path = '...';  % Define the root of the project
current_project_path = fullfile(matlab_project_path, 'PAC_network');
addpath(current_project_path);

% Toolbox setup (edit to match actual toolbox root paths)
toolbox_paths = {...
    fullfile('...', 'eeglab2021.0'),...
    fullfile('...', 'fieldtrip-master'),...
    fullfile('...', 'BrainNetViewer_20191031')
};
external_subfolders = {'mne', 'brewermap', 'matplotlib', 'cmocean', 'colorcet', 'gifti', 'spm12'};
fieldtrip_external_path = fullfile('...', 'fieldtrip-master', 'external');
cellfun(@addpath, toolbox_paths);
cellfun(@(subfolder) addpath(fullfile(fieldtrip_external_path, subfolder)), external_subfolders);

ft_defaults;
addpath(fullfile(current_project_path, 'Code'));
addpath(fullfile(current_project_path, 'Data'));
addpath(fullfile(current_project_path, 'Tools'));
addpath(genpath(fullfile(current_project_path, 'Result')));

% === INPUT/OUTPUT PATHS ===
data_method = 'TE';
input_root = fullfile(current_project_path, 'Result');
input_paths = {fullfile('Direct_net_result', 'Direct_net_data', data_method)};
output_paths = {'BrainNet_view'};
output_root = fullfile(current_project_path, 'Result');
f_createOutputFolders(output_root, output_paths);

source_net_path  = fullfile(input_root, input_paths{1});
graph_analysis_path = fullfile(output_root, output_paths{1});

% === ATLAS AND ROI DEFINITIONS ===
atlas = ft_read_atlas(fullfile('...', 'ROI_MNI_V4.nii'));
atlas = ft_convert_units(atlas, 'mm');
aal90 = atlas.tissuelabel(1:90);
channels_labels = load(fullfile('...', 'virtual_labels.mat')).virtual_label;

% Define regions of interest by network
% [... skip region assignment, sorting, and matching, already handled above ...]

% === LOAD BEHAVIORAL DATA AND TE NETWORK ===
group_name = 'Anxiety';
source_net_data_path = fullfile(source_net_path, group_name);
result_table_refine = readtable(fullfile('...', 'result_table_refine.xlsx'), 'VariableNamingRule', 'preserve');
refine_name = result_table_refine.subject;

% === LOAD AND MATCH SUBJECT DATA ===
folder_list_info = dir(fullfile(source_net_data_path, '*_te.mat'));
folder_list_sub_name = string({folder_list_info.name}');

% Determine network size using first subject
first_subject_data = load(fullfile(source_net_data_path, folder_list_sub_name{1}));
te_net = first_subject_data.te_mat;
[n_regions, ~] = size(te_net.data);

% Initialize containers
num_subjects = length(folder_list_sub_name);
all_networks = NaN(n_regions, n_regions, num_subjects);
SR_Total_scores = NaN(num_subjects, 1);
valid_sub_count = 0;

for sub_n = 1:num_subjects
    sub_file_name = folder_list_sub_name{sub_n};
    sub_name = sub_file_name(1:end-7);

    if ~ismember(sub_name, refine_name), continue; end

    te_net = load(fullfile(source_net_data_path, sub_file_name)).te_mat;
    if ~ismatrix(te_net.data) || size(te_net.data,1)~=n_regions, continue; end

    row_index = strcmp(result_table_refine.subject, sub_name);
    SR_Total_value = result_table_refine.SR_Total(row_index);
    large_net_order_sorted = te_net.data(sorted_indices, sorted_indices);

    valid_sub_count = valid_sub_count + 1;
    all_networks(:,:,valid_sub_count) = large_net_order_sorted;
    SR_Total_scores(valid_sub_count) = SR_Total_value;
end

all_networks = all_networks(:,:,1:valid_sub_count);
SR_Total_scores = SR_Total_scores(1:valid_sub_count);

% === CORRELATION ANALYSIS ===
missing_edges = all(isnan(all_networks) | all_networks==0, 3);
valid_edges = ~missing_edges;
r_values = NaN(n_regions);
p_values = NaN(n_regions);

for i = 1:n_regions
    for j = 1:n_regions
        if valid_edges(i,j)
            conn_values = squeeze(all_networks(i,j,:));
            valid_idx = ~isnan(conn_values) & ~isnan(SR_Total_scores);
            if sum(valid_idx) > 2
                [r_values(i,j), p_values(i,j)] = corr(conn_values(valid_idx), SR_Total_scores(valid_idx), 'Type', 'Spearman');
            end
        end
    end
end

% === FDR CORRECTION ===
valid_p_mask = ~isnan(p_values);
p_vector = p_values(valid_p_mask);
[p_fdr_vector, ~] = fdr_bh(p_vector, 0.05);
p_fdr = NaN(size(p_values));
p_fdr(valid_p_mask) = p_vector;
significant_edges = false(size(p_values));
significant_edges(valid_p_mask) = p_fdr(valid_p_mask) < 0.05;

% === EXPORT RESULTS ===
[row_idx, col_idx] = find(significant_edges);
edge_weights = r_values(significant_edges);

significant_connections = [brain_labels(row_idx)', brain_labels(col_idx)', num2cell(edge_weights)];
significant_table = cell2table(significant_connections, 'VariableNames', {'Source','Target','Correlation'});
writetable(significant_table, fullfile('...', 'significant_table.csv'));

% === PLOT GRAPH ===
G = digraph(col_idx, row_idx, edge_weights);
figure; p = plot(G, 'Layout', 'force', 'NodeLabel', brain_labels(unique([row_idx;col_idx])));
p.LineWidth = abs(edge_weights) * 10;
p.ArrowSize = max(abs(edge_weights) * 25, 12);
title('Significant Directed Network (Anxiety)');
colormap jet; colorbar; caxis([-1, 1]);

% === EXPORT EDGE MATRIX FOR BRAINNET VIEWER ===
significant_matrix = zeros(n_regions);
significant_matrix(significant_edges) = r_values(significant_edges);
dlmwrite(fullfile(graph_analysis_path, 'significant_connectivity_AAL34.edge'), significant_matrix, 'delimiter', '\t');
