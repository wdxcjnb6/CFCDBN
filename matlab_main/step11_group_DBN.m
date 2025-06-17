%% === CLEAN & PATH SETUP ===
clc; clearvars;
restoredefaultpath;
matlab_path = matlabroot;
project_root = 'YOUR_PROJECT_PATH'; % <-- Replace with your project root path
current_path = fullfile(project_root, 'PAC_network');
addpath(current_path);

% Toolbox paths
addpath(fullfile(matlab_path, 'toolbox', 'eeglab2021.0'));
addpath(fullfile(matlab_path, 'toolbox', 'fieldtrip-master'));
ft_defaults;

% Add FieldTrip external subfolders
ext_folders = {'mne','brewermap','matplotlib','cmocean','colorcet','gifti','spm12'};
cellfun(@(f)addpath(fullfile(matlab_path,'toolbox','fieldtrip-master','external',f)), ext_folders);

% Add project subfolders
addpath(fullfile(current_path, 'Code'));
addpath(fullfile(current_path, 'Data'));
addpath(fullfile(current_path, 'Tools'));
addpath(genpath(fullfile(current_path, 'Result')));

%% === DEFINE PATHS ===
data_method = 'TE';
input_root = fullfile(current_path, 'Result');
output_root = input_root;
input_paths = {fullfile('Direct_net_result', 'Direct_net_data', data_method)};
output_paths = {
    fullfile('Direct_net_result','Direct_net_figure','Data'),...
    fullfile('R_result','large_direct_net')};

f_createOutputFolders(output_root, output_paths);
source_net_path = fullfile(input_root, input_paths{1});
graph_data_path = fullfile(output_root, output_paths{1});
graph_R_path = fullfile(output_root, output_paths{2});

% Tool paths
tool_root = fullfile(current_path, 'Tools');
source_tool = fullfile(tool_root, 'souce_relevent');
sensor_tool = fullfile(tool_root, 'sensor_relevent');
excel_tool = fullfile(tool_root, 'excel_relevent');

%% === LOAD ROI DEFINITION ===
chan_labels = load(fullfile(source_tool, 'virtual_labels.mat')).virtual_label;

roi_struct = {
    {'S_Hippocampus_L', 'S_Hippocampus_R'},...
    {'S_Frontal_Sup_Medial_L', 'S_Frontal_Sup_Medial_R'},...
    {'S_Precuneus_L', 'S_Precuneus_R'},...
    {'S_Angular_L', 'S_Angular_R'},...
    {'S_Cingulum_Post_L', 'S_Cingulum_Post_R'},...
    {'S_Insula_L', 'S_Insula_R'},...
    {'S_Cingulum_Ant_L', 'S_Cingulum_Ant_R'},...
    {'S_Amygdala_L', 'S_Amygdala_R'},...
    {'S_Thalamus_L', 'S_Thalamus_R'},...
    {'S_Frontal_Mid_L', 'S_Frontal_Mid_R'},...
    {'S_Parietal_Inf_L', 'S_Parietal_Inf_R'},...
    {'S_Frontal_Inf_Tri_L', 'S_Frontal_Inf_Tri_R'},...
    {'S_Frontal_Sup_L', 'S_Frontal_Sup_R'},...
    {'S_Parietal_Sup_L', 'S_Parietal_Sup_R'},...
    {'S_Postcentral_L', 'S_Postcentral_R'},...
    {'S_Precentral_L', 'S_Precentral_R'},...
    {'S_Supp_Motor_Area_L', 'S_Supp_Motor_Area_R'}};

roi_flat = horzcat(roi_struct{:});
roi_idx = find(ismember(chan_labels, roi_flat));
roi_labels = chan_labels(roi_idx);

% Load ROI abbreviations
aal90_abbr = readtable(fullfile(source_tool, 'Node_AAL90.txt'), 'Delimiter','\t','ReadVariableNames',false);
roi_abbr = table2cell(aal90_abbr(roi_idx,6));

% Reorder according to roi_struct
sorted_idx = arrayfun(@(x)find(strcmp(roi_labels, roi_flat{x})), 1:numel(roi_flat));
sorted_idx = sorted_idx(~isnan(sorted_idx));
[~, back_idx] = sort(sorted_idx);
roi_labels_sorted = roi_labels(sorted_idx);
roi_abbr_sorted = roi_abbr(sorted_idx);

%% === LOAD GROUP NETWORK DATA ===
group = {'Anxiety', 'HC'};
demo_tbl = readtable(fullfile(current_path,'Result','Demo_result','result_table_refine.xlsx'), 'VariableNamingRule','preserve');
subject_list = demo_tbl.subject;

group_data = cell(1,numel(group));
group_avg = cell(1,numel(group));

for g = 1:numel(group)
    gname = group{g};
    gpath = fullfile(source_net_path, gname);
    sub_files = dir(fullfile(gpath, '*_te.mat'));
    sub_names = string({sub_files.name}');
    
    net_all = [];
    net_sum = 0;
    valid_n = 0;
    
    for s = 1:numel(sub_names)
        sname = extractBefore(sub_names{s}, '_te.mat');
        if ~ismember(sname, subject_list), continue; end
        mat_path = fullfile(gpath, [sname, '_te.mat']);
        data = load(mat_path).te_mat.data;
        net = data(sorted_idx, sorted_idx);
        net_all = cat(3, net_all, net);
        net_sum = net_sum + net;
        valid_n = valid_n + 1;
    end
    
    if valid_n > 0
        avg = net_sum / valid_n;
        avg_norm = (avg - min(avg(:))) / (max(avg(:)) - min(avg(:)));
        upper = triu(avg_norm,1); val = upper(upper>0);
        thresh = sort(val,'descend');
        cutoff = thresh(round(0.2*numel(thresh)));
        avg_thresh = avg_norm; avg_thresh(avg_thresh<cutoff)=0;
        avg_thresh = avg_thresh / max(avg_thresh(:));
        group_avg{g} = avg_thresh;
        group_data{g} = net_all;
        save(fullfile(graph_data_path, [gname '_group.mat']), 'avg_thresh', '-v6');
    end
end

%% === VISUALIZE GROUP AVG MATRIX ===
for g = 1:numel(group)
    if isempty(group_avg{g}), continue; end
    figure;
    imagesc(group_avg{g}); colorbar;
    title([group{g} ' Group Network']);
    xticks(1:numel(roi_abbr_sorted)); xticklabels(roi_abbr_sorted);
    yticks(1:numel(roi_abbr_sorted)); yticklabels(roi_abbr_sorted);
end

%% === STATISTICAL COMPARISON BETWEEN GROUPS ===
if numel(group)==2
    [n_roi,~,n_sub] = size(group_data{1});
    p_mat = nan(n_roi); t_mat = nan(n_roi);
    for i = 1:n_roi
        for j = 1:n_roi
            d1 = squeeze(group_data{1}(i,j,:));
            d2 = squeeze(group_data{2}(i,j,:));
            [~,p,~,s] = ttest2(d1,d2,'Vartype','unequal');
            p_mat(i,j) = p; t_mat(i,j) = s.tstat;
        end
    end
    sig_links = p_mat<0.05;
    weight_map = (1 - p_mat) .* sig_links;
    norm_map = (weight_map - min(weight_map(:))) / (max(weight_map(:)) - min(weight_map(:)));
    norm_map(isnan(norm_map)) = 0;
    
    save(fullfile(graph_data_path, 'significant_links_group.mat'), 'sig_links','-v6');
    save(fullfile(graph_data_path, 'stronger_in_group1.mat'), 'sig_links', '-v6');
    save(fullfile(graph_data_path, 'stronger_in_group2.mat'), 'sig_links', '-v6');
    
    figure;
    imagesc(sig_links); colorbar;
    title('Significant Group Differences');
    xticks(1:n_roi); xticklabels(roi_abbr_sorted);
    yticks(1:n_roi); yticklabels(roi_abbr_sorted);
    
    % Export significant edges to CSV
    [i,j] = find(sig_links);
    sig_table = table(roi_abbr_sorted(i)', roi_abbr_sorted(j)', t_mat(sig_links), p_mat(sig_links), ...
        'VariableNames', {'Source','Target','T_Value','P_Value'});
    writetable(sig_table, fullfile('YOUR_EXPORT_PATH', 'significant_links_stats.csv')); % <-- Replace this
end

%% === SAVE R DATA ===
save(fullfile(graph_data_path, 'R_data.mat'), 'sorted_idx', 'group_data', '-v7');
