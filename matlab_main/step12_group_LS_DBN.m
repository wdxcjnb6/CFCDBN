
clc; clear; restoredefaultpath;

% === PATH INITIALIZATION ===
matlab_path = matlabroot;
project_root = 'PROJECT_ROOT_PATH'; % Replace with your root
current_path = fullfile(project_root, 'PAC_network');
addpath(current_path);

% === TOOLBOXES ===
toolbox_paths = {
    fullfile(matlab_path, 'toolbox', 'eeglab2021.0');
    fullfile(matlab_path, 'toolbox', 'fieldtrip-master')
};
external_subfolders = {'mne','brewermap','matplotlib','cmocean','colorcet','gifti','spm12'};
ft_ext_path = fullfile(matlab_path, 'toolbox', 'fieldtrip-master', 'external');
cellfun(@addpath, toolbox_paths);
cellfun(@(f)addpath(fullfile(ft_ext_path,f)), external_subfolders);
ft_defaults;

% === PROJECT PATHS ===
addpath(fullfile(current_path, 'Code'));
addpath(fullfile(current_path, 'Data'));
addpath(fullfile(current_path, 'Tools'));
addpath(genpath(fullfile(current_path, 'Result')));

% === DEFINE PATHS ===
method = 'TE';
format = '';
input_root = fullfile(current_path, 'Result');
output_root = input_root;

input_data_path = fullfile(input_root, 'Direct_net_result', 'Direct_net_data', method, format);
output_data_path = fullfile(output_root, 'Direct_net_result', 'Direct_net_figure', 'Data');
output_r_path = fullfile(output_root, 'R_result', 'large_direct_net');

f_createOutputFolders(output_root, { ...
    'Direct_net_result\Direct_net_figure\Data\', ...
    'R_result\large_direct_net\' ...
});

% === TOOLS ===
tool_root = fullfile(current_path, 'Tools');
source_tool_path = fullfile(tool_root, 'souce_relevent');
excel_tool_path = fullfile(tool_root, 'excel_relevent');
sensor_tool_path = fullfile(tool_root, 'sensor_relevent');

% === LOAD ROI CHANNEL LABELS ===
labels = load(fullfile(source_tool_path, 'virtual_labels.mat')).virtual_label;

% === SELECT AAL90 REGIONS FOR FUNCTIONAL NETWORKS ===
regions_network = {
    {'S_Hippocampus_L', 'S_Hippocampus_R'}, ...
    {'S_Frontal_Sup_Medial_L', 'S_Frontal_Sup_Medial_R'}, ...
    {'S_Precuneus_L', 'S_Precuneus_R'}, ...
    {'S_Angular_L', 'S_Angular_R'}, ...
    {'S_Cingulum_Post_L', 'S_Cingulum_Post_R'}, ...
    {'S_Insula_L', 'S_Insula_R'}, ...
    {'S_Cingulum_Ant_L', 'S_Cingulum_Ant_R'}, ...
    {'S_Amygdala_L', 'S_Amygdala_R'}, ...
    {'S_Thalamus_L', 'S_Thalamus_R'}, ...
    {'S_Frontal_Mid_L', 'S_Frontal_Mid_R'}, ...
    {'S_Parietal_Inf_L', 'S_Parietal_Inf_R'}, ...
    {'S_Frontal_Inf_Tri_L', 'S_Frontal_Inf_Tri_R'}, ...
    {'S_Frontal_Sup_L', 'S_Frontal_Sup_R'}, ...
    {'S_Parietal_Sup_L', 'S_Parietal_Sup_R'}, ...
    {'S_Postcentral_L', 'S_Postcentral_R'}, ...
    {'S_Precentral_L', 'S_Precentral_R'}, ...
    {'S_Supp_Motor_Area_L', 'S_Supp_Motor_Area_R'}
};

all_regions = horzcat(regions_network{:});
roi_indices = cellfun(@(x)find(strcmp(labels,x)), all_regions);
roi_indices = unique([roi_indices{:}]);
roi_labels = labels(roi_indices);

% === GET REGION ABBREVIATIONS ===
aal90_table = readtable('Node_AAL90.txt','Delimiter','\t','ReadVariableNames',false);
roi_abbr = table2cell(aal90_table(roi_indices,6));

% === SORT ROI LABELS ===
[~, sort_order] = sort(roi_indices);
sorted_labels = roi_labels(sort_order);
sorted_abbr = roi_abbr(sort_order);

% === GROUPS & NETWORK SETUP ===
group_list = {'Anxiety', 'HC'};
network_names = {'DMN', 'SN', 'FPN', 'DAN', 'SMN'};
network_sizes = [10, 8, 6, 4, 6];

% === LOAD SUBJECTS ===
subject_table = readtable(fullfile(current_path, 'Result', 'Demo_result', 'result_table_refine.xlsx'));
subject_names = subject_table.subject;

% === INITIALIZE STORAGE ===
group_network_data = cell(1, numel(group_list));

% === BUILD INDIVIDUAL 5x5 NETWORK MATRICES ===
for g = 1:numel(group_list)
    gname = group_list{g};
    gpath = fullfile(input_data_path, gname);
    files = dir(fullfile(gpath, '*_te.mat'));
    files = string({files.name});
    g_subnets = [];

    for f = 1:numel(files)
        subj = extractBefore(files(f), '_te.mat');
        if ~ismember(subj, subject_names), continue; end
        mat = load(fullfile(gpath, files(f))).te_mat.data;
        mat = mat(roi_indices, roi_indices);

        % === BUILD 5x5 FUNCTIONAL NETWORK AVERAGE ===
        net5x5 = zeros(5);
        for i = 1:5
            si = sum(network_sizes(1:i-1))+(1:network_sizes(i));
            for j = 1:5
                ti = sum(network_sizes(1:j-1))+(1:network_sizes(j));
                sub_net = mat(si, ti);
                net5x5(i,j) = mean(sub_net(sub_net>0));
            end
        end
        g_subnets = cat(3, g_subnets, net5x5);
    end
    g_subnets(isnan(g_subnets)) = 0;
    group_network_data{g} = g_subnets;
end

% === GROUP-LEVEL DIFFERENCE TEST ===
diff_avg = zeros(5); tvals = zeros(5); pvals = zeros(5);
for i = 1:5
    for j = 1:5
        data1 = group_network_data{1}(:,i,j);
        data2 = group_network_data{2}(:,i,j);
        [~,p,~,stats] = ttest2(data1, data2);
        diff_avg(i,j) = mean(data1 - data2);
        tvals(i,j) = stats.tstat;
        pvals(i,j) = p;
    end
end

% === SAVE RESULTS ===
r_out = fullfile(current_path, 'R_code', 'net5x5');
if ~exist(r_out, 'dir'), mkdir(r_out); end
csvwrite(fullfile(r_out, 'group_diff_avg.csv'), diff_avg);
csvwrite(fullfile(r_out, 't_stats.csv'), tvals);
csvwrite(fullfile(r_out, 'p_values.csv'), pvals);

% === PLOT DIFFERENCE MATRIX ===
figure;
imagesc(diff_avg); colorbar;
xticks(1:5); yticks(1:5);
xticklabels(network_names); yticklabels(network_names);
xlabel('Target'); ylabel('Source');
title('Group Difference (Source → Target)');
saveas(gcf, fullfile(r_out, 'group_diff_avg.png'));

% === PLOT EACH GROUP AVERAGE ===
for g = 1:numel(group_list)
    avg_mat = squeeze(mean(group_network_data{g},3));
    figure;
    imagesc(avg_mat); colorbar;
    xticks(1:5); yticks(1:5);
    xticklabels(network_names); yticklabels(network_names);
    title([group_list{g}, ' Average Network']);
    saveas(gcf, fullfile(r_out, [group_list{g}, '_avg_5x5_network.png']));
    writematrix(avg_mat, fullfile(r_out, [group_list{g}, '_avg_5x5_network.csv']));
end
