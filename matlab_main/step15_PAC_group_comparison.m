%% CLEAN AND SET PATHS
clc; clear variables;
restoredefaultpath;

% Set base paths
matlab_path = matlabroot;
project_path = 'E:\PAC_network';
addpath(project_path);

% Toolbox paths
toolbox_paths = {
    fullfile(matlab_path, 'toolbox', 'eeglab2021.0'),...
    fullfile(matlab_path, 'toolbox', 'fieldtrip-master')
};
external_folders = {'mne', 'brewermap', 'matplotlib', 'cmocean', 'colorcet', 'gifti', 'spm12'};
cellfun(@addpath, toolbox_paths);
cellfun(@(f) addpath(fullfile(matlab_path, 'toolbox', 'fieldtrip-master', 'external', f)), external_folders);
ft_defaults;

% Project structure
addpath(fullfile(project_path, 'Code'));
addpath(fullfile(project_path, 'Data'));
addpath(fullfile(project_path, 'Tools'));
addpath(genpath(fullfile(project_path, 'Result')));

%% DEFINE INPUT / OUTPUT PATHS
data_method = 'TE';
input_root = fullfile(project_path, 'Result');
output_root = fullfile(project_path, 'Result');

input_paths = {fullfile('Direct_net_result', 'Direct_net_data', data_method)};
output_paths = {...
    fullfile('PAC_result', 'Data'),...
    fullfile('PAC_result', 'Fig')};

f_createOutputFolders(output_root, output_paths);

source_net_path  = fullfile(input_root, input_paths{1});
pac_data_path = fullfile(output_root, output_paths{1});
pac_fig_path = fullfile(output_root, output_paths{2});

% Tool paths
tools_root = fullfile(project_path, 'Tools');
source_tool_path = fullfile(tools_root, 'souce_relevent');
excel_tool_path = fullfile(tools_root, 'excel_relevent');
sensor_tool_path = fullfile(tools_root, 'sensor_relevent');

%% LOAD LABELS AND ROI DEFINITIONS
channels_labels = load(fullfile(source_tool_path, 'virtual_labels.mat')).virtual_label;

% Define selected AAL90 ROIs for multiple networks
regions_AAL90_part = {
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

% Flatten and match to channel labels
all_regions_part = horzcat(regions_AAL90_part{:});
roi_indices = unique(cell2mat(cellfun(@(r)find(strcmp(channels_labels,r)), all_regions_part, 'UniformOutput', false)));
roi_names = channels_labels(roi_indices);

% Load AAL90 abbreviations
aal90_abbr = readtable(fullfile(source_tool_path, 'Node_AAL90.txt'), 'Delimiter', '\t', 'ReadVariableNames', false);
roi_names_abbr = table2cell(aal90_abbr(roi_indices, 6));

% Reorder to match the defined structure
sorted_indices = cellfun(@(r)find(strcmp(roi_names,r)), all_regions_part);
sorted_roi_names = roi_names(sorted_indices);
sorted_roi_names_abbr = roi_names_abbr(sorted_indices);

%% LOAD GROUP PAC DATA
group = {'Anxiety', 'HC'};
group_pac_data = cell(1, numel(group));
input_dir = fullfile(output_root, 'PAC_result', 'Fig');

for g = 1:numel(group)
    T = readtable(fullfile(input_dir, [group{g} '_pac_data.xlsx']), 'ReadRowNames', true);
    group_pac_data{g} = table2array(T);
    fprintf('%s group data loaded. Size: %d x %d\n', group{g}, size(group_pac_data{g},1), size(group_pac_data{g},2));
end

%% DRAW DEFAULT-ORDER BAR PLOTS FOR R
save_dir = fullfile(project_path, 'R_code', 'PAC_strength', 'Data');
group_names = {'Anxiety', 'HC'};
bar_color = [0.2 0.6 0.8];

for g = 1:2
    data = group_pac_data{g};
    mean_vals = mean(data, 2);
    std_vals = std(data, 0, 2);
    p = polyfit(1:34, mean_vals', 3);
    y_fit = polyval(p, 1:34);

    figure('Color','w','Position',[100,100,1200,500]);
    bar(mean_vals,'FaceColor',bar_color,'EdgeColor','none'); hold on;
    errorbar(1:34,mean_vals,std_vals,'k.','LineWidth',1.2);
    plot(1:34,y_fit,'r--','LineWidth',2);

    xticks(1:34);
    xticklabels(sorted_roi_names_abbr);
    xtickangle(45); xlim([0.5 34.5]);
    ylim([min(mean_vals-std_vals)-0.002, max(mean_vals+std_vals)+0.002]);
    ylabel('PAC Mean');
    title([group_names{g} ' Group: Mean PAC (Default Order)']);
    legend({'Mean','Std','Fit'},'Location','northeast');
    grid on; box on; set(gca,'FontSize',11,'LineWidth',1.2);

    % Export for R
    roi_names = string(sorted_roi_names_abbr(:));
    export_table = table(roi_names, mean_vals, std_vals, 'VariableNames', {'ROI','Mean','Std'});
    writetable(export_table, fullfile(save_dir, [group_names{g} '_PAC_default_for_R.csv']));
end

%% RANKSUM TEST FOR GROUP DIFFERENCE (FOR RAINPLOT)
group_avg_pac = cellfun(@(x) mean(x,2,'omitnan'), group_pac_data, 'UniformOutput', false);
p_values = arrayfun(@(i) ranksum(group_pac_data{1}(i,:), group_pac_data{2}(i,:)), 1:34);
sig_indices = find(p_values < 0.05);

% Plot group difference on significant ROI
figure('Position',[100,100,800,400]);
colors = {[0 0.447 0.741], [0.85 0.325 0.098]};
dot_size = 50; jitter_amount = 0.1;

sig_pac_values = [group_avg_pac{1}(sig_indices), group_avg_pac{2}(sig_indices)];
sig_labels = sorted_roi_names_abbr(sig_indices);

b = bar(sig_pac_values,'grouped','BarWidth',0.4); hold on;
for g = 1:2, b(g).FaceColor = colors{g}; end

for g = 1:2
    for i = 1:length(sig_indices)
        scatter(i + (g-1)*0.2 + jitter_amount*randn(size(group_pac_data{g},2),1), ...
                group_pac_data{g}(sig_indices(i),:)', ...
                dot_size, colors{g}, 'filled');
    end
end

xticks(1:length(sig_labels)); xticklabels(sig_labels); xtickangle(45);
ylabel('Average PAC Value');
title('Significant ROI PAC Comparison'); grid on;
set(gca,'FontSize',12,'LineWidth',1.5);

for i = 1:length(sig_labels)
    text(i, max(sig_pac_values(i,:)) + 0.08, '*', 'FontSize', 20, 'HorizontalAlignment', 'center');
end

saveas(gcf, fullfile(pac_fig_path, 'Significant_PAC_Group_Comparison.png'));

%% EXPORT CSV FOR RAINPLOT
% Flatten PAC data for Rainplot
for g = 1:2
    data = group_pac_data{g}(sig_indices,:);
    subjects = repelem(1:size(data,2), size(data,1))';
    roi_labels = repmat(sorted_roi_names_abbr(sig_indices)', 1, size(data,2)); roi_labels = roi_labels(:);
    group_labels = repmat({group{g}}, numel(subjects), 1);
    pac_vals = data(:);
    T = table(subjects, roi_labels, group_labels, pac_vals, 'VariableNames', {'Subject','ROIs','Group','PAC_Value'});
    writetable(T, fullfile(pac_fig_path, ['subject_pac_values_' group{g} '.csv']));
end
