%% CLEAN AND SET PATHS
clc; clear variables;
restoredefaultpath;

% Base path settings
matlab_path = matlabroot;
matlab_project_path = 'E:\';
current_project_path = fullfile(matlab_project_path, 'PAC_network');
addpath(current_project_path);

% Toolbox paths
toolbox_paths = {
    fullfile(matlab_path, 'toolbox', 'eeglab2021.0'),...
    fullfile(matlab_path, 'toolbox', 'fieldtrip-master')
};
external_subfolders = {'mne', 'brewermap', 'matplotlib', 'cmocean', 'colorcet', 'gifti', 'spm12'};
fieldtrip_external_path = fullfile(matlab_path, 'toolbox', 'fieldtrip-master', 'external');
cellfun(@addpath, toolbox_paths);
cellfun(@(subfolder) addpath(fullfile(fieldtrip_external_path, subfolder)), external_subfolders);

% Initialize FieldTrip
ft_defaults;
addpath(fullfile(current_project_path, 'Code'));
addpath(fullfile(current_project_path, 'Data'));
addpath(fullfile(current_project_path, 'Tools'));
addpath(genpath(fullfile(current_project_path, 'Result')));

%% Define paths
% Input/output directories
data_method = 'TE';
input_paths = {['Direct_net_result\Direct_net_data\' data_method]};
input_root = fullfile(current_project_path, 'Result');
output_paths = {'PAC_result\Data\', 'PAC_result\Fig\'};
output_root = fullfile(current_project_path, 'Result');
f_createOutputFolders(output_root, output_paths);
source_net_path  = fullfile(input_root, input_paths{1});
pac_data_path = fullfile(output_root, output_paths[1]);
pac_fig_path = fullfile(output_root, output_paths[2]);

% Tool directories
input_tool_root = fullfile(current_project_path, 'Tools');
source_tool_path = fullfile(input_tool_root, 'souce_relevent\');
excel_tool_path = fullfile(input_tool_root, 'excel_relevent\');
sensor_tool_path = fullfile(input_tool_root, 'sensor_relevent\');

%% Load channel labels and define ROI regions
channels_labels = load(fullfile(source_tool_path, 'virtual_labels.mat')).virtual_label;

regions_AAL90_part = {
    {'S_Hippocampus_L', 'S_Hippocampus_R'}, 
    {'S_Frontal_Sup_Medial_L', 'S_Frontal_Sup_Medial_R'}, 
    {'S_Precuneus_L', 'S_Precuneus_R'}, 
    {'S_Angular_L', 'S_Angular_R'}, 
    {'S_Cingulum_Post_L', 'S_Cingulum_Post_R'},
    {'S_Insula_L', 'S_Insula_R'}, 
    {'S_Cingulum_Ant_L', 'S_Cingulum_Ant_R'}, 
    {'S_Amygdala_L', 'S_Amygdala_R'}, 
    {'S_Thalamus_L', 'S_Thalamus_R'},
    {'S_Frontal_Mid_L', 'S_Frontal_Mid_R'}, 
    {'S_Parietal_Inf_L', 'S_Parietal_Inf_R'}, 
    {'S_Frontal_Inf_Tri_L', 'S_Frontal_Inf_Tri_R'},
    {'S_Frontal_Sup_L', 'S_Frontal_Sup_R'}, 
    {'S_Parietal_Sup_L', 'S_Parietal_Sup_R'},
    {'S_Postcentral_L', 'S_Postcentral_R'}, 
    {'S_Precentral_L', 'S_Precentral_R'}, 
    {'S_Supp_Motor_Area_L', 'S_Supp_Motor_Area_R'}
};

all_regions_part = horzcat(regions_AAL90_part{:});
roi_indices = [];
for i = 1:length(all_regions_part)
    index = find(strcmp(channels_labels, all_regions_part{i}));
    if ~isempty(index)
        roi_indices = [roi_indices, index];
    end
end
roi_indices = unique(roi_indices);
roi_names = channels_labels(roi_indices);

aal90_abbr = readtable('Node_AAL90.txt', 'Delimiter', '\t', 'ReadVariableNames', false);
roi_names_abbr = table2cell(aal90_abbr(roi_indices, 6));

sorted_indices = zeros(size(all_regions_part));
for i = 1:length(all_regions_part)
    index = find(strcmp(roi_names, all_regions_part{i}));
    sorted_indices(i) = isempty(index) * NaN + ~isempty(index) * index;
end
sorted_indices = sorted_indices(~isnan(sorted_indices));
[~, sorted_indices_back] = sort(sorted_indices);
sorted_roi_names = roi_names(sorted_indices);
sorted_roi_names_abbr = roi_names_abbr(sorted_indices);

disp('Sorted ROI Names:');
disp(sorted_roi_names);
disp('Sorted ROI Names Abbreviations:');
disp(sorted_roi_names_abbr);

%% Extract left/right hemisphere indices
left_brain_indices = find(contains(sorted_roi_names_abbr, '.L'));
right_brain_indices = find(contains(sorted_roi_names_abbr, '.R'));

%% Define input file path and load PAC data for both groups
input_dir = 'E:\PAC_network\Result\PAC_result\Fig\';
group = {'Anxiety', 'HC'};
group_pac_data = cell(1, numel(group));

for group_n = 1:numel(group)
    group_name = group{group_n};
    input_file = fullfile(input_dir, [group_name, '_pac_data.xlsx']);
    T = readtable(input_file, 'ReadRowNames', true);
    data_matrix = table2array(T);
    group_pac_data{group_n} = data_matrix;
    fprintf('%s group data loaded, size: %d×%d\n', group_name, size(data_matrix,1), size(data_matrix,2));
end

%% Compute average PAC for each group (across subjects)
group_avg_pac = cellfun(@(x) mean(x, 2, 'omitnan'), group_pac_data, 'UniformOutput', false);

%% --- Plot: Sorted PAC Bar Chart with Polynomial Fit ---
group_names = {'Anxiety', 'HC'};
bar_color = [0.2 0.6 0.8];

for group_i = 1:2
    data_group = squeeze(group_pac_data{group_i});
    mean_vals = mean(data_group, 2);
    std_vals  = std(data_group, 0, 2);
    [mean_sorted, idx_sort] = sort(mean_vals, 'descend');
    std_sorted = std_vals(idx_sort);
    roi_sorted = sorted_roi_names_abbr(idx_sort);

    p = polyfit(1:34, mean_sorted', 3);
    y_fit = polyval(p, 1:34);

    figure('Color', 'w', 'Position', [100, 100, 1200, 500]);
    bar(mean_sorted, 'FaceColor', bar_color, 'EdgeColor', 'none'); hold on;
    errorbar(1:34, mean_sorted, std_sorted, 'k.', 'LineWidth', 1.2);
    plot(1:34, y_fit, 'r--', 'LineWidth', 2);

    xticks(1:34);
    xticklabels(roi_sorted);
    xtickangle(45);
    xlim([0.5 34.5]);
    ylim([min(mean_sorted - std_sorted)-0.002, max(mean_sorted + std_sorted)+0.002]);
    ylabel('Mean PAC');
    title([group_names{group_i} ' group: PAC (sorted) and fit curve'], 'FontSize', 14);
    legend({'Mean', 'Std', '3rd-order Fit'}, 'Location', 'northeast');
    grid on; box on;
    set(gca, 'FontSize', 11, 'LineWidth', 1.2);

    % Export to CSV for R analysis
    save_dir = 'E:\PAC_network\R_code\PAC_strength\Data';
    roi_sorted = string(roi_sorted);
    export_table = table(roi_sorted(:), mean_sorted(:), std_sorted(:), ...
        'VariableNames', {'ROI', 'Mean', 'Std'});
    save_path = fullfile(save_dir, [group_names{group_i} '_PAC_sorted_for_R.csv']);
    writetable(export_table, save_path);
end

%% --- Plot: Unsorted PAC Bar Chart with Polynomial Fit ---
for group_i = 1:2
    data_group = squeeze(group_pac_data{group_i});
    mean_vals = mean(data_group, 2);
    std_vals  = std(data_group, 0, 2);

    p = polyfit(1:34, mean_vals', 3);
    y_fit = polyval(p, 1:34);

    figure('Color', 'w', 'Position', [100, 100, 1200, 500]);
    bar(mean_vals, 'FaceColor', bar_color, 'EdgeColor', 'none'); hold on;
    errorbar(1:34, mean_vals, std_vals, 'k.', 'LineWidth', 1.2);
    plot(1:34, y_fit, 'r--', 'LineWidth', 2);

    xticks(1:34);
    xticklabels(sorted_roi_names_abbr);
    xtickangle(45);
    xlim([0.5 34.5]);
    ylim([min(mean_vals - std_vals)-0.002, max(mean_vals + std_vals)+0.002]);
    ylabel('Mean PAC');
    title([group_names{group_i} ' group: PAC (unsorted) and fit curve'], 'FontSize', 14);
    legend({'Mean', 'Std', '3rd-order Fit'}, 'Location', 'northeast');
    grid on; box on;
    set(gca, 'FontSize', 11, 'LineWidth', 1.2);

    % Export to CSV for R
    roi_names = string(sorted_roi_names_abbr(:));
    export_table = table(roi_names, mean_vals(:), std_vals(:), ...
        'VariableNames', {'ROI', 'Mean', 'Std'});
    save_path = fullfile(save_dir, [group_names{group_i} '_PAC_default_for_R.csv']);
    writetable(export_table, save_path);
end

%% Save data to .mat
data_filename = fullfile(pac_data_path, 'PAC_Data.mat');
save(data_filename, 'group_pac_data', 'sorted_roi_names_abbr');
fprintf('PAC data saved to %s\n', data_filename);



%% --- Plot: Left vs. Right Hemisphere PAC Comparison ---
% Recompute group-level PAC (across subjects, if not already done)
group_avg_pac = cellfun(@(x) mean(x, 3, 'omitnan'), group_pac_data, 'UniformOutput', false);

% Plot settings
colors = {[0 0.447 0.741], [0.85 0.325 0.098]};  % Blue (Anxiety) & Orange (HC)
dot_size = 30;
jitter_amount = 0.1;
bar_width = 0.4;

% Create figure
figure('Position', [100, 100, 1200, 500]);

% --- Left Hemisphere ---
subplot(1, 2, 1);
bar_data_left = [group_avg_pac{1}(left_brain_indices), group_avg_pac{2}(left_brain_indices)];
b = bar(bar_data_left, 'grouped', 'BarWidth', bar_width); hold on;

% Set bar colors
for g = 1:2
    b(g).FaceColor = colors{g};
end

% Overlay scatter dots
for g = 1:2
    for i = 1:length(left_brain_indices)
        scatter(repmat(i + (g-1)*0.2, size(group_pac_data{g}, 3), 1), ...
                squeeze(group_pac_data{g}(left_brain_indices(i), :)), ...
                dot_size, colors{g}, 'filled', ...
                'jitter', 'on', 'jitterAmount', jitter_amount);
    end
end

hold off;
legend({'Anxiety', 'HC'}, 'Location', 'northeast');
xticks(1:length(left_brain_indices));
xticklabels(sorted_roi_names_abbr(left_brain_indices));
xtickangle(45);
ylabel('Average PAC Value');
title('Left Hemisphere PAC Comparison');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

% --- Right Hemisphere ---
subplot(1, 2, 2);
bar_data_right = [group_avg_pac{1}(right_brain_indices), group_avg_pac{2}(right_brain_indices)];
b = bar(bar_data_right, 'grouped', 'BarWidth', bar_width); hold on;

% Set bar colors
for g = 1:2
    b(g).FaceColor = colors{g};
end

% Overlay scatter dots
for g = 1:2
    for i = 1:length(right_brain_indices)
        scatter(repmat(i + (g-1)*0.2, size(group_pac_data{g}, 3), 1), ...
                squeeze(group_pac_data{g}(right_brain_indices(i), :)), ...
                dot_size, colors{g}, 'filled', ...
                'jitter', 'on', 'jitterAmount', jitter_amount);
    end
end

hold off;
legend({'Anxiety', 'HC'}, 'Location', 'northeast');
xticks(1:length(right_brain_indices));
xticklabels(sorted_roi_names_abbr(right_brain_indices));
xtickangle(45);
ylabel('Average PAC Value');
title('Right Hemisphere PAC Comparison');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

% Save figure
fig_filename = fullfile(pac_fig_path, 'PAC_Group_Comparison_LeftRight.png');
saveas(gcf, fig_filename);
fprintf(' PAC group comparison figure saved: %s\n', fig_filename);


%% --- Compute Significant ROIs Based on Rank-Sum Test (Uncorrected) ---
p_values = zeros(1, length(sorted_roi_names_abbr));
for i = 1:length(sorted_roi_names_abbr)
    p_values(i) = ranksum(squeeze(group_pac_data{1}(i, :)), squeeze(group_pac_data{2}(i, :)));
end

sig_indices = find(p_values < 0.05);
if isempty(sig_indices)
    fprintf('No significant ROIs found (p < 0.05). Analysis stopped.\n');
    return;
end

%% --- Plot Group Comparison for Significant ROIs ---
figure('Position', [100, 100, 800, 400]);
colors = {[0 0.447 0.741], [0.85 0.325 0.098]};
dot_size = 50;
jitter_amount = 0.1;

sig_pac_values = [group_avg_pac{1}(sig_indices), group_avg_pac{2}(sig_indices)];
sig_labels = sorted_roi_names_abbr(sig_indices);

bar_data = [sig_pac_values(:,1), sig_pac_values(:,2)];
b = bar(bar_data, 'grouped', 'BarWidth', 0.4); hold on;

for g = 1:2
    b(g).FaceColor = colors{g};
end

num_sig_roi = length(sig_labels);
for g = 1:2
    for i = 1:num_sig_roi
        scatter(repmat(i + (g-1)*0.2, size(group_pac_data{g}, 3), 1), ...
                squeeze(group_pac_data{g}(sig_indices(i), :)), ...
                dot_size, colors{g}, 'filled', 'jitter', 'on', 'jitterAmount', jitter_amount);
    end
end

xticks(1:num_sig_roi);
xticklabels(sig_labels);
xtickangle(45);
ylabel('Average PAC Value');
title('Group Comparison: Significant ROIs');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

% Add significance markers
for i = 1:num_sig_roi
    text(i, max(bar_data(i,:)) + 0.08, '*', 'FontSize', 20, 'Color', 'k', 'HorizontalAlignment', 'center');
end

fig_filename = fullfile(pac_fig_path, 'Significant_PAC_Group_Comparison.png');
saveas(gcf, fig_filename);
fprintf('Significant ROI group comparison figure saved: %s\n', fig_filename);

%% --- Export Subject-Level PAC Values for Raincloud Plot in R ---
% Reshape PAC data into long-format tables for HC and Anxiety
hc_data = squeeze(group_pac_data{1}(sig_indices, 1, :));
anxiety_data = squeeze(group_pac_data{2}(sig_indices, 1, :));

hc_subjects = repmat(1:size(hc_data, 2), length(sig_indices), 1);
hc_labels = repmat(sorted_roi_names_abbr(sig_indices)', size(hc_data, 2), 1);
hc_group_labels = repmat({'HC'}, numel(hc_data), 1);
hc_pac_values = hc_data(:);

hc_table = table(hc_subjects(:), hc_labels(:), hc_group_labels, hc_pac_values, ...
    'VariableNames', {'Subject', 'ROIs', 'Group', 'PAC_Value'});
writetable(hc_table, fullfile(pac_fig_path, 'subject_pac_values_HC.csv'));
fprintf('HC subject-level PAC data exported for raincloud plot.\n');

anxiety_subjects = repmat(1:size(anxiety_data, 2), length(sig_indices), 1);
anxiety_labels = repmat(sorted_roi_names_abbr(sig_indices)', size(anxiety_data, 2), 1);
anxiety_group_labels = repmat({'Anxiety'}, numel(anxiety_data), 1);
anxiety_pac_values = anxiety_data(:);

anxiety_table = table(anxiety_subjects(:), anxiety_labels(:), anxiety_group_labels, anxiety_pac_values, ...
    'VariableNames', {'Subject', 'ROIs', 'Group', 'PAC_Value'});
writetable(anxiety_table, fullfile(pac_fig_path, 'subject_pac_values_Anxiety.csv'));
fprintf(' Anxiety subject-level PAC data exported for raincloud plot.\n');

%% --- Plot Standard Radar Charts for Left/Right Hemispheres ---

% Define functional networks and their ROI indices
network_dict = struct('DMN', 1:10, 'SN', 11:18, 'FPN', 19:24, 'DAN', 25:28, 'SMN', 29:34);
networks = fieldnames(network_dict);
num_networks = length(networks);

% Identify left/right hemisphere ROI indices
left_brain_indices = find(contains(sorted_roi_names_abbr, '.L'));
right_brain_indices = find(contains(sorted_roi_names_abbr, '.R'));

% Initialize radar chart data
left_data_radar = NaN(num_networks, 2);  % [Anxiety, HC]
right_data_radar = NaN(num_networks, 2);
labels_radar = networks;

for i = 1:num_networks
    network_range = network_dict.(networks{i});
    left_indices = intersect(left_brain_indices, network_range);
    right_indices = intersect(right_brain_indices, network_range);
    
    if ~isempty(left_indices)
        left_data_radar(i, 1) = mean(group_avg_pac{1}(left_indices), 'omitnan');
        left_data_radar(i, 2) = mean(group_avg_pac{2}(left_indices), 'omitnan');
    end
    if ~isempty(right_indices)
        right_data_radar(i, 1) = mean(group_avg_pac{1}(right_indices), 'omitnan');
        right_data_radar(i, 2) = mean(group_avg_pac{2}(right_indices), 'omitnan');
    end
end

% Determine max scale for radar
max_radar = max([max(left_data_radar, [], 'all'), max(right_data_radar, [], 'all')]);
r_ticks = linspace(0, max_radar, 5);

% --- Left Hemisphere Radar ---
figure;
theta = linspace(0, 2*pi, num_networks + 1);
left_data_radar = [left_data_radar; left_data_radar(1,:)];  % close the loop
p1 = polarplot(theta, left_data_radar(:, 1), '-o', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0 0.447 0.741]);
hold on;
p2 = polarplot(theta, left_data_radar(:, 2), '-o', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0.85 0.325 0.098]);
hold off;

legend([p1, p2], {'Anxiety', 'HC'}, 'Location', 'northeast');
title('Left Hemisphere PAC Radar Plot');
ax = gca;
ax.ThetaTick = rad2deg(theta(1:end-1));
ax.ThetaTickLabel = labels_radar;
ax.RLim = [0 max_radar];
ax.RTick = r_ticks;
grid on;

% Save figure
fig_filename_left = fullfile(pac_fig_path, 'PAC_Left_Hemisphere_Radar_Plot.png');
saveas(gcf, fig_filename_left);
fprintf('✅ Left hemisphere radar chart saved: %s\n', fig_filename_left);

% --- Right Hemisphere Radar ---
figure;
right_data_radar = [right_data_radar; right_data_radar(1,:)];  % close the loop
p1 = polarplot(theta, right_data_radar(:, 1), '-o', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0 0.447 0.741]);
hold on;
p2 = polarplot(theta, right_data_radar(:, 2), '-o', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0.85 0.325 0.098]);
hold off;

legend([p1, p2], {'Anxiety', 'HC'}, 'Location', 'northeast');
title('Right Hemisphere PAC Radar Plot');
ax = gca;
ax.ThetaTick = rad2deg(theta(1:end-1));
ax.ThetaTickLabel = labels_radar;
ax.RLim = [0 max_radar];
ax.RTick = r_ticks;
grid on;

fig_filename_right = fullfile(pac_fig_path, 'PAC_Right_Hemisphere_Radar_Plot.png');
saveas(gcf, fig_filename_right);
fprintf('✅ Right hemisphere radar chart saved: %s\n', fig_filename_right);

%% --- Export Mean & Std for Radar Plots to CSV ---

% Extract mean and std PAC values per hemisphere per network
network_names = labels_radar(:);
left_mean_anxiety  = left_data_radar(1:end-1, 1);
left_mean_hc       = left_data_radar(1:end-1, 2);
right_mean_anxiety = right_data_radar(1:end-1, 1);
right_mean_hc      = right_data_radar(1:end-1, 2);

left_std_anxiety = NaN(num_networks, 1);
left_std_hc = NaN(num_networks, 1);
right_std_anxiety = NaN(num_networks, 1);
right_std_hc = NaN(num_networks, 1);

for i = 1:num_networks
    net_idx = network_dict.(networks{i});
    left_roi = intersect(left_brain_indices, net_idx);
    right_roi = intersect(right_brain_indices, net_idx);
    
    if ~isempty(left_roi)
        left_std_anxiety(i) = std(group_avg_pac{1}(left_roi), 'omitnan');
        left_std_hc(i) = std(group_avg_pac{2}(left_roi), 'omitnan');
    end
    if ~isempty(right_roi)
        right_std_anxiety(i) = std(group_avg_pac{1}(right_roi), 'omitnan');
        right_std_hc(i) = std(group_avg_pac{2}(right_roi), 'omitnan');
    end
end

radar_table = table(network_names, ...
    left_mean_anxiety, left_std_anxiety, ...
    left_mean_hc, left_std_hc, ...
    right_mean_anxiety, right_std_anxiety, ...
    right_mean_hc, right_std_hc);

csv_path = 'E:\PAC_network\R_code\leidatu\left_right_radar_data.csv';
writetable(radar_table, csv_path);
fprintf(' Radar chart data saved to: %s\n', csv_path);

%% --- Advanced Radar Charts (using radarChart object) ---

% Functional network definition
network_dict = struct('DMN', 1:10, 'SN', 11:18, 'FPN', 19:24, 'DAN', 25:28, 'SMN', 29:34);
networks = fieldnames(network_dict);
num_networks = length(networks);

% Hemisphere index assignment
left_brain_indices = find(contains(sorted_roi_names_abbr, '.L'));
right_brain_indices = find(contains(sorted_roi_names_abbr, '.R'));

% Initialize radar chart data
left_data_radar = NaN(num_networks, 2);
right_data_radar = NaN(num_networks, 2);
labels_radar = networks;

for i = 1:num_networks
    network_range = network_dict.(networks{i});
    left_indices = intersect(left_brain_indices, network_range);
    right_indices = intersect(right_brain_indices, network_range);
    
    if ~isempty(left_indices)
        left_data_radar(i, 1) = mean(group_avg_pac{1}(left_indices), 'omitnan');
        left_data_radar(i, 2) = mean(group_avg_pac{2}(left_indices), 'omitnan');
    end
    if ~isempty(right_indices)
        right_data_radar(i, 1) = mean(group_avg_pac{1}(right_indices), 'omitnan');
        right_data_radar(i, 2) = mean(group_avg_pac{2}(right_indices), 'omitnan');
    end
end

% Color and transparency settings
colorList = [78 101 155; 138 140 191] ./ 255;
alpha_val = 0.3;

%% --- Left Hemisphere radarChart drawing ---
X_left = left_data_radar';  % rows = group, columns = networks
figure;
RC_left = radarChart(X_left, 'Type', 'Line');
RC_left.PropName = networks;
RC_left.ClassName = {'HC', 'Anxiety'};
RC_left = RC_left.draw();
RC_left.legend();

for n = 1:RC_left.ClassNum
    if strcmp(RC_left.Type, 'Line')
        RC_left.setPatchN(n, ...
            'Color', colorList(n,:), ...
            'LineWidth', 2, ...
            'Marker', 'o', ...
            'MarkerFaceColor', colorList(n,:), ...
            'MarkerSize', 6);
    elseif strcmp(RC_left.Type, 'Patch') || strcmp(RC_left.Type, 'Both')
        RC_left.setPatchN(n, ...
            'EdgeColor', colorList(n,:), ...
            'LineWidth', 1.8, ...
            'FaceAlpha', alpha_val, ...
            'MarkerFaceColor', colorList(n,:), ...
            'MarkerSize', 6);
    end
end

% Styling
RC_left.setThetaTick('LineWidth', 2, 'Color', [0.6, 0.6, 0.8]);
RC_left.setRTick('LineWidth', 1.5, 'Color', [0.8, 0.6, 0.6]);
RC_left.setPropLabel('FontSize', 15, 'FontName', 'Times New Roman', 'Color', [0 0 0]);
RC_left.setRLabel('FontSize', 15, 'FontName', 'Times New Roman', 'Color', [0.8 0 0]);
RC_left.setBkg('FaceColor', 'none');

% Save figure
filename_left = fullfile(pac_fig_path, "PAC_Left_Hemisphere_radarChart.tiff");
print(gcf, '-dtiff', '-r600', filename_left);
fprintf(" Advanced left hemisphere radar chart saved: %s\n", filename_left);

%% --- Right Hemisphere radarChart drawing ---
X_right = right_data_radar';
figure;
RC_right = radarChart(X_right, 'Type', 'Line');
RC_right.PropName = networks;
RC_right.ClassName = {'HC', 'Anxiety'};
RC_right = RC_right.draw();
RC_right.legend();

for n = 1:RC_right.ClassNum
    if strcmp(RC_right.Type, 'Line')
        RC_right.setPatchN(n, ...
            'Color', colorList(n,:), ...
            'LineWidth', 2, ...
            'Marker', 'o', ...
            'MarkerFaceColor', colorList(n,:), ...
            'MarkerSize', 6);
    elseif strcmp(RC_right.Type, 'Patch') || strcmp(RC_right.Type, 'Both')
        RC_right.setPatchN(n, ...
            'EdgeColor', colorList(n,:), ...
            'LineWidth', 1.8, ...
            'FaceAlpha', alpha_val, ...
            'MarkerFaceColor', colorList(n,:), ...
            'MarkerSize', 6);
    end
end

RC_right.setThetaTick('LineWidth', 2, 'Color', [0.6, 0.6, 0.8]);
RC_right.setRTick('LineWidth', 1.5, 'Color', [0.8, 0.6, 0.6]);
RC_right.setPropLabel('FontSize', 15, 'FontName', 'Times New Roman', 'Color', [0 0 0]);
RC_right.setRLabel('FontSize', 15, 'FontName', 'Times New Roman', 'Color', [0.8 0 0]);
RC_right.setBkg('FaceColor', 'none');

% Save figure
filename_right = fullfile(pac_fig_path, "PAC_Right_Hemisphere_radarChart.tiff");
print(gcf, '-dtiff', '-r600', filename_right);
fprintf(" Advanced right hemisphere radar chart saved: %s\n", filename_right);

%% --- Whole-Brain PAC Radar Chart (Combined Hemispheres) ---

% Initialize data matrix: rows = networks, columns = [Anxiety, HC]
overall_data_radar = NaN(num_networks, 2);

for i = 1:num_networks
    network_range = network_dict.(networks{i});
    overall_data_radar(i, 1) = mean(group_avg_pac{1}(network_range), 'omitnan');
    overall_data_radar(i, 2) = mean(group_avg_pac{2}(network_range), 'omitnan');
end

% Color and transparency settings
colorList = [78 101 155; 138 140 191] ./ 255;
alpha_val = 0.3;

% Draw radar chart
X_all = overall_data_radar';  % rows = group (HC/Anxiety), columns = networks
figure;
RC_all = radarChart(X_all, 'Type', 'Line');
RC_all.PropName = networks;
RC_all.ClassName = {'HC', 'Anxiety'};
RC_all = RC_all.draw();
RC_all.legend();

% Style and fill settings
for n = 1:RC_all.ClassNum
    if strcmp(RC_all.Type, 'Line')
        RC_all.setPatchN(n, ...
            'Color', colorList(n,:), ...
            'LineWidth', 2, ...
            'Marker', 'o', ...
            'MarkerFaceColor', colorList(n,:), ...
            'MarkerSize', 6);
    elseif strcmp(RC_all.Type, 'Patch') || strcmp(RC_all.Type, 'Both')
        RC_all.setPatchN(n, ...
            'EdgeColor', colorList(n,:), ...
            'LineWidth', 1.8, ...
            'FaceAlpha', alpha_val, ...
            'MarkerFaceColor', colorList(n,:), ...
            'MarkerSize', 6);
    end
end

% Aesthetic settings
RC_all.setThetaTick('LineWidth', 2, 'Color', [0.6, 0.6, 0.8]);
RC_all.setRTick('LineWidth', 1.5, 'Color', [0.8, 0.6, 0.6]);
RC_all.setPropLabel('FontSize', 15, 'FontName', 'Times New Roman', 'Color', [0 0 0]);
RC_all.setRLabel('FontSize', 15, 'FontName', 'Times New Roman', 'Color', [0.8 0 0]);
RC_all.setBkg('FaceColor', 'none');

% Save figure
filename_all = fullfile(pac_fig_path, "PAC_WholeBrain_radarChart.tiff");
print(gcf, '-dtiff', '-r600', filename_all);
fprintf(" Whole-brain radar chart saved: %s\n", filename_all);
