
clc; clear; restoredefaultpath;

%% === PATH SETUP ===
project_root = '<PROJECT_ROOT>';
toolbox_root = '<TOOLBOX_ROOT>';
addpath(fullfile(toolbox_root, 'fieldtrip-master')); ft_defaults;
addpath(fullfile(project_root,'Code'));
addpath(fullfile(project_root,'Data'));
addpath(genpath(fullfile(project_root,'Result')));

output_root = fullfile(project_root,'Result');
source_net_path = fullfile(output_root, 'Direct_net_result', 'Direct_net_data', 'TE');
pac_path_root = fullfile(output_root, 'Source_result', 'pac_sequence');
graph_save_path = fullfile(output_root, 'Graph_for_python');

%% === SUBJECTS AND ROI INFO ===
file = readtable(fullfile(output_root,'Demo_result','result_table_refine.xlsx'));
subject_list = string(file{:,1});
labels = load(fullfile(project_root,'Tools','souce_relevent','virtual_labels.mat')).virtual_label;
aal_table = readtable(fullfile(project_root,'Tools','souce_relevent','Node_AAL90.txt'),'Delimiter','\t','ReadVariableNames',false);

% Define ROI list
regions = {
    {'S_Hippocampus_L','S_Hippocampus_R'},...
    {'S_Frontal_Sup_Medial_L','S_Frontal_Sup_Medial_R'},...
    {'S_Precuneus_L','S_Precuneus_R'},...
    {'S_Angular_L','S_Angular_R'}
... 
    };
all_rois = horzcat(regions{:});
roi_idx = find(ismember(labels, all_rois));
roi_idx = unique(roi_idx);
[~, reorder_idx] = sort(roi_idx);

%% === GRAPH CONSTRUCTION ===
groups = {'Anxiety','HC'};
graph_id = 0;

for g = 1:numel(groups)
    grp = groups{g};
    pac_path = fullfile(pac_path_root, grp);
    files = dir(fullfile(pac_path, '*.mat'));

    for i = 1:numel(files)
        sub_name = files(i).name(1:end-4);
        if ~ismember(sub_name, subject_list), continue; end

        % Load PAC sequence
        pac_data = load(fullfile(pac_path, files(i).name)).sub_pac_sequence;
        flds = fieldnames(pac_data);
        N = numel(flds);
        pac_mat = zeros(N, numel(pac_data.(flds{1})));
        for j = 1:N
            pac_mat(j,:) = pac_data.(flds{j})';
        end
        pac_mat = pac_mat(roi_idx,:);
        pac_mat = pac_mat(reorder_idx,:);

        % Node feature: average PAC over all segments
        node_feature = mean(pac_mat, 2); % size = [num_ROI x 1]

        % Load graph structure (TE)
        te_path = fullfile(source_net_path, grp, [sub_name '_te.mat']);
        if ~exist(te_path, 'file'), continue; end
        te_data = load(te_path).te_mat.data;
        [row, col, val] = find(te_data);
        edge_index = [row'; col'] - 1;  % for PyG

        % Create graph struct
        graph = struct();
        graph.edge_adj = edge_index';
        graph.edge_weight = val';
        graph.x = node_feature'; % 1 x num_ROI
        graph.y = double(strcmp(grp, 'Anxiety')); % 1 for Anxiety, 0 for HC

        save(fullfile(graph_save_path, ['graph_' num2str(graph_id) '.mat']), 'graph');
        graph_id = graph_id + 1;
    end
end
