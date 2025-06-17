
clc; clear variables;
restoredefaultpath;

% === BASIC PATH SETUP ===
matlab_path = matlabroot;
project_root = '...';  % Replace with general project path
current_project_path = fullfile(project_root, 'PAC_network');

% === TOOLBOX SETUP ===
toolbox_paths = {
    fullfile(matlab_path, 'toolbox', 'eeglab2021.0');
    fullfile(matlab_path, 'toolbox', 'fieldtrip-master')
};
external_subfolders = {'mne', 'brewermap', 'matplotlib', 'cmocean', 'colorcet', 'gifti', 'spm12'};
fieldtrip_external_path = fullfile(matlab_path, 'toolbox', 'fieldtrip-master', 'external');
cellfun(@addpath, toolbox_paths);
cellfun(@(sf) addpath(fullfile(fieldtrip_external_path, sf)), external_subfolders);
ft_defaults;

% === MUTE TOOL ADDITION ===
mute_path = fullfile(current_project_path, 'Tools', 'direction_relevent', 'mute_neurips-master');
addpath(mute_path);
addpath(genpath(mute_path));

% === PATH STRUCTURE ===
input_root = fullfile(current_project_path, 'Result');
output_root = fullfile(current_project_path, 'Result');
source_pac_sequence_path = fullfile(input_root, 'Source_result', 'pac_sequence');
source_net_path = fullfile(output_root, 'Direct_net_result', 'TE');
source_tool_path = fullfile(current_project_path, 'Tools', 'souce_relevent');
sensor_tool_path = fullfile(current_project_path, 'Tools', 'sensor_relevent');
channels_labels = load(fullfile(source_tool_path, 'virtual_labels.mat')).virtual_label;

% === PARAMETERS ===
group = {'Anxiety','HC'};
seg_name = 'seg';
numProcessors = 6;
comp_channel = 90;
channel_loc_file = fullfile(sensor_tool_path, 'HBN.sfp');

% === ROI MAPPING ===
regions_AAL90 = {... % AAL90 regions grouped by network
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

all_regions = horzcat(regions_AAL90{:});
roi_indices = find(ismember(channels_labels, all_regions));
roi_names = channels_labels(roi_indices);

% === LOAD SUBJECT LIST ===
subject_table = readtable(fullfile(current_project_path, 'Result', 'Demo_result', 'result_table_refine.xlsx'));
subject_names = string(subject_table{:, 1});

% === GROUP LOOP ===
for group_n = 1:numel(group)
    group_name = group{group_n};
    group_path = fullfile(source_pac_sequence_path, group_name);
    subject_files = dir(fullfile(group_path, '*.mat'));

    for sub_n = 1:numel(subject_files)
        sub_file = subject_files(sub_n).name;
        sub_name = sub_file(1:end-4);
        if ~ismember(sub_name, subject_names), continue; end

        % === LOAD PAC SEQUENCE ===
        pac_struct = load(fullfile(group_path, sub_file));
        fields = fieldnames(pac_struct.sub_pac_sequence);
        segment_count = numel(pac_struct.sub_pac_sequence.(fields{1}));

        pac_matrix = zeros(numel(fields), segment_count);
        for f = 1:numel(fields)
            pac_matrix(f, :) = pac_struct.sub_pac_sequence.(fields{f});
        end
        pac_matrix = pac_matrix(roi_indices,:);
        pac_tensor = reshape(pac_matrix, 34, [], 10);

        te_all = {}; mute_all = {};
        for seg = 1:size(pac_tensor,3)
            data = pac_tensor(:,:,seg);
            sub_sub_name = sprintf('%s_%d', sub_name, seg);

            res_path = fullfile(source_net_path, seg_name, group_name, sub_sub_name);
            copy_path = fullfile(res_path, 'entropyMatrices');
            result_dir = fullfile(res_path, 'results');
            cellfun(@(d) mkdir(d), {res_path, copy_path, result_dir});

            save(fullfile(res_path, [sub_sub_name, '.mat']), 'data');

            mute_file_path = fullfile(source_net_path, [sub_name, '_mute.mat']);
            if exist(mute_file_path, 'file'), continue; end

            % === MUTE PARAMETERS ===
            mute_params = struct(
                'samplingRate', 5,
                'pointsToDiscard', 0,
                'channels', 1:numel(roi_indices),
                'modelOrder', 5,
                'numQuantLevels', 6,
                'numSurrogates', 100,
                'alphaPercentile', 0.05,
                'multi_bivAnalysis', 'multiv'
            );

            [mute_output, mute_cfg] = f_MyparametersAndMethods(...
                dir(fullfile(res_path, [sub_sub_name, '*'])), ...
                mute_params.samplingRate, mute_params.pointsToDiscard, mute_params.channels, ...
                ones(1,8), zeros(1,8), result_dir, res_path, copy_path, numProcessors, ...
                'binnue', [], [], [], mute_params.modelOrder, mute_params.multi_bivAnalysis, ...
                mute_params.numQuantLevels, @evaluateNonUniformEntropy, @quantization, [1 1], ...
                mute_params.numSurrogates, mute_params.alphaPercentile, @generateConditionalTerm, 0, 0);

            mute_struct.output = mute_output;
            mute_struct.params = mute_cfg;
            te_matrix = load(fullfile(res_path, 'results', 'binnue_meanReshapeMtx.mat'));

            te_data.labels = roi_names;
            te_data.data = te_matrix.meanRes';

            te_all{seg} = te_data;
            mute_all{seg} = mute_struct;
        end

        save(fullfile(source_net_path, [sub_name, '_mute.mat']), 'mute_all');
        save(fullfile(source_net_path, [sub_name, '_te.mat']), 'te_all');
    end
end

delete(gcp('nocreate'));
