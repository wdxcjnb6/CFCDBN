%% CLEAN AND SET PATHS
% This script visualizes brain networks with color-coded ROIs
clc; clear variables;
restoredefaultpath; % Restore default MATLAB path

% === Basic project and toolbox path setup ===
project_root = 'YOUR_PROJECT_ROOT'; % <- Replace with your project root path
addpath(project_root);

% Toolbox paths (modify if needed)
toolbox_paths = {
    'YOUR_EEGLAB_PATH', ...      % Replace with EEGLAB toolbox path
    'YOUR_FIELDTRIP_PATH' ...    % Replace with FieldTrip toolbox path
};
external_subfolders = {'mne', 'brewermap', 'matplotlib', 'cmocean', 'colorcet', 'gifti', 'spm12'};
fieldtrip_external_path = fullfile('YOUR_FIELDTRIP_PATH', 'external');
cellfun(@addpath, toolbox_paths);
cellfun(@(sf) addpath(fullfile(fieldtrip_external_path, sf)), external_subfolders);

ft_defaults;

addpath(fullfile(project_root, 'Code'));
addpath(fullfile(project_root, 'Data'));
addpath(fullfile(project_root, 'Tools'));
addpath(genpath(fullfile(project_root, 'Result')));

%% Define output folders
output_paths = {'Direct_net_result/Direct_net_figure/Figure'};
output_root = fullfile(project_root, 'Result');
f_createOutputFolders(output_root, output_paths);

% Tool paths (replace with your organization structure if needed)
input_tool_root = fullfile(project_root, 'Tools');
source_tool_path = fullfile(input_tool_root, 'souce_relevent');
excel_tool_path = fullfile(input_tool_root, 'excel_relevent');
sensor_tool_path = fullfile(input_tool_root, 'sensor_relevent');

addpath(genpath('YOUR_SOURCEMESH_PATH')); % <- Replace with path to SourceMesh if used

%% Load channel and ROI definitions
channels_labels = load(fullfile(source_tool_path, 'virtual_labels.mat')).virtual_label;

% Define brain networks (AAL regions)
DMN_roi = {'Hippocampus_L','Hippocampus_R','Frontal_Sup_Medial_L','Frontal_Sup_Medial_R','Precuneus_L','Precuneus_R','Angular_L','Angular_R','Cingulum_Post_L','Cingulum_Post_R'};
SN_roi = {'Insula_L','Insula_R','Cingulum_Ant_L','Cingulum_Ant_R','Amygdala_L','Amygdala_R','Thalamus_L','Thalamus_R'};
FPN_roi = {'Frontal_Mid_L','Frontal_Mid_R','Parietal_Inf_L','Parietal_Inf_R','Frontal_Inf_Tri_L','Frontal_Inf_Tri_R'};
DAN_roi = {'Frontal_Sup_L','Frontal_Sup_R','Parietal_Sup_L','Parietal_Sup_R'};
SMN_roi = {'Postcentral_L','Postcentral_R','Precentral_L','Precentral_R','Supp_Motor_Area_L','Supp_Motor_Area_R'};

% Load brain template
load('New_AALROI_6mm.mat');
positions = template_sourcemodel.pos;
roi_tissue = all_roi_tissueindex;
atlas = 'aal90';

% Prepare network info
network_names = {'Default Mode Network (DMN)','Salience Network (SN)','Frontoparietal Network (FPN)','Dorsal Attention Network (DAN)','Somatosensory-Motor Network (SMN)'};
roi_list = {DMN_roi, SN_roi, FPN_roi, DAN_roi, SMN_roi};

%% Visualization: different colors within same network
for net_idx = 1:length(roi_list)
    overlay = zeros(numel(AAL_Labels), 1);
    roi = sort(roi_list{net_idx}(:));
    color = 1;
    for i = 1:length(roi)
        idx = find(strcmp(AAL_Labels, roi{i}));
        if ~isempty(idx)
            overlay(idx) = color;
            if mod(i,2) == 0, color = color+1; end
        else
            fprintf('ROI %s not found in AAL_Labels\n', roi{i});
        end
    end
    
    figure('position', [481 176 1204 987]);
    atemplate('mesh','def4','overlay',overlay,'method',{'aal_super','aal_light'});
    colormap([repmat([0.8 0.8 0.8],1,1); parula(256)]);
    caxis([min(overlay) max(overlay)]);
    alpha 0.5;
    sgtitle(['Network: ' network_names{net_idx}]);
    slice3();
end

%% Visualization: same color for entire network
for net_idx = 1:length(roi_list)
    overlay = zeros(numel(AAL_Labels), 1);
    roi = sort(roi_list{net_idx}(:));
    for i = 1:length(roi)
        idx = find(strcmp(AAL_Labels, roi{i}));
        if ~isempty(idx)
            overlay(idx) = 1;
        else
            fprintf('ROI %s not found in AAL_Labels\n', roi{i});
        end
    end
    
    figure('position', [481 176 1204 987]);
    atemplate('mesh','def4','overlay',overlay,'method',{'aal_super',[]});
    f_brain_slice2();
    map = parula(256);
    selected_color = map(50*net_idx,:);
    colormap([repmat([0.9 0.9 0.9],1,1); selected_color]);
    caxis([min(overlay) max(overlay)]);
    alpha 0.7;
    sgtitle(['Network: ' network_names{net_idx}]);
end
