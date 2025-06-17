
clc; clear variables; restoredefaultpath;

%% === Initialize FieldTrip ===
matlab_path = matlabroot;
project_root = 'E:\PAC_network';
fieldtrip_path = fullfile(matlab_path, 'toolbox', 'fieldtrip-master');

% Add FieldTrip and external dependencies
ext_subfolders = {'mne','brewermap','matplotlib','cmocean','colorcet','gifti','spm12','freesurfer'};
cellfun(@(f) addpath(fullfile(fieldtrip_path, 'external', f)), ext_subfolders);
addpath(fieldtrip_path); ft_defaults;

% Add project-specific paths
addpath(fullfile(project_root, 'Code'));
addpath(fullfile(project_root, 'Data'));
addpath(fullfile(project_root, 'Tools'));
addpath(fullfile(project_root, 'Result'));

cd(fullfile(project_root, 'Tools', 'brainnetwork_relevent'));

%% === Load FreeSurfer Mesh (Left Hemisphere) ===
% The input should include 'lh.aparc.annot' and 'lh.pial'
fs_lh_atlas = ft_read_atlas({'lh.aparc.annot', 'lh.pial'});

%% === Create High-Resolution Plot Window ===
figure;
set(gcf, 'Position', [100, 100, 1000, 800]);  % Window size in pixels

%% === Plot the Brain Mesh ===
ft_plot_mesh(fs_lh_atlas, 'vertexcolor', fs_lh_atlas.aparc);
view([-1 0 0]);  % Left hemisphere view
material dull;
camlight headlight;  % Add headlight source

%% === Optional: Further Lighting Enhancements ===
% lighting gouraud;
% material shiny;
% camlight right;
% set(gca, 'Color', 'w');
% axis off;

%% === Export as High-Resolution Image (600 DPI) ===
save_dir = fullfile(project_root, 'Figure', 'Braintemplate');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

exportgraphics(gcf, fullfile(save_dir, 'brain_network_high_res.png'), 'Resolution', 600);

disp('✅ Brain surface mesh rendered and exported successfully.');
