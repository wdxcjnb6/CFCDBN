%% === STEP 0: Visualize MNI Template as Glass Brain ===
clc; clear variables; restoredefaultpath;

%% === Initialize FieldTrip ===
matlab_path = matlabroot;
project_root = 'E:\PAC_network';
fieldtrip_path = fullfile(matlab_path, 'toolbox', 'fieldtrip-master');

% Add FieldTrip and external toolboxes
ext_subfolders = {'mne','brewermap','matplotlib','cmocean','colorcet','gifti','spm12','freesurfer'};
cellfun(@(f) addpath(fullfile(fieldtrip_path, 'external', f)), ext_subfolders);
addpath(fieldtrip_path); ft_defaults;

% Add project directories
addpath(fullfile(project_root, 'Code'));
addpath(fullfile(project_root, 'Data'));
addpath(fullfile(project_root, 'Tools'));
addpath(fullfile(project_root, 'Result'));

%% === Load MNI T1 Template MRI ===
mni_template = ft_read_mri('single_subj_T1_1mm.nii');

%% === Plot MNI Brain Using Glassbrain View ===
cfg = [];
cfg.method = 'glassbrain';  % 3D transparent overview
ft_sourceplot(cfg, mni_template);

disp('MNI glass brain plot rendered successfully.');
