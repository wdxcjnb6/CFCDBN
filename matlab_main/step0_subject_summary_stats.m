%% CLEAN AND SET PATHS
clc; clear variables;
restoredefaultpath; % Restore default MATLAB path

% === Project and toolbox path configuration ===
project_root = 'PROJECT_ROOT';
addpath(project_root);

toolbox_paths = {
    'EEGLAB_PATH', ...
    'FIELDTRIP_PATH'
};
external_subfolders = {'mne', 'brewermap', 'matplotlib', 'cmocean', 'colorcet', 'gifti', 'spm12'};
fieldtrip_external_path = fullfile('FIELDTRIP_PATH', 'external');

cellfun(@addpath, toolbox_paths);
cellfun(@(subfolder) addpath(fullfile(fieldtrip_external_path, subfolder)), external_subfolders);

% Initialize FieldTrip
ft_defaults;

% Add relevant project subfolders
addpath(fullfile(project_root, 'Code'));
addpath(fullfile(project_root, 'Data'));
addpath(fullfile(project_root, 'Tools'));
addpath(genpath(fullfile(project_root, 'Result')));

%% Define relevant inputoutput directories
fooof_path = fullfile(project_root, 'Result', 'Fooof_result', 'fooof_sub_select_periodic');
demo_path = fullfile(project_root, 'Code', 'Code', 'scale_sum', 'Demo_summary');
anxiety_path = fullfile(project_root, 'Code', 'Code', 'scale_sum', 'Anxiety');
save_path = fullfile(project_root, 'Result', 'Demo_result');

%% Read required files
f_file_info = dir(fullfile(fooof_path, '.xlsx'));
d_file_info = dir(fullfile(demo_path, 'HC_ADdemo'));
ch_file_info = dir(fullfile(anxiety_path, 'Scared_SR'));
pa_file_info = dir(fullfile(anxiety_path, 'Scared_P'));

if isempty(f_file_info)  isempty(d_file_info)  isempty(ch_file_info)  isempty(pa_file_info)
    error('One or more required files are missing.');
end

f_file_name = f_file_info(1).name;
d_file_name = d_file_info(1).name;
ch_file_name = ch_file_info(1).name;
pa_file_name = pa_file_info(1).name;

try
    f_data = readtable(fullfile(fooof_path, f_file_name), 'Sheet', '_Beta_', 'VariableNamingRule', 'preserve');
    d_data = readtable(fullfile(demo_path, d_file_name), 'VariableNamingRule', 'preserve');
    ch_data = readtable(fullfile(anxiety_path, ch_file_name), 'VariableNamingRule', 'preserve');
    pa_data = readtable(fullfile(anxiety_path, pa_file_name), 'VariableNamingRule', 'preserve');
catch ME
    error('Error reading files %s', ME.message);
end

%% Match subject-related data
sub_name = unique(f_data.Subject);
sub_group = f_data.Group;

demo_vars = {'Age', 'Sex'};
[matched_ages, matched_sexes] = f_match_table_data(sub_name, d_data, 'EID', demo_vars);

fooof_vars = {'Freq Peak Fr', 'Freq Peak Bw'};
[matched_Beta_Peak_Fr, matched_Beta_Peak_Bw] = f_match_table_data(sub_name, f_data, 'Subject', fooof_vars);

scared_ch_vars = {'SR_GD', 'SR_PN', 'SR_SC', 'SR_SH', 'SR_SP', 'SR_Total'};
[matched_SR_GD, matched_SR_PN, matched_SR_SC, matched_SR_SH, matched_SR_SP, matched_SR_Total] = ...
    f_match_table_data(sub_name, ch_data, 'EID', scared_ch_vars);

scared_pa_vars = {'P_GD', 'P_PN', 'P_SC', 'P_SH', 'P_SP', 'P_Total'};
[matched_P_GD, matched_P_PN, matched_P_SC, matched_P_SH, matched_P_SP, matched_P_Total] = ...
    f_match_table_data(sub_name, pa_data, 'EID', scared_pa_vars);

%% Construct and save result table
result_table = table(sub_name, sub_group, matched_ages, matched_sexes, matched_Beta_Peak_Fr, matched_Beta_Peak_Bw, ...
    matched_SR_GD, matched_SR_PN, matched_SR_SC, matched_SR_SH, matched_SR_SP, matched_SR_Total, ...
    matched_P_GD, matched_P_PN, matched_P_SC, matched_P_SH, matched_P_SP, matched_P_Total, ...
    'VariableNames', {'subject', 'group', 'age', 'sex', 'Beta_Fr','Beta_Bw', ...
    'SR_GD', 'SR_PN', 'SR_SC', 'SR_SH', 'SR_SP', 'SR_Total', ...
    'P_GD', 'P_PN', 'P_SC', 'P_SH', 'P_SP', 'P_Total'});

if ~exist(save_path, 'dir'), mkdir(save_path); end
writetable(result_table, fullfile(save_path, 'result_table.xlsx'));

%% Load refined result table
result_table_refine = readtable(fullfile(save_path, 'result_table_refine.xlsx'), 'VariableNamingRule', 'preserve');
if isempty(result_table_refine)
    error('Failed to read refined result table.');
end

%% Define group indices
AD_idx = strcmp(result_table_refine.group, 'Anxiety');
HC_idx = strcmp(result_table_refine.group, 'HC');
if ~any(AD_idx)  ~any(HC_idx)
    error('Invalid group labels.');
end

%% Extract variables for each group
variables = {'age', 'SR_Total', 'P_Total', 'sex', 'Beta_Fr', 'Beta_Bw'};
data = struct();
for i = 1numel(variables)
    var = variables{i};
    data.(var).AD = result_table_refine.(var)(AD_idx);
    data.(var).HC = result_table_refine.(var)(HC_idx);
end

%% Perform t-tests and descriptive stats
[h_age, p_age, m1, s1, m2, s2] = f_perform_t_test(data.age.AD, data.age.HC);
fprintf('Age comparison p = %.4f  AD %.2f±%.2f  HC %.2f±%.2fn', p_age, m1, s1, m2, s2);

[h_sr, p_sr, m1, s1, m2, s2] = f_perform_t_test(data.SR_Total.AD, data.SR_Total.HC);
fprintf('SR_Total comparison p = %.4f  AD %.2f±%.2f  HC %.2f±%.2fn', p_sr, m1, s1, m2, s2);

[h_pt, p_pt, m1, s1, m2, s2] = f_perform_t_test(data.P_Total.AD, data.P_Total.HC);
fprintf('P_Total comparison p = %.4f  AD %.2f±%.2f  HC %.2f±%.2fn', p_pt, m1, s1, m2, s2);

%% Perform chi-square test on sex
sex_table = [sum(data.sex.AD == 0), sum(data.sex.AD == 1); sum(data.sex.HC == 0), sum(data.sex.HC == 1)];
[h_sex, p_sex] = chi2test(sex_table);
fprintf('Sex comparison p = %.4f  AD Male %.2f  HC Male %.2fn', p_sex, ...
    sum(data.sex.AD == 1)length(data.sex.AD), sum(data.sex.HC == 1)length(data.sex.HC));

%% Perform FOOOF feature comparison
[h_fr, p_fr, m1, s1, m2, s2] = f_perform_t_test(data.Beta_Fr.AD, data.Beta_Fr.HC);
fprintf('Beta Fr p = %.4f  AD %.2f±%.2f  HC %.2f±%.2fn', p_fr, m1, s1, m2, s2);

[h_bw, p_bw, m1, s1, m2, s2] = f_perform_t_test(data.Beta_Bw.AD, data.Beta_Bw.HC);
fprintf('Beta Bw p = %.4f  AD %.2f±%.2f  HC %.2f±%.2fn', p_bw, m1, s1, m2, s2);
