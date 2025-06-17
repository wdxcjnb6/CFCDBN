
clc; clear; restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED;

% === Initialize Paths ===
project_root = '<PROJECT_ROOT>';
toolbox_root = fullfile(project_root, 'toolbox');

% Toolboxes
addpath(fullfile(toolbox_root, 'eeglab2021.0'));
addpath(fullfile(toolbox_root, 'fieldtrip-20240110'));
ft_defaults;

external_paths = {'mne','brewermap','matplotlib','cmocean','colorcet','gifti','spm12'};
cellfun(@(d)addpath(fullfile(toolbox_root, 'fieldtrip-master', 'external', d)), external_paths);

addpath(genpath(fullfile(toolbox_root, 'fooof_mat-main')));

% Project structure
code_path = fullfile(project_root, 'Code');
data_path = fullfile(project_root, 'Data');
tools_path = fullfile(project_root, 'Tools');
result_path = fullfile(project_root, 'Result');

addpath(code_path);
addpath(data_path);
addpath(tools_path);
addpath(genpath(result_path));
cd(code_path);

% Define data paths
source_path = fullfile(result_path, 'Source_result', 'source_location_data');
group_plot_path = fullfile(result_path, 'Source_result', 'source_location_group_plot');
atlas_path = fullfile(tools_path, 'souce_relevent', 'ROI_MNI_V4.nii');
mri_template_path = fullfile(tools_path, 'souce_relevent', 'single_subj_T1.nii');
surf_file = fullfile(toolbox_root, 'fieldtrip-20240110', 'template', 'sourcemodel', 'surface_inflated_both.mat');

% === Load Atlas & Electrode Template ===
elec_path = fullfile(tools_path, 'sensor_relevent', 'HBN.sfp');
elec = ft_read_sens(elec_path); elec = ft_convert_units(elec,'mm');
atlas = ft_read_atlas(atlas_path); atlas = ft_convert_units(atlas, 'mm');

% === Load Subject List ===
T = readtable(fullfile(result_path, 'Demo_result', 'result_table_refine.xlsx'));
subnames = string(T{:,1});

% === Load Source Data ===
groups = {'HC','Anxiety'};
bands = {'Delta','Beta'};
source_all = struct();

for b = 1:numel(bands)
    band = bands{b};
    for g = 1:numel(groups)
        grp = groups{g};
        data_dir = fullfile(source_path, grp);
        files = dir(fullfile(data_dir, ['*' band '*']));
        names = {files.name};
        names = names(contains(names, subnames));
        tmp = cell(1,numel(names));
        for i = 1:numel(names)
            try
                d = load(fullfile(data_dir, names{i}));
                s = d.source;
                if isfield(s, 'cfg'), s = rmfield(s,'cfg'); end
                if isfield(s, 'method'), s = rmfield(s,'method'); end
                if isfield(s, 'time'), s = rmfield(s,'time'); end
                if isfield(s, 'avg')
                    fn = fieldnames(s.avg);
                    s.avg = rmfield(s.avg, setdiff(fn, {'pow'}));
                end
                tmp{i} = s;
            catch
                fprintf('Failed: %s\n', names{i});
            end
        end
        source_all.(grp).(band) = tmp;
    end
end

%% === Cluster-Based Permutation Test ===
band = 'Beta';
group1 = source_all.Anxiety.(band);
group2 = source_all.HC.(band);

cfg = [];
cfg.dim = group1{1}.dim;
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_indepsamplesT';
cfg.parameter = 'avg.pow';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail = 1;
cfg.clustertail = 1;
cfg.alpha = 0.05;
cfg.numrandomization = 4000;
cfg.computecritval = 'yes';
cfg.design = [ones(1,numel(group2)), 2*ones(1,numel(group1))];
cfg.ivar = 1;

fprintf('Running cluster-permutation for %s...\n', band);
stat = ft_sourcestatistics(cfg, group2{:}, group1{:});
fprintf('Minimum p-value: %.5f\n', min(stat.prob(:)));

%% === Visualization ===
stat_vals = stat.stat(:); stat_vals = stat_vals(~isnan(stat_vals));
if isempty(stat_vals), warning('Empty t-map'); return; end
vmax = max(abs(stat_vals)); if vmax==0 || isnan(vmax), vmax=1; end
use_mask = isfield(stat, 'mask') && any(stat.mask(:));

cfg_plot = [];
cfg_plot.funparameter = 'stat';
cfg_plot.zlim = [-vmax vmax];
cfg_plot.location = 'max';
cfg_plot.maskparameter = ternary(use_mask, 'mask', '');

figure; ft_sourceplot(cfg_plot, stat);
colormap(ft_hastoolbox('brewermap',1)*flipud(brewermap(64,'RdBu')) + ~ft_hastoolbox('brewermap',1)*parula);
title(['T-map (raw) - ', band]);
print(['stat_raw_', band], '-dpng', '-r300');

%% === Interpolation to MRI ===
mri = ft_read_mri(mri_template_path);
cfg_interp = []; cfg_interp.parameter = 'stat'; cfg_interp.interpmethod = 'nearest'; cfg_interp.voxelcoord = 'no';
statint = ft_sourceinterpolate(cfg_interp, stat, mri);

if use_mask
    cfg_interp.parameter = 'mask';
    maskint = ft_sourceinterpolate(cfg_interp, stat, mri);
    statint.mask = maskint.mask;
end

cfgp = [];
cfgp.method = 'ortho';
cfgp.funparameter = 'stat';
cfgp.zlim = [-vmax vmax];
cfgp.location = 'max';
cfgp.maskparameter = ternary(use_mask, 'mask', '');

figure; ft_sourceplot(cfgp, statint);
colormap(ft_hastoolbox('brewermap',1)*flipud(brewermap(64,'RdBu')) + ~ft_hastoolbox('brewermap',1)*parula);
title(['T-map (MRI) - ', band]);