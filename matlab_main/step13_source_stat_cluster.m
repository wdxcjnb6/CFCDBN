%% SETUP AND TOOLBOX PATHS
clc; clear variables;
restoredefaultpath; clear RESTOREDEFAULTPATH_EXECUTED;

% Define base project path
project_root = '[PROJECT_ROOT]'; % <-- Replace with relative path
code_path = fullfile(project_root, 'Code');
data_path = fullfile(project_root, 'Data');
tool_path = fullfile(project_root, 'Tools');
result_path = fullfile(project_root, 'Result');

% Add necessary toolboxes
addpath(fullfile('[MATLAB_TOOLBOX_PATH]', 'eeglab2021.0'));
addpath(fullfile('[MATLAB_TOOLBOX_PATH]', 'fieldtrip-20240110'));
addpath(genpath(fullfile('[MATLAB_TOOLBOX_PATH]', 'fooof_mat-main')));
external_subfolders = {'mne','brewermap','matplotlib','cmocean','colorcet','gifti','spm12'};
cellfun(@(s)addpath(fullfile('[FIELDTRIP_EXTERNAL]', s)), external_subfolders);
ft_defaults;

% Add project subfolders
addpath(code_path);
addpath(data_path);
addpath(tool_path);
addpath(genpath(result_path));
cd(code_path);

%% DEFINE PATHS
input_source_path = fullfile(result_path, 'Source_result', 'source_location_data');
output_stat_path = fullfile(result_path, 'Source_result', 'source_location_ana');
output_plot_path = fullfile(result_path, 'Source_result', 'source_location_group_plot');
f_createOutputFolders(result_path, {output_stat_path, output_plot_path});

elec_file = fullfile(tool_path, 'sensor_relevent', 'HBN.sfp');
elec = ft_read_sens(elec_file); elec = ft_convert_units(elec, 'mm');

atlas = ft_read_atlas(fullfile(tool_path, 'souce_relevent', 'ROI_MNI_V4.nii'));
atlas = ft_convert_units(atlas, 'mm');
ROI = atlas.tissuelabel(1:90);

%% LOAD SUBJECT INFO AND INITIALIZE
group_list = {'HC', 'Anxiety'};
band_list = {'Delta', 'Beta'};
summary_data = struct();
subject_table = readtable(fullfile(result_path, 'Demo_result', 'result_table_refine.xlsx'));
subject_names = string(subject_table{:,1});

if ~exist(output_plot_path, 'dir'), mkdir(output_plot_path); end

%% LOAD AND CLEAN SOURCE STRUCTS
for b = 1:numel(band_list)
    band = band_list{b};
    for g = 1:numel(group_list)
        group = group_list{g};
        group_folder = fullfile(input_source_path, group);
        file_list = dir(fullfile(group_folder, ['*' band '*']));
        file_names = {file_list.name};
        file_names = file_names(~ismember(file_names,{'.','..'}));
        valid_files = file_names(contains(file_names, subject_names));
        group_sources = cell(1, numel(valid_files));
        for i = 1:numel(valid_files)
            try
                d = load(fullfile(group_folder, valid_files{i}));
                source = d.source;
                rmf = intersect(fieldnames(source), {'cfg','method','time'});
                source = rmfield(source, rmf);
                if isfield(source, 'avg')
                    fn = fieldnames(source.avg);
                    source.avg = rmfield(source.avg, setdiff(fn, {'pow'}));
                end
                group_sources{i} = source;
            catch
                fprintf('Error loading: %s\n', valid_files{i});
            end
        end
        summary_data.(group).(band) = group_sources;
        fprintf('Loaded %d source files for %s - %s\n', numel(group_sources), group, band);
    end
end

%% PERFORM CLUSTER-BASED STATISTICS
save_base = fullfile(output_plot_path, 'figures');
for b = 1:numel(band_list)
    band = band_list{b};
    src1 = summary_data.Anxiety.(band);
    src2 = summary_data.HC.(band);

    cfg = [];
    cfg.dim = src1{1}.dim;
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
    cfg.design = [ones(1,numel(src2)), 2*ones(1,numel(src1))];
    cfg.ivar = 1;

    stat = ft_sourcestatistics(cfg, src2{:}, src1{:});
    fprintf('Cluster test for %s completed. Min p: %.5f\n', band, min(stat.prob(:)));

    %% VISUALIZATION
    save_path = fullfile(save_base, band);
    if ~exist(save_path, 'dir'), mkdir(save_path); end
    cfg_plot = [];
    cfg_plot.method = 'surface';
    cfg_plot.funparameter = 'stat';
    cfg_plot.projmethod = 'sphere_avg';
    cfg_plot.sphereradius = 10;

    views = {[90 0], [-90 0]};
    hemis = {'left','right'};
    view_labels = {'rightside','leftside'};

    for h = 1:numel(hemis)
        hemi = hemis{h};
        cfg_plot.surffile = sprintf('surface_white_%s.mat', hemi);
        cfg_plot.surfinflated = sprintf('surface_inflated_%s_caret.mat', hemi);
        for v = 1:numel(views)
            fig = figure('Visible', 'off'); set(fig, 'Renderer', 'painters');
            ft_sourceplot(cfg_plot, stat);
            view(views{v}); lighting gouraud; camlight('headlight');
            material dull; shading interp; axis vis3d off;
            if ft_hastoolbox('brewermap',1)
                colormap(flipud(brewermap(64, 'RdBu')));
            else
                colormap(parula);
            end
            title(sprintf('T-map | %s | %s | %s', view_labels{v}, hemi, band));
            figname = sprintf('Tmap_%s_%s_%s.png', band, view_labels{v}, hemi);
            exportgraphics(fig, fullfile(save_path, figname), 'Resolution', 600);
            close(fig);
        end
    end
end
