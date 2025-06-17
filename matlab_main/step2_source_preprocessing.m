
clc; clear variables;
restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED;

project_root = '<PROJECT_ROOT>';
toolbox_root = '<MATLAB_TOOLBOX>';

% Toolbox paths
toolbox_paths = {
    fullfile(toolbox_root, 'eeglab2021.0'), ...
    fullfile(toolbox_root, 'fieldtrip-master')
};
external_subfolders = {'mne', 'brewermap', 'matplotlib', 'cmocean', 'colorcet', 'gifti', 'spm12'};
fieldtrip_external_path = fullfile(toolbox_root, 'fieldtrip-master', 'external');

cellfun(@addpath, toolbox_paths);
cellfun(@(s) addpath(fullfile(fieldtrip_external_path, s)), external_subfolders);
addpath(genpath(fullfile(toolbox_root, 'fooof_mat-main')));

ft_defaults;

addpath(fullfile(project_root, 'Code'));
addpath(fullfile(project_root, 'Data'));
addpath(fullfile(project_root, 'Tools'));
addpath(genpath(fullfile(project_root, 'Result')));
cd(fullfile(project_root, 'Code'));

%% DEFINE I/O PATHS
input_root = fullfile(project_root, 'Result');
input_paths = {
    'Preprocess_result/EEG_pure', ...
    'Fooof_result/fooof_sub_select_periodic'
};

input_tool_root = fullfile(project_root, 'Tools');
input_tool_paths = {
    'souce_relevent', ...
    'excel_relevent', ...
    'sensor_relevent'
};

output_root = fullfile(project_root, 'Result');
output_paths = {
    'Source_result/pac_sequence', ...
    'Source_result/vitual_signal', ...
    'Source_result/source_location'
};
f_createOutputFolders(output_root, output_paths);

EEG_pure_path = fullfile(input_root, input_paths{1});
EEG_fooof_data_path = fullfile(input_root, input_paths{2});
pac_sequence_path = fullfile(output_root, output_paths{1});
vitual_signal_path = fullfile(output_root, output_paths{2});
source_location_path = fullfile(output_root, output_paths{3});

Source_tool_path = fullfile(input_tool_root, input_tool_paths{1});
Excel_tool_path = fullfile(input_tool_root, input_tool_paths{2});
Channel_location_tool_path = fullfile(input_tool_root, input_tool_paths{3});

%% INITIALIZE TOOLS
numProcessors = 5;
elec_file = fullfile(Channel_location_tool_path, 'HBN.sfp');
elec = ft_read_sens(elec_file);
elec = ft_convert_units(elec, 'mm');

atlas_file = fullfile(Source_tool_path, 'ROI_MNI_V4.nii');
atlas = ft_read_atlas(atlas_file);
atlas = ft_convert_units(atlas, 'mm');
ROI = atlas.tissuelabel(1:90);

percell_type = 'centroid';
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui');

%% PROCESS GROUPS
groups = {'HC','Anxiety','ASD','Depression','LD','ADHD'};
for g = 1:numel(groups)
    group_name = groups{g};
    group_eeg_path = fullfile(EEG_pure_path, group_name);
    fooof_group_path = fullfile(EEG_fooof_data_path, group_name);

    % Create output folders
    mkdir(fullfile(pac_sequence_path, group_name));
    mkdir(fullfile(vitual_signal_path, group_name));
    mkdir(fullfile(source_location_path, group_name));

    % Load template data
    folder_info = dir(fullfile(group_eeg_path, '*.set'));
    if isempty(folder_info), continue; end
    folder_list = string({folder_info.name}');
    
    headmodel = load(fullfile(Source_tool_path, 'standard_bem.mat')).vol;
    headmodel = ft_convert_units(headmodel, 'mm');
    mri = ft_read_mri(fullfile(Source_tool_path, 'standard_mri.mat'));
    mri = ft_convert_units(mri, 'mm');
    elec_aligned = load(fullfile(Source_tool_path, 'HBN')).elec_aligned;
    elec_aligned = ft_convert_units(elec_aligned, 'mm');
    
    cfg = [];
    cfg.headmodel = headmodel;
    cfg.channel = 'eeg';
    cfg.elec = elec_aligned;
    cfg.resolution = 6;
    cfg.sourcemodel.unit = 'mm';
    cfg.warpmni = 'yes';
    sourcemodel = ft_prepare_sourcemodel(cfg);

    cfg = [];
    cfg.interpmethod = 'nearest';
    cfg.parameter = 'tissue';
    grid_aal = ft_sourceinterpolate(cfg, atlas, sourcemodel);

    ndipole = find(sourcemodel.inside);
    roi_indices = cell(numel(ROI), 1);
    for r = 1:numel(ROI)
        roi_name = ROI{r};
        roi_indices{r} = find(grid_aal.tissue(ndipole) == find(strcmp(grid_aal.tissuelabel, roi_name)));
    end

    cfg = [];
    cfg.sourcemodel = sourcemodel;
    cfg.headmodel = headmodel;
    cfg.reducerank = 3;
    cfg.elec = elec_aligned;
    cfg.tight = 'yes';
    leadfield = ft_prepare_leadfield(cfg);

    f_plot_all_alige(atlas, headmodel, elec_aligned, leadfield);
    close all;

    fprintf('Number of subjects: %d\n', numel(folder_list));

    for s = 1:numel(folder_list)
        sub_name = folder_list{s};
        save_file = fullfile(vitual_signal_path, group_name, sprintf('%s.mat', sub_name(1:min(12,end))));
        if isfile(save_file)
            fprintf('Virtual signal for subject %s exists. Skipping.\n', sub_name);
            continue;
        end

        % === Load and convert EEG ===
        sub_eeg = f_safe_loadset(sub_name, group_eeg_path);
        sub_ft = eeglab2fieldtrip(sub_eeg, 'raw');

        % === Downsample for speed ===
        cfg = []; cfg.resamplefs = 250;
        sub_ft = ft_resampledata(cfg, sub_ft);

        % === Demean ===
        cfg = []; cfg.demean = 'yes';
        sub_ft = ft_preprocessing(cfg, sub_ft);

        % === Compute covariance ===
        cfg = []; cfg.covariance = 'yes'; cfg.keeptrials = 'yes';
        sub_timelock = ft_timelockanalysis(cfg, sub_ft);

        % === Source localization (eLORETA) ===
        cfg = [];
        cfg.method = 'eloreta';
        cfg.sourcemodel = leadfield;
        cfg.headmodel = headmodel;
        cfg.elec = elec_aligned;
        cfg.eloreta.keepfilter = 'yes';
        cfg.eloreta.normalize = 'yes';
        cfg.eloreta.lambda = 0.05;
        cfg.eloreta.projectnoise = 'yes';
        sub_source = ft_sourceanalysis(cfg, sub_timelock);

        % === Interpolate source for visualization ===
        cfg = [];
        cfg.parameter = 'avg.pow';
        cfg.interpmethod = 'spline';
        sub_source_int = ft_sourceinterpolate(cfg, sub_source, mri);

        cfg = [];
        cfg.method = 'ortho';
        cfg.funparameter = 'avg.pow';
        cfg.maskparameter = 'mask';
        cfg.maskstyle = 'colormap';
        cfg.maskcolor = [1 0 0];
        cfg.funcolormap = parula;
        cfg.atlas = atlas;
        ft_sourceplot(cfg, sub_source_int);
        saveas(gcf, fullfile(source_location_path, group_name, sprintf('%s.png', sub_name(1:min(12,end)))));
        close all;

        % === Extract virtual signals ===
        ntrial = numel(sub_ft.trial);
        ntime = size(sub_ft.time{1}, 2);
        signal_mom_1d = cell(1, ntrial);
        f_setupParallelPool(numProcessors);

        parfor t = 1:ntrial
            signal_mom_1d{t} = NaN(length(ndipole), ntime);
            for v = 1:length(ndipole)
                f = sub_source.avg.filter{ndipole(v)};
                if isempty(f), continue; end
                moments = f * squeeze(sub_ft.trial{t});
                if any(moments(:) ~= 0)
                    u = svd(moments, 'econ');
                    m = u(:,1)' * moments;
                else
                    m = zeros(1, ntime);
                end
                signal_mom_1d{t}(v,:) = m;
            end
        end

        virtual_signal = struct;
        virtual_signal.label = arrayfun(@(n) ['S_' ROI{n}], 1:numel(ROI), 'UniformOutput', false)';
        virtual_signal.sampleinfo = sub_ft.sampleinfo;
        virtual_signal.time = sub_ft.time;
        virtual_signal.fsample = sub_ft.fsample;
        virtual_signal.ntrials = ntrial;
        virtual_signal.trial = cell(ntrial,1);

        for t = 1:ntrial
            V = zeros(numel(ROI), ntime);
            for r = 1:numel(ROI)
                idx = roi_indices{r};
                if isempty(idx), continue; end
                switch percell_type
                    case 'centroid'
                        pos = sourcemodel.pos(sourcemodel.inside,:);
                        center = round(mean(pos(idx,:), 1));
                        box = all(pos(idx,:) >= center - 10 & pos(idx,:) <= center + 10, 2);
                        vox = idx(box);
                    case 'PCA'
                        vox = idx;
                end
                data_roi = signal_mom_1d{t}(vox, :);
                if isempty(data_roi)
                    V(r,:) = zeros(1, ntime);
                else
                    if strcmp(percell_type, 'PCA')
                        [~, score] = pca(data_roi');
                        V(r,:) = score(:,1)';
                    else
                        V(r,:) = mean(data_roi, 1);
                    end
                end
            end
            virtual_signal.trial{t} = V;
        end

        save(save_file, 'virtual_signal');
    end
end
