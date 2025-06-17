%% STEP 1 - EEG PREPROCESSING
clc; clear variables;
restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED;

project_root = '<PROJECT_ROOT>';
toolbox_root = '<MATLAB_TOOLBOX>';

% Add EEGLAB
addpath(fullfile(toolbox_root, 'eeglab2021.0'));

% Add FieldTrip and external toolboxes
fieldtrip_path = fullfile(toolbox_root, 'fieldtrip-master');
external_subfolders = {'mne', 'brewermap', 'matplotlib', 'cmocean', 'colorcet', 'gifti', 'spm12'};
cellfun(@(f) addpath(fullfile(fieldtrip_path, 'external', f)), external_subfolders);
addpath(fieldtrip_path); ft_defaults;

% Add project paths
addpath(fullfile(project_root, 'Code'));
addpath(fullfile(project_root, 'Data'));
addpath(fullfile(project_root, 'Tools'));
addpath(fullfile(project_root, 'Result'));
cd(fullfile(project_root, 'Code'));

% Define groups
group_names = {'HC','Anxiety','ASD','Depression','LD','ADHD'};
for g = 1:numel(group_names)
    group = group_names{g};
    
    % Input paths
    input_root = fullfile(project_root, 'Data');
    input_paths = {fullfile('EEG', group)};
    f_createOutputFolders(input_root, input_paths);
    
    % Tool paths
    input_tool_root = fullfile(project_root, 'Tools');
    input_tool_paths = {'souce_relevent', 'excel_relevent', 'sensor_relevent'};
    
    % Output paths
    output_root = fullfile(project_root, 'Result', 'Preprocess_result');
    output_paths = {
        fullfile('EEG_pure', group), ...
        fullfile('ICA', group), ...
        fullfile('Bad_channel_fig', group), ...
        fullfile('Adjust_report_file', group), ...
        fullfile('Process_Log', group)
    };
    f_createOutputFolders(output_root, output_paths);
    
    % Set full paths
    EEG_raw_path = fullfile(input_root, input_paths{1});
    EEG_pure_path = fullfile(output_root, output_paths{1});
    ICA_path = fullfile(output_root, output_paths{2});
    Bad_channel_fig_path = fullfile(output_root, output_paths{3});
    Adjust_report_file_path = fullfile(output_root, output_paths{4});
    Processing_log_path = fullfile(output_root, output_paths{5});
    Processing_log_bad_seg_path = fullfile(Processing_log_path, 'bad_seg');
    mkdir(Processing_log_bad_seg_path);
    
    % Tools
    Channel_location = fullfile(input_tool_root, input_tool_paths{3}, 'HBN.sfp');
    
    % List subjects
    sub_dirs = dir(EEG_raw_path); sub_dirs(1:2) = [];
    sub_names = string({sub_dirs.name}');
    fprintf('Found %d subjects in group %s.\n', numel(sub_names), group);
    
    % Log structure
    bdataLog = table('Size', [0 2], 'VariableTypes', {'string', 'cell'}, ...
        'VariableNames', {'SubjectName', 'BadChannels'});
    
    %% Begin subject loop
    eeglab nogui;
    for s = 1:numel(sub_names)
        sub_name = sub_names{s};
        EEG_pure_name = [sub_name '_pure.set'];
        EEG_pure_file = fullfile(EEG_pure_path, EEG_pure_name);
        
        if exist(EEG_pure_file, 'file') == 2
            fprintf('Subject %s already processed. Skipping.\n', sub_name);
            continue;
        end
        
        data_file = fullfile(EEG_raw_path, sub_name, 'RestingState.mat');
        if ~isfile(data_file)
            f_log_bad_data([], sub_name, s, Processing_log_path, 'empty_data');
            continue;
        end
        
        try
            EEG = load(data_file).EEG;
        catch ME
            fprintf('Error loading %s: %s\n', data_file, ME.message);
            f_log_bad_data([], sub_name, s, Processing_log_path, 'no_EEG');
            continue;
        end
        EEG = eeg_checkset(EEG);

        % Fix latency if missing
        if ~isfield(EEG.event, 'latency') && isfield(EEG.event, 'sample')
            for e = 1:length(EEG.event)
                EEG.event(e).latency = EEG.event(e).sample;
            end
        end

        % Find event markers
        event_types = {EEG.event.type};
        event_samples = [EEG.event.sample];
        idx_30 = find(strcmp(event_types, '30  ')); % Eyes closed
        idx_20 = find(strcmp(event_types, '20  ')); % Eyes open
        idx_90 = find(strcmp(event_types, '90  ')); % Marker
        
        if isempty(idx_30) || isempty(idx_20) || isempty(idx_90)
            f_log_bad_data(EEG, sub_name, s, Processing_log_path, 'bad_event');
            continue;
        end

        start_closed = idx_30(idx_30 > max(idx_90));
        start_open = idx_20(idx_20 > max(idx_90));

        if numel(start_closed) < 4 || numel(start_open) < 5
            f_log_bad_data(EEG, sub_name, s, Processing_log_path, 'insufficient_events');
            continue;
        end

        closed_segments = {};
        for i = 1:min(length(start_closed), length(start_open) - 1)
            s1 = event_samples(start_closed(i));
            s2 = event_samples(start_open(i+1));
            if s1 < s2
                closed_segments{end+1} = pop_select(EEG, 'point', [s1, s2]);
            end
        end

        if isempty(closed_segments)
            f_log_bad_data(EEG, sub_name, s, Processing_log_path, 'no_closed_segment');
            continue;
        end

        EEG = closed_segments{1};
        for i = 2:length(closed_segments)
            EEG = pop_mergeset(EEG, closed_segments{i}, 'noeventwarning');
        end
        EEG = pop_editeventvals(EEG, 'delete', find(strcmp({EEG.event.type}, 'boundary')));

        % Set channel location
        EEG.event = [];
        EEG = pop_chanedit(EEG, 'load', {Channel_location, 'filetype', 'sfp'});
        EEG = eeg_checkset(EEG);

        % Filtering
        EEG = pop_eegfiltnew(EEG, 'locutoff', 0.5, 'hicutoff', 45);
        EEG = pop_eegfiltnew(EEG, 'locutoff', 59, 'hicutoff', 61, 'revfilt', 1);

        % Trim edges
        EEG = pop_select(EEG, 'time', [2 EEG.xmax - 2]);

        % Remove bad segments
        ori_size = size(EEG.data, 2);
        [~, bad_interval] = pop_rejcont(EEG, 'freqlimit', [0.5 45], ...
            'threshold', 10, 'epochlength', 0.5, 'contiguous', 4, ...
            'addlength', 0.25, 'taper', 'hamming', 'eegplot', 'off');
        EEG = eeg_eegrej(EEG, bad_interval);

        if size(EEG.data, 2) < 0.6 * ori_size
            f_log_bad_data(EEG, sub_name, s, Processing_log_path, 'excessive_rejection');
            continue;
        end

        EEG.event = [];
        EEG = pop_select(EEG, 'nochannel', {'E1','E8','E14','E17','E21','E25','E32','E48','E49','E56','E63',...
                                            'E68','E73','E81','E88','E94','E99','E107','E113','E119','E125','E126','E127','E128'});
        
        % Clean channels
        [EEG_Clean, removed_channels] = f_myclean_channels(EEG);
        bad_idx = find(removed_channels);
        if numel(bad_idx) > round(0.2 * EEG.nbchan)
            f_log_bad_data(EEG, sub_name, s, Processing_log_path, 'too_many_bad_channels');
            continue;
        end
        EEG_ori_locs = EEG.chanlocs;
        if ~isempty(bad_idx)
            EEG = pop_select(EEG, 'nochannel', bad_idx);
        end

        if EEG.srate < 500
            f_log_bad_data(EEG, sub_name, s, Processing_log_path, 'low_sampling_rate');
            continue;
        end

        % ICA rank
        EEG_rank = f_getrank(EEG.data);
        if EEG_rank < EEG.nbchan - 5
            f_log_bad_data(EEG, sub_name, s, Processing_log_path, 'low_rank');
            continue;
        end

        % ICA decomposition
        ica_file = fullfile(ICA_path, [sub_name '_ICA.set']);
        if ~isfile(ica_file)
            if isempty(bad_idx)
                EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1);
            else
                EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, 'pca', EEG_rank);
            end
            pop_saveset(EEG, 'filename', [sub_name '_ICA.set'], 'filepath', ICA_path);
        end

        % ICA artifact rejection
        EEG = pop_loadset('filename', [sub_name '_ICA.set'], 'filepath', ICA_path);
        EEG = interface_ADJ(EEG, fullfile(Adjust_report_file_path, sub_name), 'report.txt');
        EEG = iclabel(EEG, 'default');
        EEG = pop_icflag(EEG, [NaN NaN;0.5 1;0.5 1;0.5 1;0.5 1;0.5 1;NaN NaN]);
        comp_auto = find(EEG.reject.gcompreject == 1);
        comp_manual = [];
        if isfile(fullfile(Adjust_report_file_path, sub_name, 'List_Filtering.mat'))
            tmp = load(fullfile(Adjust_report_file_path, sub_name, 'List_Filtering.mat'));
            comp_manual = unique([tmp.blink, tmp.vert, tmp.horiz, tmp.disc])';
        end
        bad_comp = unique([comp_manual; comp_auto]);
        if numel(bad_comp) > 0.5 * EEG_rank
            f_log_bad_data(EEG, sub_name, s, Processing_log_path, 'excessive_ica');
            continue;
        end
        EEG = pop_subcomp(EEG, bad_comp, 0);

        % Interpolate
        if ~isempty(bad_idx)
            EEG = pop_interp(EEG, EEG_ori_locs, 'spherical');
        end

        % Re-reference
        EEG = pop_reref(EEG, []);

        % Epoch
        if strcmpi(data_type, 'rest')
            EEG = eeg_regepochs(EEG, 'recurrence', 2, 'limits', [0 2], ...
                                'extractepochs', 'on', 'rmbase', NaN);
        else
            EEG = pop_epoch(EEG, {'mark'}, [-1 2]);
        end

        % Reject bad epochs
        [EEG, rej] = pop_eegthresh(EEG, 1, 1:EEG.nbchan, -100, 100, EEG.xmin, EEG.xmax, 0, 0);
        if ~isempty(rej) && numel(rej) > 0.8 * size(EEG.data, 3)
            f_log_bad_data(EEG, sub_name, s, Processing_log_path, 'excessive_bad_epochs');
            continue;
        end
        EEG = pop_rejepoch(EEG, rej, 0);

        if numel(EEG.epoch) < 70
            f_log_bad_data(EEG, sub_name, s, Processing_log_path, 'insufficient_epochs');
            continue;
        end

        % Select best 60 continuous epochs
        all_trials = setdiff(1:EEG.trials, rej);
        if numel(all_trials) >= 60
            best = 1; min_diff = Inf;
            for i = 1:(numel(all_trials)-59)
                seq = all_trials(i:i+59);
                d = sum(diff(seq));
                if d < min_diff
                    min_diff = d;
                    best = i;
                end
            end
            selected = all_trials(best:best+59);
        else
            selected = all_trials;
        end
        EEG = pop_select(EEG, 'trial', selected);
        EEG = eeg_checkset(EEG);

        % Save result
        EEG = pop_saveset(EEG, 'filename', EEG_pure_name, 'filepath', EEG_pure_path);
        f_log_bad_data(EEG, sub_name, s, Processing_log_path, ...
                       ['bad_cha_ ' num2str(bad_idx)]);
        close all;
    end
end
