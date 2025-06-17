
% Loop through defined groups
group = {'HC','Anxiety'};
vitual_signal_path = '...';  % Define correct virtual signal storage path
pac_sequence_path = '...';   % Define correct output path for PAC results
numProcessors = 5;
percell_type = 'centroid';
Delta_ranges = [1 4];
coi_name = 'Beta';

for group_n = 1:numel(group)
    group_name = group{group_n};
    virtual_group_path = fullfile(vitual_signal_path, group_name);
    pac_group_path = fullfile(pac_sequence_path, group_name);
    mkdir(pac_group_path);

    virtual_files = dir(fullfile(virtual_group_path, '*.mat'));
    for f = 1:numel(virtual_files)
        sub_file = virtual_files(f).name;
        sub_id = sub_file(1:min(12,length(sub_file)));
        pac_out_file = fullfile(pac_group_path, [sub_id '.mat']);
        if isfile(pac_out_file), continue; end

        virtual_signal = load(fullfile(virtual_group_path, sub_file));
        virtual_signal = virtual_signal.virtual_signal;

        % Merge trials into continuous data, then segment
        cfg = []; cfg.trials = 10:30; cfg.continuous = 'yes';
        virtual_signal_cont = ft_redefinetrial(cfg, virtual_signal);
        cfg = []; cfg.length = 5;
        virtual_signal_seg = ft_redefinetrial(cfg, virtual_signal_cont);

        % PAC computation for each channel pair
        sub_pac_sequence = struct();
        compute_band = @(c,bw)[c-0.5*bw,c+0.5*bw];
        low_fre_band = Delta_ranges;
        high_fre_band = compute_band(12, 6);  % Replace 12,6 with real FOOOF-derived center and bandwidth if needed

        tic;
        for low_ch = 1:numel(virtual_signal_seg.label)
            high_ch = low_ch;
            pair_pac_sequence = [];
            f_setupParallelPool(numProcessors);

            for t = 1:virtual_signal_seg.ntrials
                cfg = []; cfg.trials = t;
                trial_data = ft_selectdata(cfg, virtual_signal_seg);
                trial_eeglab = f_fieldtrip_to_eeglab(trial_data);

                part_virtial_pac = pop_pac(trial_eeglab, 'Channels', low_fre_band, high_fre_band, low_ch, high_ch, ...
                    'nfreqs1', 3, 'nfreqs2', 3, 'freqscale', 'log', 'method', 'instmipac', 'nboot', 200, 'alpha', [], ...
                    'bonfcorr', 0, 'tlimits', [0, trial_eeglab.xmax]);

                avg_pacval = squeeze(mean(mean(part_virtial_pac.etc.eegpac.instmipac.pacval, 1), 2));
                pair_pac_sequence = [pair_pac_sequence; avg_pacval];
            end

            sub_pac_sequence.(['LowCh_', num2str(low_ch), '_HighCh_', num2str(high_ch)]) = pair_pac_sequence;
        end
        toc;
        save(pac_out_file, 'sub_pac_sequence');
        fprintf('Saved PAC sequence: %s\n', sub_id);
    end
end
