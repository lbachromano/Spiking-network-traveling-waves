function Concatenate_trials(folder_name,Nr_trials)
    addpath(fullfile(pwd, '..', '..', 'functions'));
    
    skip_times = 100; 
    Prep_T = 1200; 
    freqrange = [5 250]; 
    fs = 1000;
    order = 4;

    if exist('Noisy_LFP_anis.mat', 'file')
        delete('Noisy_LFP_anis.mat');
    end

    % --- PRE-ALLOCATION STEP ---
    % Load the first file just to determine the number of channels
    first_file = sprintf('%s/LFP_absSum_i1_J1.mat', folder_name);

    temp_data = load(first_file);
    n_channel = size(temp_data.data', 1); 
    
    % Pre-allocate the full matrix: [Channels x (Prep_T * Total Trials)]
    all_prep = zeros(n_channel, Prep_T * Nr_trials);
    % ---------------------------

    for index_J = 1:Nr_trials
        name = fullfile(folder_name, sprintf('LFP_absSum_i%d_J1.mat', index_J));
        fprintf('Processing trial %d...\n', index_J); % Cleaner logging
        
        LFP_from_brian = load(name);
        lfp_zscored = (LFP_from_brian.data)';
        
        % Remove skip_times
        lfp_zscored(:, 1:skip_times) = [];
        
        % Filter data
        % Pre-allocate filtered_lfp for this specific trial to avoid growth
        filtered_lfp = zeros(n_channel, size(lfp_zscored, 2));
        for u = 1:n_channel
            filtered_lfp(u, :) = buttfilt(lfp_zscored(u, :), freqrange, fs, 'bandpass', order);
        end
        
        % Determine where to insert in the pre-allocated matrix
        start_col = ((index_J - 1) * Prep_T) + 1;
        end_col = index_J * Prep_T;
        
        % Store the prep period directly into the pre-allocated array
        all_prep(:, start_col:end_col) = filtered_lfp(:, 1:Prep_T);
    end

    prep_LFP = all_prep;
    save('Noisy_LFP_anis.mat', 'prep_LFP', '-v7.3'); % -v7.3 handles larger files
end