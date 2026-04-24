function Concatenate_trials(folder,number_of_trials_per_simulations, number_independent_simulations)

    if exist('Noisy_LFP_anis.mat', 'file')
        delete('Noisy_LFP_anis.mat'); 
    end
    
    fs = 1000;
    single_trial_duration = 1000; % ms
    total_time_points = single_trial_duration * number_of_trials_per_simulations * number_independent_simulations;

    % --- PRE-ALLOCATION STEP ---
    first_file = fullfile(folder, sprintf('LFP_i%d_J%d.mat', 1,1));
    temp_data = load(first_file);
    n_channel = size(temp_data.data, 2); % Assuming [Time x Channels]
    
    % Preallocate the final master matrix
    all_Z = zeros(n_channel, total_time_points); 
    current_idx = 1;

    for k = 1:number_of_trials_per_simulations
        % Preallocate the batch matrix for the current trial
        % Dimensions: [Channels x (Time * Sims)]
        batch_width = single_trial_duration * number_independent_simulations;
        all_prep = zeros(n_channel, batch_width);
        
        for index_J = 1:number_independent_simulations
            name=fullfile(folder, sprintf('LFP_i%d_J%d.mat',index_J, k));
            
            if exist(name, 'file')
                temp = load(name);
                % Slot the data directly into its specific window
                start_col = (index_J - 1) * single_trial_duration + 1;
                end_col = index_J * single_trial_duration;
                
                all_prep(:, start_col:end_col) = temp.data'; 
            end
        end
        
        % Z-score calculation
        mu = mean(all_prep, 2, 'omitnan');
        sd = std(all_prep, 0, 2, 'omitnan');
        sd(sd == 0) = 1;
        
        Z = (all_prep - mu) ./ sd;
        
        % --- CORRECT PREALLOCATION FILL ---
        % Slot the Z-scored batch into the master matrix
        all_Z(:, current_idx : current_idx + batch_width - 1) = Z;
        current_idx = current_idx + batch_width;
    end
 
    prep_LFP = all_Z;
    save('Noisy_LFP_anis.mat', 'prep_LFP', '-v7.3');
end