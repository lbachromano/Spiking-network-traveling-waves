% =========================================================
% Computes beta-band suppression onset timing
%   across 96 Utah array channels for an example trial,
%   and visualizes the spatial distribution across the array.
% Dependencies: Signal Processing Toolbox, Deep Learning
%   Toolbox if the autoencoder smoothing is used (autoencoder_usage=1).
%
% Data available on Zenodo: 10.5281/zenodo.19420479
% The script will prompt you to select the directory where data is stored.
% =========================================================

% LOAD THE FILES 
addpath(fullfile(pwd, '..', 'functions'));

dataDir = uigetdir(pwd, 'Select the folder containing your .mat files');
if ischar(dataDir)
    chanlfp = load(fullfile(dataDir, 'rs1050225_MI_clean_LFP.mat'));
    load(fullfile(dataDir, 'Beh_matrix.mat'));  
else
    error('Folder not selected. Script stopped.');
end


fields_lfp   = fieldnames(chanlfp);
units_lfp    = fields_lfp(startsWith(fields_lfp,'lfp'));
num_channels = length(units_lfp);
size_lfp=length(chanlfp.(units_lfp{1}));

%%

%PARAMETERS

% Threshold for detecting beta suppression onset (normalized units).
% This value is used for illustrative purposes only — the goal of this
% script is to visualize examples of suppression timing across the array.
threshold=0.8;   
fs=1000;
bias=round(0.05 * fs);  %guarding against filter edge effects 
order=2;
width_filter=3; %Hz
% Time axis in milliseconds relative to movement onset.
% Used for axis labeling and mapping results back to real time only —
% does not affect the data. Window chosen to contain the beta suppression
% period identified visually in Macaque_lfp_spectral_envelopes_move_aligned.m
time_to_GO = [-800:1:200];

% Sub-window (in array indices) restricting threshold-crossing detection
% to the peri-movement period. Chosen by hand based on visual inspection
% of the envelope plots — see Macaque_lfp_spectral_envelopes_move_aligned.m
% for justification. Corresponds to time_to_GO([620:900]) = [~-181 to ~99] ms relative to movement onset.
selected_window_in_steps = [620:900];

% Convert selected window indices to milliseconds for interpretable output.
% This bridges array indices (selected_window_in_steps) and real time (time_to_GO).
Selected_times = time_to_GO(selected_window_in_steps);

autoencoder_usage=1; % Set to 1 to use autoencoder denoising, 0 to use raw normalized envelopes.

%%
 
all_temp = zeros(size_lfp, num_channels);
for u = 1:num_channels
    all_temp(:,u) = chanlfp.(units_lfp{u});
end
% Frequency band
ave_LFP=mean(all_temp,2);
[pxx,f] = pwelch(ave_LFP,[],[],[],fs);
fp_m=gaussfilt(f,pxx,1);
[p_fr peak_p]=findpeaks(fp_m);
zp=f(peak_p);
freqrange=[zp(1)-width_filter zp(1)+width_filter]; %(most prominent frequency +- 3 Hz)
 
for n_tr=49 %example used in the paper 
    envelopes=NaN(num_channels,length(time_to_GO));
    for u=1:num_channels
        lfp_temp = all_temp(:,u)';   
        initial_time = round(ceil(beh(n_tr,5)*fs) - bias);
        final_time = round(initial_time + length(time_to_GO) - 1);       
        filtered_lfp = (buttfilt(lfp_temp(initial_time:final_time), freqrange, fs, 'bandpass', order));
        envelopes(u,:) = abs(hilbert(filtered_lfp));
    end
    figure
    plot(time_to_GO,envelopes(:,:)')
    xlabel('Time (ms) relative to movement')
    set(gca,'fontsize',18)
end


%%
selected_envelopes=envelopes(:,selected_window_in_steps);
if autoencoder_usage==1
    autoenc = trainAutoencoder(selected_envelopes,floor(size(selected_envelopes,1)/2));
    XReconstructed  = predict(autoenc,selected_envelopes);
    Normalized_reconstructed_envelopes=XReconstructed./max(XReconstructed,[],2);
    considered_envelopes=Normalized_reconstructed_envelopes;
else
    considered_envelopes=selected_envelopes./max(selected_envelopes,[],2);
end
figure
plot(time_to_GO(selected_window_in_steps)',selected_envelopes./max(selected_envelopes,[],2));
xlabel('Time (ms) relative to movement')
set(gca,'fontsize',18)
axis square
Times = NaN(num_channels,1);
for u = 1:num_channels
    sig = considered_envelopes(u, :);    
    % Find indices where the signal is ABOVE threshold AND the next point is BELOW
    % This logic ensures we only catch downward movement
    down_cross_indices = find(sig(1:end-1) >= threshold & sig(2:end) < threshold);    
    if ~isempty(down_cross_indices)
        % Grab the very first occurrence
        first_down_idx = down_cross_indices(1);
        
        % We use +1 to get the time when it is actually BELOW the threshold
        Times(u) = Selected_times(first_down_idx + 1);
    else
        % If it never starts above and goes below, return NaN
        Times(u) = NaN;
    end
end
 
figure
histogram(Times)
outlierIdx = isoutlier(Times);
Times(outlierIdx) = NaN;
LFP_phases_nan=NaN(10);
location_channels=chanlfp.MIchan2rc;
for u=1:num_channels
        r = location_channels(u,2);
        c = location_channels(u,1);
        LFP_phases_nan(r, c) = Times(u);
end

figure 
imagesc(LFP_phases_nan)
set(gca, 'YDir', 'normal')
set(gca,'fontsize',16)
axis square
colorbar()
xticks([])
yticks([])

%% Plot a subset of the array to compare it to simulations
 
interval_row=6:10;
interval_col=5:9; 
Sub_matrix_plot=LFP_phases_nan(interval_row,interval_col);

figure
Sub_matrix_plot_filled = inpaint_nans(Sub_matrix_plot, 4);  % 4 = smooth interpolation
imagesc(Sub_matrix_plot_filled)
set(gca, 'YDir', 'normal')
set(gca,'fontsize',16)
axis square
colorbar()
xticks([])
yticks([])