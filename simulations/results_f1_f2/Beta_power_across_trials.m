

%% LFP Beta Envelope Extraction Pipeline
% Processes raw LFP data to extract and normalize the amplitude 
% envelope of Beta oscillations.
%
% Methodology: 
%   1. Data Concatenation: Loads N trials to estimate dominant frequency via PSD.
%   2. Filtering: Applies a bandpass Butterworth filter based on dominant frequency.
%   3. Envelope Extraction: Uses Hilbert transform for instantaneous amplitude.
%   4. Normalization: Per-channel, per-trial min-max scaling to [0,1] to
%      account for impedance variations across channels.
%
% Inputs: 
%   - folder_name: Directory containing LFP_absSum_i[index]_J1.mat files
%   - N: Number of trials
%
% Outputs:
%   - Plots: Average envelope across channels and trials and Gaussian-smoothed visualization
% -------------------------------------------------------------------------


clear
addpath(fullfile(pwd, '..', '..', 'functions'));
folder_name='examples_current_paper';   %used as example here, replace 

%% Parameters
N=10;
skip_initial_window=50; %skip the first 100 ms containing artifacts

Tpreparation=1000; %baseline external input; see Plot_external_rate.m
fs=1000;
order=2; %for filtering
beta_range=[13,30]; %range where to look for the dominant peak
width_filter=3;
s2=0.08; %for Gaussian filtering for second figure


%% FILTER
% Concatenate trials at baseline external input to determine the dominant
% frequency of oscillation

name = fullfile(folder_name, sprintf('LFP_absSum_i%d_J1.mat',1));
lfp_temp=(load(name).data)';
n_channel=size(lfp_temp,1);
T=size(lfp_temp,2)-skip_initial_window; %lengths are guaranteed uniform

lfp_prep_temp = NaN(n_channel, Tpreparation - skip_initial_window, N);  % 800
for index_J=1:N
    name = fullfile(folder_name, sprintf('LFP_absSum_i%d_J1.mat', index_J));
    lfp_temp=(load(name).data)';
    lfp_prep_temp(:,:,index_J) = lfp_temp(:, skip_initial_window+1 : Tpreparation);
    
end
% Reshape the 3D array (n_channel x time_points x N) 
% into a 2D array (n_channel x new_T)
T_concatenation = (Tpreparation - skip_initial_window) * N;
lfp_preparation = reshape(lfp_prep_temp, [n_channel, T_concatenation]);

%Use a smaller nperseg if duration is short to allow more averaging
nperseg  = min(floor(T_concatenation/5), 512); 
noverlap = floor(nperseg/2);
nfft     = 2^nextpow2(max(1024, nperseg));
% 1. Detrend and compute PSD
x_mean = detrend(mean(lfp_preparation, 1, 'omitnan'));
[Pm, f_psd] = pwelch(x_mean, hamming(nperseg), noverlap, nfft, fs);
% 2. Log-transform or Flatten (to handle 1/f)
Pm_log = 10*log10(Pm); 
% 3. Define Beta Search Window
beta_idx = (f_psd >= beta_range(1)) & (f_psd <= beta_range(2));
f_beta = f_psd(beta_idx);
p_beta = Pm_log(beta_idx);
% 4. Find the local peak within that window
[pks, locs] = findpeaks(p_beta, f_beta, 'SortStr', 'descend');
if ~isempty(pks)
    domFreq_LFP = locs(1); % The highest local peak in the Beta band
else
    % Fallback: if no local peak, take the maximum in range
    [~, k] = max(p_beta);
    domFreq_LFP = f_beta(k);
end
freqrange = [max(0.5, domFreq_LFP - width_filter), min(fs/2 - 1e-6, domFreq_LFP + width_filter)];
fprintf('LFP-only band centered at %.2f Hz.\n', domFreq_LFP);



%% Compute envelopes

mean_envelope = zeros(N, T);
for index_J=1:N

    name = fullfile(folder_name, sprintf('LFP_absSum_i%d_J1.mat', index_J));
    lfp_temp=(load(name).data)';
    lfp_temp(:,1:skip_initial_window)=[]; 
    envelopes=zeros(n_channel,T);
    for u=1:n_channel
        filtered_lfp=(buttfilt(lfp_temp(u,:),freqrange,fs,'bandpass',order));
        x =filtered_lfp(1:T);
        envelopes(u,:)=abs(hilbert(x));
    end
    normalized_full_envelopes = envelopes ./ max(envelopes, [], 2);
    mean_envelope(index_J,:)=mean(normalized_full_envelopes,1);

end
 

%% Plots
% Plot original data
figure
plot(mean(mean_envelope,1),'k','linewidth',3)
xlabel('Time (ms)')
ylabel('Averaged Beta Envelope')
set(gca,'fontsize',15)
pbaspect([1 1 1])
 
%%
% Apply Gaussian filter to smooth the signal
t = (0:T-1)./fs;
t_ms = (0:T-1);
smoothed_data = gaussfilt(t,mean(mean_envelope,1),s2);
figure
plot(t_ms,smoothed_data, 'k', 'linewidth', 3)  % Smoothed data
xlabel('Time (ms)')
ylabel('Averaged Beta Envelope')
set(gca, 'fontsize', 15)
xlim([0 1800])
ylim([0 0.75])