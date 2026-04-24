% =========================================================================
% BETA POWER ENVELOPE ANALYSIS
% =========================================================================
% Description:
%   Computes trial-averaged beta-band power envelopes from multi-channel
%   LFP recordings, time-locked to movement onset. The dominant beta
%   frequency is estimated from the average LFP power spectrum, and the
%   signal is bandpass filtered (+/- 3 Hz around the peak) before
%   extracting the Hilbert envelope. Envelopes are normalized per channel
%   and averaged across trials.
% Data available on Zenodo: 10.5281/zenodo.19420479
% The script will prompt you to select the directory where data is stored.
% Output:
%   - Figures of single-trial and trial-averaged normalized beta envelopes
% =========================================================================

% LOAD THE FILES 
addpath(fullfile(pwd, '..',  'functions'));

dataDir = uigetdir(pwd, 'Select the folder containing your .mat files');
if ischar(dataDir)
    chanlfp = load(fullfile(dataDir, 'rs1050225_MI_clean_LFP.mat'));
    load(fullfile(dataDir, 'Beh_matrix.mat'));
else
    error('Folder not selected. Script stopped.');
end

%% extract channels
fields_lfp   = fieldnames(chanlfp);
units_lfp    = fields_lfp(startsWith(fields_lfp,'lfp'));
num_channels = length(units_lfp);
size_lfp=length(chanlfp.(units_lfp{1}));

%%
% FILTERED
fs=1000;
order=2;
width_filter=3; %Hz
n_tr=200; % number of trial for averaging
% For visualization purposes, I select this window with respect to movment
%onset
step_before=1000; %from 1000 ms before movement onset
step_after=450; %to 450 ms after movement onset
time_to_GO=-step_before:1:step_after;

%% Filter
% Frequency band

all_temp = zeros(size_lfp, num_channels);
for u = 1:num_channels
    all_temp(:,u) = chanlfp.(units_lfp{u});
end
ave_LFP=mean(all_temp,2);
[pxx,f] = pwelch(ave_LFP,[],[],[],fs);
fp_m=gaussfilt(f,pxx,1);
[p_fr peak_p]=findpeaks(fp_m);
zp=f(peak_p);
freqrange=[zp(1)-width_filter zp(1)+width_filter]; %(most prominent frequency +- 3 Hz)

%% INITIALIZATION

% Pre-allocate the matrix for normalized envelopes (Trials x Time)
% Note: This is redefined inside the channel loop for each channel

single_trial_env=NaN(n_tr,length(time_to_GO));
for k=1:n_tr  
        initial_time=ceil(beh(k,5)*fs) -step_before;
        final_time=ceil(beh(k,5)*fs)+step_after;
        del_t=final_time-initial_time;
     if (floor(beh(k,6)*fs)-floor(beh(k,5)*fs) < step_after)
         continue
     else 
        envelopes=[];
        norm_env=NaN(num_channels,length(time_to_GO));
        for u=1:num_channels
            lfp_temp = all_temp(:,u)';      
            lfp_filt = buttfilt(lfp_temp(initial_time-step_before:final_time+step_before), freqrange, fs, 'bandpass', order);
            x = lfp_filt(step_before:end-step_before);
            envelopes =abs(hilbert(x));
            norm_env(u,:)=envelopes(1:end-1)./max(envelopes(1:end-1),[],2);
        end
        single_trial_env(k,:)=mean(norm_env,1);
     end
end

%%

window_size = 50; 
filtered_envelope=nanmean(single_trial_env,1);
smoothed_envelope = smoothdata(filtered_envelope, 'gaussian', window_size);

figure
plot(time_to_GO,smoothed_envelope,'k','linewidth',3)
xlabel('Time (ms) w.r. to movement')
ylabel('Normalized beta power')
set(gca,'fontsize',18)
xlim([-step_before step_after])
 