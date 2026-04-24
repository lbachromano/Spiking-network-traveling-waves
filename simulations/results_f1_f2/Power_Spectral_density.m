
% Computes and plots the z-scored power spectral density (PSD) from
% simulations (when external rate is at baseline) and compares it to the one from data 
% during the preparation period. Wavelet spectra are averaged
% across channels and runs, normalized to the beta band (15–35 Hz), and
% displayed with ±1 SD shading for comparison.

clear

addpath(fullfile(pwd, '..', '..', 'functions'));
N=10; % Number of trials to consider
folder_name='examples_paper'; %Folder with examples of LFPs from simulations

%% params
% We normalize both the power spectral density from data and simulation
% with a z-score normalization (zero mean, unit variance) computed over
% the beta band, using the frequencies in the range defined by th_low and thr_h


skip_initial_window=100; %skip the first 100 ms in line with the other processing code
T_end_prep=900; %time when the external input rate starts to increase; see Plot_external_rate.m

th_low = 15; 
thr_h  = 35;
v_low  = 1; %for visualization
v_high = 50; % for visualization
fs     = 1000;
order  = 2;
fwave  = logspace(log10(1), log10(100), 200);  
buffer = 200;

%% collect channel-averaged spectra per run

prep=[];
for index_J = 1:N
    % --- load & slice ---
    name = fullfile(folder_name, sprintf('LFP_absSum_i%d_J1.mat', index_J));
  
    S    = load(name);
    lfp  = S.data.';                       % [nChan x T]
    one_prep = lfp(:,skip_initial_window:T_end_prep);             
    prep=[prep one_prep];
end

%%    


% --- wavelet power per channel (log10), then average over time ---
    nChan = size(prep,1);
    chan_spec = zeros(nChan, numel(fwave));  % [nChan x nFreq]
    for i = 1:nChan
        x = prep(i,:);
        [~, pow] = multiphasevec3(fwave, x, fs, 6);   % [nFreq x 1 x T]
        pow = log10(pow(:,:,buffer+1:end-buffer));    % trim edges
        chan_spec(i,:) = mean(real(squeeze(pow)), 2).';  % time-avg
    end

    % --- average across channels for this run ---
    run_means = mean(chan_spec, 1);          % 1 x nFreq
     


%% average across runs (equal weight), then normalize once
sim_pw_p = mean(run_means, 1);    
std_dev_pw_p=std(run_means, 1);% 1 x nFreq
idx_sim  = fwave >= th_low & fwave <= thr_h;

% z-score (band-based)
mu_sims = mean(sim_pw_p(idx_sim));
sd_sims = std( sim_pw_p(idx_sim) );
final_z_prep = (sim_pw_p - mu_sims) ./ sd_sims;
 
%% load & normalize data (analogous)
ld    = load('LFP_prep_power_spect_1_50.mat');  % has fields f, ff
fr_p  = ld.f; 
pw_p  = ld.ff;
std_d=ld.std_freq;
idx_d = fr_p >= th_low & fr_p <= thr_h;
mu_d  = mean(pw_p(idx_d));
sd_d  = std( pw_p(idx_d) );
z_pw_p = (pw_p - mu_d) / sd_d;
z_std_p= std_d./ sd_d;


%%


df    = median(diff(fr_p));
FWHM  = 2;                            % Hz
sigma = 2*FWHM/2.355;
R     = ceil(6*sigma/df);             % ~±3σ
k     = (-R:R)*df;
g     = exp(-0.5*(k/sigma).^2); g = g/sum(g);
pw_p_smooth = conv(z_pw_p, g, 'same');



% --- std band in z-score units ---
z_std = (std_dev_pw_p ./ sd_sims)./sqrt(nChan);      % std transforms by division (shift doesn't affect std)
yl = final_z_prep - z_std;
yu = final_z_prep + z_std;
dyl = pw_p_smooth - z_std_p;
dyu = pw_p_smooth + z_std_p;
figure
fill([fwave fliplr(fwave)], [yl fliplr(yu)], 'b', ...
     'FaceAlpha', 0.15, 'EdgeColor', 'none'); 
hold on
plot(fwave, final_z_prep, 'b', 'LineWidth', 2);
hold on
fill([fr_p fliplr(fr_p)], [dyl fliplr(dyu)], 'k', ...
     'FaceAlpha', 0.15, 'EdgeColor', 'none');  
hold on
plot(fr_p,  pw_p_smooth,       'k', 'LineWidth', 2);
legend('sim ±1 SD','sim mean','data','Location','best')
xlim([v_low v_high]); xlabel('Frequency (Hz)'); ylabel('Z-scored PSD');
set(gca,'FontSize',14);  

