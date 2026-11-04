%% =========================================================
% FINITE-SIZE ANALYSIS
%
% For each network size:
%   1. Load and concatenate the 4 consecutive simulation segments.
%   2. Read and plot the LFP.
%   3. Read E and I population rates.
%   4. Compute:
%
%       r_pop(t) = 0.8*r_E(t) + 0.2*r_I(t)
%
%   5. Compute PSD of r_pop(t).
%   6. Compute oscillation strength as RMS band-limited amplitude:
%
%       A_osc = sqrt( integral_band PSD(f) df )
%
%   7. Save oscillation strength for each network size.
%   8. Plot oscillation strength vs N.
%   9. Plot oscillation strength vs 1/sqrt(N).
%
% Filenames correspond exactly to the Python script:
%
%   LFP_i0_J0.mat ... LFP_i3_J0.mat
%   RATE_E_i0_J0.mat ...
%   RATE_I_i0_J0.mat ...
%
% J0 corresponds to Size_x = 500 um
% ...
% J8 corresponds to Size_x = 2500 um
% =========================================================

clear;
close all;
clc;

%% ---------------------------------------------------------
% PARAMETERS
% ----------------------------------------------------------

Size_list = [500, 750, 1000, 1250, 1500, ...
             1750, 2000, 2250, 2500];       % um

n_sizes = length(Size_list);

% Four 500-ms runs are saved for each network
n_iterations = 4;
seed_list = 1:10;
n_seeds = length(seed_list);
segment_duration_s = 0.5;

% LFP StateMonitor has dt = 1 ms in the Python code
fs_LFP = 1000;                               % Hz

% ----------------------------------------------------------
% SET THIS:
% Frequency interval used to quantify oscillation strength
% Example syntax:
% osc_band = [f_low f_high];
% ----------------------------------------------------------

osc_band = [15 100];

if isempty(osc_band)
    error(['Set osc_band at the beginning of the script, e.g. ', ...
           'osc_band = [f_low f_high];']);
end

% Welch window length.
% With 2 s of total rate data, 1 s gives useful frequency resolution.
welch_window_s = 1.0;


%% ---------------------------------------------------------
% NUMBER OF NEURONS
% ----------------------------------------------------------
% Python:
%
% Total_nr_neurons = 28000 * (Volume/mm^2)
%
% and Size_x = Size_y.

N_total = 28000 .* (Size_list ./ 1000).^2;


%% ---------------------------------------------------------
% STORAGE
% ----------------------------------------------------------

% Rows = network sizes
% Columns = independent random initializations
osc_strength = nan(n_sizes,n_seeds);
osc_power    = nan(n_sizes,n_seeds);
peak_freq    = nan(n_sizes,n_seeds);
fs_rate_all  = nan(n_sizes,n_seeds);

LFP_osc_strength = nan(n_sizes,n_seeds);
LFP_osc_power    = nan(n_sizes,n_seeds);
LFP_peak_freq    = nan(n_sizes,n_seeds);

%% =========================================================
% LOOP OVER NETWORK SIZES
% =========================================================

for isize = 1:n_sizes

    Size_x = Size_list(isize);
    Jindex = isize - 1;

    for iseed = 1:n_seeds

        seed_id = seed_list(iseed);

        rate_E = [];
        rate_I = [];
        LFP=[];
        samples_per_rate_segment = [];

        for iteration = 0:n_iterations-1

            % E population rate
            filename_E = sprintf('RATE_E_i%d_J%d_%d.mat', ...
                                 iteration, Jindex, seed_id);

            tmp = load(filename_E, 'data');
            current_E = double(tmp.data(:));

            % I population rate
            filename_I = sprintf('RATE_I_i%d_J%d_%d.mat', ...
                                 iteration, Jindex, seed_id);

            tmp = load(filename_I, 'data');
            current_I = double(tmp.data(:));

            if length(current_E) ~= length(current_I)
                error(['E and I rate files have different lengths ', ...
                       'for J%d, seed %d, iteration %d.'], ...
                       Jindex, seed_id, iteration);
            end

            if isempty(samples_per_rate_segment)
                samples_per_rate_segment = length(current_E);
            elseif length(current_E) ~= samples_per_rate_segment
                error('Rate segments have inconsistent lengths.');
            end

            filename_LFP = sprintf('LFP_i%d_J%d_%d.mat', ...
                                 iteration, Jindex, seed_id);

            tmp = load(filename_LFP, 'data');
            current_LFP = tmp.data;
            rate_E = [rate_E; current_E];
            rate_I = [rate_I; current_I];
            LFP=[LFP; current_LFP];

        end


        % Sampling frequency
        fs_rate = samples_per_rate_segment / segment_duration_s;
        fs_rate_all(isize,iseed) = fs_rate;
        rate_pop = 0.8 .* rate_E + 0.2 .* rate_I;
        
        %% PSD rate
        rate_pop_zero_mean = rate_pop - mean(rate_pop);
        window_samples = round(welch_window_s * fs_rate);
        window_samples = min(window_samples, ...
                             length(rate_pop_zero_mean));
        overlap_samples = floor(window_samples / 2);
        nfft = window_samples;
        [Pxx, f] = pwelch(rate_pop_zero_mean, ...
                          hamming(window_samples), ...
                          overlap_samples, ...
                          nfft, ...
                          fs_rate);
        band_idx = f >= osc_band(1) & f <= osc_band(2);
        if sum(band_idx) < 2
            error('Not enough frequency bins inside osc_band.');
        end

        osc_power(isize,iseed) = trapz(f(band_idx), Pxx(band_idx));
        osc_strength(isize,iseed) = sqrt(osc_power(isize,iseed));
        f_band = f(band_idx);
        P_band = Pxx(band_idx);

        [~, max_idx] = max(P_band);
        peak_freq(isize,iseed) = f_band(max_idx);

        %% PSD LFP

        LFP = LFP - mean(LFP,1);
         
        %% PSD LFP

        % LFP is [time x channels]
        % Each channel has already been demeaned over the full 2-s recording.
        
        window_samples_LFP = round(welch_window_s * fs_LFP);
        window_samples_LFP = min(window_samples_LFP, size(LFP,1));        
        overlap_samples_LFP = floor(window_samples_LFP / 2);
        nfft_LFP = window_samples_LFP;
        % Pxx_LFP has dimensions:
        %       frequency x channel
        [Pxx_LFP, f_LFP] = pwelch(LFP, ...
                                  hamming(window_samples_LFP), ...
                                  overlap_samples_LFP, ...
                                  nfft_LFP, ...
                                  fs_LFP);
        
        % Average PSD across the 36 channels
        Pxx_LFP_mean = mean(Pxx_LFP, 2);        
        % Frequency band of interest
        band_idx_LFP = f_LFP >= osc_band(1) & ...
                       f_LFP <= osc_band(2);
        
        if sum(band_idx_LFP) < 2
            error('Not enough frequency bins inside osc_band for LFP.');
        end
        
        % Oscillatory power from channel-averaged PSD
        LFP_osc_power(isize,iseed) = ...
            trapz(f_LFP(band_idx_LFP), ...
                  Pxx_LFP_mean(band_idx_LFP));
        
        LFP_osc_strength(isize,iseed) = ...
            sqrt(LFP_osc_power(isize,iseed));
        
        % Peak frequency from channel-averaged PSD
        f_band_LFP = f_LFP(band_idx_LFP);
        P_band_LFP = Pxx_LFP_mean(band_idx_LFP);
        [~, max_idx] = max(P_band_LFP);
        LFP_peak_freq(isize,iseed) = ...
            f_band_LFP(max_idx);

    end
    figure
    plot(LFP)
end

mean_osc_strength = mean(osc_strength, 2, 'omitnan');
std_osc_strength  = std(osc_strength, 0, 2, 'omitnan');
n_valid = sum(~isnan(osc_strength), 2);
sem_osc_strength = std_osc_strength ./ sqrt(n_valid);
mean_osc_power = mean(osc_power, 2, 'omitnan');
mean_peak_freq = mean(peak_freq, 2, 'omitnan');

LFP_mean_osc_strength = mean(LFP_osc_strength, 2, 'omitnan');
LFP_std_osc_strength  = std(LFP_osc_strength, 0, 2, 'omitnan');
LFP_n_valid = sum(~isnan(LFP_osc_strength), 2);
LFP_sem_osc_strength = LFP_std_osc_strength ./ sqrt(n_valid);
LFP_mean_osc_power = mean(LFP_osc_power, 2, 'omitnan');
LFP_mean_peak_freq = mean(LFP_peak_freq, 2, 'omitnan');

%% =========================================================
% SAVE RESULTS
% =========================================================


save('FiniteSize_nu1.3.mat','N_total', 'mean_osc_strength', 'std_osc_strength', 'LFP_mean_osc_strength', 'LFP_std_osc_strength')

%% =========================================================
% FINITE-SIZE SUMMARY
% =========================================================

figure('Color','w');
subplot(1,2,1)
errorbar(N_total, mean_osc_strength, std_osc_strength, ...
         'o-', 'LineWidth', 1.3, 'MarkerSize', 7);
xlabel('Number of neurons, N');
ylabel('Oscillation strength (Hz)');
% title('Oscillation strength vs network size');
grid on; 
% Linear fit:
% A_osc = a / sqrt(N) + b
subplot(1,2,2)
x_finite = 1 ./ sqrt(N_total);
% Sort for plotting
[x_sorted, idx] = sort(x_finite);
y_sorted = mean_osc_strength(idx);
yerr_sorted = std_osc_strength(idx);
errorbar(x_sorted, y_sorted, yerr_sorted, ...
         'o', 'MarkerSize', 7);
hold on;
% Linear fit:
% A_osc = a / sqrt(N) + b
p = polyfit(x_finite, mean_osc_strength', 1);
x_fit = linspace(min(x_finite), max(x_finite), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, '-', 'LineWidth', 1.3);
xlabel('1 / sqrt(N)');
ylabel('Oscillation strength (Hz)');
% title('Test of finite-size scaling');
grid on;
legend('Simulation', ...
       sprintf('Linear fit: A = %.3g / sqrt(N) + %.3g', ...
               p(1), p(2)), ...
       'Location', 'best');
title('rates')

figure('Color','w');
subplot(1,2,1)
errorbar(N_total, LFP_mean_osc_strength, LFP_std_osc_strength, ...
         'o-', 'LineWidth', 1.3, 'MarkerSize', 7);
xlabel('Number of neurons, N');
ylabel('Oscillation strength (Hz)');
% title('Oscillation strength vs network size');
grid on; 

grid on; 
% Linear fit:
% A_osc = a / sqrt(N) + b
subplot(1,2,2)
x_finite = 1 ./ sqrt(N_total);
% Sort for plotting
[x_sorted, idx] = sort(x_finite);
y_sorted = LFP_mean_osc_strength(idx);
yerr_sorted =  LFP_std_osc_strength(idx);
errorbar(x_sorted, y_sorted, yerr_sorted, ...
         'o', 'MarkerSize', 7);
hold on;
% Linear fit:
% A_osc = a / sqrt(N) + b
p = polyfit(x_finite,  LFP_mean_osc_strength', 1);
x_fit = linspace(min(x_finite), max(x_finite), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, '-', 'LineWidth', 1.3);
xlabel('1 / sqrt(N)');
ylabel('Oscillation strength (Hz)');
% title('Test of finite-size scaling');
grid on;
legend('Simulation', ...
       sprintf('Linear fit: A = %.3g / sqrt(N) + %.3g', ...
               p(1), p(2)), ...
       'Location', 'best');
title('LFP')