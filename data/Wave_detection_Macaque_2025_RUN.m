% =========================================================================
% WAVE DETECTION AND STATISTICS - MACAQUE LFP DATA
% =========================================================================
%
% DESCRIPTION:
%   This script detects traveling waves in macaque local field potential
%   (LFP) recordings and computes their kinematic statistics (speed,
%   wavelength, direction). Option to visualize a few detected waves.
%
%   Detection results are saved to disk and subsequently loaded for
%   visualization of wave velocity, wavelength, and direction distributions.
%
% PIPELINE:
%   1. Load LFP, spike, and behavioral data
%   2. Run wave_detection_Macaque_2025_test_R() across parameter grid
%   3. Load saved results and plot wave statistics
%
% INPUT DATA:
% Data available on Zenodo: 10.5281/zenodo.19420479
% The script will prompt you to select the directory where data is stored.
%
% KEY PARAMETERS:
%   number_trials   : number of trials to analyse  (default: 100)
%   duration        : trial window length in ms from trial onset (default: 1200)
%   width_filters   : spatial smoothing kernel width in electrodes (default: 3)
%   threshol_inputs : PGD threshold percentile for wave detection (default: 99.99)
%   th_es           : residual error threshold percentile (default: 20)
%   block_sz        : spatial block size for surrogate computation (default: 1)
%   K_surr_per_frame: number of surrogates per frame (default: 2)
%
% OUTPUT:
%   .mat file named: Macaque_R_trials_<N>_DWfilter_<wf>_thPGD_<th>_the_<the>.mat
%   containing a structure with fields:
%       .all_speeds       : wave propagation speeds (cm/s)
%       .wave_wavelengths : spatial wavelengths (cm)
%       .wave_directions  : propagation directions (radians)
%
% USAGE:
%   Set parameters in the "Define parameters" section, then run the script.
%   Particularly the number of waves to plot: max_nr_waves_to_plot
% =========================================================================

%% Define the parameters that will be used for wave detection

number_trials = 10; 
number_trials=max(2,number_trials);%for reliable frequency computation, use at least 2
duration = 1200; %ms from start trial
max_nr_waves_to_plot=3;
use_modes = {'spatialblock'};  % or {'phase'} etc.
K_surr_per_frame=2;
block_sz = 1; % spatial block size for surrogate shuffling (in electrodes).
              % 1 = fully random shuffle (least conservative, breaks all structure).
              % 2 = shuffles 2x2 tiles (preserves local within-block gradients).
              % 5 = shuffles 5x5 tiles (very conservative; grid must be divisible).
              % Larger blocks = more conservative null; default 1 for 10x10 arrays.
              % NOTE: grid size must be divisible by block_sz.
width_filters    = 3;
threshol_inputs  = 99.99;   % PGD percentile
th_es            = 20;         % residual percentile



% LOAD THE FILES 
addpath(fullfile(pwd, '..',  'functions'));
dataDir = uigetdir(pwd, 'Select the folder containing your .mat files');
if ischar(dataDir)
    chanlfp = load(fullfile(dataDir, 'rs1050225_MI_clean_LFP.mat'));
    load(fullfile(dataDir, 'Beh_matrix.mat'));
else
    error('Folder not selected. Script stopped.');
end

  
for wf = width_filters
    for thw = threshol_inputs
        for the = th_es
            try
 
            wave_detection_Macaque_2025_test_R(chanlfp, beh, number_trials, duration, wf, thw, the, max_nr_waves_to_plot, use_modes, K_surr_per_frame, block_sz);

                tag = sprintf('WF_%g_PGD_%g_E_%g', wf, thw, the);
                tag = strrep(tag, '.', 'p');  % e.g., 99.9 -> 99p9
                fprintf('Done: %s\n', tag);

            catch ME
                warning('Failed for wf=%g, thPGD=%g, thE=%g: %s', wf, thw, the, ME.message);
            end
        end
    end
end

%% Load and Plot the resultss scomputed from 100 trials
number_trials=100; 
width_filter    = 3;          %[3 6];  
threshold_input  = 99.99;         % [95 97 99 99.9]  
th_es            = 20;         % [5 15 20]  

titlefile=sprintf('Macaque_R_trials_%d_DWfilter_%d_thPGD_%.2f_the_%d.mat',number_trials,width_filter,threshold_input,th_es);
ld=load(titlefile); 
Macaque_final=ld.Macaque_final;
wave_selected_speed=Macaque_final.all_speeds;

p_low  = prctile(wave_selected_speed, 1);
p_high = prctile(wave_selected_speed, 95);
clean_speeds = wave_selected_speed(wave_selected_speed >= p_low & wave_selected_speed <= p_high);


figure;
histogram(wave_selected_speed, 100, 'Normalization','probability')
xlabel('Wave Velocity (cm/s)'); ylabel('Probability'); title('Distribution of Wave Velocities'); grid on

figure;
histogram(clean_speeds, 100, 'Normalization','probability')
xlabel('Wave Velocity [outliers removed] (cm/s)'); ylabel('Probability'); title('Distribution of Wave Velocities'); grid on

figure;
histogram(Macaque_final.wave_wavelengths, 100, 'Normalization','probability')
xlabel('Wave length (cm)'); ylabel('Probability')

figure;
polarhistogram(Macaque_final.wave_directions, 15)
title('Polar Histogram of Wave Directions')
  