# Supplementary Code — Beta Traveling Waves in Motor Cortex

This repository contains the analysis and simulation code accompanying the paper:

> **A Spatially Structured Spiking Network Model of Beta Traveling Waves and Their Attenuation in Motor Cortex**


## Dependencies

- **Python 3** with [Brian2] (network simulations)
- **MATLAB** (data analysis and plots)
                                                        

### Data files

| File | Description |
|------|-------------|
| `Beh_matrix.mat` | Behavioral timing data from a center-out reaching task. Rows = trials. Columns (from 2nd onward): start of trial, target on, go cue, start of movement, end of movement, reward, target identity. (First two columns are identical.) |
| `rs1050225_MI_clean_LFP.mat` | Low-pass filtered LFP from 96 channels of a Utah Array. The field `Michan2rc` contains channel coordinates. |


## Repository structure

```
├── data/                    # MATLAB scripts for experimental data analysis (Figs 1–2 and Fig 5)
├── simulations/             # Python scripts for Brian2 network simulations
│   ├── results_f1_f2/       # Analysis of spatially-extended network simulations with ramping external input (Figs 1–2)
│   ├── results_f3/          # Analysis of homogeneous network simulations at constant input (Fig 3)
│   └── results_f4_f5/       # Analysis of spatially-extended network simulations at constant input (Figs 4–5)
└── theory/
    ├── figure3/             # Plots results from linear stability analysis compared to simulations  (Fig 3)
    ├── figure4/             # Plots results from linear stability analysis compared to simulations  (Fig 4)
    └── figure5/             # Plots results from linear stability analysis (Fig 5)

├── functions/              # shared helper functions used across the data analysis, simulation, and theory scripts.
                                                                         
```


## 1. `data/` — Experimental data analysis (MATLAB)

### `Wave_detection_Macaque_2025_RUN.m`
Detects traveling waves in macaque LFP recordings and computes kinematic statistics (speed, wavelength, direction). Plots histograms of these statistics (Figs. 1c, 1d), examples of detected waves (Fig 1b), and the distribution of wave directions (Fig 5d).

### `Macaque_lfp_spectral_envelopes_move_aligned.m`
Computes trial-averaged beta-band power envelopes from multi-channel LFP, time-locked to movement onset. Estimates the dominant beta frequency from the average power spectrum, bandpass filters around the peak (±3 Hz), extracts the Hilbert envelope, normalizes per channel, and averages across trials (Fig 2a).

### `Autoencoder_Analysis_mov.m`
Computes beta-band suppression onset timing across all 96 Utah Array channels for an example trial and visualizes the spatial distribution across the array (Figs 2c–d).

### `Plot_LFP_and_ACF.m`
Visualizes LFP amplitude during movement preparation and execution for a selected trial. Plots z-scored LFP traces and autocorrelation functions (ACF) separately for pre- and post-movement periods (Fig 2g).

 

## 2. `simulations/` — Network simulations (Python / Brian2)

The two scripts simulate the network dynamics using Brian2.

### `Simulations_sigmoidal_external_rate.py`
Spatially extended network - uses a ramping external input field. Inter-site distance: `grid_dist_r = 400 µm`. Results are written to `results_f1_f2/` and reproduce figures 1–2.
To run:
./run_sigmoidal_external_rate.sh
The shell script allows you to specify the number of simulations to run.

### `Simulations_const_ext_rate.py`
Uses constant external input. Two modes are available:
- **Homogeneous connectivity** → results in `results_f3/`, reproduces Fig 3.
- **Spatially extended network** (inter-site distance `grid_dist_r = 141 µm`) → results in `results_f4_f5/`, reproduces Figs 4–5.
To run:
./run_const_external_rate.sh
In the shell script, select options for: homogeneous vs. spatially extended connectivity, value of the external field, simulation duration, and anisotropy of EE connectivity (`anisotropy_EE`).
> Note: For reliable wave direction histograms, we recommend running at least 50 simulations.

 
## 3. Results folders — MATLAB post-processing scripts

### `results_f1_f2/` (ramping external input)

| Script | Description | Figure |
|--------|-------------|--------|
| `Read_multiple_runs_step.m` | Loads and preprocesses simulated LFP. Plots LFP amplitude and compares auto-covariance functions during preparation vs. execution phases. | Fig 2h |
| `Beta_attenuation_map.m` | Generates spatial maps of beta attenuation latencies via Hilbert envelopes. Identifies dominant beta peaks, filters LFP, and computes threshold-crossing latencies. | Figs 2e–f |
| `Plot_waves.m` | Detects and characterizes traveling waves in simulated LFP at baseline external rate. Plots wave examples and velocity distributions compared to data. | Figs 1b, 1d |
| `Power_Spectral_density.m` | Computes and plots z-scored PSD from simulations at baseline, compared to data during the preparation period. | Fig 1c |
| `Beta_power_across_trials.m` | Computes and plots the beta power envelope across trials, averaging normalized envelopes across channels and trials. | Fig 2b |
| `Plot_external_rate.m` | Plots the rate of the external input population. | Fig 2b |
                                                                            
### `results_f3/` (homogeneous network)

| Script | Description | Figure |
|--------|-------------|--------|
| `Read_all.m` | Loads simulation output for a single trial. Plots: z-scored LFP traces; ACF with estimated oscillation frequency; raster plot for E (800) and I (200) neurons. | Fig 3c |

### `results_f4_f5/` (spatially extended network, constant input)

| Script | Description | Figure |
|--------|-------------|--------|
| `Snapshots.m` | Shows spike activity on the 2D array. | Figs 4d–e |
| `Plot_waves.m` | Detects and characterizes planar traveling waves in simulated LFP. Plots wave examples and wave direction distributions. | Figs 4d–e, 5e–f |

---

## 4. `theory/` — Mean-field and stability analysis (MATLAB)

### `figure3/`

| Script | Description | Figure |
|--------|-------------|--------|
| `Plot_phase_diagram.m` | Plots the phase diagram. | Fig 3a |
| `Read_rates.m` | Computes firing rates from the mean-field model at g = 3.625 for varying external field values, and compares to simulations. | Fig 3b (upper) |
| `Plot_lambda_omega.m` | Plots the imaginary part of the eigenvalue at g = 3.625 for varying external field values, compared to oscillation frequencies from simulations. | Fig 3b (lower) |

### `figure4/`

| Script | Description | Figure |
|--------|-------------|--------|
| `W_print_lambdas.m` | Stability analysis: smoothed eigenvalue curves vs. external drive. | Fig 4a |
| `Plot_omega_and_lambda.m` | Plots spatial frequency (k) and oscillation frequency (ω) vs. external drive, comparing simulations to mean-field predictions. | Fig 4b |
| `W_print_lambdas_vs_k.m` | Visualizes the maximum real part of the eigenvalue spectrum as a function of spatial wavenumber k for different propagation velocities. Shows dispersion-like curves Re(λ)(k) and selected instability peaks k*(v). | Fig 4c |

### `figure5/`

Plots Re(λ(kx, ky)) as a 2D heatmap in the space of rotated wave numbers u (direction of maximal connectivity spread) and v (direction of minimal spread).

| Script | Description | Figure |
|--------|-------------|--------|
| `Isotropic_data/W_print_peak_isotropy.m` | Isotropic case. | Fig 5b |
| `Anisotropic_data/W_print_peak_anisotropy.m` | Anisotropic case. | Fig 5c |

## Data availability

All experimental data are available on Zenodo: [10.5281/zenodo.19420479]
                                                                         

## Citation

If you use this code, please cite the paper: Bachschmid-Romano L, Hatsopoulos N, Brunel N. A Spatially Structured Spiking Network Model of Beta Traveling Waves and Their Attenuation in Motor Cortex. bioRxiv. 2026:2026-03. (full citation to be added upon publication).
                          
If you use this data, please cite the paper:  Rubino D, Robbins KA, Hatsopoulos NG. Propagating waves mediate information transfer in the motor cortex. Nature neuroscience. 2006 Dec;9(12):1549-57.
