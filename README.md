# Supplementary Code — Beta Traveling Waves in Motor Cortex

This repository contains the analysis and simulation code accompanying the paper:

> **A Spatially Structured Spiking Network Model of Beta Traveling Waves and Their Attenuation in Motor Cortex**

The repository includes Python notebooks used for experimental traveling-wave analyses during revision of the paper, the original MATLAB analysis scripts, Brian2 network simulations, and MATLAB code for mean-field and linear stability analyses.

## Dependencies

- **Experimental data analysis notebooks:** Python 3, NumPy, SciPy, Matplotlib, and JupyterLab.
- **Network simulations:** Python 3 with Brian2.
- **MATLAB analyses and plotting:** MATLAB, with additional toolboxes as required by individual scripts.

Install the Python notebook dependencies with:

```bash
pip install numpy scipy matplotlib jupyterlab
```

Brian2 is required for the network simulations but not for the experimental data analysis notebooks.

## Experimental data

Download the experimental data separately from [Zenodo](https://doi.org/10.5281/zenodo.19420479).

| File | Description |
|------|-------------|
| `Beh_matrix.mat` | Behavioral timing data from a center-out reaching task. Rows correspond to trials. Columns from the second onward contain: start of trial, target onset, go cue, movement onset, movement end, reward, and target identity. The first two columns are identical. |
| `rs1050225_MI_clean_LFP.mat` | Low-pass filtered LFP recordings from 96 channels of a Utah array. The field `MIchan2rc` contains channel coordinates. |

For the Python notebooks, place both files in `data/Macaque_data/`. This folder is created locally by the user; the input files are not included in the repository.

## Repository structure

| Folder | Contents |
|--------|----------|
| `data/` | Python notebooks used during revision and original MATLAB experimental data analyses |
| `simulations/` | Brian2 network simulation scripts |
| `simulations/results_f1_f2/` | Analysis of spatially extended network simulations with ramping external input (Figs. 1–2) |
| `simulations/results_f3/` | Analysis of homogeneous network simulations with constant input (Fig. 3) |
| `simulations/results_f4_f5/` | Analysis of spatially extended network simulations with constant input (Figs. 4–5) |
| `theory/figure3/` | Mean-field and stability results compared with simulations (Fig. 3) |
| `theory/figure4/` | Linear stability results compared with simulations (Fig. 4) |
| `theory/figure5/` | Linear stability results for isotropic and anisotropic connectivity (Fig. 5) |
| `functions/` | Shared MATLAB helper functions |

## 1. `data/` — Experimental data analysis

### Python notebooks used during revision

The following notebooks analyse macaque LFP recordings, detect traveling-wave episodes, and estimate propagation direction, speed, and wavelength. They implement two different measures of spatial phase-gradient coherence.

| Notebook | Detection method |
|----------|------------------|
| `Macaque_wave_detection_fixed_r.ipynb` | Uses the resultant length of normalized spatial phase gradients, giving each valid grid location equal weight. The default criterion is `r > 0.4` for at least 10 ms. |
| `Macaque_wave_detection_as_Rubino.ipynb` | Uses phase-gradient directionality (PGD): the magnitude of the spatially averaged gradient divided by the spatial average of gradient magnitudes. The default criterion is `PGD > 0.5` for at least 10 ms. |

Detailed descriptions of preprocessing, detection, and wave-property estimation are included in each notebook.

Both notebooks export MATLAB-compatible `.mat` files containing event-level measurements, trial-level diagnostics, and analysis settings. They also generate summary plots.

### Running the Python notebooks

1. Download this repository and the experimental data from [Zenodo](https://doi.org/10.5281/zenodo.19420479).
2. Create a folder named `Macaque_data` inside `data/`, alongside the notebooks.
3. Place `Beh_matrix.mat` and `rs1050225_MI_clean_LFP.mat` inside `data/Macaque_data/`.
4. Install the Python notebook dependencies listed above.
5. From the repository root, launch JupyterLab:

   ```bash
   cd data
   jupyter lab
   ```

6. Open either notebook, review its configuration cell, and run the cells in order.

Input paths, trial counts, analysis windows, frequency-band settings, and detection thresholds can be adjusted in the configuration cell. If the input files are stored elsewhere, update `DATA_DIR`.

With the working directory specified above, exported results are saved in `data/results/`, which is created automatically.

### Original MATLAB experimental analyses

The original MATLAB scripts are retained alongside the Python notebooks for reference and access to the original analysis workflows.

#### `Wave_detection_Macaque_2025_RUN.m`

Detects traveling waves in macaque LFP recordings and computes speed, wavelength, and direction. Plots histograms of wave properties (Figs. 1c, 1d), examples of detected waves (Fig. 1b), and the distribution of wave directions (Fig. 5d).

#### `Macaque_lfp_spectral_envelopes_move_aligned.m`

Computes trial-averaged beta-band power envelopes from multichannel LFP, time-locked to movement onset. Estimates the dominant beta frequency from the average power spectrum, bandpass filters around the peak (±3 Hz), extracts the Hilbert envelope, normalizes per channel, and averages across trials (Fig. 2a).

#### `Autoencoder_Analysis_mov.m`

Computes beta-band suppression onset timing across the Utah array for an example trial and visualizes its spatial distribution (Figs. 2c–d).

#### `Plot_LFP_and_ACF.m`

Visualizes LFP amplitude during movement preparation and execution for a selected trial. Plots z-scored LFP traces and autocorrelation functions separately for pre- and post-movement periods (Fig. 2g).

## 2. `simulations/` — Network simulations (Python / Brian2)

The two simulation scripts generate network activity using Brian2. Run the corresponding shell scripts from the `simulations/` directory.

### `Simulations_sigmoidal_external_rate.py`

Simulates a spatially extended network with ramping external input. The recording-site spacing is `grid_dist_r = 400 µm`. Results are written to `results_f1_f2/` for the analyses associated with Figs. 1–2.

To run:

```bash
./run_sigmoidal_external_rate.sh
```

The shell script allows you to specify the number of simulations.

### `Simulations_const_ext_rate.py`

Simulates a network with constant external input. Two modes are available:

- **Homogeneous connectivity:** results are written to `results_f3/` for Fig. 3.
- **Spatially extended network:** recording-site spacing is `grid_dist_r = 141 µm`; results are written to `results_f4_f5/` for Figs. 4–5.

To run:

```bash
./run_const_external_rate.sh
```

In the shell script, select the connectivity mode, external input strength, simulation duration, and EE connectivity anisotropy (`anisotropy_EE`).

For reliable wave-direction histograms, we recommend running at least 50 simulations.

## 3. Results folders — MATLAB post-processing scripts

The following folders are located inside `simulations/`.

### `results_f1_f2/` — Ramping external input

| Script | Description | Figure |
|--------|-------------|--------|
| `Read_multiple_runs_step.m` | Loads and preprocesses simulated LFP. Plots LFP amplitude and compares autocovariance functions during preparation and execution. | Fig. 2h |
| `Beta_attenuation_map.m` | Generates spatial maps of beta attenuation latencies using Hilbert envelopes. Identifies dominant beta peaks, filters LFP, and computes threshold-crossing latencies. | Figs. 2e–f |
| `Plot_waves.m` | Detects and characterizes traveling waves in simulated LFP at baseline external input. Plots wave examples and velocity distributions compared with data. | Figs. 1b, 1d |
| `Power_Spectral_density.m` | Computes and plots z-scored power spectra from simulations at baseline, compared with experimental data during preparation. | Fig. 1c |
| `Beta_power_across_trials.m` | Computes beta-band power envelopes across trials, averaging normalized envelopes across channels and trials. | Fig. 2b |
| `Plot_external_rate.m` | Plots the rate of the external input population. | Fig. 2b |

### `results_f3/` — Homogeneous network

| Script | Description | Figure |
|--------|-------------|--------|
| `Read_all.m` | Loads simulation output for a single trial. Plots z-scored LFP traces, autocorrelation functions with estimated oscillation frequency, and a raster for 800 excitatory and 200 inhibitory neurons. | Fig. 3c |

### `results_f4_f5/` — Spatially extended network, constant input

| Script | Description | Figure |
|--------|-------------|--------|
| `Snapshots.m` | Shows spike activity on the two-dimensional array. | Figs. 4d–e |
| `Plot_waves.m` | Detects and characterizes planar traveling waves in simulated LFP. Plots wave examples and wave-direction distributions. | Figs. 4d–e, 5e–f |

## 4. `theory/` — Mean-field and stability analysis (MATLAB)

### `figure3/`

| Script | Description | Figure |
|--------|-------------|--------|
| `Plot_phase_diagram.m` | Plots the phase diagram. | Fig. 3a |
| `Read_rates.m` | Computes mean-field firing rates at `g = 3.625` for varying external input strengths and compares them with simulations. | Fig. 3b, upper |
| `Plot_lambda_omega.m` | Plots the imaginary part of the eigenvalue at `g = 3.625` for varying external input strengths and compares it with oscillation frequencies from simulations. | Fig. 3b, lower |

### `figure4/`

| Script | Description | Figure |
|--------|-------------|--------|
| `W_print_lambdas.m` | Plots smoothed eigenvalue curves as a function of external input strength. | Fig. 4a |
| `Plot_omega_and_lambda.m` | Plots spatial frequency (`k`) and oscillation frequency (`ω`) as a function of external input strength, comparing simulations with mean-field predictions. | Fig. 4b |
| `W_print_lambdas_vs_k.m` | Plots the maximal real part of the eigenvalue spectrum as a function of spatial wavenumber for different propagation velocities, together with selected instability peaks. | Fig. 4c |

### `figure5/`

Plots `Re(λ(kx, ky))` as a two-dimensional heatmap in rotated wavenumber coordinates: `u` follows the direction of maximal connectivity spread, and `v` follows the direction of minimal spread.

| Script | Description | Figure |
|--------|-------------|--------|
| `Isotropic_data/W_print_peak_isotropy.m` | Isotropic connectivity. | Fig. 5b |
| `Anisotropic_data/W_print_peak_anisotropy.m` | Anisotropic connectivity. | Fig. 5c |

## Data availability

All experimental data are available on [Zenodo](https://doi.org/10.5281/zenodo.19420479). Download the input MAT files separately and follow the setup instructions above.

## Citation

If you use this code, please cite:

Bachschmid-Romano L, Hatsopoulos N, Brunel N. *A Spatially Structured Spiking Network Model of Beta Traveling Waves and Their Attenuation in Motor Cortex.* bioRxiv. 2026. Full citation to be added upon publication.

If you use the experimental data, please also cite:

Rubino D, Robbins KA, Hatsopoulos NG. Propagating waves mediate information transfer in the motor cortex. *Nature Neuroscience.* 2006;9(12):1549–1557.

