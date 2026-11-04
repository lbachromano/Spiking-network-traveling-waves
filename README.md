# Beta Traveling Waves in Motor Cortex

Code accompanying *A Spatially Structured Spiking Network Model of Beta Traveling Waves and Their Attenuation in Motor Cortex*.

The `theory/` directory contains the mean-field and linear stability calculations. The `simulations/` directory contains the Brian2 spiking network simulations and scripts for analyzing their output. 

## Getting started

Install the Python packages used for simulation and analysis, then start JupyterLab from the repository root.

```bash
python -m pip install numpy scipy pandas matplotlib brian2 jupyterlab
python -m jupyterlab
```

Open each notebook in its original folder so its relative paths resolve. Precomputed theoretical results are included. To regenerate them, follow the instructions in the README for [Figure 3](theory/Fig3/README.md), [Figure 4](theory/Fig4/README.md), or [Figure 5](theory/Fig5/README.md). To plot results, use the notebooks listed below.

## Figures 1–2

The MATLAB post-processing scripts are in [simulations/results_f1_f2/](simulations/results_f1_f2/).

### `results_f1_f2/` — Ramping external input

| Script | Description | Figure |
|--------|-------------|--------|
| `Read_multiple_runs_step.m` | Loads and preprocesses simulated LFP. Plots LFP amplitude and compares autocovariance functions during preparation and execution. | Fig. 2h |
| `Beta_attenuation_map.m` | Generates spatial maps of beta attenuation latencies using Hilbert envelopes. Identifies dominant beta peaks, filters LFP, and computes threshold-crossing latencies. | Figs. 2e–f |
| `Plot_waves.m` | Detects and characterizes traveling waves in simulated LFP at baseline external input. Plots wave examples and velocity distributions compared with data. | Figs. 1b, 1d |
| `Power_Spectral_density.m` | Computes and plots z-scored power spectra from simulations at baseline, compared with experimental data during preparation. | Fig. 1c |
| `Beta_power_across_trials.m` | Computes beta-band power envelopes across trials, averaging normalized envelopes across channels and trials. | Fig. 2b |
| `Plot_external_rate.m` | Plots the rate of the external input population. | Fig. 2b |

Precomputed theoretical results are included. To regenerate them, follow the instructions for [Figure 3](theory/Fig3/README.md), [Figure 4](theory/Fig4/README.md), or [Figure 5](theory/Fig5/README.md).

## Figure 3

| Panel or analysis | Files |
| --- | --- |
| 3a · Phase diagram | Plot the phase diagram with [Phase_diagram.ipynb](theory/Fig3/Phase_diagram.ipynb). |
| 3b · Dominant frequency | Plot the precomputed simulation results with [Plot_dominant_frequency.ipynb](simulations/results_f3/Fig3b/Plot_dominant_frequency.ipynb). To rerun the simulations, use [Simulate_dynamics.py](simulations/results_f3/Fig3b/Simulate_dynamics.py). |
| 3c · Example simulation | Plot the precomputed example with [Fig_3c.ipynb](simulations/results_f3/Examples_3c_g3.159_nu_1.4/Fig_3c.ipynb). To rerun the simulation, use [Simulate_dynamics.py](simulations/results_f3/Examples_3c_g3.159_nu_1.4/Simulate_dynamics.py). |
| Finite-size analysis | In [Results_g3.16](simulations/results_f3/Finite_size/Results_g3.16/), each input condition has its own subfolder containing `FiniteSize_analysis_peak_centered.ipynb`. Plot the results from both conditions with [Plot_finite_size_last.ipynb](simulations/results_f3/Finite_size/Results_g3.16/Plot_finite_size_last.ipynb). |

## Figure 4

| Panel | Files |
| --- | --- |
| 4a | [Fig4_a_b.ipynb](theory/Fig4/Fig4_a_b.ipynb) plots panel 4a and prepares the mean-field prediction used in panel 4b. |
| 4b | Use [Fig4b_fixed_r.ipynb](simulations/results_f4/Fig4b_fixed_r.ipynb) to plot the precomputed simulation results together with the mean-field prediction. To rerun the simulations, use [Simulate_dynamics.py](simulations/results_f4/Isotropic_g3.3/Simulate_dynamics.py). |
| 4c | Plot the theoretical results with [Fig4_c_plot.ipynb](theory/Fig4/Fig4_c_plot.ipynb). |
| 4d–e | In [Print_waves_spikes](simulations/results_f4/Print_waves_spikes/), use `LFP_spike_snapshots.ipynb` to plot the precomputed LFP and spike snapshots. To rerun the simulations, use `Simulations.py`.|

## Figure 5 and supplementary figures

| Panel | Files |
| --- | --- |
| 5b–c | Plot the dominant eigenvalue surfaces with [Dominant_eigenvalue_surfaces.ipynb](theory/Fig5/Dominant_eigenvalue_surfaces.ipynb). |
| 5d | Plot the parameter sweep results with [plot_5d.ipynb](theory/Fig5/plot_5d.ipynb). To rerun the sweeps, use [sweep_g.py](theory/Fig5/sweep_g.py) and [sweep_nu.py](theory/Fig5/sweep_nu.py). |
| 5f | Run the isotropic network simulations with [Simulate_dynamics_isotropic.py](simulations/results_f5/Simulate_dynamics_isotropic.py). |
| 5g–h | Analyze wave propagation direction in the precomputed simulations with [Wave_propagation_direction.ipynb](simulations/results_f5/Wave_propagation_direction.ipynb). To rerun the simulations, use [Simulate_dynamics.py](simulations/results_f5/Simulate_dynamics.py). |
| S3 | Generate the simulation data by following the instructions in the [Supplementary_S3 folder](simulations/Supplementary_S3/). Then run [Supplementary_beta_attenuation_planarity.ipynb](simulations/Supplementary_S3/Supplementary_beta_attenuation_planarity.ipynb) to compute beta attenuation latency maps and assess their planarity. |
| S9 | Analyze wave propagation direction in the precomputed simulations with [Wave_propagation_direction.ipynb](simulations/results_f5/Wave_propagation_direction.ipynb). To rerun the simulations for a given panel, use `Simulate_dynamics.py` in the corresponding subfolder of [Supplementary_results](simulations/results_f5/Supplementary_results/). Each of the four subfolders corresponds to one panel. |

## Simulations with time-varying input

From `simulations/`, use `bash run_sigmoidal_external_rate.sh` to launch `Simulations_sigmoidal_external_rate.py`. The shell script sets the run indices, seed, external input and inhibitory coupling.

