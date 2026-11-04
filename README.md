# Beta Traveling Waves in Motor Cortex

Code accompanying *A Spatially Structured Spiking Network Model of Beta Traveling Waves and Their Attenuation in Motor Cortex*.

The `theory/` directory contains the mean-field and linear stability calculations. The `simulations/` directory contains the Brian2 spiking network simulations and scripts for analyzing their output. The tables below identify the files used for each figure.

## Getting started

Install the Python packages used for simulation and analysis, then start JupyterLab from the repository root.

```bash
python -m pip install numpy scipy pandas matplotlib brian2 jupyterlab
python -m jupyterlab
```

Open each notebook in its original folder so its relative paths resolve. To generate theoretical results, follow the README in the corresponding figure folder. To plot results, use the notebooks listed below.

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

## Figure 3

| Panel or analysis | Files |
| --- | --- |
| 3a · Phase diagram | Follow the [stability analysis instructions](theory/Fig3/README.md), then open [Phase_diagram.ipynb](theory/Fig3/Phase_diagram.ipynb). |
| 3b · Dominant frequency | In [Fig3b](simulations/results_f3/Fig3b/), generate simulation data with `Simulate_dynamics.py` and plot them with `Plot_dominant_frequency.ipynb`. |
| 3c · Example simulation | In [Examples_3c_g3.159_nu_1.4](simulations/results_f3/Examples_3c_g3.159_nu_1.4/), use `Simulate_dynamics.py` to generate the example and `Fig_3c.ipynb` to plot it. |
| Finite-size analysis | Under [Results_g3.16](simulations/results_f3/Finite_size/Results_g3.16/), each input condition has a `FiniteSize_analysis_peak_centered.ipynb` notebook. Use `Plot_finite_size_last.ipynb` in their parent folder for the combined plot. |

## Figure 4

See the [Figure 4 theory instructions](theory/Fig4/README.md) for the stability calculations.

| Panel | Files |
| --- | --- |
| 4a | [Fig4_a_b.ipynb](theory/Fig4/Fig4_a_b.ipynb) plots panel 4a and prepares the mean-field prediction used in 4b. |
| 4b | Generate the isotropic simulation data with [Simulate_dynamics.py](simulations/results_f4/Isotropic_g3.3/Simulate_dynamics.py). Open [Fig4b_fixed_r.ipynb](simulations/results_f4/Fig4b_fixed_r.ipynb) to plot the simulation results together with the mean-field prediction. |
| 4c | [Fig4_c_plot.ipynb](theory/Fig4/Fig4_c_plot.ipynb) plots the theoretical results. |
| 4d–e | The simulation is in [Print_waves_spikes](simulations/results_f4/Print_waves_spikes/). The plotting notebook is `LFP_spike_snapshots.ipynb`. |

## Figure 5 and supplementary figures

| Panel | Files |
| --- | --- |
| 5b–c | Follow the [Figure 5 theory instructions](theory/Fig5/README.md) and plot the dominant eigenvalue surfaces with [Dominant_eigenvalue_surfaces.ipynb](theory/Fig5/Dominant_eigenvalue_surfaces.ipynb). |
| 5d | The parameter sweeps are in `theory/Fig5/sweep_g.py` and `theory/Fig5/sweep_nu.py`. Use [plot_5d.ipynb](theory/Fig5/plot_5d.ipynb) for the panel. |
| 5f | [Simulate_dynamics_isotropic.py](simulations/results_f5/Simulate_dynamics_isotropic.py) generates the isotropic network simulations. |
| 5g–h and S9 | [Simulate_dynamics.py](simulations/results_f5/Simulate_dynamics.py) generates the simulations for 5g–h. [Wave_propagation_direction.ipynb](simulations/results_f5/Wave_propagation_direction.ipynb) analyzes propagation direction for 5g–h and S9. |
| S3 | [Supplementary_beta_attenuation_planarity.ipynb](simulations/Supplementary_S3/Supplementary_beta_attenuation_planarity.ipynb) computes beta attenuation latency maps and assesses their planarity. |

## Simulations with time-varying input

From `simulations/`, use `bash run_sigmoidal_external_rate.sh` to launch `Simulations_sigmoidal_external_rate.py`. The shell script sets the run indices, seed, external input and inhibitory coupling.

