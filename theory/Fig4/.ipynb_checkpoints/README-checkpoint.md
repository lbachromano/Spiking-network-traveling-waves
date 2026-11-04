# Spatial stability analysis of a spiking E–I network

Codes spatial_stability.py and peak_isotropy_sweep.py implement the linear stability analysis of a spatially structured network of excitatory (E) and inhibitory (I) leaky integrate-and-fire neurons. The analysis computes the homogeneous mean-field operating point and numerically searches for minima of the absolute value of the characteristic determinant as a function of the complex growth rate, external drive, and spatial wavenumber.

The calculation incorporates spatially dependent connectivity, propagation delays, synaptic filtering, and an effective frequency-dependent correction to the neuronal response functions.  

## Files

- **`spatial_stability.py`** — Main analysis module and executable. Contains the model parameters, mean-field rate solver, response functions, characteristic determinant, and minimization routines. Its command-line sweep varies external drive for each spatial wavenumber.
- **`peak_isotropy_sweep.py`** — Uses functions from `spatial_stability_last.py` to follow determinant minima while varying the y-component of the spatial wavenumber at a fixed external drive. It can repeat the calculation for several propagation velocities.

Both files must be in the same directory because the second script imports the first.

## Requirements

- Python 3.9 or later
- NumPy
- SciPy

Install the required packages with:

```bash
python -m pip install numpy scipy
```

No external data files are required for the calculations.

## Running the analyses

### 1. External-drive and wavenumber sweep

```bash
python spatial_stability.py --output-dir drive_sweep_results
```

By default, this script uses `g = 3.3`, evaluates 51 values of `k` from 0 to 1.5, and sweeps 100 values of the normalized external drive from 0.7 to 1.2. It sets `kx = ky = k`, so the squared spatial wavenumber entering the determinant is `q² = 2k²`. The propagation velocity is set to 300 mm/s in `Set_parameters()`.
 

Some useful options are `--g`, `--nk`, `--k-min`, `--k-max`, `--ninput`, `--input-min`, `--input-max`, and `--output-dir`. See `python spatial_stability_last.py --help` for all options.

The script writes one `.mat` file and one `.npz` file for each sampled `k` value. Both contain the accepted determinant minima and the associated mean-field rates; the output also records the sweep parameters and other metadata.

### 2. Wavenumber sweep at fixed external drive

```bash
python "peak_isotropy_sweep.py" --output-dir isotropy_sweep_results
```

By default, this script uses `g = 3.3`, `nu_ext / nt_alpha = 0.85`, a propagation velocity of 600 mm/s, `kx = 0`, `rho = 0`, and 20 values of `ky` from 0 to 1.5. Its effective squared wavenumber is

```text
q² = (1 + rho) kx² + (1 - rho) ky²
```

The script follows numerical minima as `ky` increases, using the minima found at one step as starting points for the next. If no acceptable minimum is found, the next step restarts from the initial search grid.

To run several velocities, for example:

```bash
python "peak_isotropy_sweep.py" --velocities 300 600 --output-dir velocity_sweep_results
```

Use `--nu-ext-ratio` to change the external drive and `--drive-reference nutheta` to normalize it by `nutheta` instead of `nt_alpha`. Other options include `--points`, `--ky-min`, `--ky-max`, `--kx`, and `--rho`. See `python "peak_isotropy_sweep(1).py" --help`.

The script saves one MATLAB `.mat` file per propagation velocity.

## Output format

Each accepted minimum is stored as a row with three columns:

1. `Re(lambda)` — real part of the candidate eigenvalue, in s⁻¹.
2. `Im(lambda) / (2π)` — oscillation frequency, in Hz.
3. `log(|det| + 1e-7)` — minimized objective value.

For `spatial_stability.py`, `local_minima_all` stores a list of minima for each external-drive value (`list_varying_param`); `rate_history` contains the corresponding E and I firing rates, in Hz. The compressed `.npz` file stores the minima in `roots`, with `offsets` marking the beginning and end of the results for each drive value.

For `peak_isotropy_sweep.py`, the `.mat` file contains `List_ky` and `local_minima_all`, a MATLAB cell array with one entry per `ky` value. Empty entries mean that no minimum passed the acceptance criterion at that point.

## Numerical and interpretation notes

- The scripts minimize `log(|det| + 1e-7)` over the real and imaginary parts of the complex growth rate using a multistart Nelder–Mead search. The default determinant tolerance is `5e-4` (the first step of the drive sweep uses `1e-3`).
- The mean-field rate solver constrains firing rates to be at least `1e-6` Hz. Very low-drive regimes should therefore be interpreted with care, particularly if the physical network is effectively silent.
- Existing output files are protected against accidental replacement; pass `--overwrite` only when intentional.

Model parameters and response-function approximations are defined in `spatial_stability.py`; consult that file for the implementation used in a particular run.
