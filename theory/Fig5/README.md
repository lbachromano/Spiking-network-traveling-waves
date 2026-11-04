# Spatial linear stability of an excitatory–inhibitory network

Use this folder to find oscillatory modes of the linearized network at different spatial wavevectors. For each sampled wavevector `(u, v)`, the code searches for complex roots of the characteristic determinant and saves the root with the largest real part among those found. Its real part gives the growth rate in s⁻¹. A positive value means the mode grows. Its imaginary part divided by 2π gives the oscillation frequency in Hz. The calculation requires no external data files.

## Choose an anisotropy condition

Each kernel has its own anisotropy parameter. Supply `--rho-ee` for EE. The optional arguments `--rho-ii`, `--rho-ie`, and `--rho-ei` control II, IE, and EI respectively and each defaults to zero. All four values must lie strictly between −1 and 1. Zero gives an isotropic kernel. For example:

| Condition | Command-line arguments |
|---|---|
| All kernels isotropic | `--rho-ee 0.0` |
| EE anisotropy only | `--rho-ee 0.6` |
| EI and II anisotropy only | `--rho-ee 0.0 --rho-ei 0.6 --rho-ii 0.6` |

## First run

Put all seven Python files in the same directory and use the filenames in the table below. If your copies have a parenthesized number or timestamp before `.py`, remove that suffix so the imports resolve. From that directory, run the isotropic example below.

```bash
python -m pip install numpy scipy
python run_spatial_stability.py --help
python run_spatial_stability.py --rho-ee 0.0 --rho-ei 0.0 --rho-ii 0.0 --nu-ext 2000 --g 3.3 \
    --velocity 300 --n-u 2 --n-v 2 --n-initial-grid 4 --n-audit-grid 4 \
    --output-dir Results_test
```

This small run checks the workflow and writes one MATLAB-compatible file per completed `u` row. Its coarse search can miss roots. For the default resolution, omit the four grid-size options and use a separate output directory, such as `Results`. Runs with the same model parameters can overwrite row files because the filenames omit the grid settings and round `u` to two decimal places.

| Option | Meaning |
|---|---|
| `--nu-ext` | Absolute external input rate in Hz. Convert any threshold-normalized drive to Hz before calling the script. |
| `--g` | Inhibitory coupling parameter. The wrapper sets `gi = g` and `ge = 1.07*g`. |
| `--velocity` | Propagation speed in mm/s. |
| `--n-u`, `--n-v` | Number of wavevector samples on each axis. Both default to 21, from zero through `--u-max` or `--v-max` (default 1.5 mm⁻¹). |
| `--n-initial-grid`, `--n-audit-grid` | Samples per axis of the complex eigenvalue search. Defaults are 50 for the initial search and 18 for later global checks. |

## Sweeps over coupling or external drive

The two sweep scripts call `Wrap_peak_anisotropy_EE.py` repeatedly. For each tested parameter value, they identify the **largest real part of an accepted characteristic root across the sampled `(u, v)` grid**, together with the wavevector and frequency at which that maximum was found. This provides a one-dimensional summary of the dominant growth rate versus `g` or versus the absolute external rate `nu_ext` (in Hz). Both sweeps set the II, IE, and EI anisotropies to zero; `rho_ee` controls EE anisotropy.

### Sweep the inhibitory coupling `g`: `sweep_g.py`

`sweep_g.py` defines `sweep_g(...)` but has no command-line entry point. From the directory containing the scripts, run:

```bash
python - <<'PY'
import numpy as np
from sweep_g import sweep_g

sweep_g(
    np.linspace(3.3, 3.6, 4),
    rho_ee=0.6,
    nu_ext=2000,
    velocity=300,
    n_u=50, n_v=1, u_max=1.7, v_max=0.0,
    progress=False,
    output_file="Results_5d_summary/max_real_vs_g.npz",
)
PY
```

This varies `g` while holding `nu_ext`, velocity, and anisotropy fixed. The `n_u`, `n_v`, `u_max`, and `v_max` arguments control the spatial grid passed to the wrapper. Without these overrides, the wrapper uses a 21 × 21 grid spanning 0–1.5 mm⁻¹ on each axis.

### Sweep the external rate `nu_ext`: `sweep_nu.py`

Run directly from the terminal:

```bash
python sweep_nu.py --nu-min 2000 --nu-max 2050 --n-nu 4 \
    --g 3.6 --rho-ee 0.6 --velocity 300 \
    --n-u 50 --n-v 1 --u-max 1.7 --v-max 0 \
    --output Results_5d_summary/max_real_vs_nu_ext_g3.6_rho0.6.npz
```

The command samples four evenly spaced **absolute** external input rates between 2000 and 2050 Hz, including the endpoints. By default this script scans a one-dimensional wavevector slice with `n_u=50`, `n_v=1`, `u_max=1.7`, and `v_max=0`; these options can be changed on the command line. Use `python sweep_nu.py --help` to see its available options.

### Where the outputs go

The commands above explicitly place the **summary** `.npz` files in `Results_5d_summary/`, regardless of the default output paths in the local copies of the scripts. This folder is created automatically. Both scripts save a checkpoint after **each completed parameter value**, updating the same summary file. An interrupted run therefore leaves a partial summary, but the scripts do not automatically resume the remaining points.

The wrapper independently writes **detailed** MATLAB-compatible `.mat` files (one per `u` value for each parameter setting) under its default `Results_5d_test/` directory, in subfolders identified by the model parameters. Passing `--output` to `sweep_nu.py` or `output_file` to `sweep_g(...)` changes only the summary destination, not these `.mat` files. To make `Results_5d_summary/` the *default* without specifying paths, update the output defaults in both sweep scripts, including the command-line default constructed in `sweep_nu.py`.

Each summary contains the swept parameter array (`g` or `nu_ext`), `max_real`, `u_star`, `v_star`, `frequency_hz`, `n_success`, `n_total`, and `completed`, along with fixed model parameters and the sampled `List_u` and `List_v`. `max_real` is in s⁻¹, and `frequency_hz` is in Hz. `n_success` counts sampled wavevector points where a finite accepted root was found; `n_total` counts all sampled wavevector points. Entries that have not been computed, or for which no accepted root was found, have `max_real = NaN`; use `completed` to distinguish the former from the latter. As with a single spatial scan, **failure to find a root is not proof of stability**, and results depend on the search resolution and sampled wavevector domain.

To inspect a summary:

```python
import numpy as np

with np.load("Results_5d_summary/max_real_vs_g.npz") as summary:
    for g, growth, done in zip(summary["g"], summary["max_real"], summary["completed"]):
        print(f"g={g:.3f}: max Re(lambda)={growth:.6g}; completed={done}")
```

## Files and calculation

| File | Role |
|---|---|
| `run_spatial_stability.py` | Terminal entry point. Parses arguments and starts one scan. |
| `sweep_g.py` | Runs the stability calculation at a sequence of inhibitory coupling values `g`, and saves the maximum growth rate over the sampled wavevector grid for each value. Called as a Python function (no command-line parser). |
| `sweep_nu.py` | Runs the calculation at a sequence of absolute external rates `nu_ext`, and saves the maximum growth rate over the sampled wavevector grid for each rate. Can be run from the terminal. |
| `Wrap_peak_anisotropy_EE_grid_dominant.py` | Sets the run condition, obtains stationary rates, searches for roots, and writes results. Also callable from Python or Jupyter. |
| `Set_parameters_spatial.py` | Supplies the spatial model parameters used by the wrapper. Lengths are in mm and times in seconds. |
| `stability.py` | Supplies the stationary-rate solver and response functions used in the scan. |
| `Function_determinant_abs_nov_EE.py` | Evaluates the determinant, including response derivatives, synaptic filters, and distance-dependent propagation delays. The finite-square spatial integrals use cached Gauss–Legendre quadrature (64 nodes per axis). |

Each `u` row starts with a search at `v = 0`. At subsequent `v` values the wrapper continues roots from the preceding point, with global checks at the middle and last samples and whenever continuation finds no roots. It retains up to three roots. The search covers growth rates from −50 to 80 s⁻¹ and frequencies from 1 to 50 Hz. A missing root means the search did not accept a solution there. It does not establish stability.

## Read a result

For the example above, the first row is stored under `Results_test/rhoEE0.0_rhoII0.0_rhoIE0.0_rhoEI0.0_g3.30_nuext2000.00_v300/`. The filename ends with `kx0.00.mat`. To read that row from the directory where you ran the command:

```python
from pathlib import Path
import numpy as np
from scipy.io import loadmat

path = next(Path("Results_test").glob("*/Stability_*_kx0.00.mat"))
row = loadmat(path, squeeze_me=True)
v = row["List_v"]
found = np.asarray(row["dominant_success"], dtype=bool)
growth, frequency, log_residual = np.vstack(row["local_minima_all"]).T
print(np.column_stack((v[found], growth[found], frequency[found])))
```

Each entry of `local_minima_all` corresponds to one `v` and contains the dominant growth rate, frequency, and search cost, in that order. The cost is `log(abs(det) + 1e-7)`. If no root was accepted, growth and frequency are NaN. Use `all_roots` to inspect every retained root, `number_roots_found` for their counts, and `audit_was_run` to identify additional global searches. Each file also stores both wavevector axes and the root-search configuration.
