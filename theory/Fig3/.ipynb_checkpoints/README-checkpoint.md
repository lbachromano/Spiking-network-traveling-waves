# Linear stability analysis

The Python scripts perform a linear stability analysis of the excitatory–inhibitory network model. For each inhibitory coupling parameter `g`, they sweep the external input, compute the corresponding mean-field firing rates, and search for minima of the absolute value of the characteristic determinant in the complex growth-rate plane.

## Requirements

- Python 3
- NumPy
- SciPy
  
Install the dependencies with:

```bash
python -m pip install numpy scipy
```

## Running the analysis

Keep both scripts in the same directory, named:

- `stability.py` — model parameters, mean-field calculations, determinant, and minimization routines.
- `Wrap_minimization.py` — command-line script for running the parameter sweep.

From that directory, run:

```bash
python Wrap_minimization.py
```

By default, this scans **21 values of `g` from 2.8 to 3.5** and **100 external-input ratios** (`nu_ext / nutheta`) from 0.7 to 1.5. This can be computationally intensive.

For a shorter test run:

```bash
python Wrap_minimization.py --g 3.625 --points 5 --output-dir test_results
```

To change the parameter range, use `--g-start`, `--g-stop`, `--g-points`, or supply specific values with `--g`. Use `--points` for the number of external-input values and `--output-dir` for the results folder.

## Output

The code saves one MATLAB-compatible `.mat` file per `g` value in `stability_results/` (or the specified output directory). Each contains:

- `list_varying_param`: external-input ratios (`nu_ext / nutheta`).
- `local_minima_all`: minima found at each input ratio, with columns **Re(λ) [s⁻¹]**, **Im(λ)/(2π) [Hz]**, and **log(|det| + 10⁻⁷)**.

The numerical search identifies candidate roots through local minimization; results depend on the initial guesses and should not be interpreted as an exhaustive list of eigenvalues.
