"""Trace stability minima while increasing ``k_y`` in an isotropic network.

Keep this file in the same directory as ``spatial_stability.py``; it reuses that module's
parameter setup, mean-field solver, determinant, minimizer, and tolerance
based de-duplication.

Example
-------
python peak_isotropy_sweep.py --output-dir stability_results

The default external drive uses ``nt_alpha * nu_ext_ratio``, matching the
current convention in ``spatial_stability.py``.  Pass ``--drive-reference
nutheta``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.io import savemat

from spatial_stability_last import Set_parameters, Set_rates_IC, find_minima


def starting_points(args: argparse.Namespace) -> list[tuple[float, float]]:
    """Return (Re(lambda), Im(lambda)) starting points for Nelder--Mead.

    ``find_minima`` expects the imaginary component in rad/s.  The four
    default points reproduce the active starting grid documented in the
    attached ``spatial_stability.py``: Re(lambda) = 0 or 20 and f = 10 or
    20 Hz.
    """
    return [
        (real_part, 2.0 * np.pi * frequency)
        for real_part in np.linspace(args.start_re_min, args.start_re_max,
                                     args.start_re_count)
        for frequency in np.linspace(args.start_freq_min, args.start_freq_max,
                                     args.start_freq_count)
    ]


def set_external_drive(param, ratio: float, reference: str) -> None:
    """Set the external input using either the current or legacy reference."""
    if reference == "nt-alpha":
        param.nu_ext = param.nt_alpha * ratio
    else:
        param.nu_ext = param.nutheta * ratio


def objective_threshold(args: argparse.Namespace) -> float:
    """Return the acceptance threshold for Function_determinant_abs_nov.

    The imported function returns log(abs(det) + 1e-7), not abs(det).  A
    direct log threshold remains available for an exact legacy setting whose
    value is already known.
    """
    if args.log_threshold is not None:
        return args.log_threshold
    return float(np.log(args.det_tol + 1e-7))


def output_path(args: argparse.Namespace, velocity: float) -> Path:
    name = ("Stability_peak_isotropy_2_"
            f"nu_ext{args.nu_ext_ratio:.2f}_g{args.g:.2f}_v{velocity:.1f}.mat")
    return Path(args.output_dir) / name


def run(args: argparse.Namespace) -> None:
    param = Set_parameters()
    param.gi = args.g
    param.ge = 1.07 * param.gi

    list_ky = np.linspace(args.ky_min, args.ky_max, args.points)
    starts = starting_points(args)
    threshold = objective_threshold(args)
    initial_nu = np.asarray(args.initial_rates, dtype=float)

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    print(f"nutheta={param.nutheta:.9g}, nt_alpha={param.nt_alpha:.9g}")
    print(f"accepting log(abs(det) + 1e-7) < {threshold:.9g}")

    for velocity in args.velocities:
        param.velocity = velocity
        # MATLAB cell(1, N_pt_v), including an explicit empty matrix when no
        # acceptable minimum is found at a k_y value.
        local_minima_all = np.empty((1, list_ky.size), dtype=object)
        previous_local_minima = np.empty((0, 3), dtype=float)

        for index, ky in enumerate(list_ky):
            # rho=0 and kx=0 by default, exactly as in the MATLAB wrapper.
            param.qksquared = ((1.0 + args.rho) * args.kx**2
                               + (1.0 - args.rho) * ky**2)

            set_external_drive(param, args.nu_ext_ratio, args.drive_reference)
            rates = Set_rates_IC(param, initial_nu)
            initial_nu = np.array([rates.nu_E, rates.nu_I])

            # Compute_starting_points_isanis is not defined in the supplied
            # Python module.  On the first k_y we use its documented
            # multistart grid; thereafter we continue from the prior roots.
            if previous_local_minima.size:
                current_starts = [
                    (root[0], root[1] * 2.0 * np.pi)
                    for root in previous_local_minima
                ]
            else:
                current_starts = starts

            local_minima = find_minima(param, rates, current_starts, threshold)
            local_minima_all[0, index] = local_minima

            # Unlike the general sweep in spatial_stability.py, the MATLAB
            # code explicitly resets continuation when a step has no roots.
            previous_local_minima = local_minima

            print(
                f"k_y {index + 1}/{list_ky.size} = {ky:.6g}; "
                f"rates = ({rates.nu_E:.6g}, {rates.nu_I:.6g}) Hz; "
                f"accepted minima = {len(local_minima)}",
                flush=True,
            )

        path = output_path(args, velocity)
        if path.exists() and not args.overwrite:
            raise FileExistsError(
                f"{path} already exists. Choose another --output-dir or pass --overwrite."
            )

        # Preserve the two variables and cell-array orientation saved by
        # ``save(myname, 'List_ky', 'local_minima_all')`` in MATLAB.
        savemat(path, {
            "List_ky": list_ky.reshape(1, -1),
            "local_minima_all": local_minima_all,
        })
        print(f"Saved {path}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", default=".")
    parser.add_argument("--g", type=float, default=3.3)
    parser.add_argument("--nu-ext-ratio", type=float, default=0.85)
    parser.add_argument("--drive-reference", choices=("nt-alpha", "nutheta"),
                        default="nt-alpha",
                        help="Use nt-alpha (current convention) or nutheta (legacy MATLAB convention).")
    parser.add_argument("--velocities", type=float, nargs="+", default=[600.0])

    parser.add_argument("--points", type=int, default=20)
    parser.add_argument("--ky-min", type=float, default=0.0)
    parser.add_argument("--ky-max", type=float, default=1.5)
    parser.add_argument("--kx", type=float, default=0.0)
    parser.add_argument("--rho", type=float, default=0.0)

    parser.add_argument("--initial-rates", type=float, nargs=2, metavar=("NU_E", "NU_I"),
                        default=[0.1, 0.2])
    parser.add_argument("--start-re-min", type=float, default=0.0)
    parser.add_argument("--start-re-max", type=float, default=20.0)
    parser.add_argument("--start-re-count", type=int, default=2)
    parser.add_argument("--start-freq-min", type=float, default=10.0)
    parser.add_argument("--start-freq-max", type=float, default=20.0)
    parser.add_argument("--start-freq-count", type=int, default=2)

    parser.add_argument("--det-tol", type=float, default=5e-4,
                        help="Accept |det| below this value (converted to the log objective).")
    parser.add_argument("--log-threshold", type=float, default=None,
                        help="Use this direct threshold on log(abs(det) + 1e-7).")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    values = np.array([
        args.g, args.nu_ext_ratio, args.ky_min, args.ky_max, args.kx, args.rho,
        *args.velocities, *args.initial_rates, args.start_re_min, args.start_re_max,
        args.start_freq_min, args.start_freq_max, args.det_tol,
    ])
    if not np.all(np.isfinite(values)):
        parser.error("All numeric arguments must be finite.")
    if (args.points < 1 or args.start_re_count < 1 or args.start_freq_count < 1
            or args.ky_max < args.ky_min or args.det_tol <= 0
            or args.nu_ext_ratio <= 0 or any(v <= 0 for v in args.velocities)):
        parser.error("Require positive counts, drive ratio, determinant tolerance, and velocities; "
                     "also require --ky-max >= --ky-min.")
    if args.log_threshold is not None and not np.isfinite(args.log_threshold):
        parser.error("--log-threshold must be finite.")
    return args


if __name__ == "__main__":
    run(parse_args())
