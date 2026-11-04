from pathlib import Path
import argparse
import numpy as np

from Wrap_peak_anisotropy_EE import (
    Wrap_peak_anisotropy_EE,
)


def sweep_nu_ext(nu_ext_values, *, rho_ee, g, velocity,
                 output_file="Results_5d_summary/max_real_vs_nu_ext.npz",
                 **grid_options):

    # Default to a 50-point slice, including both endpoints.
    grid_options = dict(n_u=50, n_v=1, u_max=1.7, v_max=0.0) | grid_options
    nu_ext_values = np.asarray(nu_ext_values, dtype=float)
    if nu_ext_values.ndim != 1 or not nu_ext_values.size:
        raise ValueError("nu_ext_values must be a nonempty 1D sequence.")

    n = len(nu_ext_values)
    max_real = np.full(n, np.nan)
    u_star = np.full(n, np.nan)
    v_star = np.full(n, np.nan)
    frequency_hz = np.full(n, np.nan)
    n_success = np.zeros(n, dtype=int)
    n_total = np.zeros(n, dtype=int)
    completed = np.zeros(n, dtype=bool)

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    for inu, nu_ext in enumerate(nu_ext_values):
        rows = Wrap_peak_anisotropy_EE(
            rho_EE_input=rho_ee,
            rho_EI_input=0,
            rho_IE_input=0,
            rho_II_input=0,
            nu_from_input=nu_ext,
            g_from_input=g,
            v_fom_input=velocity,
            **grid_options,
        )

        for iu, row in enumerate(rows):
            u = row["List_u"].ravel()[iu]
            v_values = row["List_v"].ravel()
            success = row["dominant_success"].ravel()
            n_total[inu] += len(v_values)

            for iv in np.flatnonzero(success):
                # Columns: real part, frequency in Hz, cost.
                real, frequency, _ = row["local_minima_all"][0, iv][0]
                if not np.isfinite(real):
                    continue

                n_success[inu] += 1
                if np.isnan(max_real[inu]) or real > max_real[inu]:
                    max_real[inu] = real
                    u_star[inu] = u
                    v_star[inu] = v_values[iv]
                    frequency_hz[inu] = frequency

        completed[inu] = True

        # Save a checkpoint after each completed external rate.
        np.savez(
            output_file,
            nu_ext=nu_ext_values,
            max_real=max_real,
            u_star=u_star,
            v_star=v_star,
            frequency_hz=frequency_hz,
            n_success=n_success,
            n_total=n_total,
            completed=completed,
            rho_ee=rho_ee,
            g=g,
            velocity=velocity,
            List_u=rows[0]["List_u"].ravel(),
            List_v=rows[0]["List_v"].ravel(),
        )

        print(
            f"nu_ext={nu_ext:g} Hz: "
            f"max Re(lambda)={max_real[inu]:.6g}; "
            f"roots found at {n_success[inu]}/{n_total[inu]} points",
            flush=True,
        )

    return max_real


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--nu-min", type=float, required=True)
    parser.add_argument("--nu-max", type=float, required=True)
    parser.add_argument("--n-nu", type=int, required=True)
    parser.add_argument("--g", type=float, required=True)
    parser.add_argument("--rho-ee", type=float, required=True)
    parser.add_argument("--velocity", type=float, required=True)
    parser.add_argument("--n-u", type=int, default=50)
    parser.add_argument("--n-v", type=int, default=1)
    parser.add_argument("--u-max", type=float, default=1.7)
    parser.add_argument("--v-max", type=float, default=0.0)
    
    # 1. Set default to None so we can handle it dynamically
    parser.add_argument("--output", default=None, help="Path to output npz file")
    
    args = parser.parse_args()

    # 2. Generate the dynamic default if the user didn't provide --output
    if args.output is None:
        args.output = f"Results_5d_summary/max_real_vs_nu_ext_g{args.g:.1f}_rho{args.rho_ee:.1f}.npz"

    if args.n_nu < 1 or args.nu_max < args.nu_min:
        parser.error("Require n-nu >= 1 and nu-max >= nu-min.")
    if args.nu_min < 0:
        parser.error("External rates must be nonnegative.")
    if min(args.n_u, args.n_v) < 1:
        parser.error("n-u and n-v must be at least 1.")
    if args.velocity <= 0:
        parser.error("velocity must be positive.")

    sweep_nu_ext(
        np.linspace(args.nu_min, args.nu_max, args.n_nu),
        rho_ee=args.rho_ee,
        g=args.g,
        velocity=args.velocity,
        output_file=args.output,
        n_u=args.n_u,
        n_v=args.n_v,
        u_max=args.u_max,
        v_max=args.v_max,
        progress=False,
    )
