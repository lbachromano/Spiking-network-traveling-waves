from pathlib import Path
import numpy as np

from Wrap_peak_anisotropy_EE import (
    Wrap_peak_anisotropy_EE,
)


def sweep_g(g_values, *, rho_ee, nu_ext, velocity,
            output_file="Results_5d_summary/max_real_vs_g.npz", **grid_options):

    g_values = np.asarray(g_values, dtype=float)
    if g_values.ndim != 1 or not g_values.size:
        raise ValueError("g_values must be a nonempty 1D sequence.")

    n = len(g_values)
    max_real = np.full(n, np.nan)
    u_star = np.full(n, np.nan)
    v_star = np.full(n, np.nan)
    frequency_hz = np.full(n, np.nan)
    n_success = np.zeros(n, dtype=int)
    n_total = np.zeros(n, dtype=int)
    completed = np.zeros(n, dtype=bool)

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    for ig, g in enumerate(g_values):
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
            n_total[ig] += len(v_values)

            for iv in np.flatnonzero(success):
                # Stored columns: real part, frequency in Hz, cost.
                real, frequency, _ = row["local_minima_all"][0, iv][0]
                if not np.isfinite(real):
                    continue

                n_success[ig] += 1
                if np.isnan(max_real[ig]) or real > max_real[ig]:
                    max_real[ig] = real
                    u_star[ig] = u
                    v_star[ig] = v_values[iv]
                    frequency_hz[ig] = frequency

        completed[ig] = True
        np.savez(
            output_file,
            g=g_values,
            max_real=max_real,
            u_star=u_star,
            v_star=v_star,
            frequency_hz=frequency_hz,
            n_success=n_success,
            n_total=n_total,
            completed=completed,
            rho_ee=rho_ee,
            nu_ext=nu_ext,
            velocity=velocity,
            List_u=rows[0]["List_u"].ravel(),
            List_v=rows[0]["List_v"].ravel(),
        )

        print(
            f"g={g:g}: max Re(lambda)={max_real[ig]:.6g}; "
            f"roots found at {n_success[ig]}/{n_total[ig]} points",
            flush=True,
        )

    return max_real
