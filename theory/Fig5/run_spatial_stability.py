
from __future__ import annotations

import argparse
from pathlib import Path

from Wrap_peak_anisotropy_EE import (
    Wrap_peak_anisotropy_EE,
)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Find dominant roots of the spatial EE-anisotropy determinant.")
    parser.add_argument("--rho-ee", type=float, required=True,
                        help="EE anisotropy rho_EE.")
    parser.add_argument("--rho-ii", type=float, default=0.0,
                        help="Anisotropy of the II kernels (default: 0.0).")
    parser.add_argument("--rho-ie", type=float, default=0.0,
                        help="Anisotropy of the IE kernels (default: 0.0).")
    parser.add_argument("--rho-ei", type=float, default=0.0,
                        help="Anisotropy of the EI kernels (default: 0.0).")
    parser.add_argument("--nu-ext", type=float, required=True,
                        help="Absolute external input rate (Hz).")
    parser.add_argument("--g", type=float, required=True,
                        help="Inhibitory coupling parameter g.")
    parser.add_argument("--velocity", type=float, required=True,
                        help="Propagation velocity in mm/s (e.g. 300).")
    parser.add_argument("--n-u", type=int, default=21,
                        help="Number of u wavevector values (default: 21).")
    parser.add_argument("--n-v", type=int, default=21,
                        help="Number of v wavevector values (default: 21).")
    parser.add_argument("--u-max", type=float, default=1.5,
                        help="Largest u value in mm^-1 (default: 1.5).")
    parser.add_argument("--v-max", type=float, default=1.5,
                        help="Largest v value in mm^-1 (default: 1.5).")
    parser.add_argument("--n-initial-grid", type=int, default=50,
                        help="Initial eigenvalue-search grid size per dimension (default: 50).")
    parser.add_argument("--n-audit-grid", type=int, default=18,
                        help="Audit eigenvalue-search grid size per dimension (default: 18).")
    parser.add_argument("--output-dir", type=Path, default=Path("Results"),
                        help="Directory for the output subfolder (default: Results).")
    parser.add_argument("--quiet", action="store_true",
                        help="Suppress per-grid-point progress messages.")
    return parser


def main():
    args = build_parser().parse_args()
    if min(args.n_u, args.n_v, args.n_initial_grid, args.n_audit_grid) < 2:
        raise SystemExit("--n-u, --n-v, --n-initial-grid and --n-audit-grid must be at least 2.")
    if args.velocity <= 0:
        raise SystemExit("--velocity must be positive.")

    results = Wrap_peak_anisotropy_EE(
        rho_EE_input=args.rho_ee,
        rho_II_input=args.rho_ii,
        rho_IE_input=args.rho_ie,
        rho_EI_input=args.rho_ei,
        nu_from_input=args.nu_ext,
        g_from_input=args.g,
        v_fom_input=args.velocity,
        n_u=args.n_u,
        n_v=args.n_v,
        u_max=args.u_max,
        v_max=args.v_max,
        n_initial_grid=args.n_initial_grid,
        n_audit_grid=args.n_audit_grid,
        output_dir=args.output_dir,
        progress=not args.quiet,
    )
    print(f"Finished. Saved {len(results)} u rows under {args.output_dir.resolve()}.")


if __name__ == "__main__":
    main()
