"""Command-line entry point for the stability forward sweep."""

import argparse
import numpy as np

from stability import Wrap_minimization_forw_dec_OK


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)

    # Number of external-input ratios in each forward sweep
    parser.add_argument("--points", type=int, default=100)

    # By default, run the requested 100-value g sweep.
    # Providing --g overrides this, e.g. --g 3.625
    parser.add_argument("--g", type=float, nargs="+", default=None)
    parser.add_argument("--g-start", type=float, default=2.8)
    parser.add_argument("--g-stop", type=float, default=3.5)
    parser.add_argument("--g-points", type=int, default=21)

    parser.add_argument("--output-dir", default="stability_results")
    args = parser.parse_args()

    allg = (
        args.g
        if args.g is not None
        else np.linspace(args.g_start, args.g_stop, args.g_points)
    )

    Wrap_minimization_forw_dec_OK(
        allg=allg,
        N_pt=args.points,
        output_dir=args.output_dir,
    )
