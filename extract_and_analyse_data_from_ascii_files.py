import argparse

import numpy as np

from plot_evolution_at_depth import evolution_at_depth


DEFAULT_TSTAR_YEARS = 13190.0


def extract_and_analyse(
    ascii_file,
    *,
    make_plots=True,
    print_diagnostics=True,
    allow_missing_diagnostics=False,
    tstar_years=DEFAULT_TSTAR_YEARS,
):
    """Load evolution data from an ASCII file and return its diagnostics."""
    data = np.loadtxt(ascii_file)

    return evolution_at_depth(
        data[:, 0] * tstar_years / 1_000,
        data[:, 1],
        data[:, 2],
        data[:, 3],
        data[:, 4],
        data[:, 5],
        data[:, 6],
        make_plots=make_plots,
        print_diagnostics=print_diagnostics,
        allow_missing_diagnostics=allow_missing_diagnostics,
    )


def main():
    """Load evolution data from an ASCII file and analyse it."""
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("ascii_file", help="Rhythmite solution_x_*.ascii file")
    parser.add_argument(
        "--tstar-years", type=float, default=DEFAULT_TSTAR_YEARS,
        help=f"Rhythmite time scale in years (default: {DEFAULT_TSTAR_YEARS:g})",
    )
    args = parser.parse_args()
    extract_and_analyse(args.ascii_file, tstar_years=args.tstar_years)


if __name__ == "__main__":
    main()
