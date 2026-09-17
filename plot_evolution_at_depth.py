import argparse

import numpy as np
import pylab as pyl

from utility import first_zero, oscillation_onset, time_first_minimum


FLOAT_FMT = ".3f"
QUANTITIES = ("CA", "CC", "cCa", "cCO3", "Phi")


def _plot_event_times(t_onset, t_first_zero_U, t_first_minimum_U):
    """Add the calculated event times to the current axes."""
    pyl.axvline(
        t_onset,
        color="k",
        linestyle="--",
        label=f"t_onset = {t_onset:{FLOAT_FMT}}",
    )
    pyl.axvline(
        t_first_zero_U,
        color="g",
        linestyle="-.",
        label=f"time of first zero U = {t_first_zero_U:{FLOAT_FMT}}",
    )
    pyl.axvline(
        t_first_minimum_U,
        color="r",
        linestyle=":",
        label=f"time of first minimum U = {t_first_minimum_U:{FLOAT_FMT}}",
    )


def plot_evolution(
    times,
    CA,
    CC,
    cCa,
    cCO3,
    Phi,
    U_values,
    extrema,
    t_onset,
    t_first_zero_U,
    t_first_minimum_U,
    *,
    U_times=None,
    quantity="cCO3",
):
    """Make the three plots used to inspect the evolution at one depth."""
    if U_times is None:
        U_times = times

    pyl.plot(times, CA, "ro")
    pyl.plot(times, CC, "gv")
    pyl.plot(times, cCa, "k>")
    pyl.plot(times, cCO3, "b<")
    pyl.plot(times, Phi, "y+")
    pyl.legend(["CA", "CC", "cCa", "cCO3", "Phi"])
    pyl.show()
    pyl.clf()

    selected_values = dict(zip(QUANTITIES, (CA, CC, cCa, cCO3, Phi)))[quantity]
    pyl.plot(times, selected_values, "rx", label=quantity)
    pyl.plot(
        times[extrema], selected_values[extrema], "bo",
        label=f"selected extrema of {quantity}",
    )
    _plot_event_times(t_onset, t_first_zero_U, t_first_minimum_U)
    pyl.legend()
    pyl.show()
    pyl.clf()

    pyl.plot(U_times, U_values, "rx")
    _plot_event_times(t_onset, t_first_zero_U, t_first_minimum_U)
    pyl.xlabel("time [kyr]")
    pyl.ylabel("U at the bottom of the grid (unitless)")
    pyl.legend()
    pyl.show()


def evolution_at_depth(
    times,
    CA,
    CC,
    cCa,
    cCO3,
    Phi,
    U_values,
    *,
    U_times=None,
    quantity="cCO3",
    make_plots=True,
    print_diagnostics=True,
    allow_missing_diagnostics=False,
):
    """Calculate diagnostics and plot the evolution at a selected depth.

    ``U_times`` can be supplied when U was sampled at different times from the
    other quantities. If omitted, ``times`` is also used for U. Plotting and
    printing can be disabled for batch processing. Set
    ``allow_missing_diagnostics`` to return NaNs for diagnostics that cannot be
    determined because a time series is too short.

    Returns a dictionary containing all calculated diagnostics.
    """
    times, CA, CC, cCa, cCO3, Phi, U_values = map(
        np.asarray, (times, CA, CC, cCa, cCO3, Phi, U_values)
    )
    U_times = times if U_times is None else np.asarray(U_times)
    if quantity not in QUANTITIES:
        raise ValueError(f"quantity must be one of: {', '.join(QUANTITIES)}")
    selected_values = dict(zip(QUANTITIES, (CA, CC, cCa, cCO3, Phi)))[quantity]

    try:
        t_onset, period, extrema = oscillation_onset(
            times, selected_values, prominence=0.01, min_extrema=6
        )
    except ValueError:
        if not allow_missing_diagnostics:
            raise
        t_onset = period = np.nan
        extrema = np.array([], dtype=int)
    try:
        t_first_zero_U = first_zero(U_times, U_values)
    except ValueError:
        if not allow_missing_diagnostics:
            raise
        t_first_zero_U = np.nan
    try:
        t_first_minimum_U = time_first_minimum(U_times, U_values)
    except ValueError:
        if not allow_missing_diagnostics:
            raise
        t_first_minimum_U = np.nan

    extremum_times = np.full(4, np.nan)
    selected_extrema = times[extrema[:4]]
    extremum_times[: len(selected_extrema)] = selected_extrema
    if len(selected_extrema) < 4 and not allow_missing_diagnostics:
        raise ValueError("At least four extrema are required for the diagnostics.")
    first, second, third, fourth = extremum_times

    results = {
        "t_onset": t_onset,
        "period": period,
        "extrema": extrema,
        "t_first_zero_U": t_first_zero_U,
        "t_first_minimum_U": t_first_minimum_U,
        "t_first_extremum": first,
        "t_second_extremum": second,
        "t_third_extremum": third,
        "t_fourth_extremum": fourth,
        "time_between_first_and_third_extremum": third - first,
        "time_between_second_and_fourth_extremum": fourth - second,
        "time_of_first_extremum_minus_quarter_period": first - period / 4.0,
        "time_of_second_extremum_minus_half_period": second - period / 2.0,
    }

    if print_diagnostics:
        _print_diagnostics(results)
    if make_plots:
        plot_evolution(
            times,
            CA,
            CC,
            cCa,
            cCO3,
            Phi,
            U_values,
            extrema,
            t_onset,
            t_first_zero_U,
            t_first_minimum_U,
            U_times=U_times,
            quantity=quantity,
        )
    return results


def _print_diagnostics(results):
    labels = (
        ("first minimum time U", "t_first_minimum_U"),
        ("time of first zero U", "t_first_zero_U"),
        ("t_onset", "t_onset"),
        ("period", "period"),
        ("first extremum", "t_first_extremum"),
        ("second extremum", "t_second_extremum"),
        ("third extremum", "t_third_extremum"),
        ("fourth extremum", "t_fourth_extremum"),
        (
            "time between first and third extremum",
            "time_between_first_and_third_extremum",
        ),
        (
            "time between second and fourth extremum",
            "time_between_second_and_fourth_extremum",
        ),
        (
            "time of first extremum minus period/4",
            "time_of_first_extremum_minus_quarter_period",
        ),
        (
            "time of second extremum minus period/2",
            "time_of_second_extremum_minus_half_period",
        ),
    )
    for label, key in labels:
        print(f"{label} = {results[key]:{FLOAT_FMT}}")


def main():
    """Load a result file and analyse the selected depth."""
    import h5py

    parser = argparse.ArgumentParser(
        description=(
            "Analyse the evolution at a grid node in a marlpde HDF5 result "
            "file. By default, use the deepest node. Print the oscillation "
            "onset, period, and characteristic times for U, then show three "
            "diagnostic plots. U is always measured at the bottom of the grid."
        ),
        epilog=(
            "Example: python plot_evolution_at_depth.py "
            "../Results/L_500/LMAHeureuxPorosityDiff.hdf5"
        ),
    )
    parser.add_argument(
        "file",
        help="path to a marlpde LMAHeureuxPorosityDiff.hdf5 result file",
    )
    parser.add_argument(
        "--depth-index",
        type=int,
        metavar="INDEX",
        help="zero-based grid node index (default: deepest node)",
    )
    parser.add_argument(
        "--quantity",
        choices=QUANTITIES,
        default="cCO3",
        help=(
            "quantity used for oscillation analysis and the second plot "
            "(default: cCO3; unavailable diagnostics are shown as nan)"
        ),
    )
    args = parser.parse_args()

    with h5py.File(args.file, "r") as hf:
        tstar_years = float(hf.attrs["Tstar"])
        times = np.asarray(hf["times"]) * tstar_years / 1_000
        data = np.asarray(hf["data"])
        U_arr = np.asarray(hf["U"]["U_at_bottom"])

    depth_index = data.shape[2] - 1 if args.depth_index is None else args.depth_index
    if not 0 <= depth_index < data.shape[2]:
        parser.error(
            f"--depth-index must be between 0 and {data.shape[2] - 1} "
            f"for this result file"
        )

    evolution_at_depth(
        times,
        data[:, 0, depth_index],
        data[:, 1, depth_index],
        data[:, 2, depth_index],
        data[:, 3, depth_index],
        data[:, 4, depth_index],
        U_arr[:, 1],
        U_times=U_arr[:, 0] * tstar_years / 1_000,
        quantity=args.quantity,
        allow_missing_diagnostics=args.quantity != "cCO3",
    )


if __name__ == "__main__":
    main()
