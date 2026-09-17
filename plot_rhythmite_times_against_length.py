import argparse
from pathlib import Path

import h5py
import numpy as np
import pylab as pyl
from scipy.optimize import curve_fit

from extract_and_analyse_data_from_ascii_files import (
    DEFAULT_TSTAR_YEARS,
    extract_and_analyse,
)
from utility import oscillation_onset


# L=1625 is deliberately excluded: its apparent oscillations are numerical
# artefacts rather than physical model behaviour.
MODEL_LENGTHS = np.array([500, 625, 750, 875, 1000, 1125, 1250])
PREDICTION_LENGTHS = np.array([1125, 1250])
MARLPDE_LENGTHS = np.array([500, 625, 750, 875, 1000, 1125])


def power_law(length, time_at_1000, exponent):
    """Use a scaled length for numerically stable fitting; time is in kyr."""
    return time_at_1000 * (np.asarray(length) / 1000.0)**exponent


def ascii_file_for_length(data_directory, reference_500_file, length):
    """Return the rhythmite output file at the bottom of a model."""
    if length == 500:
        return reference_500_file
    depth_index = int(length * 200 / 500 - 1)
    directory_name = f"L_{length}"
    if length in (1125, 1250):
        directory_name += "_longer_integration"
    return data_directory / directory_name / f"solution_x_{depth_index:06d}.ascii"


def collect_times(
    data_directory, reference_500_file, model_lengths=MODEL_LENGTHS,
    *, tstar_years=DEFAULT_TSTAR_YEARS,
):
    """Collect the three characteristic times for each rhythmite model size."""
    collected = {
        "t_onset": [],
        "t_first_zero_U": [],
        "t_first_minimum_U": [],
    }

    for length in model_lengths:
        ascii_file = ascii_file_for_length(data_directory, reference_500_file, length)
        results = extract_and_analyse(
            ascii_file,
            make_plots=False,
            print_diagnostics=False,
            allow_missing_diagnostics=True,
            tstar_years=tstar_years,
        )
        for name in collected:
            collected[name].append(results[name])

    return {name: np.asarray(values) for name, values in collected.items()}


def collect_marlpde_onsets(
    results_directory,
    model_lengths=MARLPDE_LENGTHS,
):
    """Find the cCO3 oscillation onset at the deepest node of each run."""
    onsets = []
    for length in model_lengths:
        result_file = (
            results_directory / f"L_{length}" / "LMAHeureuxPorosityDiff.hdf5"
        )
        with h5py.File(result_file, "r") as result:
            times = np.asarray(result["times"]) * float(result.attrs["Tstar"]) / 1_000
            cCO3_at_bottom = np.asarray(result["data"][:, 3, -1])
        try:
            onset, _, _ = oscillation_onset(
                times, cCO3_at_bottom, prominence=0.01, min_extrema=6
            )
        except ValueError:
            onset = np.nan
        onsets.append(onset)
        print(f"Marlpde L={length} cm: cCO3 t_onset={onset:.3f} kyr")
    return np.asarray(onsets)


def fit_scores(observed, predicted, parameter_count):
    """Return centred R² and its residual-degrees-of-freedom adjustment.

    parameter_count counts ALL fitted parameters, including any intercept.
    For a nonlinear fit the n - parameter_count adjustment is approximate.
    """
    n = observed.size
    total_sum = np.sum((observed - observed.mean())**2)
    if total_sum <= 0:
        return np.nan, np.nan
    residual_sum = np.sum((observed - predicted)**2)
    r_squared = 1 - residual_sum / total_sum
    adjusted = (
        1 - (residual_sum / (n - parameter_count)) / (total_sum / (n - 1))
        if n > parameter_count else np.nan
    )
    return r_squared, adjusted


def plot_times(
    model_lengths, times, marlpde_lengths, marlpde_onsets,
    *, output_stem="rhythmite_marlpde_onset", show=True,
):
    """Compare measured onsets, show paired differences, and export the figure."""
    lengths = np.asarray(model_lengths, dtype=float)
    rhythmite_onsets = np.asarray(times["t_onset"], dtype=float)
    marlpde_lengths = np.asarray(marlpde_lengths, dtype=float)
    marlpde_onsets = np.asarray(marlpde_onsets, dtype=float)

    fig, (ax, difference_ax) = pyl.subplots(
        2, 1, sharex=True, figsize=(8.5, 7), constrained_layout=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08},
    )
    ax.plot(
        lengths, rhythmite_onsets, "o", color="tab:blue",
        markerfacecolor="none", markersize=10, markeredgewidth=1.5,
        label="Rhythmite: $t_\\mathrm{onset}$", zorder=3,
    )
    ax.plot(
        marlpde_lengths, marlpde_onsets, "x", color="tab:orange",
        markersize=7, markeredgewidth=1.8,
        label="Marlpde: $t_\\mathrm{onset}$ (cCO3)", zorder=4,
    )
    ax.plot(
        lengths, times["t_first_zero_U"], "s", color="0.35",
        markerfacecolor="none", markersize=7, markeredgewidth=1.3,
        label="Rhythmite: first zero of U", zorder=5,
    )
    ax.plot(
        lengths, times["t_first_minimum_U"], "^", color="0.55",
        markerfacecolor="none", markersize=10, markeredgewidth=1.3,
        label="Rhythmite: first minimum of U", zorder=6,
    )

    valid = np.isfinite(lengths) & np.isfinite(rhythmite_onsets) & (lengths > 0)
    x, y = lengths[valid], rhythmite_onsets[valid]
    if np.unique(x).size >= 2:
        parameters, covariance = curve_fit(
            power_law, x, y, p0=(y[np.argmax(x)], 2.0),
            bounds=([0, 0], [np.inf, np.inf]),
        )
        time_at_1000, b = parameters
        fit_lengths = np.linspace(x.min(), x.max(), 200)
        ax.plot(
            fit_lengths, power_law(fit_lengths, *parameters), "--",
            color="tab:blue", linewidth=1.3,
            label="Rhythmite: power-law fit", zorder=2,
        )
        r_squared, adjusted = fit_scores(y, power_law(x, *parameters), 2)
        ax.text(
            0.02, 0.98,
            rf"Rhythmite fit: $t={time_at_1000:.1f}(L/1000\,\mathrm{{cm}})^{{{b:.3f}}}$ kyr"
            + "\n" + rf"$R^2_\mathrm{{adj}}={adjusted:.5f}$",
            transform=ax.transAxes, va="top", fontsize=9,
        )
        print(f"Power-law fit: t = A (L / 1000 cm)^b, A={time_at_1000:.8f} kyr, b={b:.8f}")
        print(f"Power law: n={y.size}, k=2, R²={r_squared:.8f}, adjusted R²={adjusted:.8f}")
        x_squared = x**2
        parabola_a = np.dot(x_squared, y) / np.dot(x_squared, x_squared)
        parabola_r2, parabola_adjusted = fit_scores(y, parabola_a * x_squared, 1)
        print(f"Parabola: n={y.size}, k=1, R²={parabola_r2:.8f}, adjusted R²={parabola_adjusted:.8f}")
        for length, prediction in zip(PREDICTION_LENGTHS, power_law(PREDICTION_LENGTHS, *parameters)):
            print(f"L={length} cm: fitted t_onset={prediction:.2f} kyr")

    rhythmite_by_length = dict(zip(lengths, rhythmite_onsets))
    paired = [
        (length, 100 * (onset / rhythmite_by_length[length] - 1))
        for length, onset in zip(marlpde_lengths, marlpde_onsets)
        if length in rhythmite_by_length
        and np.isfinite(onset)
        and np.isfinite(rhythmite_by_length[length])
        and rhythmite_by_length[length] > 0
    ]
    difference_ax.axhline(0, color="0.5", linewidth=0.8)
    if paired:
        paired_lengths, differences = np.asarray(paired).T
        difference_ax.plot(
            paired_lengths, differences, "D-", color="tab:orange",
            markersize=4,
        )
    difference_ax.set_ylabel(
        "Relative difference in $t_\\mathrm{onset}$ [%]\n"
        "100 × (Marlpde - Rhythmite) / Rhythmite"
    )
    difference_ax.set_xlabel("model length L [cm]")
    difference_ax.grid(axis="y", color="0.9")
    ax.set_ylabel("time [kyr]")
    ax.legend(loc="upper left", bbox_to_anchor=(0, 0.85), frameon=False, fontsize=9)
    fig.align_ylabels()
    output_stem = Path(output_stem)
    for extension in ("pdf", "svg"):
        fig.savefig(output_stem.with_suffix(f".{extension}"), bbox_inches="tight")
    if show:
        pyl.show()
    pyl.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Plot Rhythmite and Marlpde onset times and their paired differences."
    )
    parser.add_argument(
        "--output-stem", default="rhythmite_marlpde_onset",
        help="path without extension for the PDF and SVG figures",
    )
    parser.add_argument(
        "--rhythmite-directory", type=Path, required=True,
        help="directory containing the Rhythmite L_* output directories",
    )
    parser.add_argument(
        "--reference-500-file", type=Path, required=True,
        help="Rhythmite 500 cm solution_x_000199.ascii file",
    )
    parser.add_argument(
        "--marlpde-results-directory", type=Path, required=True,
        help="directory containing the Marlpde L_* output directories",
    )
    parser.add_argument(
        "--rhythmite-tstar-years", type=float, default=DEFAULT_TSTAR_YEARS,
        help=f"Rhythmite time scale in years (default: {DEFAULT_TSTAR_YEARS:g})",
    )
    parser.add_argument(
        "--no-show", action="store_true", help="save figures without opening a window"
    )
    args = parser.parse_args()
    times = collect_times(
        args.rhythmite_directory, args.reference_500_file,
        tstar_years=args.rhythmite_tstar_years,
    )
    marlpde_onsets = collect_marlpde_onsets(args.marlpde_results_directory)
    plot_times(
        MODEL_LENGTHS, times, MARLPDE_LENGTHS, marlpde_onsets,
        output_stem=args.output_stem, show=not args.no_show,
    )


if __name__ == "__main__":
    main()
