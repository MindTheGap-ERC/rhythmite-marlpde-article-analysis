# Rhythmite and Marlpde onset analysis

Analysis code and figure for the manuscript comparing oscillation onset in
Rhythmite and Marlpde. The four Python scripts here are copies of the working
scripts in the Marlpde checkout; ongoing model runs do not depend on this copy.

## Inputs

- Rhythmite ASCII output: `L_*` directories plus the separate 500 cm reference
  file. The ASCII files have no stored time scale, so the analysis uses
  `Tstar = 13190` years by default; change it with
  `--rhythmite-tstar-years` if the Rhythmite runs used another value.
- Marlpde HDF5 output: `L_*/LMAHeureuxPorosityDiff.hdf5`. The time scale is
  read from each file's `Tstar` attribute.
- Current Marlpde lengths: 500, 625, 750, 875, 1000, 1125, 1250, 1625 cm. Add new lengths
  to `MARLPDE_LENGTHS` when their completed output is available.

The input data are not copied into this repository. Record the model commit,
run parameters, file locations, and archive identifiers with the manuscript.

## Environment

Use Python 3.11 and [Pipenv](https://pipenv.pypa.io/). Install the locked
dependencies with:

```bash
pipenv sync
```

`Pipfile.lock` records the complete environment used for the analysis. Use
`pipenv install` when intentionally changing a dependency and commit the
updated `Pipfile.lock`.

## Make the figure

```bash
pipenv run python plot_rhythmite_times_against_length.py \
  --rhythmite-directory /path/to/rhythmite_manuscript_runs/Larger_system \
  --reference-500-file /path/to/rhythmite_manuscript_runs/resolution/200nnx/solution_x_000199.ascii \
  --marlpde-results-directory /path/to/Results \
  --output-stem figures/rhythmite_marlpde_onset --no-show
```

The script writes PDF and SVG to `figures/` by default, creating the
output directory if needed. Use `--output-stem` to select another path.
The upper panel shows onset and velocity diagnostics; the lower panel shows
`100 * (Marlpde - Rhythmite) / Rhythmite` for the onset times.

The power-law fit uses all eight Rhythmite onsets over 500–1625 cm,
including the fine-sampled 1625 cm result as an ordinary observation.
Weighted least squares minimizes `sum(((fit - onset) / onset)**2)`,
using `sigma=onset` in `curve_fit` (weights proportional to `1/onset**2`).
This weights equal percentage errors equally; the weights are not measured
uncertainties. This is not a logarithmic fit. The figure reports RMS and
maximum absolute percentage errors instead of R²; signed per-point residuals
are printed to the console.
The logarithmic display axis does not change the fitting objective.
The fit is an empirical description, not evidence of a physical scaling law.

Onset detection is unchanged: cCO3 at the deepest grid point, prominence
0.01, and six alternating extrema. Numerical convergence remains an open
question for the full series, without singling out 1625 cm.

`plot_evolution_at_depth.py` can inspect an individual Marlpde HDF5 file;
`extract_and_analyse_data_from_ascii_files.py` can inspect an individual
Rhythmite ASCII file. Their `-h` output describes the available arguments.

## Diagnose an LSODA run

Run `scripts/diagnose_lsoda.py` with the **model's** Python environment (it requires
the model's py-pde dependencies). From the Marlpde checkout:

```bash
pipenv run python ../rhythmite-marlpde-article-analysis/scripts/diagnose_lsoda.py \
  --model-directory marlpde \
  --output-directory ../Results/LSODA_diagnostic_01
```

The output directory must be new. The wrapper uses the current scenario and
solver settings, requires LSODA with `lband=uband=None`, and leaves model files
unchanged. It replaces the trial-time progress bar with periodic JSON records:
`accepted_kyr` is accepted model time, `last_h` is the last accepted step in
dimensionless time, and `rhs_calls`, `njev`, and `nlu` count solver work.
Reports occur about every ten seconds when Python callbacks are reached;
`stacks.log` independently records Python stacks every minute, including while
native solver work blocks callbacks. Those stacks cannot resolve native Fortran
or BLAS internals. Counters for Jacobian work can lag during a solver step.

`checkpoint.npz` is atomically replaced at reports with an accepted state;
it can therefore lag the latest accepted step. It contains `t`, flat `y`, and
`Tstar`, not LSODA's internal history, and is not a seamless solver restart.
Non-finite trial states or derivatives stop the run immediately and are saved
in `invalid_evaluation.npz`. This diagnostic stop is stricter than normal
integration, since trial states need not be accepted states. `parameters.json`
records settings; `sampled_solution.npz` is written if solve_ivp returns.
These are diagnostic artifacts, not the normal manuscript HDF5 input.
Use a new output directory for every run. For a brief setup check, optionally
add `--end-time 0.000002` (dimensionless time).
