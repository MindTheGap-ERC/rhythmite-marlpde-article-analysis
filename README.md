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
- Current Marlpde lengths: 500, 625, 750, 875, 1000, 1125, 1250 cm. Add new lengths
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

The script writes PDF and SVG. Its upper panel shows the measured onset times,
the Rhythmite power-law fit, and two Rhythmite velocity diagnostics. The lower
panel shows `100 * (Marlpde onset / Rhythmite onset - 1)` at matching lengths.
Onset is estimated from cCO3 at the deepest grid node with prominence 0.01 and
six alternating extrema. The fit uses only measured Rhythmite onsets; values
printed for 1125 and 1250 cm are fitted values within the measured range.

`plot_evolution_at_depth.py` can inspect an individual Marlpde HDF5 file;
`extract_and_analyse_data_from_ascii_files.py` can inspect an individual
Rhythmite ASCII file. Their `-h` output describes the available arguments.
