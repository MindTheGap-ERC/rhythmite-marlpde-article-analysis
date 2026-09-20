# Project context for Codex

## Purpose

This repository contains the reproducible analysis for a manuscript comparing
oscillation onset in the Rhythmite and Marlpde implementations of the
L'Heureux limestone-marl model. Keep analysis code here rather than in either
model implementation repository.

## Environment and workflow

- Use Python 3.11 with Pipenv.
- Install the locked environment with `pipenv sync`.
- Run commands with `pipenv run python ...`.
- Keep `Pipfile.lock` committed when dependencies change.
- Input model output is external and must not be committed to this repository.
- Generate the manuscript figure as both PDF and SVG under `figures/`.
- After changing the analysis, run the full figure command from `README.md` and
  inspect the PDF visually before committing.

## Analysis definition

- Determine onset from cCO3 at the deepest grid node.
- Use peak prominence 0.01 and require six alternating extrema.
- `oscillation_onset` defines onset as half an estimated period before the
  first extremum in the sustained oscillatory regime.
- Read Marlpde `Tstar` from each HDF5 file.
- Rhythmite ASCII output does not store `Tstar`; use 13190 years unless the run
  metadata establish another value.
- The main figure uses a logarithmic time axis because onset spans more than
  two orders of magnitude.
- The lower panel is
  `100 * (Marlpde t_onset - Rhythmite t_onset) / Rhythmite t_onset`.

## Current data and decisions

Marlpde cCO3 onset times at the deepest node, in kyr:

- 500 cm: 153.334
- 625 cm: 372.354
- 750 cm: 781.507
- 875 cm: 1509.595
- 1000 cm: 2773.197
- 1125 cm: 4948.888
- 1250 cm: 8653.629
- 1625 cm: 44252.450

The fine-sampled Rhythmite run at 1625 cm gives 44280.808 kyr. A longer,
coarser Rhythmite output gives 44275.533 kyr. Include the fine-sampled value as
a provisional point and in the paired relative difference. Exclude it from the
power-law fit until its numerical robustness has been investigated further.

Fit the Rhythmite power law only over 500--1250 cm. The current fit is
`t = 2824.35599703 * (L / 1000 cm)^4.99201811 kyr`, with adjusted R-squared
0.99888093. Its extrapolation to 1625 cm is 31878.89 kyr. Show the fit as a
dashed line over its fitted range and its extrapolation as a dotted line.

At 1625 cm, Marlpde is about 28.36 kyr earlier than fine-sampled Rhythmite, a
relative difference of about -0.064 percent. The detected onset is unchanged
for tested prominence values from 0.001 to 0.05 in Rhythmite and from 0.005 to
0.03 in Marlpde. Both implementations used the same spatial resolution, so a
shared resolution effect has not yet been excluded.

## Figure conventions

- Rhythmite onset: large open blue circle.
- Marlpde onset: orange cross.
- Rhythmite first zero of U: open dark-gray square.
- Rhythmite first minimum of U: open gray triangle.
- Label the Rhythmite 1625 cm point as provisional.
- Do not horizontally offset overlapping data markers.

## Open scientific work

- Investigate whether the 1625 cm oscillations persist at another spatial
  resolution before treating the Rhythmite point as definitive.
- Decide how the 1625 cm result affects the interpretation of the power-law
  extrapolation in the manuscript.
- Insert the final PDF figure and an updated caption into the Overleaf article.
- Record model commit identifiers, parameters, input locations, and archive
  identifiers before publication.
