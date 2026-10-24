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

## Fixed-depth branch — 24 October 2026

Explicit user decision: replace bottom-node diagnostics with measurements at
500 cm for Figure 5. Interpolation/extrapolation is authorized. Work on
analysis/fixed-depth-500cm; do not modify the manuscript or merge/commit without
an instruction. The historical bottom analysis below remains a validation
reference. Use plot_fixed_depth.py for the fixed-depth figure and record the
actual bracketing nodes, source files, and VTU/ASCII validation in its JSON.
Do not rerun simulations unless missing data make this necessary and discussed.

## Analysis definition (historical bottom-node figure)

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

Decision updated 14 October 2026 by explicit user instruction: treat the
fine-sampled Rhythmite result at 1625 cm (44280.808 kyr) as an ordinary
observation. Include all eight lengths in the relative-error least-squares
power-law fit, minimizing sum(((fit - onset) / onset)**2) with sigma=onset.
Report RMS and maximum absolute percentage errors instead of R². Do not mark it provisional or draw an extrapolation segment.
The coarser output gives 44275.533 kyr; use the fine-sampled result.

Numerical convergence remains to be investigated for the whole series.

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
- Do not horizontally offset overlapping data markers.

## Open scientific work

- Investigate spatial and temporal convergence for the full series.
- Assess whether the empirical power law adequately describes the data.
- Insert the final PDF figure and an updated caption into the Overleaf article.
- Record model commit identifiers, parameters, input locations, and archive
  identifiers before publication.
## AI contribution records in commit messages

Explicit user instruction (14 October 2026): every commit with any substantive
scientific/manuscript/analysis contribution from an AI assistant must include
an AI contribution section in its commit message. This includes analysis code,
method changes, figures and captions. Do not describe such work merely as
language polishing. This record will later support an author-reviewed AI-use
disclosure in Materials and Methods under AGU's AI Tools Policy:
https://www.agu.org/publications/authors/policies

Use a section such as:

AI assistance:
- Tool: OpenAI Codex (record the actual tool; include model/version only when
  reliably known, never infer it).
- Contribution: concrete files, code, analysis choices, generated text or figures
  created or modified with AI assistance.
- Provenance: for imported figures, identify source repository, commit and script;
  record uncommitted source changes if applicable. Distinguish plotting model
  output with code from generative image creation.
- Validation: checks actually performed by the assistant and checks actually
  confirmed by a human. Do not imply co-author/scientific approval from a
  successful build or an approval of layout; state outstanding review if relevant.

Include this information in proposed commit messages as well as commits made
by the assistant. Keep the title concise and put details in the body. Do not
rewrite existing commit history to add disclosures without explicit permission.

Manuscript text changes require Hanno's prior approval so he can consult
co-authors. An authorization to edit a figure or code is not blanket permission
to revise manuscript prose or insert the final AI-use statement.
