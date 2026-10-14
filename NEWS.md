# ConfSets 1.1.1

## Changes

- Set the package citation to the 2024 Journal of Biomechanics article.
- Remove the unbalanced-data caution from the README.
- Clarify the balanced-data definition and negative ICC estimates in the vignette.

# ConfSets 1.1.0

## New features

- Add the `resample_trials` option to `Conf_sets()`. The default (`FALSE`) resamples individuals while keeping their observed trial rows; `TRUE` also resamples trial rows within each individual and session.
- Document the option in the README, function help, and package vignette.
