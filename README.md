# ConfSets

**Confidence sets for functional intraclass correlation coefficients (ICC) in test–retest studies.**

ConfSets estimates reliability from repeated curve measurements and supports confidence intervals for integrated ICC values and confidence bands for ICC curves. It includes simulated test and retest datasets to help you get started.

## Install

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("mr-seydi/ConfSets")
```

Then load the package:

```r
library(ConfSets)
```

## Quick start

The package includes `TEST_data` and `RETEST_data`, matrices with 30 simulated curves (rows) measured at 100 evaluation points (columns). Load them explicitly with `data()`:

```r
data("TEST_data", package = "ConfSets")
data("RETEST_data", package = "ConfSets")

icc <- functional_ICC(
  test = TEST_data,
  retest = RETEST_data,
  designation = "ICC(A,1)",
  alpha = 0.05
)

icc$integrated_ICC
```

To construct a confidence interval for the integrated ICC and simultaneous confidence bands for the ICC curve, use `Conf_sets()`:

```r
set.seed(123)

results <- Conf_sets(
  test = TEST_data,
  retest = RETEST_data,
  type = c("Rank", "Normal"),
  n_boot = 1000,
  designation = "ICC(A,1)",
  alpha = 0.05
)

results$CSets$Rank    # simultaneous bands for the ICC curve
results$CSets$Normal  # interval for the integrated ICC
plot_ICC(results)
```

`n_boot` controls the number of bootstrap samples. Increase it for final analyses when runtime allows.

## Main functions

- `functional_ICC()` calculates pointwise and integrated ICC estimates.
- `Conf_sets()` computes confidence sets: `"Rank"` for simultaneous ICC curve bands and `"Normal"` for an integrated ICC confidence interval.
- `plot_ICC()` plots the estimated ICC curve with Rank confidence bands.
- `pointwise_F_bands_ICC()` calculates pointwise F-based confidence bands. These are pointwise bands and do not provide simultaneous coverage.

The functions support `"ICC(A,1)"` (single measurement) and `"ICC(A,k)"` (average measurement) designations where documented. Supply matrices with evaluation points in columns. For balanced test and retest data, each row should correspond to the same individual in both matrices; leave `individuals = NULL`. For unbalanced `"ICC(A,1)"` data, provide the combined individual-ID vector for the test rows followed by the retest rows.

**Confidence sets for unbalanced data have not been validated.** Interpret those results with care. The package documentation also describes the data format and methods in more detail.

## Documentation

See the package vignette for the data-generation details, ICC designations, examples, and confidence-set methods:

```r
vignette("ConfSets", package = "ConfSets")
```

Function help is available in R with, for example, `?functional_ICC` and `?Conf_sets`.

## Citation

If you use ConfSets, please cite:

Seydi MR, Pini A, Pataky TC, Schelin L. Confidence sets for intraclass correlation coefficients in test–retest curve measurements. *Journal of Biomechanics*. 2024;173:112232.

```r
citation("ConfSets")
```

## Authors

- Mohammad Reza Seydi — [ORCID](https://orcid.org/0009-0002-3278-5554)
- Lina Schelin — [ORCID](https://orcid.org/0000-0001-7917-5687)

## License

MIT. See [LICENSE](LICENSE) for details.
