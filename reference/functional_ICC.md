# Functional intraclass correlation coefficient, ICC, in Test-Retest Reliability

This function calculates the pointwise and integrated intraclass
correlation coefficient, ICC, for a functional data set in Test-Retest
setting.

## Usage

``` r
functional_ICC(
  test,
  retest,
  individuals = NULL,
  designation = "ICC(A,1)",
  alpha = 0.05
)
```

## Arguments

- test:

  A matrix of test measurements with rows representing individuals and
  columns representing evaluation points.

- retest:

  A matrix of retest measurements with rows representing individuals and
  columns representing evaluation points.

- individuals:

  A vector of individual identifiers. It must be the combined vector of
  individual identifiers for test and retest. If NULL, the function
  assumes that the individuals in the test and retest matrices are the
  same.

- designation:

  A character string specifying the designation of the ICC. Possible
  values are "ICC(A,1)" and "ICC(A,k)" Possible designation is
  "ICC(A,1)", "ICC(A,k)" based on the paper: McGraw KO, Wong SP. Forming
  inferences about some intraclass correlation coefficients.
  Psychological methods. 1996 Mar;1(1):30. The default value is
  "ICC(A,1)".

- alpha:

  A numeric value specifying the significance level for the ICC
  calculation. The default value is 0.05.

## Value

A list containing the pointwise and integrated ICC values.

## Examples

``` r
# The package contains the datasets, TEST_data, and RETEST_data which are the
#test and retest data, respectively.

output_balanced_ICCA1 <- functional_ICC(test = TEST_data, retest = RETEST_data,
individuals = NULL, designation = "ICC(A,1)", alpha = 0.05)

# The output_balanced_ICCA1 contains the pointwise and integrated ICC values.
output_balanced_ICCA1$pointwise_ICC
#>   [1] 0.7232935 0.7173577 0.7166171 0.7214663 0.7314508 0.7453887 0.7617302
#>   [8] 0.7789596 0.7958703 0.8116618 0.8259045 0.8384457 0.8493108 0.8586240
#>  [15] 0.8665537 0.8732795 0.8789742 0.8837947 0.8878794 0.8913482 0.8943037
#>  [22] 0.8968338 0.8990134 0.9009069 0.9025696 0.9040495 0.9053884 0.9066230
#>  [29] 0.9077854 0.9089038 0.9100028 0.9111038 0.9122247 0.9133800 0.9145810
#>  [36] 0.9158353 0.9171468 0.9185158 0.9199389 0.9214091 0.9229161 0.9244465
#>  [43] 0.9259844 0.9275120 0.9290100 0.9304587 0.9318378 0.9331277 0.9343099
#>  [50] 0.9353667 0.9362829 0.9370528 0.9376803 0.9381713 0.9385330 0.9387735
#>  [57] 0.9389016 0.9389269 0.9388591 0.9387084 0.9384843 0.9381966 0.9378542
#>  [64] 0.9374653 0.9370374 0.9365768 0.9360888 0.9355773 0.9350451 0.9344935
#>  [71] 0.9339223 0.9333299 0.9327130 0.9320663 0.9313830 0.9306537 0.9298667
#>  [78] 0.9290075 0.9280582 0.9269967 0.9257964 0.9244249 0.9228432 0.9210042
#>  [85] 0.9188517 0.9163195 0.9133311 0.9098007 0.9056368 0.9007507 0.8950715
#>  [92] 0.8885714 0.8813006 0.8734309 0.8652935 0.8573884 0.8503387 0.8447760
#>  [99] 0.8411871 0.8397857
output_balanced_ICCA1$integrated_ICC
#> [1] 0.8948601

#Consider the scenario with unbalanced case, in total we have 28 individuals,
#and for example individual 2 has 1 test and 2 retest trials, individual 3 has
#3 test and 2 retest trials, and the other individuals have one trial per session.
#Then the individuals vector should be like
# c(1:2, rep(3, 3), 4:28, 1, rep(2,2), rep(3, 2), 4:28).

output_unbalanced_ICCA1 <- functional_ICC(test = TEST_data, retest = RETEST_data,
individuals = c(1:2, rep(3, 3), 4:28, 1, rep(2,2), rep(3, 2), 4:28),
 designation = "ICC(A,1)", alpha = 0.05)

output_unbalanced_ICCA1$pointwise_ICC
#>   [1] 0.4048062 0.3875587 0.3837391 0.3945786 0.4190687 0.4542392 0.4960855
#>   [8] 0.5406308 0.5846673 0.6260292 0.6635164 0.6966576 0.7254548 0.7501767
#>  [15] 0.7712163 0.7890021 0.8039490 0.8164347 0.8267905 0.8353005 0.8422051
#>  [22] 0.8477061 0.8519725 0.8551458 0.8573451 0.8586716 0.8592127 0.8590450
#>  [29] 0.8582379 0.8568557 0.8549599 0.8526113 0.8498710 0.8468017 0.8434686
#>  [36] 0.8399391 0.8362828 0.8325711 0.8288756 0.8252675 0.8218152 0.8185832
#>  [43] 0.8156301 0.8130073 0.8107577 0.8089147 0.8075021 0.8065338 0.8060143
#>  [50] 0.8059397 0.8062978 0.8070673 0.8082191 0.8097204 0.8115355 0.8136263
#>  [57] 0.8159529 0.8184743 0.8211490 0.8239360 0.8267954 0.8296890 0.8325810
#>  [64] 0.8354381 0.8382304 0.8409310 0.8435165 0.8459668 0.8482647 0.8503960
#>  [71] 0.8523490 0.8541142 0.8556837 0.8570510 0.8582104 0.8591565 0.8598839
#>  [78] 0.8603866 0.8606577 0.8606886 0.8604691 0.8599867 0.8592263 0.8581704
#>  [85] 0.8567994 0.8550923 0.8530288 0.8505923 0.8477762 0.8445922 0.8410835
#>  [92] 0.8373412 0.8335224 0.8298648 0.8266852 0.8243493 0.8232061 0.8234960
#>  [99] 0.8252710 0.8283677
output_unbalanced_ICCA1$integrated_ICC
#> [1] 0.7938654

```
