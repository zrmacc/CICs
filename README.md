 # Inference for Cumulative Incidence Curves

Zachary McCaw <br>
Updated: 20-12-24



### Description

This package provides functions for inference on the difference in AUCs, event rates, and quantiles comparing two cumulative incidence curves. Also see [MCC](https://github.com/zrmacc/MCC) for inference on mean cumulative count curves. 

## Installation


```r
devtools::install_github(repo = "zrmacc/CICs")
```

## Examples

### Data

The function `GenTwoSampleData` simulates example data in the format expected by this package. The censoring, event, and death times are drawn from independent exponential distributions. Note that 'rate' refers to the arrival rate of the corresponding exponential rather than the proportion of the sample. 


```r
library(CICs)

# Generate data.
data <- GenTwoSampleData(
  n1 = 200,
  n0 = 200,
  censor_rate1 = 0.25,
  censor_rate0 = 0.25,
  event_rate1 = 0.50,
  event_rate0 = 0.75,
  death_rate1 = 0.25,
  death_rate0 = 0.25
)

# Add strata.
strata <- rmultinom(n = 400, size = 1, prob = c(0.3, 0.4, 0.3))
data$strata <- apply(strata, 2, which.max)
head(data)
```

```
##        time status arm strata
## 1 0.4669380      1   1      2
## 2 2.0153491      1   1      3
## 3 0.9062862      0   1      3
## 4 1.9460783      0   1      1
## 5 0.6581468      2   1      2
## 6 0.5718448      1   1      1
```

In these data, `arm` is the treatment arm, 0 for reference, 1 for treatment; `time` is the observation time; and `status` is the event type, 0 for censoring, 1 for an event, 2 for death. For analysing other data sets, `arm` should likewise be coded as 0/1, and status as 0/1/2, with status 1 identifying the event of interest.

### Compare AUCs

To find a confidence interval and p-vaue for the difference and ratio in areas under the cumulative incidence curve at time $\tau = 2$:

```r
aucs <- CompareCICs(
  time = data$time,
  status = data$status,
  arm = data$arm,
  strata = data$strata,
  sum_stat = "AUC",
  param = 2,
  reps = 100,
  alpha = 0.05
)
show(aucs)
```

```
## Marginal Stats:
## # A tibble: 2 x 4
##     arm stat      n   est
##   <dbl> <chr> <int> <dbl>
## 1     0 AUC     200 0.845
## 2     1 AUC     200 0.683
## 
## 
## CIs:
##            contrast        est      lower       upper
## boot_diff     A1-A0 -0.1624119 -0.2684778 -0.04445343
## boot_ratio    A1/A0  0.8078946  0.6961344  0.94427490
## 
## 
## P-values:
##   contrast        est     perm_p     boot_p
## 1    A1-A0 -0.1624119 0.05940594 0.05940594
## 2    A1/A0  0.8078946 0.07920792 0.05940594
## 
## 
## Stratum Stats:
## # A tibble: 3 x 9
##   strata weight stat     n0    n1  est0  est1   diff ratio
##    <int>  <dbl> <chr> <int> <int> <dbl> <dbl>  <dbl> <dbl>
## 1      1   0.3  AUC      48    72 0.830 0.729 -0.101 0.878
## 2      2   0.44 AUC      96    80 0.867 0.651 -0.216 0.751
## 3      3   0.26 AUC      56    48 0.827 0.685 -0.142 0.828
```

Replace "AUC" by "AOC" for area over the cumulative incidence curve.

### Compare Event Rates

For inference on the difference and ratio of event rates at $\tau = 2$:

```r
rates <- CompareCICs(
  time = data$time,
  status = data$status,
  arm = data$arm,
  strata = data$strata,
  sum_stat = 'Rate',
  param = 2,
  reps = 100,
  alpha = 0.05
)
show(rates)
```

```
## Marginal Stats:
## # A tibble: 2 x 4
##     arm stat      n   est
##   <dbl> <chr> <int> <dbl>
## 1     0 Rate    200 0.660
## 2     1 Rate    200 0.505
## 
## 
## CIs:
##            contrast        est      lower       upper
## boot_diff     A1-A0 -0.1548808 -0.2309364 -0.05149196
## boot_ratio    A1/A0  0.7651985  0.6618048  0.91587069
## 
## 
## P-values:
##   contrast        est     perm_p     boot_p
## 1    A1-A0 -0.1548808 0.01980198 0.01980198
## 2    A1/A0  0.7651985 0.01980198 0.01980198
## 
## 
## Stratum Stats:
## # A tibble: 3 x 9
##   strata weight stat     n0    n1  est0  est1    diff ratio
##    <int>  <dbl> <chr> <int> <int> <dbl> <dbl>   <dbl> <dbl>
## 1      1   0.3  Rate     48    72 0.666 0.583 -0.0829 0.876
## 2      2   0.44 Rate     96    80 0.638 0.449 -0.189  0.704
## 3      3   0.26 Rate     56    48 0.689 0.508 -0.180  0.738
```

### Compare Medians

To compare the difference and ratio of medians:

```r
quants <- CompareCICs(
  time = data$time,
  status = data$status,
  arm = data$arm,
  strata = data$strata,
  sum_stat = 'Quantile',
  param = 0.5,
  reps = 100,
  alpha = 0.05
)
show(quants)
```

```
## Marginal Stats:
## # A tibble: 2 x 4
##     arm stat         n   est
##   <dbl> <chr>    <int> <dbl>
## 1     0 Quantile   200  1.06
## 2     1 Quantile   200  1.85
## 
## 
## CIs:
##            contrast       est      lower    upper
## boot_diff     A1-A0 0.7963396 0.08973706 2.416764
## boot_ratio    A1/A0 1.7527690 1.07271141 3.495009
## 
## 
## P-values:
##   contrast       est     perm_p     boot_p
## 1    A1-A0 0.7963396 0.05940594 0.04301075
## 2    A1/A0 1.7527690 0.03960396 0.04301075
## 
## 
## Stratum Stats:
## # A tibble: 3 x 9
##   strata weight stat        n0    n1  est0  est1   diff ratio
##    <int>  <dbl> <chr>    <int> <int> <dbl> <dbl>  <dbl> <dbl>
## 1      1   0.3  Quantile    48    72  1.07  1.29 0.213   1.20
## 2      2   0.44 Quantile    96    80  1.00  2.61 1.61    2.61
## 3      3   0.26 Quantile    56    48  1.14  1.23 0.0937  1.08
```
