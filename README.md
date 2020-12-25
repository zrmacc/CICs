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
## 1 0.2343015      1   1      3
## 2 2.4710715      2   1      3
## 3 1.4048439      1   1      1
## 4 0.4910617      0   1      3
## 5 0.4647632      1   1      1
## 6 0.2193612      1   1      2
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
## 1     0 AUC     200 0.833
## 2     1 AUC     200 0.555
## 
## 
## CIs:
##            contrast        est      lower      upper
## boot_diff     A1-A0 -0.2776488 -0.3984006 -0.1596547
## boot_ratio    A1/A0  0.6667226  0.5294923  0.7933343
## 
## 
## P-values:
##   contrast        est     perm_p     boot_p
## 1    A1-A0 -0.2776488 0.01980198 0.01980198
## 2    A1/A0  0.6667226 0.01980198 0.01980198
## 
## 
## Stratum Stats:
## # A tibble: 3 x 9
##   strata weight stat     n0    n1  est0  est1    diff ratio
##    <int>  <dbl> <chr> <int> <int> <dbl> <dbl>   <dbl> <dbl>
## 1      1  0.265 AUC      56    50 0.927 0.457 -0.470  0.493
## 2      2  0.365 AUC      71    75 0.688 0.650 -0.0374 0.946
## 3      3  0.37  AUC      73    75 0.909 0.532 -0.377  0.585
```

Replace "AUC" by "AOC" for area over the cumulative incidence curve.

### Compare Event Rates

To find a confidence interval and p-vaue for the difference and ratio in event rates at $\tau = 28$ days:

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
## 1     0 Rate    200 0.668
## 2     1 Rate    200 0.495
## 
## 
## CIs:
##            contrast        est      lower      upper
## boot_diff     A1-A0 -0.1729386 -0.2583795 -0.0452592
## boot_ratio    A1/A0  0.7410808  0.6350871  0.9243429
## 
## 
## P-values:
##   contrast        est     perm_p     boot_p
## 1    A1-A0 -0.1729386 0.01980198 0.01980198
## 2    A1/A0  0.7410808 0.01980198 0.01980198
## 
## 
## Stratum Stats:
## # A tibble: 3 x 9
##   strata weight stat     n0    n1  est0  est1    diff ratio
##    <int>  <dbl> <chr> <int> <int> <dbl> <dbl>   <dbl> <dbl>
## 1      1  0.265 Rate     56    50 0.687 0.428 -0.259  0.623
## 2      2  0.365 Rate     71    75 0.587 0.564 -0.0226 0.961
## 3      3  0.37  Rate     73    75 0.735 0.475 -0.260  0.647
```

### Compare Medians

To find a confidence interval and p-vaue for the difference and ratio in medians $q = 0.5$:

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
## 1     0 Quantile   200  1.18
## 2     1 Quantile   200  2.48
## 
## 
## CIs:
##            contrast      est     lower    upper
## boot_diff     A1-A0 1.305356 0.2967083 1.890002
## boot_ratio    A1/A0 2.108569 1.1752127 2.680476
## 
## 
## P-values:
##   contrast      est     perm_p     boot_p
## 1    A1-A0 1.305356 0.02083333 0.03846154
## 2    A1/A0 2.108569 0.02083333 0.03846154
## 
## 
## Stratum Stats:
## # A tibble: 3 x 9
##   strata weight stat        n0    n1  est0  est1  diff ratio
##    <int>  <dbl> <chr>    <int> <int> <dbl> <dbl> <dbl> <dbl>
## 1      1  0.265 Quantile    56    50 0.860  3.87 3.01   4.50
## 2      2  0.365 Quantile    71    75 1.54   1.68 0.140  1.09
## 3      3  0.37  Quantile    73    75 1.04   2.28 1.23   2.18
```
