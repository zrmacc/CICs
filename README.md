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
## 1 1.5342347      0   1      1
## 2 1.8967266      0   1      2
## 3 0.2310048      1   1      2
## 4 1.0370315      0   1      3
## 5 1.0091773      1   1      2
## 6 0.3223376      0   1      3
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
## 1     0 AUC     200 0.795
## 2     1 AUC     200 0.554
## 
## 
## CIs:
##            contrast        est      lower      upper
## boot_diff     A1-A0 -0.2413093 -0.3842677 -0.1174986
## boot_ratio    A1/A0  0.6964183  0.5420529  0.8408382
## 
## 
## P-values:
##   contrast        est     perm_p     boot_p
## 1    A1-A0 -0.2413093 0.01980198 0.01980198
## 2    A1/A0  0.6964183 0.01980198 0.01980198
## 
## 
## Stratum Stats:
## # A tibble: 3 x 9
##   strata weight stat     n0    n1  est0  est1    diff ratio
##    <int>  <dbl> <chr> <int> <int> <dbl> <dbl>   <dbl> <dbl>
## 1      1  0.295 AUC      52    66 0.869 0.533 -0.337  0.613
## 2      2  0.442 AUC      91    86 0.776 0.475 -0.301  0.612
## 3      3  0.262 AUC      57    48 0.742 0.709 -0.0334 0.955
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
## 1     0 Rate    200 0.637
## 2     1 Rate    200 0.428
## 
## 
## CIs:
##            contrast        est      lower      upper
## boot_diff     A1-A0 -0.2088132 -0.2936807 -0.1178695
## boot_ratio    A1/A0  0.6721585  0.5598392  0.8010745
## 
## 
## P-values:
##   contrast        est     perm_p     boot_p
## 1    A1-A0 -0.2088132 0.01980198 0.01980198
## 2    A1/A0  0.6721585 0.01980198 0.01980198
## 
## 
## Stratum Stats:
## # A tibble: 3 x 9
##   strata weight stat     n0    n1  est0  est1    diff ratio
##    <int>  <dbl> <chr> <int> <int> <dbl> <dbl>   <dbl> <dbl>
## 1      1  0.295 Rate     52    66 0.703 0.442 -0.260  0.630
## 2      2  0.442 Rate     91    86 0.616 0.365 -0.250  0.593
## 3      3  0.262 Rate     57    48 0.599 0.518 -0.0808 0.865
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
## 1     0 Quantile   200  1.24
## 2     1 Quantile   200  2.17
## 
## 
## CIs:
##            contrast       est     lower upper
## boot_diff     A1-A0 0.9271617 0.3817628   Inf
## boot_ratio    A1/A0 1.7462949 1.2497082   Inf
## 
## 
## P-values:
##   contrast       est     perm_p     boot_p
## 1    A1-A0 0.9271617 0.07920792 0.01980198
## 2    A1/A0 1.7462949 0.05940594 0.01980198
## 
## 
## Stratum Stats:
## # A tibble: 3 x 9
##   strata weight stat        n0    n1  est0  est1   diff ratio
##    <int>  <dbl> <chr>    <int> <int> <dbl> <dbl>  <dbl> <dbl>
## 1      1  0.295 Quantile    52    66 0.952  2.58 1.62    2.71
## 2      2  0.442 Quantile    91    86 1.34   2.33 0.993   1.74
## 3      3  0.262 Quantile    57    48 1.40   1.44 0.0336  1.02
```
