---
title: "README"
author: "Zachary McCaw"
date: "2026-02-22"
output: 
  html_document: 
    keep_md: TRUE
--- 

# Inference for Cumulative Incidence Curves

[![R-CMD-check](https://github.com/zrmacc/CICs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zrmacc/CICs/actions/workflows/R-CMD-check.yaml)

Zachary McCaw <br>
Updated: 2026-02-22



### Description

This package provides functions for inference on the difference in AUCs, event rates, and quantiles comparing two cumulative incidence curves. Also see [MCC](https://github.com/zrmacc/MCC) for inference on mean cumulative count curves. 

## Installation


``` r
devtools::install_github(repo = "zrmacc/CICs")
```

## Examples

### Data

The function `GenTwoSampleData` simulates example data in the format expected by this package. The censoring, event, and death times are drawn from independent exponential distributions. Note that 'rate' refers to the arrival rate of the corresponding exponential rather than the proportion of the sample. 


``` r
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
## 1 0.2727276      2   1      1
## 2 1.6937283      2   1      2
## 3 0.1898692      1   1      3
## 4 2.6857049      1   1      3
## 5 0.5788937      2   1      1
## 6 3.7567580      2   1      2
```

In these data, `arm` is the treatment arm, 0 for reference, 1 for treatment; `time` is the observation time; and `status` is the event type, 0 for censoring, 1 for an event, 2 for death. For analysing other data sets, `arm` should likewise be coded as 0/1, and status as 0/1/2, with status 1 identifying the event of interest.

### Compare AUCs

To find a confidence interval and p-vaue for the difference and ratio in areas under the cumulative incidence curve at time $\tau = 2$:

``` r
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
## # A tibble: 2 × 4
##     arm stat      n   est
##   <dbl> <chr> <int> <dbl>
## 1     0 AUC     200 0.860
## 2     1 AUC     200 0.719
## 
## 
## CIs:
##            contrast        est      lower      upper
## boot_diff     A1-A0 -0.1409598 -0.2928194 0.01627149
## boot_ratio    A1/A0  0.8360528  0.6809977 1.02110234
## 
## 
## P-values:
##   contrast        est    perm_p    boot_p
## 1    A1-A0 -0.1409598 0.1386139 0.1188119
## 2    A1/A0  0.8360528 0.1584158 0.1188119
## 
## 
## Stratum Stats:
## # A tibble: 3 × 9
##   strata weight stat     n0    n1  est0  est1    diff ratio
##    <int>  <dbl> <chr> <int> <int> <dbl> <dbl>   <dbl> <dbl>
## 1      1  0.295 AUC      54    64 0.944 0.757 -0.187  0.802
## 2      2  0.425 AUC      96    74 0.761 0.780  0.0193 1.03 
## 3      3  0.28  AUC      50    62 0.922 0.586 -0.336  0.636
```

Replace "AUC" by "AOC" for area over the cumulative incidence curve.

### Compare Event Rates

For inference on the difference and ratio of event rates at $\tau = 2$:

``` r
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
## # A tibble: 2 × 4
##     arm stat      n   est
##   <dbl> <chr> <int> <dbl>
## 1     0 Rate    200 0.683
## 2     1 Rate    200 0.535
## 
## 
## CIs:
##            contrast        est      lower      upper
## boot_diff     A1-A0 -0.1482881 -0.2664158 -0.0564455
## boot_ratio    A1/A0  0.7829757  0.6163373  0.9091545
## 
## 
## P-values:
##   contrast        est     perm_p     boot_p
## 1    A1-A0 -0.1482881 0.03960396 0.01980198
## 2    A1/A0  0.7829757 0.03960396 0.01980198
## 
## 
## Stratum Stats:
## # A tibble: 3 × 9
##   strata weight stat     n0    n1  est0  est1    diff ratio
##    <int>  <dbl> <chr> <int> <int> <dbl> <dbl>   <dbl> <dbl>
## 1      1  0.295 Rate     54    64 0.735 0.527 -0.208  0.717
## 2      2  0.425 Rate     96    74 0.611 0.575 -0.0362 0.941
## 3      3  0.28  Rate     50    62 0.738 0.482 -0.256  0.654
```

### Compare Medians

To compare the difference and ratio of medians:

``` r
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
## # A tibble: 2 × 4
##     arm stat         n   est
##   <dbl> <chr>    <int> <dbl>
## 1     0 Quantile   200  1.11
## 2     1 Quantile   200  1.53
## 
## 
## CIs:
##            contrast       est       lower   upper
## boot_diff     A1-A0 0.4220736 -0.04348134 1.42147
## boot_ratio    A1/A0 1.3817362  0.96480772 2.31713
## 
## 
## P-values:
##   contrast       est     perm_p    boot_p
## 1    A1-A0 0.4220736 0.03960396 0.0990099
## 2    A1/A0 1.3817362 0.03960396 0.0990099
## 
## 
## Stratum Stats:
## # A tibble: 3 × 9
##   strata weight stat        n0    n1  est0  est1   diff ratio
##    <int>  <dbl> <chr>    <int> <int> <dbl> <dbl>  <dbl> <dbl>
## 1      1  0.295 Quantile    54    64 0.973  1.14 0.165   1.17
## 2      2  0.425 Quantile    96    74 1.38   1.45 0.0754  1.05
## 3      3  0.28  Quantile    50    62 0.834  2.05 1.22    2.46
```
