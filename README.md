# Description

This package provides functions for inference on the difference in event rates, quantiles, and areas under the cumulative incidence curves between two treatment arms. 

# Installation


```r
devtools::install_github(repo = 'zrmacc/CICs')
```

# Examples

Synthetic example data in the format expected by this package may be loaded via:


```r
library(CICs)
data(cic_data)
head(cic_data)
```

```
##   arm time status
## 1   0    1      0
## 2   0    1      1
## 3   0    1      1
## 4   0    1      1
## 5   0    1      1
## 6   0    1      1
```

In these data, `arm` is the treatment arm, 0 for reference, 1 for treatment; `time` is the observation time in days; and `status` is the event type, 0 for censoring, 1 for recovery, 2 for death. For analysing other data sets, `arm` should likewise be coded 0/1, and the event with status 1 is assumed to be the event of interest. 

## Compare Event Rates

To find a confidence interval and p-vaue for the difference and ratio in event rates at $\tau = 28$ days:

```r
event_rates <- CompareEventRates(
  time = cic_data$time,
  status = cic_data$status,
  arm = cic_data$arm,
  tau = 28,
  reps = 2000,
  alpha = 0.05
)
show(event_rates)
```

```
##   Time      Arm0      Arm1 Contrast   Estimate           L         U          P
## 1   28 0.6931547 0.7492792    A1-A0 0.05612454 -0.01251013 0.1247889 0.09995002
## 2   28 0.6931547 0.7492792    A1/A0 1.08096972  0.98282756 1.1925455 0.10044978
```

## Compare Medians

To find a confidence interval and p-vaue for the difference and ratio in medians $q = 0.5$:

```r
medians <- CompareQuants(
  time = cic_data$time,
  status = cic_data$status,
  arm = cic_data$arm,
  q = 0.5,
  reps = 2000,
  alpha = 0.05
)
show(medians)
```

```
##   Prob Arm0 Arm1 Contrast Estimate    L          U            P
## 1  0.5   15    9    A1-A0     -6.0 -8.0 -3.0000000 0.0004997501
## 2  0.5   15    9    A1/A0      0.6  0.5  0.7857143 0.0004997501
```

## Compare AUCs

To find a confidence interval and p-vaue for the difference and ratio in areas under the curve at $\tau = 28$ days:

```r
aucs <- CompareAUCs(
  time = cic_data$time,
  status = cic_data$status,
  arm = cic_data$arm,
  tau = 28,
  reps = 2000,
  alpha = 0.05
)
show(aucs)
```

```
##   Time     Arm0     Arm1 Contrast Estimate        L        U            P
## 1   28 11.81682 14.60295    A1-A0 2.786127 1.542519 3.940925 0.0004997501
## 2   28 11.81682 14.60295    A1/A0 1.235776 1.122065 1.353300 0.0004997501
```
