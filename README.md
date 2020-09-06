# Inference for Cumulative Incidence Curves

Zachary McCaw <br>
Updated: 20-09-06



### Description

This package provides functions for inference on the difference in AUCs, event rates, and quantiles comparing two cumulative incidence curves. Also see [MCC](https://github.com/zrmacc/MCC) for inference on mean cumulative count curves. 

## Installation


```r
devtools::install_github(repo = 'zrmacc/CICs')
```

## Examples

### Data

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

### Compare AUCs

To find a confidence interval and p-vaue for the difference and ratio in areas under the curve at $\tau = 28$ days:

```r
aucs <- CompareCICs(
  time = cic_data$time,
  status = cic_data$status,
  arm = cic_data$arm,
  sum_stat = 'AUC',
  param = 28,
  reps = 100,
  alpha = 0.05
)
show(aucs)
```

```
## Stats:
##   Arm   N Stat Observed
## 1   0 600  AUC 11.81682
## 2   1 600  AUC 14.60295
## 
## 
## CIs:
##   Contrast Observed Alpha        L        U
## 1    A1-A0 2.786127  0.05 1.677072 3.968014
## 2    A1/A0 1.235776  0.05 1.132008 1.353295
## 
## 
## P-values:
##   Contrast Method          P
## 1    A1-A0   Boot 0.01980198
## 2    A1/A0   Boot 0.01980198
## 3    A1-A0   Perm 0.01980198
## 4    A1/A0   Perm 0.01980198
```

### Compare Event Rates

To find a confidence interval and p-vaue for the difference and ratio in event rates at $\tau = 28$ days:

```r
rates <- CompareCICs(
  time = cic_data$time,
  status = cic_data$status,
  arm = cic_data$arm,
  sum_stat = 'Rate',
  param = 28,
  reps = 100,
  alpha = 0.05
)
show(rates)
```

```
## Stats:
##   Arm   N Stat  Observed
## 1   0 600 Rate 0.6931547
## 2   1 600 Rate 0.7492792
## 
## 
## CIs:
##   Contrast   Observed Alpha            L         U
## 1    A1-A0 0.05612454  0.05 -0.007629904 0.1200614
## 2    A1/A0 1.08096972  0.05  0.989639079 1.1841944
## 
## 
## P-values:
##   Contrast Method         P
## 1    A1-A0   Boot 0.1386139
## 2    A1/A0   Boot 0.1386139
## 3    A1-A0   Perm 0.1188119
## 4    A1/A0   Perm 0.1188119
```

### Compare Medians

To find a confidence interval and p-vaue for the difference and ratio in medians $q = 0.5$:

```r
quants <- CompareCICs(
  time = cic_data$time,
  status = cic_data$status,
  arm = cic_data$arm,
  sum_stat = 'Quantile',
  param = 0.5,
  reps = 100,
  alpha = 0.05
)
show(quants)
```

```
## Stats:
##   Arm   N     Stat Observed
## 1   0 600 Quantile       15
## 2   1 600 Quantile        9
## 
## 
## CIs:
##   Contrast Observed Alpha    L          U
## 1    A1-A0     -6.0  0.05 -8.0 -3.0000000
## 2    A1/A0      0.6  0.05  0.5  0.7857143
## 
## 
## P-values:
##   Contrast Method          P
## 1    A1-A0   Boot 0.01980198
## 2    A1/A0   Boot 0.01980198
## 3    A1-A0   Perm 0.01980198
## 4    A1/A0   Perm 0.01980198
```
