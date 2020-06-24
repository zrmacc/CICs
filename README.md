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

In these data, `arm` is the treatment arm, 0 for reference, 1 for treatment; `time` is the observation time in days; and `status` is the event type, 0 for censoring, 1 for recovery, 2 for death. 

## Difference in event rates

To find a confidence interval and p-vaue for the difference in event rates at $\tau = 28$ days:

```r
ci <- EventRateDiffCI(cic_data$time, cic_data$status, cic_data$arm, tau = 28)
round(ci, digits = 2)
```

```
##  Time  Arm1  Arm0 Delta     L     U 
## 28.00  0.75  0.69  0.06 -0.01  0.12
```

```r
pval <- EventRateDiffP(cic_data$time, cic_data$status, cic_data$arm, tau = 28)
round(pval, digits = 2)
```

```
##  Time  Arm1  Arm0 Delta     P 
## 28.00  0.75  0.69  0.06  0.09
```

## Difference in medians

To find a confidence interval and p-vaue for the difference in medians $q = 0.5$:

```r
ci <- QuantDiffCI(cic_data$time, cic_data$status, cic_data$arm, q = 0.5)
round(ci, digits = 2)
```

```
## Quantile     Arm1     Arm0    Delta        L        U 
##      0.5      9.0     15.0     -6.0     -8.0     -3.0
```

```r
pval <- QuantDiffP(cic_data$time, cic_data$status, cic_data$arm, q = 0.5)
round(pval, digits = 2)
```

```
## Quantile     Arm1     Arm0    Delta        P 
##      0.5      9.0     15.0     -6.0      0.0
```

## Difference in areas under the curve

To find a confidence interval and p-vaue for the difference in areas under the curve at $\tau = 28$ days:

```r
ci <- AUCDiffCI(cic_data$time, cic_data$status, cic_data$arm, tau = 28)
round(ci, digits = 2)
```

```
##  Time  Arm1  Arm0 Delta     L     U 
## 28.00 14.60 11.82  2.79  1.54  3.97
```

```r
pval <- AUCDiffP(cic_data$time, cic_data$status, cic_data$arm, tau = 28)
round(pval, digits = 2)
```

```
##  Time  Arm1  Arm0 Delta     P 
## 28.00 14.60 11.82  2.79  0.00
```
