CDC Statistics
---
Jin Zhu

## Installation

You can install the released version of cdcsis from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("cdcsis")
```

## Example

### Feature screening

This is a basic example which shows you how to pick out the important feature from high-dimensional dataset:

``` r
library(cdcsis)
set.seed(1)
num <- 100
p <- 1000
x <- matrix(rnorm(num * p), nrow = num, ncol = p)
z <- rnorm(num)
y <- 3*x[, 1] + 1.5*x[, 2] + 4*z*x[, 5] + rnorm(num)
res <- cdcsis(x, y, z)
head(res[["ix"]], n = 10)
```

**cdcsis** function successfully selects the informative variables from 1000 features pool.

``` r
[1]   1   5   2 628  17  87 912 903 395 630
```

License
----------
GPL (>= 2)