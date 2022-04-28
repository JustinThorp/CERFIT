
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CERFIT

<!-- badges: start -->
<!-- badges: end -->

The CERFIT R package is an implementation of Random Forest of
Interaction Trees (RFIT) in R. Which is a modification of the random
forest algorithm to estimate individualized treatment effect in
randomized control trials and observational studies. It does this by
modifying the split rule used in the tree growing process. It chooses a
split by maximizing subgroup treatment heterogeneity. It can handle
binary, multiple, ordered, and conintous treatment.

## Installation

You can install the development version of CERFIT from
[GitHub](https://github.com/JustinThorp/CERFIT) with:

``` r
# install.packages("devtools")
devtools::install_github("JustinThorp/CERFIT")
```

## Example

Below is a small example showing how to fit a Random Forest of
Interaction Trees, get predictions from the tree, and calculate variable
importance.

``` r
library(CERFIT)
data <- data.frame(y = rnorm(100),x = rnorm(100),t = rbinom(1000,1,.5))
fit <- CERFIT(y ~ x | t,method = "RCT",data = data,ntrees = 10)
#> [1] "continous response"
#> [1] 1 2
#> [1] "binary Treatment"
#> [1] "Tree Number: 10"
ite <- predict(fit,type = "ITE")
#> [1] 0 1
importance <- MinDepth(fit)
```
