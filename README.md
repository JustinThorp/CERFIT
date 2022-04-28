
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CERFIT

<!-- badges: start -->
<!-- badges: end -->

The CERFIT R package is an implementation of Random Forest of
Interaction Trees (RFIT) in R. Which is a modification of the random
forest algorithm to estimate individualized treatment effect in
randomized control trials and observational studies.

## Installation

You can install the development version of CERFIT from
[GitHub](https://github.com/JustinThorp/CERFIT) with:

``` r
# install.packages("devtools")
devtools::install_github("JustinThorp/CERFIT")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(CERFIT)
fit <- CERFIT(y + x_1 + x_2 | t,
              data = data)
ite <- predict(fit,type = "ITE")
importance <- minDepth(fit)
```
