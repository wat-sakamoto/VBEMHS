
<!-- README.md is generated from README.Rmd. Please edit that file -->

# VBEMHS

<!-- badges: start -->

<!-- badges: end -->

Variational Bayesian Variable Selection and Inference in Bayesian Sparse
Linear Models

## Installation

You can install the development version of VBEMHS like so:

``` r
#install.packages(pak)
pak::pkg_install("wat-sakamoto/VBEMHS")
```

## Example

Below is an example of application to `diabetes` dataset in `lars`
package. The default is using the independent prior.

``` r
library(VBEMHS)
library(lars); data(diabetes)
diabetes.lm <- VBEMHS_lm(y ~ x, data=diabetes)
```

If you use correction for underestimation of the error variance, use the
option `s2_correct=TRUE`. This is recommended when $n<p$.

``` r
diabetes.VBEMHS_lm <- VBEMHS_lm(y ~ x, data=diabetes, s2_correct=TRUE)
```

If you obtain a trace along the global shrinkage coefficients on a grid
sequence, run the following code:

``` r
diabetes.VBEMHS_lm.trace <- VBEMHS_lm(y ~ x, data=diabetes, gs_trace=TRUE)
```

If you use the conjugate prior, with sigma2 estimated by posterior mode,
run the following code:

``` r
diabetes.VBEMHSc_lm <- VBEMHS_lm(y ~ x, data=diabetes,
      prior=VBEMHS_lm.prior(independent=FALSE),
      control=VBEMHS_lm.control(s2_estim="mode"))
```

<!-- 
What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: 
&#10;
``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```
&#10;You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this.
&#10;You can also embed plots, for example:
&#10;<img src="man/figures/README-pressure-1.png" alt="" width="100%" />
&#10;In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN.
-->
