
<!-- README.md is generated from README.Rmd. Please edit that file -->

# VBEMHS

<!-- badges: start -->

<!-- badges: end -->

Variational Bayesian EM algorithm with (independent) horseshoe prior for
variable selection in sparse linear model

## Installation

You can install the development version of VBEMHS like so:

``` r
#install.packages(devtools)
devtools::install_github("wat-sakamoto/VBEMHS")
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo wat-sakamoto/VBEMHS@HEAD
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>          checking for file 'C:\Users\wat-sakamoto\AppData\Local\Temp\Rtmp0kBz9x\remotes503472e84cd7\wat-sakamoto-VBEMHS-c655d53/DESCRIPTION' ...     checking for file 'C:\Users\wat-sakamoto\AppData\Local\Temp\Rtmp0kBz9x\remotes503472e84cd7\wat-sakamoto-VBEMHS-c655d53/DESCRIPTION' ...   ✔  checking for file 'C:\Users\wat-sakamoto\AppData\Local\Temp\Rtmp0kBz9x\remotes503472e84cd7\wat-sakamoto-VBEMHS-c655d53/DESCRIPTION' (782ms)
#>       ─  preparing 'VBEMHS':
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
#>       ─  checking for LF line-endings in source and make files and shell scripts (635ms)
#>   ─  checking for empty or unneeded directories
#>    Omitted 'LazyData' from DESCRIPTION
#>       ─  building 'VBEMHS_0.1.0.tar.gz'
#>      
#> 
#> Installing package into 'C:/Users/wat-sakamoto/AppData/Local/R/win-library/4.5'
#> (as 'lib' is unspecified)
```

## Example

This is an example of application to diabetes data in lars package:

``` r
library(VBEMHS)
library(lars); data(diabetes)
#> Loaded lars 1.3
diabetes.lm <- VBEMHS_lm(y ~ x, data=diabetes)
```

If you want to see the trace for fixed $\tau^2$’s on a grid sequence,
run the following code:

``` r
diabetes.trace <- VBEMHS_lm(y ~ x, data=diabetes, tau_trace=TRUE)
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
&#10;<img src="man/figures/README-pressure-1.png" width="100%" />
&#10;In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN.
-->
