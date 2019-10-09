# FastDMinR
FastDMinR is an R interface for fast-dm.

## Installation
------------

To install FastDMinR, run:

``` r
install.packages("devtools")
library(devtools)
install_github("AdrianJusepeitis/FastDMinR")
```

## How to use
------------

Simply use the function 

``` r
fast-dm()
```

to specify the input data and the diffusion model. The function returns a list of parameter estimates and cdf values for evaluating model fit. 

### Example
------------

You may try using FastDMinR starting with simple simulated data.

``` r
data = data.frame(
  sub = 1,
  cnd = 1,
  RESPONSE = sample (c(0,1), 100, p=c(0.1,0.9), replace = TRUE),
  TIME = round((rnorm(100,400,30) + rexp(100,0.01))/1000, 2))
```

Load the package and save the output of fast-dm() in an object.
``` r
library(FastDMinR)

dif_data <- fast_dm(data, 
                    Subject = "sub",
                    Conditions = "cnd",
                    TIME = "TIME",
                    RESPONSE = "RESPONSE",
                    precision = 5.0,
                    method = "ks",
                    fix_to = list(p = 0, d = 0, sv = 0, st0 = 0, szr = 0),
                    depend_on_condition = list(a = "cnd"),
                    invariant = c("zr", "v", "t0"))
```

Parameter estimates are now stored in 
```r
dif_data$indiv_estimates
```
and 
```r
dif_data$aggr_estimates
```

Cdf values are stored in 
```r
dif_data$cdf$indiv_cdf
```
and
```r
dif_data$cdf$aggr_cdf
```

It is formated in long format making it easy producing the necessary plots with ggplot2:
```r
library(ggplot2)

ggplot(dif_data$cdf$aggr_cdf, aes(x = RT, y = CDF)) + 
  geom_line(aes(lty = cdf_Type), lwd = 1) + 
  facet_grid(. ~ cnd)
```


