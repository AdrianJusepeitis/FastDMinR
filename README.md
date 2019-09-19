# FastDMinR
FastDMinR is an R interface for fast-dm.

Installation
------------

To install FastDMinR, run:

``` r
install.packages("devtools")
library(devtools)
install_github("AdrianJusepeitis/FastDMinR")
```

How to use
------------

Simply use the function 

``` r
fast-dm()
```

to specify the input data and the diffusion model. The function returns a list of parameter estimates and cdf values for evaluating model fit. 
