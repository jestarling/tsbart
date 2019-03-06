# tsbart R Package

## About
Implements the BART with Targeted Smoothing method (tsBART) detailed in Starling et al. (2018).  Also implements the Targeted Smooth Bayesian Causal Forest method (tsbcf) detailed in Starling et al. (2019).

BART with Targeted Smoothing, or tsBART, is an extension of the original BART model which induces smoothness over a single targeted covariate. 

Targeted Smooth BCF (tsbcf) is an extension of Bayesian Causal FOrests and tsBART, which allows for heterogeneous treatment effects which vary smoothly across a target covariate.

## Documentation

Function documentation is available at https://jestarling.github.io/tsbart/.

Vignette is available at https://jestarling.github.io/tsbart/articles/myvignette.html.

## Installation

The tsbart package can be installed as follows.
```
library(devtools)

install_github("jestarling/tsbart")
```

## Author

Jennifer E. Starling  
jstarling@utexas.edu  
https://jestarling.github.io  

