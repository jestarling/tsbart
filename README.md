# tsbart R Package

## About
Implements the BART with Targeted Smoothing method (tsBART) detailed in Starling et al. (2019).  BART with Targeted Smoothing, or tsBART, is an extension of the original BART model which induces smoothness over a single targeted covariate.  

Includes Projective Smooth BART options (for monotonicity and rounded responses), as detailed in Starling et. al (2019b), *Monotone function estimation in the presence of extreme data coarsening: Analysis of preeclampsia and birth weight in urban Uganda*.

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

