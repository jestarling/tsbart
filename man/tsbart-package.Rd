\name{tsbart-package}
\alias{tsbart-package}
\title{tsbart: BART with Targeted Smoothing R Package}
\description{
Implements BART with Targeted Smoothing (tsBART) as detailed in Starling et al. (2018).

Targeted Smooth BART introduces a new BART prior, which extends BART to allow for smoothing over a selected covariate (t), where the response also depends on covariates (x) which are not smoothed.  Smoothness is induced via Gaussian Process priors on the trees terminal nodes, which are now functions of t instead of scalars.
}
\author{
Jennifer E. Starling,  jstarling@utexas.edu

}



